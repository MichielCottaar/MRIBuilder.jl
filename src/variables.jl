"""
Defines the functions that can be called on parts of an MRI sequence to query or constrain any variables.

In addition this defines:
- [`variables`](@ref): dictionary containing all variables.
- [`VariableType`](@ref): parent type for any variables (whether number or JuMP variable).
- [`get_free_variable`](@ref): helper function to create new JuMP variables.
- [`VariableNotAvailable`](@ref): error raised if variable is not defined for specific [`AbstractBlock`](@ref).
- [`set_simple_constraints!`](@ref): call [`apply_simple_constraint!`](@ref) for each keyword argument.
- [`apply_simple_constraint!`](@ref): set a simple equality constraint.
- [`get_pulse`](@ref)/[`get_gradient`](@ref)/[`get_readout`](@ref): Used to get the pulse/gradient/readout part of a building block
- [`gradient_orientation`](@ref): returns the gradient orientation of a waveform if fixed.
"""
module Variables
import JuMP: @constraint, @variable, Model, @objective, objective_function, AbstractJuMPScalar
import StaticArrays: SVector
import MacroTools
import ..Scanners: gradient_strength, slew_rate, Scanner
import ..BuildSequences: global_model, global_scanner, fixed

"""
Parent type of all components, building block, and sequences that form an MRI sequence.
"""
abstract type AbstractBlock end

function fixed(ab::AbstractBlock)
    params = []
    for prop_name in propertynames(ab)
        push!(params, fixed(getproperty(ab, prop_name)))
    end
    return typeof(ab)(params...)
end


"""
    adjust_internal(block, names_used; kwargs...)

Returns the adjusted blocks and add any keywords used in the process to `names_used`.

This is a helper function used by `adjust`. 
"""
function adjust_internal end

"""
    adjustable(block)

Returns whether a sequence, building block, or component can be adjusted

Can return one of:
- `:false`: not adjustable
- `:gradient`: expects gradient adjustment parameters
- `:pulse`: expects RF pulse adjustment parameters
"""
adjustable(::AbstractBlock) = :false


abstract type AnyVariable end

"""
A sequence property that can be constrained and/or optimised.

It acts as a function, so you can call it on a sequence or building block to get the actual values (e.g., `v(sequence)`).
It can return one of the following:
- a number
- a vector of number
- a NamedTuple with the values for individual sequence components
"""
mutable struct Variable <: AnyVariable
    name :: Symbol
    f :: Function
    getter :: Union{Nothing, Function}
end

struct AlternateVariable <: AnyVariable
    name :: Symbol
    other_var :: Symbol
    to_other :: Function
    from_other :: Function
    inverse :: Bool
end


"""
    variable_defined_for(var, Val(type))

Check whether variable is defined for a specific sub-type.
"""
variable_defined_for(var::Variable, ::Val{T}) where {T <: AbstractBlock} = hasmethod(var.f, (T, ))

struct _Variables
    variables :: Dict{Symbol, AnyVariable}
end

variables = _Variables(Dict{Symbol, AnyVariable}())

Base.getindex(v::_Variables, i::Symbol) = getfield(v, :variables)[i]
Base.keys(v::_Variables) = keys(getfield(v, :variables))

Base.propertynames(v::_Variables) = Tuple(keys(getfield(v, :variables)))

Base.getproperty(v::_Variables, s::Symbol) = v[s]


macro defvar(func_def) 
    return _defvar(func_def, nothing)
end

macro defvar(getter, func_def) 
   return _defvar(func_def, getter)
end

function _defvar(func_def, getter=nothing)
    func_names = []


    if getter isa Symbol
        getter_dict = Dict(
            :pulse => get_pulse,
            :gradient => get_gradient,
            :pathway => get_pathway,
            :readout => get_readout,
        )
        if !(getter in keys(getter_dict))
            error("label in `@defvar <label> <statement>` should be one of `pulse`/`gradient`/`pathway`/`readout`, not `$getter`")
        end
        getter = getter_dict[getter]
    end

    function adjust_function(ex)
        if ex isa Expr && ex.head == :block
            return Expr(:block, adjust_function.(ex.args)...)
        end
        if ex isa Expr && ex.head == :function && length(ex.args) == 1
            push!(func_names, ex.args[1])
            return :nothing
        end
        try
            fn_def = MacroTools.splitdef(ex)
            push!(func_names, fn_def[:name])
            new_def = Dict{Symbol, Any}()
            new_def[:name] = Expr(:., fn_def[:name], QuoteNode(:f))
            new_def[:args] = esc.(fn_def[:args])
            new_def[:kwargs] = esc.(fn_def[:kwargs])
            new_def[:body] = esc(fn_def[:body])
            new_def[:whereparams] = esc.(fn_def[:whereparams])
            return MacroTools.combinedef(new_def)
        catch e
            if e isa AssertionError
                return ex
            end
            rethrow()
        end
    end
    new_func_def = adjust_function(func_def)

    function fix_function_name(ex)
        if ex in func_names
            return esc(ex)
        else
            return ex
        end
    end
    new_func_def = MacroTools.postwalk(fix_function_name, new_func_def)

    expressions = Expr[]
    for func_name in func_names
        push!(expressions, quote
            $(esc(func_name)) = if $(QuoteNode(func_name)) in keys(variables)
                variables[$(QuoteNode(func_name))]
            else
                function $(func_name) end
                getfield(variables, :variables)[$(QuoteNode(func_name))] = Variable($(QuoteNode(func_name)), $(func_name), $getter)
            end
            if $(esc(func_name)) isa AlternateVariable
                error("$($(esc(func_name)).name) is defined through $($(esc(func_name)).other_var). Please define that variable instead.")
            end
            if !isnothing($getter) && $(esc(func_name)).getter != $getter
                if isnothing($(esc(func_name)).getter)
                    $(esc(func_name)).getter = $getter
                else
                    name = $(esc(func_name)).name
                    error("$(name) is already defined as a variable for $($(esc(func_name)).getter). Cannot switch to $($getter).")
                end
            end
        end
        )
    end
    args = vcat([e.args for e in expressions]...)
    return Expr(
        :block,
        args...,
        new_func_def
    )
end

@defvar function duration end
"""
    duration(block)

Duration of the sequence or building block in ms.
""" 
duration


function def_alternate_variable!(name::Symbol, other_var::Symbol, to_other::Function, from_other::Function, inverse::Bool)
    getfield(variables, :variables)[name] = AlternateVariable(name, other_var, to_other, from_other, inverse)
end

def_alternate_variable!(:spoiler_scale, :spoiler_scale, q->1e-3 * 2π/q, l->1e-3 * 2π/l, true)

for name in [:slice_thickness, :bandwidth, :fov, :voxel_size]
    inv_name = Symbol("inverse_" * string(name))
    def_alternate_variable!(name, inv_name, inv, inv, true)
end

for (name, alt_name) in [
    (:TE, :echo_time),
    (:TR, :repetition_time),
    (:Δ, :diffusion_time),
]
    def_alternate_variable!(name, alt_name, identity, identity, false)
end


"""
Parent type for any variable in the MRI sequence.

Each variable can be one of:
- a new JuMP variable
- an expression linking this variable to other JuMP variable
- a number

Create these using [`get_free_variable`](@ref).
"""
const VariableType = Union{Number, AbstractJuMPScalar}


"""
    get_free_variable(value; integer=false, start=0.01)

Get a representation of a given `variable` given a user-defined constraint.

The result is guaranteed to be a [`VariableType`](@ref).
"""
get_free_variable(value::Number; integer=false, kwargs...) = integer ? Int(value) : Float64(value)
get_free_variable(value::VariableType; kwargs...) = value
get_free_variable(::Nothing; integer=false, start=0.01) = @variable(global_model(), start=start, integer=integer)
get_free_variable(value::Symbol; integer=false, kwargs...) = integer ? error("Cannot maximise or minimise an integer variable") : get_free_variable(Val(value); kwargs...)
function get_free_variable(::Val{:min}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    model = global_model()
    @objective model Min objective_function(model) + var
    return var
end
function get_free_variable(::Val{:max}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    model = global_model()
    @objective model Min objective_function(model) - var
    return var
end

"""
    get_pulse(building_block)

Get the pulse played out during the building block.

Any `pulse` variables not explicitly defined for this building block will be passed on to the pulse.
"""
function get_pulse end

"""
    get_gradient(building_block)

Get the gradient played out during the building block.

Any `gradient` variables not explicitly defined for this building block will be passed on to the gradient.
"""
function get_gradient end

"""
    get_readout(building_block)

Get the readout played out during the building block.

Any `readout` variables not explicitly defined for this building block will be passed on to the readout.
"""
function get_readout end

"""
    get_pathway(sequence)

Get the default spin pathway for the sequence.

Any `pathway` variables not explicitly defined for this building block will be passed on to the pathway.
"""
function get_pathway end


"""
    gradient_orientation(building_block)

Returns the gradient orientation.
"""
function gradient_orientation end


@defvar function effective_time end


function (var::Variable)(block::AbstractBlock, args...; kwargs...)
    if !applicable(var.f, block, args...) && !isnothing(var.getter)
        apply_to = var.getter(block)
        if apply_to isa AbstractBlock
            return var(apply_to, args...; kwargs...)
        elseif apply_to isa NamedTuple
            return NamedTuple(k => var(v, args...; kwargs...) for (k, v) in pairs(apply_to))
        elseif apply_to isa AbstractVector{<:AbstractBlock} || apply_to isa Tuple
            return var.(apply_to, args...; kwargs...)
        end
    end
    return var.f(block, args...; kwargs...)
end

function (var::AlternateVariable)(args...; kwargs...)
    other_var = variables[var.other_var]
    apply_from_other(res::Number) = var.from_other(res)
    apply_from_other(res::AbstractArray{<:Number}) = var.from_other.(res)
    apply_from_other(res::NamedTuple) = NamedTuple(k => apply_from_other(v) for (k, v) in pairs(res))
    return apply_from_other(other_var(args...; kwargs...))
end


for base_fn in [:qval, :gradient_strength, :slew_rate]
    fn3 = Symbol(String(base_fn) * "3")
    @eval function $fn3(bb::AbstractBlock, args...; kwargs...)
        if hasmethod(get_gradient, (typeof(bb), ))
            return $fn3(get_gradient(bb), args...; kwargs...)
        else
            value = $base_fn(bb, args...; kwargs...)
            if value isa Number && iszero(value)
                return zero(SVector{3, Float64})
            elseif value isa AbstractVector
                return value
            else
                return value .* gradient_orientation(bb)
            end
        end
    end
    @eval $fn3(nt::NamedTuple, args...; kwargs...) = map(v -> $fn3(v, args...; kwargs...), nt)
end


"""
    set_simple_constraints!(block, kwargs)

Add any constraints or objective functions to the variables of a [`AbstractBlock`](@ref).

Each keyword argument has to match one of the functions in [`variables`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(block::AbstractBlock, kwargs)
    real_kwargs = NamedTuple(key => value for (key, value) in kwargs if !isnothing(value))

    for (key, value) in real_kwargs
        var = variables[key]
        if var isa AlternateVariable
            if var.other_var in real_kwargs
                error("Set constraints on both $key and $(var.other_var), however they are equivalent.")
            end
            invert_value(value::VariableType) = var.to_other(value)
            invert_value(value::Symbol) = invert_value(Val(value))
            invert_value(::Val{:min}) = var.inverse ? Val(:max) : Val(:min)
            invert_value(::Val{:max}) = var.inverse ? Val(:min) : Val(:max)
            invert_value(value::AbstractVector) = invert_value.(value)
            apply_simple_constraint!(variables[var.other_var](block), invert_value(value))
        else
            apply_simple_constraint!(var(block), value)
        end
    end
    nothing
end

"""
    apply_simple_constraint!(variable, value)

Add a single constraint or objective to the `variable`.

`value` can be one of:
- `nothing`: do nothing
- `:min`: minimise the variable
- `:max`: maximise the variable
- `number`: fix variable to this value
- `equation`: fix variable to the result of this equation
"""
apply_simple_constraint!(variable::AbstractVector, value::Symbol) = apply_simple_constraint!(sum(variable), Val(value))
apply_simple_constraint!(variable::VariableType, value::Symbol) = apply_simple_constraint!(variable, Val(value))
apply_simple_constraint!(variable::VariableType, ::Val{:min}) = @objective global_model() Min objective_function(global_model()) + variable
apply_simple_constraint!(variable::VariableType, ::Val{:max}) = @objective global_model() Min objective_function(global_model()) - variable
apply_simple_constraint!(variable::VariableType, value::VariableType) = @constraint global_model() variable == value
apply_simple_constraint!(variable::AbstractVector, value::AbstractVector) = [apply_simple_constraint!(v1, v2) for (v1, v2) in zip(variable, value)]
apply_simple_constraint!(variable::AbstractVector, value::VariableType) = [apply_simple_constraint!(v1, value) for v1 in variable]
apply_simple_constraint!(variable::Number, value::Number) = @assert variable ≈ value "Variable set to multiple incompatible values."
function apply_simple_constraint!(variable::NamedTuple, value)
    for sub_var in variable
        apply_simple_constraint!(sub_var, value)
    end
end
function apply_simple_constraint!(variable::NamedTuple, value::NamedTuple)
    for key in keys(value)
        apply_simple_constraint!(variable[key], value[key])
    end
end


"""
    make_generic(sequence/building_block/component)

Returns a generic version of the `BaseSequence`, `BaseBuildingBlock`, or `BaseComponent`

- Sequences are all flattened and returned as a single `Sequence` containing only `BuildingBlock` objects.
- Any `BaseBuildingBlock` is converted into a `BuildingBlock`.
- Pulses are replaced with `GenericPulse` (except for instant pulses).
- Instant readouts are replaced with `ADC`.
"""
function make_generic end


"""
    scanner_constraints!(block)

Constraints [`gradient_strength`](@ref) and [`slew_rate`](@ref) to be less than the [`global_scanner`](@ref) maximum.
"""
function scanner_constraints!(bb::AbstractBlock)
    try
        global_scanner()
    catch e
        return
    end
    for f in (slew_rate, gradient_strength)
        value = nothing
        try
            value = f(bb)
        catch e
                continue
                continue
            else
                rethrow()
            end
            continue
            else
                rethrow()
            end
        end
        if value isa AbstractVector
            for v in value
                @constraint global_model() v <= f(global_scanner())
                @constraint global_model() v >= -f(global_scanner())
            end
        else
            @constraint global_model() value <= f(global_scanner())
            @constraint global_model() value >= -f(global_scanner())
        end
    end
end

end