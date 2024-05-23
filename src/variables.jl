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
struct Variable <: AnyVariable
    name :: Symbol
    f :: Function
end

"""
    variable_defined_for(var, Val(type))

Check whether variable is defined for a specific sub-type.
"""
variable_defined_for(var::Variable, ::Val{T}) where {T <: AbstractBlock} = hasmethod(var.f, (T, ))

all_variables = Dict{Symbol, AnyVariable}()

macro defvar(func_def)
    func_names = []

    function adjust_function(ex)
        try
            fn_def = MacroTools.splitdef(ex)
            push!(func_names, fn_def[:name])
            fn_def[:name] = Expr(:., fn_def[:name], QuoteNode(:f))
            return MacroTools.combinedef(fn_def)
        catch e
            if e isa AssertionError
                return ex
            end
            rethrow()
        end
    end
    new_func_def = MacroTools.prewalk(adjust_function, func_def)

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
            $(esc(func_name)) = if $(QuoteNode(func_name)) in keys(all_variables)
                all_variables[$(QuoteNode(func_name))]
            else
                function $(func_name) end
                all_variables[$(QuoteNode(func_name))] = Variable($(QuoteNode(func_name)), $(func_name))
            end
        end
        )
    end
    return Expr(
        :block,
        expressions...,
        new_func_def
    )
end

@defvar function duration end
@doc """    duration(block)

Duration of the sequence or building block in ms.
""" duration


struct AlternateVariable <: AnyVariable
    name :: Symbol
    other_var :: Symbol
    to_other :: Function
    from_other :: Function
    inverse :: Bool
end

all_variables[:qval] = AlternateVariable(:qval, :qval_square, q->q^2, sqrt, false)
all_variables[:spoiler_scale] = AlternateVariable(:spoiler_scale, :spoiler_scale, q->1e-3 * 2π/q, l->1e-3 * 2π/l, true)

for name in [:slice_thickness, :bandwidth, :fov, :voxel_size]
    inv_name = Symbol("inverse_" * string(name))
    all_variables[name] = AlternateVariable(name, inv_name, inv, inv, true)
end

for (name, alt_name) in [
    (:TE, :echo_time),
    (:TR, :repetition_time),
    (:Δ, :diffusion_time),
]
    all_variables[name] = AlternateVariable(name, alt_name, identity, identity, false)
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
    get_pulse(building_block)]

Get the pulse played out during the building block.

Any `pulse` variables not explicitly defined for this building block will be passed on to the pulse.
"""
function get_pulse end

"""
    get_gradient(building_block)]

Get the gradient played out during the building block.

Any `gradient` variables not explicitly defined for this building block will be passed on to the gradient.
"""
function get_gradient end

"""
    get_readout(building_block)]

Get the readout played out during the building block.

Any `readout` variables not explicitly defined for this building block will be passed on to the readout.
"""
function get_readout end

"""
    bmat_gradient(gradient::GradientBlock, qstart=(0, 0, 0))

Computes the diffusion-weighting matrix due to a single gradient block in rad^2 ms/um^2.

This should be defined for every `GradientBlock`, but not be called directly.
Instead, the `bmat` and `bval` should be constrained for specific `Pathways`
"""
function bmat_gradient end


"""
    gradient_orientation(building_block)

Returns the gradient orientation.
"""
function gradient_orientation end


function effective_time end


"""
    VariableNotAvailable(building_block, variable, alt_variable)

Exception raised when a variable function does not support a specific `AbstractBlock`.
"""
mutable struct VariableNotAvailable <: Exception
    bb :: Type{<:AbstractBlock}
    variable :: Function
    alt_variable :: Union{Nothing, Function}
end
VariableNotAvailable(bb::Type{<:AbstractBlock}, variable::Function) = VariableNotAvailable(bb, variable, nothing)



function Base.showerror(io::IO, e::VariableNotAvailable)
    if isnothing(e.alt_variable)
        print(io, e.variable, " is not available for block of type ", e.bb, ".")
    else
        print(io, e.variable, " is not available for block of type ", e.bb, ". Please use ", e.alt_variable, " instead to set any contsraints or objective functions.")
    end
end


for (target_name, all_vars) in all_variables_symbols
    for (variable_func, _) in all_vars
        if variable_func in [:qval3, :TR, :TE, :Δ]
            continue
        end
        get_func = Symbol("get_" * string(target_name))
        use_get_func = target_name in (:pulse, :readout, :gradient)
        @eval function Variables.$variable_func(bb::AbstractBlock)
            try
                if Variables.$variable_func in keys(alternative_variables)
                    alt_var, forward, backward, _ = alternative_variables[Variables.$variable_func]
                    try
                        value = alt_var(bb)
                        if value isa Number
                            return backward(value)
                        elseif value isa AbstractArray{<:Number}
                            return backward.(value)
                        end
                    catch e
                        if e isa VariableNotAvailable
                            throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
                        end
                        rethrow()
                    end
                    throw(VariableNotAvailable(typeof(bb), Variables.$variable_func, alt_var))
                end
                throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
            catch e
                if $use_get_func && e isa VariableNotAvailable && hasmethod($get_func, tuple(typeof(bb)))
                    apply_to = try
                        $(get_func)(bb)
                    catch
                        throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
                    end
                    if apply_to isa AbstractBlock
                        return Variables.$variable_func(apply_to)
                    elseif apply_to isa NamedTuple
                        return NamedTuple(k => Variables.$variable_func(v) for (k, v) in pairs(apply_to))
                    elseif apply_to isa AbstractVector{<:AbstractBlock} || apply_to isa Tuple
                        return Variables.$variable_func.(apply_to)
                    end
                    error("$($(use_get_func)) returned an unexpected type for $(bb).")
                end
                rethrow()
            end
        end
    end
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
    for (key, value) in kwargs
        if variables[key] in keys(alternative_variables)
            alt_var, forward, backward, to_invert = alternative_variables[variables[key]]
            invert_value(value::VariableType) = forward(value)
            invert_value(value::Symbol) = invert_value(Val(value))
            invert_value(::Val{:min}) = to_invert ? Val(:max) : Val(:min)
            invert_value(::Val{:max}) = to_invert ? Val(:min) : Val(:max)
            invert_value(value::AbstractVector) = invert_value.(value)
            invert_value(value) = value
            try
                apply_simple_constraint!(alt_var(block), invert_value(value))
                continue
            catch e
                if !(e isa VariableNotAvailable)
                    rethrow()
                end
            end
        end
        apply_simple_constraint!(variables[key](block), value)
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
apply_simple_constraint!(variable::VariableType, ::Nothing) = nothing
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
            if e isa VariableNotAvailable
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