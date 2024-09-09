"""
Defines the functions that can be called on parts of an MRI sequence to query or constrain any variables.

In addition this defines:
- [`variables`](@ref): module containing all variables.
- [`VariableType`](@ref): parent type for any variables (whether number or JuMP variable).
- [`get_free_variable`](@ref): helper function to create new JuMP variables.
- [`add_cost_function!`](@ref): add a specific term to the model cost functions.
- [`set_simple_constraints!`](@ref): call [`apply_simple_constraint!`](@ref) for each keyword argument.
- [`apply_simple_constraint!`](@ref): set a simple equality constraint.
- [`get_pulse`](@ref)/[`get_gradient`](@ref)/[`get_readout`](@ref): Used to get the pulse/gradient/readout part of a building block
- [`gradient_orientation`](@ref): returns the gradient orientation of a waveform if fixed.
"""
module Variables
import JuMP: @constraint, @variable, Model, @objective, objective_function, AbstractJuMPScalar, QuadExpr, AffExpr
import StaticArrays: SVector
import MacroTools
import ..Scanners: gradient_strength, slew_rate, Scanner
import ..BuildSequences: global_scanner, fixed, GLOBAL_MODEL

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

This is a helper function used by [`adjust`](@ref). 
It should be defined for any block that is adjustable (as defined by [`adjust_groups`](@ref)).
"""
function adjust_internal end

"""
    adjust_groups(block)

Returns an array of keywords in [`adjust`](@ref) that should affect a specfic block.

If any of these keywords are present in [`adjust`](@ref), then [`adjust_internal`](@ref) will be called.

Some standard keywords are:
- `:gradient`: expects gradient adjustment parameters
- `:pulse`: expects RF pulse adjustment parameters
"""
adjust_groups(::AbstractBlock) = Symbol[]

# further defined in post_hoc.jl
function adjust end


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
    from_other :: Function
    to_other :: Union{Nothing, Function}
    inverse :: Bool
end


"""
    base_variables([T])

Return dictionary with all [`Variable`](@ref) objects defined for a specific sequence component/block `T`.

This only returns those [`Variable`](@ref) directly defined for this component/block, not for any sub-components (through [`get_pulse`](@ref), [`get_gradient`][(@ref), etc.)

If `T` is not provided, all [`Variable`](@ref) objects are returned.
"""
function base_variables()
    all_members = (s => getproperty(variables, s) for s in names(variables, all=true))
    return Dict{Symbol, Variable}(
        s => v for (s, v) in all_members if v isa Variable
    )
end

function base_variables(T::Type{<:AbstractBlock})
    return Dict(
        s => v for (s, v) in base_variables() if hasmethod(v.f, (T, ))
    )
end


# Add Variables to the individual sequence components/blocks properties
function Base.propertynames(::T) where {T <: AbstractBlock}
    f = Base.fieldnames(T)
    var_names = [k for k in keys(base_variables(T)) if !(k in f)]
    return (f..., var_names...)
end

function Base.getproperty(block::T, v::Symbol) where T <: AbstractBlock
    if v in Base.fieldnames(T)
        return getfield(block, v)
    end
    vars = base_variables(T)
    if v in keys(vars)
        return vars[v](block)
    end
    error("Type $(T) has no field or variable $(v)")
end

function Base.setproperty!(block::T, v::Symbol, value) where T <: AbstractBlock
    if v in Base.fieldnames(T)
        return setfield!(block, v, value)
    end
    vars = base_variables(T)
    if v in keys(vars)
        orig = vars[v](block)
        return apply_simple_constraint!(orig, value)
    end
    error("Type $(T) has no field or variable $(v)")
end

"""
    variable_defined_for(var, Val(type))

Check whether variable is defined for a specific sub-type.
"""
variable_defined_for(var::Variable, ::Val{T}) where {T <: AbstractBlock} = hasmethod(var.f, (T, ))

"""
Main module containing all the MRIBuilder sequence variables.

All variables are available as members of this module, e.g.
`variables.echo_time` returns the echo time variable.
New variables can be defined using `@defvar`.

Set constraints on variables by passing them on as keywords during the sequence generation,
e.g., `seq=SpinEcho(echo_time=70)`.

After sequence generation you can get the variable values by calling
`variables.echo_time(seq)`.
For the sequence defined above this would return 70. (or a number very close to that).
"""
baremodule variables
end


"""
    @defvar([getter, ], function(s))

Defines new [`variables`](@ref).

Each variable is defined as regular Julia functions embedded within a `@defvar` macro.
For example, to define a `variables.echo_time` variable for a `SpinEcho` sequence, one can use:
```julia
@defvar echo_time(ge::SpinEcho) = 2 * (variables.effective_time(ge, :refocus) - variables.effective_time(ge, :excitation))
```

Multiple variables can be defined in a single `@defvar` by including them in a code block:
```julia
@defvar begin
    function var1(seq::SomeSequenceType)
        ...
    end
    function var2(seq::SomeSequenceType)
        ...
    end
end
```

Before the variable function definitions one can include a `getter`.
This `getter` defines the type of the sequence component for which the variables will be defined.
If the variable is not defined for the sequence, the variable will be extracted for those type of sequence components instead.
The following sequence component types are provided:
- `pulse`: use [`get_pulse`](@ref)
- `gradient`: use [`get_gradient`](@ref)
- `readout`: use [`get_readout`](@ref)
- `pathway`: use [`get_pathway`](@ref)
e.g. the following defines a `flip_angle` variable, which is marked as a property of an RF pulse.
```julia
@defvar pulse flip_angle(...) = ...
```
If after this definition, `flip_angle` is not explicitly defined for any sequence, it will be extracted for the RF pulses in that sequence instead.
"""
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
            new_def[:name] = Expr(:., Expr(:., :variables, QuoteNode(fn_def[:name])), QuoteNode(:f))
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
            if !($(QuoteNode(func_name)) in names(variables; all=true))
                function $(func_name) end
                variables.$(func_name) = Variable($(QuoteNode(func_name)), $(func_name), $getter)
            end
            if variables.$(func_name) isa AlternateVariable
                error("$($(esc(func_name)).name) is defined through $(variables.$(func_name).other_var). Please define that variable instead.")
            end
            if !isnothing($getter) && variables.$(func_name).getter != $getter
                if isnothing(variables.$(func_name).getter)
                    variables.$(func_name).getter = $getter
                else
                    name = variables.$(func_name).name
                    error("$(name) is already defined as a variable for $(variables.$(func_name).getter). Cannot switch to $($getter).")
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
variables.duration


function def_alternate_variable!(name::Symbol, other_var::Symbol, from_other::Function, to_other::Union{Nothing, Function}, inverse::Bool)
    setproperty!(variables, name, AlternateVariable(name, other_var, from_other, to_other, inverse))
end

def_alternate_variable!(:spoiler, :qval, q->1e-3 * 2π/q, l->1e-3 * 2π/l, true)
def_alternate_variable!(:qval, :qval_square, sqrt, q -> q * q, false)
def_alternate_variable!(:qval_square, :qvec, qv -> sum(q -> q * q, qv), nothing, false)

"""
    qval(gradient)

The norm of the [`variables.qvec`](@ref).
"""
variables.qval

"""
    spoiler(gradient)

Spatial scale in mm over which the spoiler gradient will dephase by 2π.

Automatically computed based on [`variables.qvec`](@ref).
"""
variables.spoiler


for vec_variable in [:gradient_strength, :slew_rate]
    vec_square = Symbol(string(vec_variable) * "_square")
    vec_norm = Symbol(string(vec_variable) * "_norm")
    def_alternate_variable!(vec_norm, vec_square, sqrt, v -> v * v, false)
    def_alternate_variable!(vec_square, vec_variable, v -> v[1] * v[1] + v[2] * v[2] + v[3] * v[3], nothing, false)
end

"""
    gradient_strength_norm(gradient)

The norm of the [`variables.gradient_strength`](@ref).
"""
variables.gradient_strength_norm

"""
    slew_rate_norm(gradient)

The norm of the [`variables.slew_rate`](@ref).
"""
variables.slew_rate_norm

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
    TE(sequence)

Alternative name to compute the [`variables.echo_time`](@ref) of a sequence in ms.
"""
variables.TE

"""
    TR(sequence)

Alternative name to compute the [`variables.repetition_time`](@ref) of a sequence in ms.
"""
variables.TR

"""
    Δ(sequence)

Alternative name to compute the [`variables.diffusion_time`](@ref) of a sequence in ms.
"""
variables.Δ


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
get_free_variable(::Nothing; integer=false, start=0.01) = @variable(GLOBAL_MODEL[][1], start=start, integer=integer)
get_free_variable(value::Symbol; integer=false, kwargs...) = integer ? error("Cannot maximise or minimise an integer variable") : get_free_variable(Val(value); kwargs...)
function get_free_variable(::Val{:min}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    add_cost_function!(var, 1)
    return var
end
function get_free_variable(::Val{:max}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    add_cost_function!(-var, 1)
    return var
end

"""
    get_pulse(sequence)

Get the main pulses played out during the sequence.

This has to be defined for individual sequences to work.

Any `pulse` variables not explicitly defined for this sequence will be passed on to the pulse.
"""
function get_pulse end

"""
    get_gradient(sequence)

Get the main gradients played out during the sequence.

This has to be defined for individual sequences to work.

Any `gradient` variables not explicitly defined for this sequence will be passed on to the gradient.
"""
function get_gradient end

"""
    get_readout(sequence)

Get the main readout events played out during the sequence.

This has to be defined for individual sequences to work.

Any `readout` variables not explicitly defined for this sequence will be passed on to the readout.
"""
function get_readout end

"""
    get_pathway(sequence)

Get the default spin pathway(s) for the sequence.

This has to be defined for individual sequences to work.

Any `pathway` variables not explicitly defined for this sequence will be passed on to the pathway.
"""
function get_pathway end


"""
    gradient_orientation(building_block)

Returns the gradient orientation.
"""
function gradient_orientation end


function (var::Variable)(block::AbstractBlock, args...; kwargs...)
    if !applicable(var.f, block, args...) && !isnothing(var.getter)
        apply_to = var.getter(block)
        if apply_to isa AbstractBlock
            return var(apply_to, args...; kwargs...)
        elseif apply_to isa NamedTuple
            return NamedTuple(k => var(v, args...; kwargs...) for (k, v) in pairs(apply_to))
        elseif apply_to isa AbstractVector{<:AbstractBlock} || apply_to isa Tuple
            return var.(apply_to, args...; kwargs...)
        else
            error("$(var.getter) returned an unexpected type: $(typeof(apply_to)).")
        end
    end
    return var.f(block, args...; kwargs...)
end

# Special case for BuildingBlock events
function (var::Variable)(event::Tuple{<:VariableType, <:AbstractBlock}, args...; kwargs...)
    if applicable(var.f, event, args...; kwargs...)
        return var.f(event, args...; kwargs...)
    end
    # falling back to just processing the `AbstractBlock`
    return var(event[2], args...; kwargs...)
end

function (var::AlternateVariable)(args...; kwargs...)
    other_var = getproperty(variables, var.other_var)
    apply_from_other(res::VariableType) = var.from_other(res)
    function apply_from_other(res::AbstractVector{<:VariableType}) 
        try
            return var.from_other(res)
        catch e
            if e isa MethodError
                return var.from_other.(res)
            end
        end
    end
    apply_from_other(res::NamedTuple) = NamedTuple(k => apply_from_other(v) for (k, v) in pairs(res))
    return apply_from_other(other_var(args...; kwargs...))
end


"""
    add_cost_function!(function, level=2)

Adds an additional term to the cost function.

This term will be minimised together with any other terms in the cost function.
Terms added at a lower level will be optimised before any terms with a higher level.

By default, the term is added to the `level=2`, which is appropriate for any cost functions added by the developer,
which will generally only be optimised after any user-defined cost functions (which are added at `level=1` by [`add_simple_constraint!`](@ref) or [`set_simple_constraints!`](@ref).

Any sequence will also have a `level=3` cost function, which minimises the total sequence duration.
"""
function add_cost_function!(func::AbstractJuMPScalar, level=2)
    push!(GLOBAL_MODEL[][2], (Float64(level), func))
end

add_cost_function!(func::Number, level=2) = nothing


"""
    set_simple_constraints!(block, kwargs)

Add any constraints or objective functions to the variables of a [`AbstractBlock`](@ref).

Each keyword argument has to match one of the functions in [`variables`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(block::AbstractBlock, kwargs)
    real_kwargs = Dict(key => value for (key, value) in kwargs if !isnothing(value))

    for (key, value) in real_kwargs
        var = getproperty(variables, key)
        if var isa AlternateVariable
            if var.other_var in keys(real_kwargs)
                error("Set constraints on both $key and $(var.other_var), however they are equivalent.")
            end
            invert_value(value::VariableType) = var.to_other(value)
            invert_value(value::Symbol) = invert_value(Val(value))
            invert_value(::Val{:min}) = var.inverse ? Val(:max) : Val(:min)
            invert_value(::Val{:max}) = var.inverse ? Val(:min) : Val(:max)
            invert_value(value::AbstractVector) = invert_value.(value)
            apply_simple_constraint!(getproperty(variables, var.other_var)(block), invert_value(value))
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


    apply_simple_constraint!(variable, :>=/:<=, value)

Set an inequality constraint to the `variable`.

`value` can be a number of an equation.
"""
apply_simple_constraint!(variable::AbstractVector, value::Symbol) = apply_simple_constraint!(sum(variable), Val(value))
apply_simple_constraint!(variable, value::Nothing) = nothing
apply_simple_constraint!(variable::NamedTuple, value::Nothing) = nothing
apply_simple_constraint!(variable::VariableType, value::Symbol) = apply_simple_constraint!(variable, Val(value))
apply_simple_constraint!(variable::VariableType, ::Val{:min}) = add_cost_function!(variable, 1)
apply_simple_constraint!(variable::VariableType, ::Val{:max}) = add_cost_function!(-variable, 1)
apply_simple_constraint!(variable::VariableType, value::VariableType) = @constraint GLOBAL_MODEL[][1] variable == value
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

apply_simple_constraint!(variable::VariableType, sign::Symbol, value) = apply_simple_constraint!(variable, Val(sign), value)
apply_simple_constraint!(variable::VariableType, ::Val{:(==)}, value) = apply_simple_constraint!(variable, value)
apply_simple_constraint!(variable::VariableType, ::Val{:(<=)}, value::VariableType) = apply_simple_constraint!(value, Val(:>=), variable)
function apply_simple_constraint!(variable::VariableType, ::Val{:(>=)}, value::VariableType)
    @constraint GLOBAL_MODEL[][1] variable >= value
    if value isa Number && iszero(value)
        @constraint GLOBAL_MODEL[][1] variable * 1e6 >= value
        @constraint GLOBAL_MODEL[][1] variable * 1e12 >= value
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

Constraints [`variables.gradient_strength`](@ref) and [`variables.slew_rate`](@ref) to be less than the [`global_scanner`](@ref) maximum.
"""
function scanner_constraints!(bb::AbstractBlock)
    for (var, max_value) in [
        (variables.slew_rate, global_scanner().slew_rate),
        (variables.gradient_strength, global_scanner().gradient),
    ]
        value = nothing
        try
            value = var(bb)
        catch
            continue
        end
        for v in value
            if v isa Number || ((v isa Union{QuadExpr, AffExpr}) && length(v.terms) == 0)
                continue
            end
            apply_simple_constraint!(v, :<=, max_value)
            apply_simple_constraint!(v, :>=, -max_value)
        end
    end
end

end