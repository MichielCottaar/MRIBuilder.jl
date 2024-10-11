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

function fixed(ab::T) where {T<:AbstractBlock}
    params = []
    for prop_name in fieldnames(T)
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


"""
    base_variables([T])

Return dictionary with all [`Variable`](@ref) objects defined for a specific sequence component/block `T`.

This only returns those [`Variable`](@ref) directly defined for this component/block, not for any sub-components (through [`get_pulse`](@ref), [`get_gradient`][(@ref), etc.)

If `T` is not provided, all [`Variable`](@ref) objects are returned.
"""
function base_variables()
    return Dict{Symbol, Function}(
        s => getproperty(variables, s) for s in names(variables) if s != :variables
    )
end

function base_variables(T::Type{<:AbstractBlock})
    return Dict(
        s => v for (s, v) in base_variables() if which(v, (T, )) !== default_generic_method[v]
    )
end


"""
    variable_defined_for(variable, args...; kwargs...)
"""
function variable_defined_for(name::Symbol, args...; kwargs...)
    variable_defined_for(getproperty(variables, name), args...; kwargs...)
end

function variable_defined_for(func::Function, args...; kwargs...)
    return which(func, typeof.(args)) !== default_generic_method[func]
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
    if v in Base.unsorted_names(variables; all=true)
        var = getproperty(variables, v)
        return var(block)
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
baremodule variables end


"""
Contains for each variable the default, generic method.

This is the method that checks for [`alternative_functions`](@ref) or uses one of the getters.
"""
const default_generic_method = IdDict{Function, Method}()

"""
Mapping of variable names to alternative ways to compute that variables.
"""
const alternative_variables = Dict{Symbol, Vector{Tuple{Symbol, Function}}}()

"""
Raised if there is no way to reach a valid function by the default, generic method ([`default_generic_method`](@ref)).
"""
struct InvalidRoute <: Exception end

"""
Assigns getters to specific variables.
"""
const getters = Dict{Symbol, Function}()

"""
    _get_variable(name, tried_names, args...; kwargs...)

Tries to find a route to get values for the variable `name` with the given `args` and `kwargs`.

The route through `tried_names` has already been attempted.

This function returns the route to get to the value and the value itself.
"""
function _get_variable(name::Symbol, tried_names::Set{Symbol}, args...; kwargs...)
    if hasproperty(variables, name)
        func = getproperty(variables, name)
        if variable_defined_for(func, args...; kwargs...)
            return Any[], func(args...; kwargs...)
        end
    end

    # alternative functions
    if name in keys(alternative_variables)
        for (new_name, converter) in alternative_variables[name]
            if !(new_name in tried_names)
                try
                    route, value = _get_variable(new_name, union(tried_names, [name]), args...; kwargs...)
                    pushfirst!(route, (:alternative, new_name))
                    return route, converter(value)
                catch e
                    if e isa InvalidRoute
                        continue
                    else
                        rethrow()
                    end
                end
            end
        end
    end

    if name in keys(getters)
        # getters
        first, others = args[1], args[2:end]
        getter = getters[name]
        route, value = _get_mult_variable(name, getter(first), others...; kwargs...)
        pushfirst!(route, (:getter, getter))
        return route, value
    end

    throw(InvalidRoute())
end

"""
Helper to call the variable for a result of a getter. Used in [`_get_variable`](@ref).
"""
_get_mult_variable(name, t, args...; kwargs...) = _get_variable(name, Set{Symbol}(), t, args...; kwargs...)
_get_mult_variable(name, t::NamedTuple, args...; kwargs...) = (Any[], NamedTuple(
    key => _get_mult_variable(name, var, args...; kwargs...) for (key, (_, var)) in pairs(t)
))
_get_mult_variable(name, t::AbstractVector, args...; kwargs...) = (Any[], [
    _get_mult_variable(name, var, args...; kwargs...) for (_, var) in t
])
_get_mult_variable(name, t::Tuple, args...; kwargs...) = (Any[], Tuple(
    _get_mult_variable(name, var, args...; kwargs...) for (_, var) in t
))

"""
    add_new_variable!(name)

Adds a new variable in [`variables`](@ref).

This is a helper function, which is called by `@defvar` for any variable that does not exist yet.
"""
function add_new_variable!(name::Symbol)
    @eval variables function $(name) end
    @eval variables $(Expr(:public, name))
    @eval function variables.$(name)(ab::AbstractBlock, args...; kwargs...) 
        try
            return _get_variable($(QuoteNode(name)), Set{Symbol}(), ab, args...; kwargs...)[2]
        catch e
            if e isa InvalidRoute
                string_name = String($(QuoteNode(name)))
                error("Variable `" * string_name * "` is not defined for block of type `" * string(typeof(ab)) * ".")
            end
            rethrow()
        end
    end

    func = getproperty(variables, name)
    default_generic_method[func] = methods(func)[1]
end


"""
    @defvar(function(s))

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
"""
macro defvar(getter, func_def)
    func_names, expr = _defvar(func_def)
    additional = :(for name in $(func_names)
        set_getter!(name, $(QuoteNode(getter)))
    end)
    return Expr(
        :block,
        expr.args...,
        additional
    )
end

macro defvar(func_def)
    _, expr = _defvar(func_def)
    return expr
end

function _defvar(func_def)
    func_names = []

    function adjust_function(ex)
        if ex isa Expr 
            if ex.head == :block
                return Expr(:block, adjust_function.(ex.args)...)
            end
            if ex.head == :function && length(ex.args) == 1
                push!(func_names, ex.args[1])
                return :nothing
            end
            if ex.head == :macrocall && ex.args[1] == GlobalRef(Core, Symbol("@doc"))
                new_expr = macroexpand(variables, ex)
                func_def, add_docs = new_expr.args
                fixed_func_def = adjust_function(func_def)
                add_docs.args[end] = esc(add_docs.args[end])
                return Expr(:block, fixed_func_def, add_docs)
            end
        end
        try
            fn_def = MacroTools.splitdef(ex)
            push!(func_names, fn_def[:name])
            new_def = Dict{Symbol, Any}()
            new_def[:name] = Expr(:., :variables, QuoteNode(fn_def[:name]))
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
                add_new_variable!($(QuoteNode(func_name)))
            end
        end
        )
    end
    args = vcat([e.args for e in expressions]...)
    return func_names, Expr(
        :block,
        args...,
        new_func_def
    )
end

@defvar begin
    """
        duration(block)

    Duration of the sequence or building block in ms.
    """ 
    function duration end
end


"""
    set_getter!(variable_name, getter)

Set the getter function for `variable_name`.

If the value for `variable` is not defined for a sequence, the value for the result of the `getter` function is returned instead.

Possible values for the `getter` function are:
- `:pulse`: [`get_pulse`](@ref)
- `:gradient`: [`get_gradient`](@ref)
- `:readout`: [`get_readout`](@ref)
- `:pathway`: [`get_pathway`](@ref)
"""
set_getter!(name::Symbol, getter::Symbol) = set_getter!(name, getter_functions[getter])

function set_getter!(name::Symbol, getter::Function)
    if !(name in keys(getters))
        getters[name] = getter
    else
        if getters[name] !== getter
            error("")
        end
    end
end

"""
    add_alternative_variable!(name, other_func, conversion)

Defines an alternative way to compute the variable with given `name`.

If the variable `name` is not defined and `other_name` is,
then the value of `name` is computed by applying `conversion` to the value of `other_name`.
"""
function add_alternative_variable!(name::Symbol, other_name::Symbol, conversion)
    if !(name in keys(alternative_variables))
        alternative_variables[name] = Tuple{Symbol, Function}[]
    end
    push!(alternative_variables[name], (other_name, conversion))
end

@defvar begin
    """
        qval(gradient)

    The norm of the [`variables.qvec`](@ref).
    """
    function qval end

    """
        qval_square(gradient)

    The square of the area under the curve ([`variables.qval`](@ref)).

    Constraining this can be more efficient than constraining [`qval`](@ref) as it avoids taking a square root.
    """
    function qval_square end

    """
        spoiler(gradient)

    Spatial scale in mm over which the spoiler gradient will dephase by 2π.

    Automatically computed based on [`variables.qvec`](@ref).
    """
    function spoiler end

    """
        qvec(gradient)

    The total integral of the area under the gradient curve as a length-3 vector.

    The norm of this vector is available as [`qval`](@ref).
    """
    function qvec end
end

add_alternative_variable!(:spoiler, :qval, q->1e-3 * 2π/q)
add_alternative_variable!(:qval, :spoiler, l->1e-3 * 2π/l)
add_alternative_variable!(:qval, :qval_square, sqrt)
add_alternative_variable!(:qval_square, :qvec, qv -> sum(q -> q * q, qv))



for vec_variable in [:gradient_strength, :slew_rate]
    vec_square = Symbol(string(vec_variable) * "_square")
    vec_norm = Symbol(string(vec_variable) * "_norm")
    add_alternative_variable!(vec_square, vec_variable, v -> v[1] * v[1] + v[2] * v[2] + v[3] * v[3])
    add_alternative_variable!(vec_norm, vec_square, sqrt)
end

@defvar begin
    """
        gradient_strength_norm(gradient)

    The norm of the [`variables.gradient_strength`](@ref).
    """
    function gradient_strength_norm end

    """
        slew_rate_norm(gradient)

    The norm of the [`variables.slew_rate`](@ref).
    """
    function slew_rate_norm end
end

for name in [:fov, :voxel_size]
    inv_name = Symbol("inverse_" * string(name))
    mult_inv(x) = inv.(x)
    add_alternative_variable!(name, inv_name, mult_inv)
    add_alternative_variable!(inv_name, name, mult_inv)
    add_new_variable!(name)
    add_new_variable!(inv_name)
end

for name in [:slice_thickness, :bandwidth]
    inv_name = Symbol("inverse_" * string(name))
    add_alternative_variable!(name, inv_name, inv)
    add_alternative_variable!(inv_name, name, inv)
    add_new_variable!(name)
    add_new_variable!(inv_name)
end


@defvar begin
    """
        echo_time(sequence)

    Computes the echo time(s) of a sequence in ms.
    """
    function echo_time end

    """
        repetition_time(sequence)

    Computes the repetition_times of a sequence in ms.
    """
    function repetition_time end

    """
        diffusion_time(sequence)

    Computes the diffusion time of a sequence in ms.
    """
    function diffusion_time end
end

@eval variables begin
    const TE = echo_time
    const TR = repetition_time
    const Δ = diffusion_time
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
Mapping of symbols to actual getter functions.

Used in [`set_getter!`](@ref).
"""
getter_functions = Dict(
    :pulse => get_pulse,
    :gradient => get_gradient,
    :pathway => get_pathway,
    :readout => get_readout,
)

"""
    gradient_orientation(building_block)

Returns the gradient orientation.
"""
function gradient_orientation end


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
        apply_simple_constraint!(var(block), value)
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