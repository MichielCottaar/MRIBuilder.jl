module BuildSequences
using JuMP
import Ipopt
import Juniper
import ..Scanners: Scanner, gradient_strength, Default_Scanner

const GLOBAL_MODEL = Ref((Model(), Tuple{Float64, AbstractJuMPScalar}[]))
const IGNORE_MODEL = GLOBAL_MODEL[]

const GLOBAL_SCANNER = Ref(Scanner())


"""
    iterate_cost()

Return a sequence of all the cost functions that should be optimised in order.
"""
function iterate_cost()
    unique_weights = sort(unique([w for (w, _) in GLOBAL_MODEL[][2]]))
    if iszero(length(unique_weights))
        return [0.]
    else
        return [sum([f for (w, f) in GLOBAL_MODEL[][2] if w <= weight]; init=0.) for weight in unique_weights]
    end
end

function total_cost_func()
    return sum([10^(-w) * f for (w, f) in GLOBAL_MODEL[][2]]; init=0.)
end


"""
Wrapper to build a sequence.

Use as 
```julia
build_sequence(scanner;) do
    ...
end
```
Within the code block you can create one or more sequences, e.g.
```
seq = Sequence(
    SincPulse(flip_angle=90, phase=0, duration=2., bandwidth=:max)
    nothing.,
    SingleReadout
)
```

You can also add any arbitrary constraints or objectives using one of:
- `set_simple_constraints!`
- `apply_simple_constraint!`
- `add_cost_function!`

As soon as the code block ends the sequence is optimised (if `optimise=true`) and returned.

## Parameters
- `scanner`: Set to a [`Scanner`](@ref) to limit the gradient strength and slew rate. When this call to `build_sequence` is embedded in another, this parameter can be set to `nothing` to indicate that the same scanner should be used. 
- `optimise`: Whether to optimise and fix the sequence as soon as it is returned. This defaults to `true` if a scanner is provided and `false` if no scanner is provided.
- `n_attempts`: How many times to restart the optimiser (default: 20)? Decrease if you want to quickly check feasibility of a simple sequence. Increase if convergence fails for a complex sequence.
- `kwargs...`: Other keywords are passed on as attributes to the `optimiser_constructor` (e.g., set `print_level=3` to make the Ipopt optimiser quieter).
"""
function build_sequence(f::Function, scanner::Union{Nothing, Scanner}, model::Tuple, optimise::Bool, n_attempts::Int)
    prev_model = GLOBAL_MODEL[]
    GLOBAL_MODEL[] = model
    prev_scanner = GLOBAL_SCANNER[]
    if !isnothing(scanner)
        GLOBAL_SCANNER[] = scanner
    elseif !isfinite(gradient_strength(GLOBAL_SCANNER[]))
        error("Scanner should be explicitly set when creating a new top-level sequence.")
    end
    try
        sequence = f()
        if optimise
            jump_model = GLOBAL_MODEL[][1]
            set_optimizer(jump_model, get_optimiser(jump_model))
            if !iszero(num_variables(jump_model))
                optimise_with_cost_func!(jump_model, total_cost_func(), n_attempts)
                prev_cost_func = nothing
                for cost_func in iterate_cost()
                    if !isnothing(prev_cost_func)
                        @constraint jump_model prev_cost_func <= value(prev_cost_func)
                    end
                    optimise_with_cost_func!(jump_model, cost_func, n_attempts)
                    prev_cost_func = cost_func
                end
            end
            return fixed(sequence)
        else
            return sequence
        end
    finally
        GLOBAL_MODEL[] = prev_model
        GLOBAL_SCANNER[] = prev_scanner
    end
end

function number_equality_constraints(model::Model)
    sum([num_constraints(model, expr, comp) for (expr, comp) in JuMP.list_of_constraint_types(model) if comp <: MOI.EqualTo])
end

"""
    model_has_integers(model)

Returns true if the model contains integer variables.
"""
model_has_integers(model::Model) = (VariableRef, MOI.Integer) in list_of_constraint_types(model)

"""
    get_optimiser(model)

Returns a JuMP solver (https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) appropriate for this model.
"""
function get_optimiser(model::Model)
    base_ipopt = optimizer_with_attributes(
        Ipopt.Optimizer, 
        "print_level"=>0, 
        "mu_strategy"=>"adaptive", 
        "max_iter"=>1000,
    )
    if model_has_integers(model)
         optimizer_with_attributes(
            Juniper.Optimizer, 
            "nl_solver" => base_ipopt,
            "log_levels" => [],
        )
    else
        return base_ipopt
    end
end


function optimise_with_cost_func!(jump_model::Model, cost_func, n_attempts)
    @objective jump_model Min cost_func
    min_objective = Inf
    min_values = nothing
    nsuccess = 0
    sample_error = Dict{TerminationStatusCode, Any}()
    errors = TerminationStatusCode[]
    for _ in num_variables(jump_model):number_equality_constraints(jump_model)
        @variable(jump_model)
    end
    for attempt in 1:n_attempts
        if attempt != 1
            new_values = rand(length(all_variables(jump_model)))
            for (var, v) in zip(all_variables(jump_model), new_values)
                set_start_value(var, v)
            end
        end
        for sub_attempt in 1:3
            if sub_attempt != 1
                old_values = value.(all_variables(jump_model))
                size_kick = 0.5 / sub_attempt
                new_values = old_values .* (2 .* size_kick .* rand(length(old_values)) .+ 1. .- size_kick)
                for (var, v) in zip(all_variables(jump_model), new_values)
                    set_start_value(var, v)
                end
            end
            optimize!(jump_model)
            #println(solution_summary(jump_model))
            if termination_status(jump_model) in (LOCALLY_SOLVED, OPTIMAL)
                nsuccess += 1
                if objective_value(jump_model) < min_objective
                    min_objective = objective_value(jump_model)
                    min_values = copy(get_inner_state(jump_model))
                end
                break
            elseif sub_attempt == 3
                push!(errors, termination_status(jump_model))
                sample_error[termination_status(jump_model)] = solution_summary(jump_model)
            end
        end
        if nsuccess > 2
            break
        end
    end
    if iszero(nsuccess)
        err_string = join(["$(String(Symbol(err))) (x$(sum([e == err for e in errors])))" for err in unique(errors)], ", ", ", and ")
        for msg in values(sample_error)
            println(msg)
            println("")
        end
        error("Optimisation failed to converge. The following errors were raised: $err_string. Example errors for each type are printed above.")
    end
    get_inner_state(jump_model) .= min_values
    @assert value(cost_func) ≈ min_objective eps(min_objective)
end

get_inner_state(backend::Juniper.JuniperProblem) = backend.solution
get_inner_state(backend::Ipopt.IpoptProblem) = backend.x
get_inner_state(model::Model) = get_inner_state(backend(model).optimizer.model.inner)

function build_sequence(f::Function, scanner::Union{Nothing, Scanner}=Default_Scanner; optimise=true, n_attempts=20)
    if optimise || GLOBAL_MODEL[] == IGNORE_MODEL
        model = (
            Model(),
            Tuple{Float64, AbstractJuMPScalar}[]
        )
    else
        model = GLOBAL_MODEL[]
    end
    build_sequence(f, scanner, model, optimise, n_attempts)
end


"""
    global_scanner()

Return the currently set [`Scanner`](@ref).

The scanner can be set using [`build_sequence`](@ref)
"""
function global_scanner()
    if !isfinite(GLOBAL_SCANNER[].gradient)
        error("No valid scanner has been set. Please provide one when calling `build_sequence`.")
    end
    return GLOBAL_SCANNER[]
end


"""
    fixed(building_block)

Return an equiavalent `BuildingBlock` with all free variables replaced by numbers.

This will only work after calling `optimize!`(@ref)([`global_model`](@ref)()).
It is used internally by [`build_sequence`](@ref).
"""
fixed(some_value) = some_value
fixed(jump_variable::AbstractJuMPScalar) = value(jump_variable)
fixed(jump_variable::GenericVariableRef) = is_integer(jump_variable) ? Int(value(jump_variable)) : value(jump_variable)
fixed(jump_variable::AbstractArray) = fixed.(jump_variable)
fixed(dict_variable::AbstractDict) = typeof(dict_variable)(k => fixed(v) for (k, v) in pairs(dict_variable))
fixed(tuple_variable::Tuple) = fixed.(tuple_variable)
fixed(pair:: Pair) = fixed(pair[1]) => fixed(pair[2])



end