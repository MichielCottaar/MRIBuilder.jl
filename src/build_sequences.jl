module BuildSequences
import JuMP: Model, optimizer_with_attributes, optimize!, AbstractJuMPScalar, value, solution_summary, termination_status, LOCALLY_SOLVED, OPTIMAL, num_variables, all_variables, set_start_value, ALMOST_LOCALLY_SOLVED, objective_value, INVALID_MODEL, @variable
import Ipopt
import Juniper
import ..Scanners: Scanner, gradient_strength, Default_Scanner

const GLOBAL_MODEL = Ref(Model())
const IGNORE_MODEL = GLOBAL_MODEL[]

const GLOBAL_SCANNER = Ref(Scanner())

"""
Wrapper to build a sequence.

Use as 
```julia
build_sequence(scanner[, optimiser_constructor];) do
    ...
end
```
Within the code block you can create one or more sequences, e.g.
```
seq = Sequence(
    SincPulse(flip_angle=90, phase=0, duration=2., bandwidth=:max)
    nothing.,
    InstantReadout
)
```

You can also add any arbitrary constraints or objectives using the same syntax as for [`JuMP`](https://jump.dev/JuMP.jl):
```
@constraint global_model() duration(seq) == 30.
```

As soon as the code block ends the sequence is optimised (if `optimise=true`) and returned.

## Parameters
- `scanner`: Set to a [`Scanner`](@ref) to limit the gradient strength and slew rate. When this call to `build_sequence` is embedded in another, this parameter can be set to `nothing` to indicate that the same scanner should be used. 
- `optimiser_constructor`: A `JuMP` solver optimiser as described in the [JuMP documentation](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/#What-is-a-solver?). Defaults to using [Ipopt](https://github.com/jump-dev/Ipopt.jl).
- `optimise`: Whether to optimise and fix the sequence as soon as it is returned. This defaults to `true` if a scanner is provided and `false` if no scanner is provided.
- `n_attempts`: How many times to restart the optimser (default: 100).
- `kwargs...`: Other keywords are passed on as attributes to the `optimiser_constructor` (e.g., set `print_level=3` to make the Ipopt optimiser quieter).
"""
function build_sequence(f::Function, scanner::Union{Nothing, Scanner}, model::Model, optimise::Bool, n_attempts::Int)
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
            if !iszero(num_variables(model))
                min_objective = Inf
                for attempt in 1:n_attempts
                    if attempt != 1
                        old_values = value.(all_variables(model))
                        size_kick = 0.5 / attempt
                        new_values = old_values .* (2 .* size_kick .* rand(length(old_values)) .+ 1. .- size_kick)
                        for (var, v) in zip(all_variables(model), new_values)
                            set_start_value(var, v)
                        end
                    end
                    optimize!(model)
                    while termination_status(model) == INVALID_MODEL
                        @variable(model, x[i=1:10])
                        optimize!(model)
                    end
                    if termination_status(model) in (LOCALLY_SOLVED, OPTIMAL)
                        if objective_value(model) < min_objective
                            min_objective = objective_value(model)
                        elseif isapprox(min_objective, objective_value(model), rtol=1e-6)
                            println("Optimisation succeeded after $(attempt-1) restarts.")
                            break
                        end
                    end
                end
                if !(termination_status(model) in (LOCALLY_SOLVED, OPTIMAL))
                    @warn "Optimisation did not report successful convergence. Please check the output sequence."
                    println(solution_summary(model))
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

function build_sequence(f::Function, scanner::Union{Nothing, Scanner}, optimiser_constructor; optimise=true, n_attempts=1000, kwargs...)
    if optimise || GLOBAL_MODEL[] == IGNORE_MODEL
        model = Model(optimizer_with_attributes(optimiser_constructor, [string(k) => v for (k, v) in kwargs]...))
    else
        model = global_model()
    end
    build_sequence(f, scanner, model, optimise, n_attempts)
end

function build_sequence(f::Function, scanner::Union{Nothing, Scanner}=Default_Scanner; print_level=0, kwargs...)
    build_sequence(f, scanner, Ipopt.Optimizer; print_level=print_level, kwargs...)
end

build_sequence(f::Function, optimiser_constructor; kwargs...) = build_sequence(f, Default_Scanner, optimiser_constructor; kwargs...)


"""
    global_model()

Return the currently set JuMP model.

The model can be set using [`build_sequence`](@ref)
"""
function global_model()
    if GLOBAL_MODEL[] == IGNORE_MODEL
        error("No global model has been set. Please embed any sequence creation within an `build_sequence` block.")
    end
    return GLOBAL_MODEL[]
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
fixed(jump_variable::AbstractArray) = fixed.(jump_variable)
fixed(dict_variable::AbstractDict) = typeof(dict_variable)(k => fixed(v) for (k, v) in pairs(dict_variable))
fixed(tuple_variable::Tuple) = fixed.(tuple_variable)
fixed(pair:: Pair) = fixed(pair[1]) => fixed(pair[2])



end