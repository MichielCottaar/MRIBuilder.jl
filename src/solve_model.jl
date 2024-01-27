module SolveModel
import JuMP: Model

const GLOBAL_MODEL = Ref(Model())
const IGNORE_MODEL = GLOBAL_MODEL[]

"""
Sets a global JuMP `Model`.

Use as 
```julia
model = set_model([model]) do
    ...
end
```
Any sequences created within the code block will be assigned the same JuMP `Model`.
If no model is provided a new one is created (using a Juniper optimizer based on the Ipopt non-linear optimizer).

The function will return the fully formed JuMP `Model`.
"""
function set_model(f::Function, model::Model)
    prev_model = global_model[]
    global_model[] = model
    try
        f()
    finally
        global_model[] = prev_model
    end
    return model
end

function set_model(f::Function)
    ipopt_opt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
    juniper_opt = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt_opt)
    model = Model(juniper_opt)
    set_model(f, model)
end


function get_global_model()
    if GLOBAL_MODEL[] == IGNORE_MODEL
        error("No global model has been set. Please explicitly set one in the constructor or set a global model using `set_model`.")
    end
    return GLOBAL_MODEL[]
end


"""
    @global_model_constructor BuildingBlockType

Add a constructor to the [`BuildingBlock`](@ref) subtype that fetches the global JuMP model (set by [`set_model`](@ref)) and assigns it to the first argument.
```
BuildingBlockType(args...; kwargs...) = BuildingBlockType(global_model::JuMP.Model, args...; kwargs...)
```
"""
macro global_model_constructor(bb)
    quote
        $(esc(bb))(args...; kwargs...) = $(esc(bb))(get_global_model(), args...; kwargs...)
    end
end


end