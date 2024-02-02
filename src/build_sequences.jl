module BuildSequences
import JuMP: Model, optimizer_with_attributes, optimize!
import Ipopt
import Juniper
import ..Scanners: Scanner, scanner_constraints!
import ..BuildingBlocks: BuildingBlock, apply_simple_constraint!, match_blocks!

const GLOBAL_MODEL = Ref(Model())
const IGNORE_MODEL = GLOBAL_MODEL[]

const GLOBAL_SCANNER = Ref(Scanner())

"""
Wrapper to build a sequence.

Use as 
```julia
build_sequence([model]) do model
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
@constraint model duration(seq) == 30.
```

As soon as the code block is the optimal sequence matching all your constraints and objectives will be returned.
"""
function build_sequence(f::Function, scanner::Scanner, model::Model)
    prev_model = GLOBAL_MODEL[]
    GLOBAL_MODEL[] = model
    prev_scanner = GLOBAL_SCANNER[]
    if !isnothing(scanner)
        GLOBAL_SCANNER[] = scanner
    end
    try
        sequence = f(model)
        optimize!(model)
        return sequence
    finally
        GLOBAL_MODEL[] = prev_model
        GLOBAL_SCANNER[] = prev_scanner
    end
end

function build_sequence(f::Function, scanner::Scanner)
    ipopt_opt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3)
    juniper_opt = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt_opt)
    model = Model(ipopt_opt)
    build_sequence(f, scanner, model)
end


function global_model()
    if GLOBAL_MODEL[] == IGNORE_MODEL
        error("No global model has been set. Please explicitly set one in the constructor or set a global model using `set_model`.")
    end
    return GLOBAL_MODEL[]
end

function global_scanner()
    if !isfinite(GLOBAL_SCANNER[].gradient)
        error("No valid scanner has been set. Please provide one when calling `build_sequence`.")
    end
    return GLOBAL_SCANNER[]
end


scanner_constraints!(bb::BuildingBlock)  = scanner_constraints!(bb, global_scanner())

apply_simple_constraint!(variable, value) = apply_simple_constraint!(global_model(), variable, value)
match_blocks!(block1::BuildingBlock, block2::BuildingBlock, property_list) = match_blocks!(global_model(), block1, block2, property_list)
scanner_constraints!(building_block::BuildingBlock, scanner::Scanner, func::Function) = scanner_constraints!(building_block, scanner, func)
scanner_constraints!(building_block::BuildingBlock) = scanner_constraints!(building_block, global_scanner())

end