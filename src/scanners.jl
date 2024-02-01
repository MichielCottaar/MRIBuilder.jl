"""
Define general [`Scanner`](@ref) type and methods as well as some concrete scanners.
"""
module Scanners
import JuMP: Model, @constraint, owner_model
import ..Variables: gradient_strength, slew_rate
import ..BuildingBlocks: BuildingBlock, get_children_blocks

const gyromagnetic_ratio = 42576.38476  # (kHz/T)

"""
    Scanner(;B0=3., gradient=Inf, slew_rate=Inf, units=:kHz)

Properties of an MRI scanner relevant for the MR signal simulations.
- [`B0`](@ref): magnetic field strength (in Tesla)
- [`gradient_strength`](@ref): maximum gradient strength long each axis.
- [`slew_rate`](@ref): maximum rate of change in the gradient strength

By default `gradient` and `slew_rate` are expected to be provided in units of, respectively, kHz/um and kHz/um/ms.
However, if the keyword `units=:Tesla` is set, the `gradient` and `slew_rate` should be provided in units of, respectively, mT/m and T/m/s.
"""
struct Scanner
    B0::Float64
    gradient::Float64
    slew_rate::Float64
    function Scanner(;B0=3., gradient=Inf, slew_rate=Inf, units=:kHz)
        if isnothing(gradient)
            gradient = Inf
        end
        if isnothing(slew_rate)
            slew_rate = Inf
        end
        if isinf(gradient) && !isinf(slew_rate)
            error("Can't have infinite gradient strength with finite slew rate.")
        end
        if units == :Tesla
            gradient *= gyromagnetic_ratio * 1e-9
            slew_rate *= gyromagnetic_ratio * 1e-9
        end
        new(Float64(B0), Float64(gradient), Float64(slew_rate))
    end
end

"""
    B0(scanner)
    B0(sequence)

Returns the magnetic field strength of the scanner in Tesla.
"""
B0(scanner::Scanner) = scanner.B0

"""
    gradient_strength(scanner[, units])

Returns the maximum magnetic field gradient of the scanner in kHz/um.
By setting `units` to :Tesla, the gradient strength can be returned in mT/m instead.
"""
gradient_strength(scanner::Scanner, units=:kHz) = units == :kHz ? scanner.gradient : scanner.gradient / (gyromagnetic_ratio * 1e-9)

"""
    slew_rate(scanner[, units])

Returns the maximum magnetic field slew rate of the scanner in kHz/um/ms.
By setting `units` to :Tesla, the slew rate can be returned in T/m/s instead.
"""
slew_rate(scanner::Scanner, units=:kHz) = units == :kHz ? scanner.slew_rate : scanner.slew_rate / (gyromagnetic_ratio * 1e-9)


"""
Siemens MAGNETOM 3T Prisma MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/3t-mri-scanner/magnetom-prisma).
"""
Siemens_Prisma = Scanner(B0=3., gradient=80, slew_rate=200, units=:Tesla)

"""
Siemens MAGNETOM 7T Terra MRI scanner (https://www.siemens-healthineers.com/en-uk/magnetic-resonance-imaging/7t-mri-scanner/magnetom-terra)
"""
Siemens_Terra = Scanner(B0=7., gradient=80, slew_rate=200, units=:Tesla)

"""
Siemens 3T Connectom MRI scanner ([fan22_MappingHumanConnectome](@cite)).
"""
Siemens_Connectom = Scanner(B0=3., gradient=300, slew_rate=200, units=:Tesla)

predefined_scanners = Dict(
    :Siemens_Prisma => Siemens_Prisma,
    :Siemens_Terra => Siemens_Terra,
    :Siemens_Connectom => Siemens_Connectom,
)

"""
    scanner_constraints!([model, ]building_block, scanner)

Adds any constraints from a specific scanner to a [`BuildingBlock`]{@ref}.
"""
function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    scanner_constraints!(owner_model(building_block), building_block, scanner)
end

function scanner_constraints!(model::Model, building_block::BuildingBlock, scanner::Scanner)
    for func in [gradient_strength, slew_rate]
        if isfinite(func(scanner))
            scanner_constraints!(model, building_block, scanner, func)
        end
    end
end

function scanner_constraints!(model::Model, building_block::BuildingBlock, scanner::Scanner, func::Function)
    if func in variables(building_block)
        # apply constraint at this level
        res_bb = func(building_block)
        if res_bb isa AbstractVector
            if isnothing(building_block.rotate)
                # no rotation; apply constraint to each dimension independently
                for expr in res_bb
                    @constraint model expr <= func(scanner)
                    @constraint model expr >= -func(scanner)
                end
            else
                # with rotation: apply constraint to total squared
                total_squared = sum(map(n->n^2, res_bb))
                @constraint model total_squared <= func(scanner)^2
            end
        else
            @constraint model res_bb <= func(scanner)
            @constraint model res_bb >= -func(scanner)
        end
    elseif building_block isa ContainerBlock
        # apply constraints at lower level
        for (_, child_block) in get_children_blocks(building_block)
            scanner_constraints!(model, child_block, scanner, func)
        end
    end
end

end