"""
Defines a set of different options for MRI gradients.
"""
module Trapezoids

import JuMP: @constraint
import StaticArrays: SVector
import LinearAlgebra: norm
import ...Variables: qval, rise_time, flat_time, slew_rate, gradient_strength, variables, duration, δ, get_free_variable, VariableType, inverse_bandwidth, effective_time, qval_square, duration, set_simple_constraints!, scanner_constraints!, inverse_slice_thickness
import ...Variables: Variables, all_variables_symbols, dwell_time, inverse_fov, inverse_voxel_size, fov, voxel_size, get_gradient, get_pulse, get_readout, gradient_orientation
import ...BuildSequences: global_model
import ...Components: ChangingGradient, ConstantGradient, RFPulseComponent, ADC
import ..BaseBuildingBlocks: BaseBuildingBlock


"""
Parent type for any [`BuildingBlock`](@ref) that has a trapezoidal gradient waveform.

Sub-types:
- [`Trapezoid`](@ref)
- [`SliceSelect`](@ref)
- [`LineReadout`](@ref)

The `N` indicates whether the gradient has a fixed orientation (N=1) or is free (N=3).
"""
abstract type BaseTrapezoid{N} <: BaseBuildingBlock end

"""
    Trapezoid(; orientation=nothing, group=nothing, variables...)

Defines a trapezoidal pulsed gradient

## Parameters
- `orientation` sets the gradient orientation (completely free by default). Can be set to a vector for a fixed orientation.
- `group`: assign the trapezoidal gradient to a specific group. This group will be used to scale or rotate the gradients after optimisation.

## Variables
Variables can be set during construction or afterwards as an attribute.
If not set, they will be determined during the sequence optimisation.
### Timing variables
- [`rise_time`](@ref): Time of the gradient to reach from 0 to maximum in ms. If explicitly set to 0, the scanner slew rate will be ignored.
- [`flat_time`](@ref): Time that the gradient stays at maximum strength in ms.
- [`δ`](@ref): effective pulse duration (`rise_time` + `flat_time`) in ms.
- [`duration`](@ref): total pulse duration (2 * `rise_time` + `flat_time`) in ms.
### Gradient variables
- [`gradient_strength`](@ref): Maximum gradient strength achieved during the pulse in kHz/um
- [`qval`](@ref): Spatial scale on which spins will be dephased due to this pulsed gradient in rad/um (given by `δ` * `gradient_strength`).

The `bvalue` can be constrained for multiple gradient pulses by creating a `Pathway`.
"""
abstract type Trapezoid{N} <: BaseTrapezoid{N} end

struct Trapezoid1D <: Trapezoid{1}
    rise_time :: VariableType
    flat_time :: VariableType
    slew_rate :: VariableType
    orientation :: SVector{3, Float64}
    group :: Union{Nothing, Symbol}
end

struct Trapezoid3D <: Trapezoid{3}
    rise_time :: VariableType
    flat_time :: VariableType
    slew_rate :: SVector{3, VariableType}
    group :: Union{Nothing, Symbol}
end

function (::Type{Trapezoid})(; orientation=nothing, rise_time=nothing, flat_time=nothing, group=nothing, slew_rate=nothing, kwargs...)
    if isnothing(orientation)
        if isnothing(slew_rate)
            slew_rate = (nothing, nothing, nothing)
        end
        res = Trapezoid3D(
            get_free_variable(rise_time),
            get_free_variable(flat_time),
            get_free_variable.(slew_rate),
            group
        )
    else
        res = Trapezoid1D(
            get_free_variable(rise_time),
            get_free_variable(flat_time),
            get_free_variable(slew_rate),
            orientation,
            group
        )
        @constraint global_model() res.slew_rate >= 0
    end

    set_simple_constraints!(res, kwargs)

    @constraint global_model() res.flat_time >= 0
    @constraint global_model() res.rise_time >= 0
    scanner_constraints!(res)
    return res
end

Base.keys(::Trapezoid) = (Val(:rise), Val(:flat), Val(:fall))

Base.getindex(pg::BaseTrapezoid{N}, ::Val{:rise}) where {N} = ChangingGradient(N == 3 ? zeros(3) : 0., slew_rate(pg), gradient_orientation(pg), rise_time(pg), get_group(pg))
Base.getindex(pg::BaseTrapezoid, ::Val{:flat}) = ConstantGradient(gradient_strength(pg), gradient_orientation(pg), flat_time(pg), get_group(pg))
Base.getindex(pg::BaseTrapezoid, ::Val{:fall}) = ChangingGradient(gradient_strength(pg), -slew_rate(pg), gradient_orientation(pg), rise_time(pg), get_group(pg))
gradient_orientation(::BaseTrapezoid{3}) = nothing
gradient_orientation(pg::BaseTrapezoid{1}) = gradient_orientation(get_gradient(pg))
gradient_orientation(pg::Trapezoid{1}) = pg.orientation

get_group(pg::Trapezoid) = pg.group
get_group(pg::BaseTrapezoid) = get_group(get_gradient(pg))

rise_time(pg::Trapezoid) = pg.rise_time
flat_time(pg::Trapezoid) = pg.flat_time
gradient_strength(g::Trapezoid) = slew_rate(g) .* rise_time(g)
slew_rate(g::Trapezoid) = g.slew_rate
δ(g::Trapezoid) = rise_time(g) + flat_time(g)
duration(g::BaseTrapezoid) = 2 * rise_time(g) + flat_time(g)

qval(g::BaseTrapezoid, ::Nothing, ::Nothing) = δ(g) .* gradient_strength(g) .* 2π

"""
    SliceSelect(pulse; orientation=nothing, group=nothing, variables...)

Defines a trapezoidal gradient with a pulse played out during the flat time.

Parameters and variables are identical as for [`TrapezoidGradient`](@ref) with the addition of:

## Parameters
- `pulse`: sub-type of [`RFPulseComponent`](@ref) that describes the RF pulse.

## Variables
- `slice_thickness`: thickness of the selected slice in mm
"""
struct SliceSelect{N} <: BaseTrapezoid{N}
    trapezoid :: Trapezoid{N}
    pulse :: RFPulseComponent
end

function SliceSelect(pulse::RFPulseComponent; orientation=nothing, rise_time=nothing, group=nothing, slew_rate=nothing, variables...)
    res = SliceSelect(
        Trapezoid(; orientation=orientation, rise_time=rise_time, flat_time=duration(pulse), group=group, slew_rate=slew_rate),
        pulse
    )
    set_simple_constraints!(res, variables)
    return res
end

Base.keys(::SliceSelect) = (Val(:rise), Val(:flat), Val(:pulse), Val(:fall))
Base.getindex(pg::SliceSelect, ::Val{:pulse}) = pg.pulse

inverse_slice_thickness(ss::SliceSelect) = 1e3 * gradient_strength(ss.trapezoid) .* inverse_bandwidth(ss.pulse)

get_pulse(ss::SliceSelect) = ss.pulse
get_gradient(ss::SliceSelect) = ss.trapezoid

"""
    LineReadout(adc; ramp_overlap=1., orientation=nothing, group=nothing, variables...)

Defines a trapezoidal gradient with an ADC readout overlaid.

Parameters and variables are identical as for [`TrapezoidGradient`](@ref) with the addition of:

## Parameters
- `adc`: [`ADC`](@ref) object that describes the readout.
- `ramp_overlap`: how much the gradient ramp should overlap with the ADC. 0 for no overlap, 1 for full overlap (default: 1). Can be set to `nothing` to become a free variable.

## Variables
- [`fov`](@ref): FOV of the output image along this single k-space line in mm.
- [`voxel_size`](@ref): size of each voxel along this single k-space line in mm.
"""
struct LineReadout{N} <: BaseTrapezoid{N}
    trapezoid :: Trapezoid{N}
    adc :: ADC
    ramp_overlap :: VariableType
end

function LineReadout(adc::ADC; ramp_overlap=nothing, orientation=nothing, group=nothing, variables...)
    res = LineReadout(
        Trapezoid(; orientation=orientation, group=group),
        adc,
        get_free_variable(ramp_overlap)
    )
    set_simple_constraints!(res, variables)
    @constraint global_model() (res.ramp_overlap * rise_time(res.trapezoid) * 2 + flat_time(res.trapezoid)) == duration(res.adc)
    return res
end

ramp_overlap(lr::LineReadout) = lr.ramp_overlap
inverse_fov(lr::LineReadout) = 1e3 * dwell_time(lr.adc) * gradient_strength(lr.trapezoid)
inverse_voxel_size(lr::LineReadout) = 1e3 * duration(lr.adc) * gradient_strength(lr.trapezoid)

get_readout(lr::LineReadout) = lr.adc
get_gradient(lr::LineReadout) = lr.trapezoid

"""
    opposite_kspace_lines(; orientation=[1, 0, 0], kwargs...)

Return a positive and negative readout of a k-space line.
"""
function opposite_kspace_lines(; orientation=[1, 0, 0], kwargs...)
    if isnothing(orientation)
        error("orientation of k-space readout should be fixed at construction.")
    end
    positive = LineReadout(ADC(); orientation=orientation, kwargs...)
    trap = positive.trapezoid
    negative = LineReadout(
        Trapezoid1D(trap.rise_time, trap.flat_time, trap.slew_rate, -orientation, trap.group),
        positive.adc,
        positive.ramp_overlap,
    )
    return (positive, negative)
end

end