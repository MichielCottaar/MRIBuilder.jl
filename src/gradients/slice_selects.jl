module SliceSelects
import JuMP: Model
import StaticArrays: SVector
import ...BuildingBlocks: RFPulseBlock, ContainerBlock, get_children_indices
import ...Variables: VariableType, slice_thickness, get_free_variable, effective_time, start_time, qvec, bmat_gradient, qval
import ...Variables: flip_angle, amplitude, phase, frequency, bandwidth, inverse_bandwidth, N_left, N_right
import ..ConstantGradientBlocks: ConstantGradientBlock

struct AbstractSliceSelect{T<:VariableType}
    model :: Model
    time_before :: T
    flat_before :: ConstantGradientBlock
    pulse :: RFPulseBlock
    time_after :: T
    flat_after :: ConstantGradientBlock
    flat_total :: ConstantGradientBlock
end

"""
The fixed equivalent of [`SliceSelect`](@ref).
"""
const FixedSliceSelect = AbstractSliceSelect{Float64}

"""
    SliceSelect(gradient_strength, time_before_pulse, pulse, time_after_pulse, rotate, scale)

Represents a slice-selective RF pulse.

This will typically be embedded within a larger gradient waveform.
Do not create it directly.

The fixed equivalent is [`FixedSliceSelect`](@ref).
"""
const SliceSelect = AbstractSliceSelect{VariableType}

function AbstractSliceSelect{VariableType}(model::Model, gradient_strength::SVector{3, VariableType}, time_before_pulse, pulse::RFPulseBlock, time_after_pulse, rotate::Symbol, scale::Symbol)
    time_before_pulse = get_free_variable(model, time_before_pulse),
    time_after_pulse = get_free_variable(model, time_after_pulse),
    return SliceSelect(
        model,
        time_before_pulse,
        ConstantGradientBlock(gradient_strength, time_before_pulse + effective_time(pulse), rotate, scale),
        pulse,
        time_after_pulse,
        ConstantGradientBlock(gradient_strength, duration(pulse) - effective_time(pulse) + time_after_pulse, rotate, scale),
        ConstantGradientBlock(gradient_strength, duration(pulse) + time_before_pulse + time_after_pulse, rotate, scale),
    )
end


get_children_indices(::SliceSelect) = (:flat_before, :pulse, :flat_after)
Base.getindex(select::SliceSelect, index::Symbol) = getproperty(select, index)

start_time(select::SliceSelect, symbol::Symbol) = start_time(select, Val(symbol))
start_time(select::SliceSelect, ::Val{:flat_before}) = 0.
start_time(select::SliceSelect, ::Val{:pulse}) = select.time_before_pulse
start_time(select::SliceSelect, ::Val{:flat_after}) = duration(select.flat_before)
start_time(select::SliceSelect, ::Val{:flat_total}) = 0.

qvec(select::SliceSelect) = qvec(select.flat_total)
bmat_gradient(select::SliceSelect) = bmat_gradient(select.flat_total)
duration(select::SliceSelect) = duration(select.flat_total)
gradient_strength(select::SliceSelect) = gradient_strength(select.flat_total)
slice_thickness(select::SliceSelect) = bandwidth(select) ./ gradient_strength(select)

for fn in (:flip_angle, :amplitude, :phase, :frequency, :bandwidth, :inverse_bandwidth, :N_left, :N_right)
    @eval $(fn)(select::SliceSelect) = $(fn)(select.pulse)
end

variables(::Type{<:SliceSelect}) = (duration, gradient_strength, slice_thickness, qvec, qval, flip_angle, amplitude, phase, frequency, bandwidth, inverse_bandwidth, N_left, N_right)

end