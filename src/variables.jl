module Variables
import JuMP: @variable, Model, @objective, objective_function, value, AbstractJuMPScalar
import ..Scanners: gradient_strength, slew_rate
import ..BuildSequences: global_model

all_variables_symbols = [
    :block => [
        :duration => "duration of the building block in ms.",
    ],
    :sequence => [
        :TR => "Time on which an MRI sequence repeats itself in ms.",
    ],

    :pulse => [
        :flip_angle => "The flip angle of the RF pulse in degrees",
        :amplitude => "The maximum amplitude of an RF pulse in kHz",
        :phase => "The angle of the phase of an RF pulse in KHz",
        :frequency => "The off-resonance frequency of an RF pulse (relative to the Larmor frequency of water) in KHz",
        :bandwidth => "Bandwidth of the RF pulse in kHz. To set constraints it is often better to use [`inverse_bandwidth`](@ref).",
        :inverse_bandwidth => "Inverse of bandwidth of the RF pulse in 1/kHz. Also see [`bandwidth`](@ref).",
        :N_left => "The number of zero crossings of the RF pulse before the main peak",
        :N_right => "The number of zero crossings of the RF pulse after the main peak",
        :slice_thickness => "Slice thickness of an RF pulse that is active during a gradient in mm. To set constraints it is often better to use [`inverse_slice_thickness`](@ref).",
        :inverse_slice_thickness => "Inverse of slice thickness of an RF pulse that is active during a gradient in 1/mm. Also, see [`slice_thickness`](@ref).",
    ],

    # gradients
    :gradient => [
        :qvec => "The spatial range and orientation on which the displacements can be detected due to this gradient in rad/um.",
        :qval => "The spatial range on which the displacements can be detected due to this gradient in rad/um (i.e., norm of [`qvec`](@ref)). To set constraints it is often better to use [`qvec`](@ref) or [`qval_square`](@ref).",
        :qval_square => "Square of [`qval`](@ref) in rad^2/um^2.",
        :δ => "Effective duration of a gradient pulse ([`rise_time`](@ref) + [`flat_time`](@ref)) in ms.",
        :rise_time => "Time for gradient pulse to reach its maximum value in ms.",
        :flat_time => "Time of gradient pulse at maximum value in ms.",
        :gradient_strength => "vector with maximum strength of a gradient along each dimension (kHz/um)",
        :slew_rate => "vector with maximum slew rate of a gradient along each dimension (kHz/um)",
    ],
    :readout => [
        :dwell_time => "Time between two samples in an `ADC` in ms.",
        :nsamples => "Number of samples during a readout. During the optimisation this might produce non-integer values.",
        :fov => "Size of the field of view in mm. To set constraints it is often better to use [`inverse_fov`](@ref).",
        :inverse_fov => "Inverse of size of the field of view in 1/mm. Also see [`fov`](@ref).",
        :voxel_size => "Size of each voxel in mm. To set constraints it is often better to use [`inverse_voxel_size`](@ref).",
        :inverse_voxel_size => "Inverse of voxel size in 1/mm. Also see [`voxel_size`](@ref).",
        :resolution => "Number of voxels in the final readout. During the optimisation this might produce non-integer values, but this will be fixed after optimsation.",
        :oversample => "How much to oversample with ([`nsamples`](@ref) / [`resolution`](@ref))",
    ]
]

variables = Dict{Symbol, Function}()



for (block_symbol, all_functions) in all_variables_symbols
    for (func_symbol, description) in all_functions
        as_string = "    $func_symbol($block_symbol)\n\n$description\n\nThis represents a variable within the sequence. Variables can be set during the construction of a [`BuildingBlock`](@ref) or used to create constraints after the fact."
        new_func = @eval begin
            function $func_symbol end
            @doc $as_string $func_symbol
            $func_symbol
        end
        variables[func_symbol] = new_func
    end
end


# helper functions
function qval_square(bb, args...; kwargs...)
    vec = qvec(bb, args...; kwargs...)
    return vec[1]^2 + vec[2]^2 + vec[3]^2
end

alternative_variables = Dict(
    qval => (qval_square, n->n^2, sqrt, false),
    slice_thickness => (inverse_slice_thickness, inv, inv, true),
    bandwidth => (inverse_bandwidth, inv, inv, true),
    fov => (inverse_fov, inv, inv, true),
    voxel_size => (inverse_voxel_size, inv, inv, true),
)


# These functions are more fully defined in building_blocks.jl
function start_time end
function end_time end
function effective_time end


const VariableType = Union{Number, AbstractJuMPScalar}


"""
    get_free_variable(value; integer=false)

Get a representation of a given `variable` given a user-defined constraint.
"""
get_free_variable(value::Number; integer=false) = integer ? Int(value) : Float64(value)
get_free_variable(value::VariableType; integer=false) = value
get_free_variable(::Nothing; integer=false) = @variable(global_model(), start=0.01, integer=integer)
get_free_variable(value::Symbol; integer=false) = integer ? error("Cannot maximise or minimise an integer variable") : get_free_variable(Val(value))
function get_free_variable(::Val{:min})
    var = get_free_variable(nothing)
    model = global_model()
    @objective model Min objective_function(model) + var
    return var
end
function get_free_variable(::Val{:max})
    var = get_free_variable(nothing)
    model = global_model()
    @objective model Min objective_function(model) - var
    return var
end

"""
    bmat_gradient(gradient::GradientBlock, qstart=(0, 0, 0))

Computes the diffusion-weighting matrix due to a single gradient block in rad^2 ms/um^2.

This should be defined for every `GradientBlock`, but not be called directly.
Instead, the `bmat` and `bval` should be constrained for specific `Pathways`
"""
function bmat_gradient end

end