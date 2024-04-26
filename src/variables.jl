"""
Defines the functions that can be called on parts of an MRI sequence to query or constrain any variables.

In addition this defines:
- [`variables`](@ref): dictionary containing all variables.
- [`VariableType`](@ref): parent type for any variables (whether number or JuMP variable).
- [`get_free_variable`](@ref): helper function to create new JuMP variables.
- [`VariableNotAvailable`](@ref): error raised if variable is not defined for specific [`AbstractBlock`](@ref).
- [`set_simple_constraints`](@ref): call [`apply_simple_constraint`](@ref) for each keyword argument.
- [`apply_simple_constraint`](@ref): set a simple equality constraint.
- [`get_pulse`](@ref)/[`get_gradient`](@ref)/[`get_readout`](@ref): Used to get the pulse/gradient/readout part of a building block
- [`gradient_orientation`](@ref): returns the gradient orientation of a waveform if fixed.
"""
module Variables
import JuMP: @constraint, @variable, Model, @objective, objective_function, AbstractJuMPScalar
import StaticArrays: SVector
import ..Scanners: gradient_strength, slew_rate, Scanner
import ..BuildSequences: global_model, global_scanner, fixed

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
"""
function adjust_internal end

"""
    adjustable(block)

Returns whether a sequence, building block, or component can be adjusted

Can return one of:
- `:false`: not adjustable
- `:gradient`: expects gradient adjustment parameters
- `:pulse`: expects RF pulse adjustment parameters
"""
adjustable(::AbstractBlock) = :false



all_variables_symbols = [
    :block => [
        :duration => "duration of the building block in ms.",
    ],
    :sequence => [
        :TR => "Time on which an MRI sequence repeats itself in ms. Defaults to the result of [`repetition_time`](@ref).",
        :repetition_time => "Time on which an MRI sequence repeats itself in ms.",
        :TE => "Echo time of the sequence in ms. Defaults to the result of [`echo_time`](@ref).",
        :echo_time => "Echo time of the sequence in ms.",
        :diffusion_time => "Diffusion time in ms (i.e., time between start of the diffusion-weighted gradients).",
        :Δ => "Diffusion time in ms (i.e., time between start of the diffusion-weighted gradients). Defaults to the result of [`diffusion_time`](@ref).",
        :qvec => "Net dephasing due to gradients in rad/um.",
        :area_under_curve => "Net dephasing due to gradients in rad/um (same as [`qvec`](@ref)).",
        :bmat => "Full 3x3 diffusion-weighting matrix in ms/um^2.",
        :bval => "Size of diffusion-weighting in ms/um^2 (trace of [`bmat`](@ref)).",
        :duration_dephase => "Net T2' dephasing experienced over a sequence (in ms).",
        :duration_transverse => "Net T2 signal loss experienced over a sequence (in ms). This is the total duration spins spent in the transverse plane.",
        :delay => "Delay between readout and spin echo in an asymmetric spin echo in ms.",
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
    :gradient => [
        :qval3 => "The spatial range with orientation on which the displacements can be detected due to this gradient in rad/um (also see [`qval`](@ref)).",
        :qval => "The spatial range on which the displacements can be detected due to this gradient in rad/um. This will be a scalar if the orientation is fixed and a scalar otherwise. If you always want a vector, use [`qvec`](@ref).",
        :qval_square => "Square of [`qval`](@ref) in rad^2/um^2.",
        :δ => "Effective duration of a gradient pulse ([`rise_time`](@ref) + [`flat_time`](@ref)) in ms.",
        :rise_time => "Time for gradient pulse to reach its maximum value in ms.",
        :flat_time => "Time of gradient pulse at maximum value in ms.",
        :gradient_strength => "Vector with maximum strength of a gradient along each dimension (kHz/um)",
        :slew_rate => "Vector with maximum slew rate of a gradient along each dimension (kHz/um/ms)",
        :spoiler_scale => "Length-scale on which spins will be dephased by exactly 2π in mm.",
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
        :ramp_overlap => "Fraction of overlap between ADC event and underlying gradient pulse ramp (between 0 and 1)."
    ],
]

"""
Collection of all functions that return variables that can be used to query or constrain their values.
"""
variables = Dict{Symbol, Function}()


for (block_symbol, all_functions) in all_variables_symbols
    for (func_symbol, description) in all_functions
        as_string = "    $func_symbol($block_symbol)\n\n$description\n\nThis represents a variable within the sequence. Variables can be set during the construction of a [`AbstractBlock`](@ref) or used to create constraints after the fact."
        new_func = @eval begin
            function $func_symbol end
            @doc $as_string $func_symbol
            $func_symbol
        end
        variables[func_symbol] = new_func
    end
end


TE(ab::AbstractBlock) = echo_time(ab)
TR(ab::AbstractBlock) = repetition_time(ab)
Δ(ab::AbstractBlock) = diffusion_time(ab)

"""
Dictionary with alternative versions of specific function.
    
Setting constraints on these alternative functions can be helpful as it avoids some operations, which the optimiser might struggle with.
"""
alternative_variables = Dict(
    qval => (qval_square, n->n^2, sqrt, false),
    slice_thickness => (inverse_slice_thickness, inv, inv, true),
    spoiler_scale => (qval, q->1e-3 * 2π/q, l->1e-3 * 2π/l, true),
    bandwidth => (inverse_bandwidth, inv, inv, true),
    fov => (inverse_fov, inv, inv, true),
    voxel_size => (inverse_voxel_size, inv, inv, true),
)


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
get_free_variable(::Nothing; integer=false, start=0.01) = @variable(global_model(), start=start, integer=integer)
get_free_variable(value::Symbol; integer=false, kwargs...) = integer ? error("Cannot maximise or minimise an integer variable") : get_free_variable(Val(value); kwargs...)
function get_free_variable(::Val{:min}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    model = global_model()
    @objective model Min objective_function(model) + var
    return var
end
function get_free_variable(::Val{:max}; kwargs...)
    var = get_free_variable(nothing; kwargs...)
    model = global_model()
    @objective model Min objective_function(model) - var
    return var
end

"""
    get_pulse(building_block)]

Get the pulse played out during the building block.

Any `pulse` variables not explicitly defined for this building block will be passed on to the pulse.
"""
function get_pulse end

"""
    get_gradient(building_block)]

Get the gradient played out during the building block.

Any `gradient` variables not explicitly defined for this building block will be passed on to the gradient.
"""
function get_gradient end

"""
    get_readout(building_block)]

Get the readout played out during the building block.

Any `readout` variables not explicitly defined for this building block will be passed on to the readout.
"""
function get_readout end

"""
    bmat_gradient(gradient::GradientBlock, qstart=(0, 0, 0))

Computes the diffusion-weighting matrix due to a single gradient block in rad^2 ms/um^2.

This should be defined for every `GradientBlock`, but not be called directly.
Instead, the `bmat` and `bval` should be constrained for specific `Pathways`
"""
function bmat_gradient end


"""
    gradient_orientation(building_block)

Returns the gradient orientation.
"""
function gradient_orientation end


function effective_time end


"""
    VariableNotAvailable(building_block, variable, alt_variable)

Exception raised when a variable function does not support a specific `AbstractBlock`.
"""
mutable struct VariableNotAvailable <: Exception
    bb :: Type{<:AbstractBlock}
    variable :: Function
    alt_variable :: Union{Nothing, Function}
end
VariableNotAvailable(bb::Type{<:AbstractBlock}, variable::Function) = VariableNotAvailable(bb, variable, nothing)

function Base.showerror(io::IO, e::VariableNotAvailable)
    if isnothing(e.alt_variable)
        print(io, e.variable, " is not available for block of type ", e.bb, ".")
    else
        print(io, e.variable, " is not available for block of type ", e.bb, ". Please use ", e.alt_variable, " instead to set any contsraints or objective functions.")
    end
end


for (target_name, all_vars) in all_variables_symbols
    for (variable_func, _) in all_vars
        if variable_func in [:qval3, :TR, :TE, :Δ]
            continue
        end
        get_func = Symbol("get_" * string(target_name))
        use_get_func = target_name in (:pulse, :readout, :gradient)
        @eval function Variables.$variable_func(bb::AbstractBlock)
            try
                if Variables.$variable_func in keys(alternative_variables)
                    alt_var, forward, backward, _ = alternative_variables[Variables.$variable_func]
                    try
                        value = alt_var(bb)
                        if value isa Number
                            return backward(value)
                        elseif value isa AbstractArray{<:Number}
                            return backward.(value)
                        end
                    catch e
                        if e isa VariableNotAvailable
                            throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
                        end
                        rethrow()
                    end
                    throw(VariableNotAvailable(typeof(bb), Variables.$variable_func, alt_var))
                end
                throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
            catch e
                if $use_get_func && e isa VariableNotAvailable && hasmethod($get_func, tuple(typeof(bb)))
                    apply_to = try
                        $(get_func)(bb)
                    catch
                        throw(VariableNotAvailable(typeof(bb), Variables.$variable_func))
                    end
                    if apply_to isa AbstractBlock
                        return Variables.$variable_func(apply_to)
                    elseif apply_to isa NamedTuple
                        return NamedTuple(k => Variables.$variable_func(v) for (k, v) in pairs(apply_to))
                    elseif apply_to isa AbstractVector{<:AbstractBlock} || apply_to isa Tuple
                        return Variables.$variable_func.(apply_to)
                    end
                    error("$($(use_get_func)) returned an unexpected type for $(bb).")
                end
                rethrow()
            end
        end
    end
end


for base_fn in [:qval, :gradient_strength, :slew_rate]
    fn3 = Symbol(String(base_fn) * "3")
    @eval function $fn3(bb::AbstractBlock, args...; kwargs...)
        if hasmethod(get_gradient, (typeof(bb), ))
            return $fn3(get_gradient(bb), args...; kwargs...)
        else
            value = $base_fn(bb, args...; kwargs...)
            if value isa Number && iszero(value)
                return zero(SVector{3, Float64})
            elseif value isa AbstractVector
                return value
            else
                return value .* gradient_orientation(bb)
            end
        end
    end
    @eval $fn3(nt::NamedTuple, args...; kwargs...) = map(v -> $fn3(v, args...; kwargs...), nt)
end


"""
    set_simple_constraints!(block, kwargs)

Add any constraints or objective functions to the variables of a [`AbstractBlock`](@ref).

Each keyword argument has to match one of the functions in [`variables`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(block::AbstractBlock, kwargs)
    for (key, value) in kwargs
        if variables[key] in keys(alternative_variables)
            alt_var, forward, backward, to_invert = alternative_variables[variables[key]]
            invert_value(value::VariableType) = forward(value)
            invert_value(value::Symbol) = invert_value(Val(value))
            invert_value(::Val{:min}) = to_invert ? Val(:max) : Val(:min)
            invert_value(::Val{:max}) = to_invert ? Val(:min) : Val(:max)
            invert_value(value::AbstractVector) = invert_value.(value)
            invert_value(value) = value
            try
                apply_simple_constraint!(alt_var(block), invert_value(value))
                continue
            catch e
                if !(e isa VariableNotAvailable)
                    rethrow()
                end
            end
        end
        apply_simple_constraint!(variables[key](block), value)
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
"""
apply_simple_constraint!(variable::VariableType, ::Nothing) = nothing
apply_simple_constraint!(variable::AbstractVector, value::Symbol) = apply_simple_constraint!(sum(variable), Val(value))
apply_simple_constraint!(variable::VariableType, value::Symbol) = apply_simple_constraint!(variable, Val(value))
apply_simple_constraint!(variable::VariableType, ::Val{:min}) = @objective global_model() Min objective_function(global_model()) + variable
apply_simple_constraint!(variable::VariableType, ::Val{:max}) = @objective global_model() Min objective_function(global_model()) - variable
apply_simple_constraint!(variable::VariableType, value::VariableType) = @constraint global_model() variable == value
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

Constraints [`gradient_strength`](@ref) and [`slew_rate`](@ref) to be less than the [`global_scanner`](@ref) maximum.
"""
function scanner_constraints!(bb::AbstractBlock)
    try
        global_scanner()
    catch e
        return
    end
    for f in (slew_rate, gradient_strength)
        value = nothing
        try
            value = f(bb)
        catch e
            if e isa VariableNotAvailable
                continue
            else
                rethrow()
            end
        end
        if value isa AbstractVector
            for v in value
                @constraint global_model() v <= f(global_scanner())
                @constraint global_model() v >= -f(global_scanner())
            end
        else
            @constraint global_model() value <= f(global_scanner())
            @constraint global_model() value >= -f(global_scanner())
        end
    end
end

end