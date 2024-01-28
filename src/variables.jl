module Variables
import JuMP: @variable, Model, @objective, objective_function, owner_model, has_values, value, AbstractJuMPScalar
import ..Scanners: gradient_strength, slew_rate

all_variables_symbols = [
    # general
    :duration => (:block, "duration of the building block in ms."),
    :TR => (:sequence, "Time on which an MRI sequence repeats itself in ms."),

    # RF pulse
    :flip_angle => (:pulse, "The flip angle of the RF pulse in degrees"),
    :amplitude => (:pulse, "The maximum amplitude of an RF pulse in kHz"),
    :phase => (:pulse, "The angle of the phase of an RF pulse in KHz"),
    :frequency => (:pulse, "The off-resonance frequency of an RF pulse (relative to the Larmor frequency of water) in KHz"),
    :bandwidth => (:pulse, "Bandwidth of the RF pulse in kHz. If you are going to divide by the bandwidth, it can be more efficient to use the [`inverse_bandwidth`](@ref)."),
    :inverse_bandwidth => (:pulse, "Inverse of the [`bandwidth`](@ref) of the RF pulse in ms"),
    :N_left => (:pulse, "The number of zero crossings of the RF pulse before the main peak"),
    :N_right => (:pulse, "The number of zero crossings of the RF pulse after the main peak"),
    :slice_thickness => (:pulse, "Slice thickness of an RF pulse that is active during a gradient.")

    # gradients
    :qval => (:gradient, "The spatial range on which the displacements can be detected due to this gradient in 1/um (equivalent to [`area_under_curve`](@ref))"),
    :area_under_curve => (:gradient, "Area under the curve of the gradient in 1/um (equivalent to [`qval`](@ref))."),
    :δ => (:gradient, "Effective duration of a gradient pulse ([`rise_time`](@ref) + [`flat_time`](@ref)) in ms."),
    :rise_time => (:gradient, "Time for gradient pulse to reach its maximum value in ms."),
    :flat_time => (:gradient, "Time of gradient pulse at maximum value in ms."),
    :slew_rate => (:gradient, "maximum increase of a gradient (kHz/um/ms)"),
    :gradient_strength => (:gradient, "maximum strength of a gradient (kHz/um)"),
]

symbol_to_func = Dict{Symbol, Function}()



for (func_symbol, (block_symbol, description)) in all_variables_symbols
    as_string = "    $func_symbol($block_symbol)\n\n$description\n\nThis represents a variable within the sequence. Variables can be set during the construction of a [`BuildingBlock`](@ref) or used to create constraints after the fact."
    new_func = @eval begin
        function $func_symbol end
        @doc $as_string $func_symbol
        $func_symbol
    end
    symbol_to_func[func_symbol] = new_func
end


"""
    variables(building_block)
    variables()

Returns all functions representing properties of a [`BuildingBlock`](@ref) object.
"""
variables() = [values(symbol_to_func)...]


# Some universal truths
area_under_curve(bb) = qval(bb)


# These functions are more fully defined in building_blocks.jl
function start_time end
function end_time end
function effective_time end


const VariableType = Union{Number, AbstractJuMPScalar}


"""
    get_free_variable(model, value; integer=false)

Get a representation of a given `variable` given a user-defined constraint.
"""
get_free_variable(::Model, value::Number; integer=false) = integer ? Int(value) : Float64(value)
function get_free_variable(model::Model, value::VariableType; integer=false)
    if owner_model(value) != model
        if has_values(value)
            return value(value)
        end
        error("Cannot set any constraints between sequences stored in different JuMP models.")
    end
    return value
end
get_free_variable(model::Model, ::Nothing; integer=false) = @variable(model, integer=integer)
get_free_variable(model::Model, value::Symbol; integer=false) = integer ? error("Cannot maximise or minimise an integer variable") : get_free_variable(model, Val(value))
function get_free_variable(model::Model, ::Val{:min})
    var = get_free_variable(model, nothing)
    @objective model Min objective_function(model) + var
    return var
end
function get_free_variable(model::Model, ::Val{:max})
    var = get_free_variable(model, nothing)
    @objective model Min objective_function(model) - var
    return var
end

function bval end

end