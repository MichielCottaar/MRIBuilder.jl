module Properties

"""
    flip_angle(pulse_block)

The flip angle of the RF pulse in a [`BuildingBlock`](@ref) in degrees.
"""
function flip_angle end

"""
    phase(pulse_block)

The angle of the phase at the start of the RF pulse in a [`BuildingBlock`](@ref) in degrees.
"""
function phase end

"""
    amplitude(pulse_block)

The maximum amplitude during the RF pulse in a [`BuildingBlock`](@ref) in kHz.
"""
function amplitude end

"""
    frequency(pulse_block)

The maximum frequency during the RF pulse in a [`BuildingBlock`](@ref) relative to the Larmor frequency in kHz.
"""
function frequency end


"""
    bandwidth(pulse_block)

FWHM of the frequency content of the RF pulse in kHz.
"""
function bandwidth end

end