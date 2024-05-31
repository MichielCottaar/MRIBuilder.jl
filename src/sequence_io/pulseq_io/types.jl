"""
Define the main types forming a [`PulseqSequence`](@ref).

Extensions and sections types are defined in their own modules.
"""
module Types

"""
    PulseqSection(:<title>)(lines)

Represents a section in the pulseq file format.
"""
struct PulseqSection{T}
    content :: Vector{String}
end

"""
    PulseqShape(samples)

Define the shape of a [`PulseqRFPulse`](@ref) or [`PulseqGradient`](@ref).
"""
struct PulseqShape
    samples :: Vector{Float64}
end

"""
Super-type for any RF pulses/gradients/ADC/extensions that can play out during a [`PulseqBlock`](@ref).
"""
abstract type AnyPulseqComponent end


"""
    PulseqRFPulse(amplitude::Number, magnitude::PulseqShape, phase::PulseqShape, time::PulseqShape, delay::Int, frequency::Number, phase_offset::Number)

An RF pulse defined in Pulseq (see [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf)).
"""
struct PulseqRFPulse <: AnyPulseqComponent
    amplitude :: Float64
    magnitude :: PulseqShape
    phase :: PulseqShape
    time :: Union{Nothing, PulseqShape}
    delay :: Int
    frequency :: Float64
    phase_offset :: Float64
end

"""
Super-type of Pulseq gradients:
- [`PulseqGradient`](@ref)
- [`PulseqTrapezoid`](@ref)
"""
abstract type AnyPulseqGradient <: AnyPulseqComponent end

"""
    PulseqGradient(amplitude::Number, shape::PulseqShape, time::PulseqShape, delay::int)

A generic gradient waveform defined in Pulseq (see [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf)).
"""
struct PulseqGradient <: AnyPulseqGradient
    amplitude :: Float64
    shape :: PulseqShape
    time :: Union{Nothing, PulseqShape}
    delay :: Int
end

"""
    PulseqTrapezoid(amplitude::Number, rise::Int, flat::Int, fall::Int, delay::Int)

A trapezoidal gradient pulse defined in Pulseq (see [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf)).
"""
struct PulseqTrapezoid <:AnyPulseqGradient
    amplitude :: Float64
    rise :: Int
    flat :: Int
    fall :: Int
    delay :: Int
end

"""
    PulseqADC(num::Int, dwell::Float64, delay::Int, frequency::Number, phase::Number)

A trapezoidal gradient pulse defined in Pulseq (see [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf)).
"""
struct PulseqADC <: AnyPulseqComponent
    num :: Int
    dwell :: Float64
    delay :: Int
    frequency :: Float64
    phase :: Float64
end

"""
    PulseqExtensionDefinition(name, content)

Abstract definition of an unknown Pulseq extension.
"""
struct PulseqExtensionDefinition{N}
    content :: Vector{String}
end

"""
    PulseqExtension(definition::PulseqExtensionDefinition, id::Int)

Reference to a specific implementation of a [`PulseqExtensionDefinition`](@ref).
"""
struct PulseqExtension{N} <: AnyPulseqComponent
    definition::PulseqExtensionDefinition{N}
    id :: Int
end

"""
    PulseqBlock(duration::Int, rf::PulseqRFPulse, gx::AnyPulseqGradient, gy::AnyPulseqGradient, gz::AnyPulseqGradient, adc::PulseqADC, ext)

Defines a Building Block with the Pulseq sequence (see [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf)).

The RF pulse, gradients, and ADC can be set to `nothing`.

The `ext` is a sequence of extension blocks that will be played out.
Set this to a sequence of zero length to not have any extensions.
"""
struct PulseqBlock
    duration :: Int
    rf :: Union{Nothing, PulseqRFPulse}
    gx :: Union{Nothing, AnyPulseqGradient}
    gy :: Union{Nothing, AnyPulseqGradient}
    gz :: Union{Nothing, AnyPulseqGradient}
    adc :: Union{Nothing, PulseqADC}
    ext :: Vector
end

"""
    PulseqSequence(version::VersionNumber, definitions::NamedTuple, blocks::Vector{PulseqBlock})

A full sequence defined according to the Pulseq [specification](https://raw.githubusercontent.com/pulseq/pulseq/master/doc/specification.pdf).
"""
struct PulseqSequence
    version:: VersionNumber
    definitions:: NamedTuple
    blocks:: Vector{PulseqBlock}
end

end