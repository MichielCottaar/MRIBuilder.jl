"""
Defines the Pulseq IO extension interface with default implementations
"""
module Extensions
import ..Types: PulseqExtension, PulseqExtensionDefinition

struct UnknownExtensionMapper
    ext::PulseqExtensionDefinition
end

Base.getindex(mapper::UnknownExtensionMapper, i::Int) = PulseqExtension(mapper.ext, i)

"""
    parse_extension(ext::PulseqExtensionDefinition{name})

Parse a Pulseq extension definition into a dictionary-like object
that maps integer reference IDs to any object describing the extension.

This can be overriden to support the reading of a specific extension type.
For example, to define a parser for an extension with the name "LABELSET":
```
function PulseqIO.parse_extension(ext::PulseqExtensionDefinition{:LABELSET})
    ...
end
```
"""
function parse_extension(ext::PulseqExtensionDefinition{N}) where {N}
    @warn "Parsing unknown extension: {N}"
    return UnknownExtensionMapper(ext)
end


"""
    get_extension_name(obj)

Get the name under which the given `obj` should be stored in a Pulseq extension.

To write an object to a Pulseq file extension, 
one needs to define both this function and [`add_extension_definition`](@ref).
"""
get_extension_name(::PulseqExtension{N}) where {N} = N

"""
    add_extension_definition!(content::Vector{String}, obj)

Add the object to the extension definition and returns the reference index.

The extension definition is passed on as a vector of strings.
This vector can be appended to, when adding a new object.

If the object is already in the `definition` the reference index of the already existing object should be returned instead.
"""
add_extension_definition!(content::Vector{String}, ext::PulseqExtension) = ext.index

end