```@meta
CollapsedDocStrings = true
```

# MRIBuilder.jl internal API
```@docs
MRIBuilder
```
## Type diagram
```@eval
import MRIBuilder.Variables: AbstractBlock
import InteractiveUtils: subtypes
using AbstractTrees
import Markdown
AbstractTrees.children(x::Type) = subtypes(x)
io = IOBuffer()
print_tree(io, AbstractBlock)
seek(io, 0)
Markdown.parse("```\n" * read(io, String) * "```")
```

## Sequence components
```@autodocs
Modules = [
    MRIBuilder.Components,
    MRIBuilder.Components.AbstractTypes,
    MRIBuilder.Components.GradientWaveforms,
    MRIBuilder.Components.GradientWaveforms.NoGradientBlocks,
    MRIBuilder.Components.GradientWaveforms.ConstantGradientBlocks,
    MRIBuilder.Components.GradientWaveforms.ChangingGradientBlocks,
    MRIBuilder.Components.InstantGradients,
    MRIBuilder.Components.Pulses,
    MRIBuilder.Components.Pulses.GenericPulses,
    MRIBuilder.Components.Pulses.InstantPulses,
    MRIBuilder.Components.Pulses.ConstantPulses,
    MRIBuilder.Components.Pulses.SincPulses,
    MRIBuilder.Components.Pulses.CompositePulses,
    MRIBuilder.Components.Readouts,
    MRIBuilder.Components.Readouts.ADCs,
    MRIBuilder.Components.Readouts.SingleReadouts,
]
```
## Containers for sequence components
```@autodocs
Modules = [
    MRIBuilder.Containers,
    MRIBuilder.Containers.Abstract,
    MRIBuilder.Containers.BuildingBlocks,
    MRIBuilder.Containers.BaseSequences,
    MRIBuilder.Containers.Alternatives,
]
```

## Pre-defined sequence parts
There are [helper functions](@ref helper_functions) available to actually add these to a sequence.
```@autodocs
Modules = [
    MRIBuilder.Parts,
    MRIBuilder.Parts.Trapezoids,
    MRIBuilder.Parts.SpoiltSliceSelects,
    MRIBuilder.Parts.SliceSelectRephases,
    MRIBuilder.Parts.EPIReadouts,
]
```
## Sequence I/O
```@autodocs
Modules = [
    MRIBuilder.SequenceIO,
    MRIBuilder.SequenceIO.Pulseq,
    MRIBuilder.SequenceIO.PulseqIO,
    MRIBuilder.SequenceIO.PulseqIO.Types,
    MRIBuilder.SequenceIO.PulseqIO.Parsers,
    MRIBuilder.SequenceIO.PulseqIO.BasicParsers,
    MRIBuilder.SequenceIO.PulseqIO.Components,
    MRIBuilder.SequenceIO.PulseqIO.SectionsIO,
    MRIBuilder.SequenceIO.PulseqIO.ParseSections,
]
```

## Plot
```@autodocs
Modules = [
    MRIBuilder.Plot,
]
```