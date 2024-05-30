# [Scanners](@id scanners)
The MRI scanner that is used during acquisition puts various constraints on the MR sequences that can be used.
These constraints include safety considerations, such as tissue heating, and hardware constraints, such as maximum gradient strength and slew rate.
Currently, MRIBuilder only considers the latter.

To define a sequence appropriate for a specific scanner, a user would define a new [`Scanner`](@ref) with the appropriate [`B0`](@ref), maximum [`variables.gradient_strength`](@ref), and maximum [`variables.slew_rate`](@ref).
This scanner would then be passed on to the [sequence optimisation](@ref sequence_optimisation).

For ease of use, the `gradient_strength` and `slew_rate` of many scanners have already been pre-defined.
These are listed below.

## Scanners API
```@meta
CollapsedDocStrings = true
```
```@autodocs
Modules = [
    MRIBuilder.Scanners,
]
```