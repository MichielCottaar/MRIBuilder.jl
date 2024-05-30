# [Sequence optimisation](@id sequence_optimisation)
In MRIBuilder an MR [`Sequence`](@ref) is defined as a sequence of [`BuildingBlock`](@ref) objects.
Most `BuildingBlock` objects will contain free parameters determining, for example, the duration of the block or the strength/orientation of the MR gradient.
In most MR sequence building software, the user will have to set all of these free parameters by computing the appropriate values given a desired echo time, b-value, etc.

In MRIBuilder the internal free parameters are not set directly. 
Instead, they are inferred using a non-linear, constrained optimisation.
For each sequence type, the developer defines how to compute various summary variables from the `BuildingBlock` free parameters, such as [`echo_time`](@ref), [`duration`](@ref), [`resolution`](@ref), [`gradient_strength`](@ref), [`diffusion_time`](@ref), etc.
A user can then create a specific instantiation of the sequence by fixing any of these summary variables to their desired values (or setting them to `:min`/`:max` to minimise/maximise them).
In addition to the user-defined constraints, this optimisation will also take into account any [scanner-defined constraints](@ref scanners).
Internally, MRIBuilder will then optimise the `BuildingBlock` free parameters to match any user-defined constraints and/or objectives.
This optimisation uses the [Ipopt](https://github.com/coin-or/Ipopt) optimiser accessed through the [JuMP.jl](https://jump.dev/JuMP.jl/stable/) library.

## [Summary variables](@id variables)
All variables are available as members of the [`variables`](@ref) structure.
```@meta
CollapsedDocStrings = true
```
```@docs
variables.
```

## Variables interface
```@autodocs
Modules = [
    MRIBuilder.Variables,
]
```