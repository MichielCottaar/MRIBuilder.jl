# [Using implemented sequences](@id implemented_sequences)

## Usage

MRIBuilder comes with several sequences pre-defined.
Each sequence can be created through a simple function call.
To get help on a specific sequence, either follow the link in the sequence list below or type `?<function_name>` in Julia.

When reading the help, you will notice that there are two type of expected inputs:
- `Parameters`: these define the type of components that will be included, e.g., the shape of the excitation pulse or the readout strategy. These parameters have to be set to certain fixed values. If not set, they will be determined by their default value as defined in the documentation.
- `Variables`: These are a special type of parameters. In addition to being set to a fixed value, they can also be set to `:min` or `:max` to minimise or maximise the variable. If they are not set, they will be determined based on the optimisation or fixing of other variables. For more details, see the section on [sequence optimisation](@ref sequence_optimisation). All variables are available in the [`variables`](@ref) module. They can also be accessed as properties of the specific pulse/gradient/readout/block (i.e., `block.<variable_name>`).

As an example, the following creates and plots a [`DiffusionSpinEcho`](@ref) that has a b-value of 1, a minimal echo time, and a slice thickness of 2 mm:
```@example
using MRIBuilder
using CairoMakie
sequence = DiffusionSpinEcho(bval=1., TE=:min, slice_thickness=2)
f = plot(sequence)
f
save("dwi_1_min_2.png", f) # hide
nothing # hide
```
![DWI sequence diagram](dwi_1_min_2.png)

If we want a specific [`variables.diffusion_time`](@ref), we can just add it to the constraints, and the rest of the sequence will adapt as needed:
```@example
using MRIBuilder # hide
using CairoMakie # hide
sequence = DiffusionSpinEcho(bval=1., diffusion_time=80, TE=:min, slice_thickness=2)
f = plot(sequence)
f
save("dwi_1_80_min_2.png", f) # hide
nothing # hide
```
![DWI sequence diagram with fixed diffusion time](dwi_1_80_min_2.png)

We can even directly set some aspect of one of the sequence components, such as slowing down the gradient [`variables.rise_time`](@ref)
and the additional constraint will just be included in the [sequence optimisation](@ref sequence_optimisation):
```@example
using MRIBuilder # hide
using CairoMakie # hide
sequence = DiffusionSpinEcho(bval=1., diffusion_time=80, TE=:min, slice_thickness=2, gradient=(rise_time=15, ))
f = plot(sequence)
f
save("dwi_1_80_min_2_15.png", f) # hide
nothing # hide
```
![DWI sequence diagram with fixed diffusion time and rise time](dwi_1_80_min_2_15.png)

Note that the previous sequences do not contain a realistic readout.
Most sequences will only include an instant readout, unless you directly set the [`variables.voxel_size`](@ref) and [`variables.resolution`](@ref).
```@example
using MRIBuilder # hide
using CairoMakie # hide
sequence = DiffusionSpinEcho(bval=1., TE=:min, voxel_size=2, resolution=(20, 20, 20))
f = plot(sequence)
f
save("dwi_1_80_min_2_15_epi.png", f) # hide
nothing # hide
```
![DWI sequence diagram with EPI readout](dwi_1_80_min_2_15_epi.png)

## Available sequences
```@meta
CollapsedDocStrings = true
```
```@autodocs
Modules = [
    MRIBuilder.Sequences,
    MRIBuilder.Sequences.GradientEchoes,
    MRIBuilder.Sequences.SpinEchoes,
    MRIBuilder.Sequences.DiffusionSpinEchoes,
]
```
