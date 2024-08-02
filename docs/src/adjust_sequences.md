# [Post-hoc adjustment of sequences](@id adjust_sequences)
Typically a sequence is repeated for multiple [`repetition_time`](@ref) ([`TR`](@ref)) to allow the acquisition of multiple k-space lines, slice selection, or other sequence parameters (e.g., varying diffusion-weighting parameters for [`DiffusionSpinEcho`](@ref)).
MRIBuilder supports this by allowing the creation of a longer sequence out of multiple repeats of a base sequence with some minor alterations.

To support post-hoc alterations, each RF pulse or gradient waveform in the sequence can be given a label. 
Some commonly-used labels in MRIBuilder are:
- `:diffusion` used for diffusion-weighted gradients produced by [`dwi_gradients`](@ref).
- `:FOV` used for gradients that should align with the field-of-view (slice-select gradients in [`excitation_pulse`](@ref) or [`refocus_pulse`](@ref) and readout gradients in [`readout_event`](@ref)).

Post-hoc alterations can be applied to gradients or RF pulses with a specific labels (or to all gradients/RF pulses) using [`adjust`](@ref).
Some example usages are:
- Reduce the RF pulse amplitude by 20% (e.g., to model the effect of transmit bias field): `adjust(sequence, pulse=(scale=0.8, ))`
- Repeat sequence 2 times with different diffusion-weighted gradient orientations (x- and y-direction) and gradient strength reduced by 30%: `adjust(sequence, diffusion=(orientation=[[1., 0., 0], [0., 1., 0.]], scale=0.7))`
- Repeat the sequence by shifting the excited slice by the given number of millimetres in the slice-select direction: `adjust(sequence, FOV=(shift=[-7.5, -2.5, 2.5, 7.5, -5., 0., 5., 10.]))`. These shifts represent an interleaved acquisition scheme, where the acquired slices/bands are 2.5 mm apart.
- Rotations defined using the [`Rotations.jl`](https://github.com/JuliaGeometry/Rotations.jl) package can be applied to gradient orientations or the field of view. For example, to rotate the field of view by 45 degrees around the y-axis:
```julia
using Rotations
rotation = Rotations.AngleAxis(deg2rad(45), 0., 1., 0.)
adjust(sequence, FOV=(rotation=rotation, ))
```
When repeating the same sequence, a spoiler gradient and/or dead time can be added in between each pair of repeats by supplying the appropriate keywords to the `merge` parameter in [`adjust`](@ref) (e.g., `merge=(wait_time=10., )`). These parameters are described in more detail in [`merge_sequences`](@ref).

## Post-hoc adjustments API
```@meta
CollapsedDocStrings = true
```
```@autodocs
Modules = [
    MRIBuilder.PostHoc,
]
```

