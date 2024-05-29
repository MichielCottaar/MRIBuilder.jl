# [Defining new sequences](@id defining_sequences)

Examples of sequence definitions can be found [here](`https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl/-/tree/main/src/sequences?ref_type=heads`).
These are the actual sequence definitions used for the [pre-implemented sequences](@ref implemented_sequences) in MRIBuilder. In this example, we will look at the file for [`DiffusionSpinEcho`](https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl/-/blob/main/src/sequences/diffusion_spin_echoes.jl?ref_type=heads) as an example. 

We can see that there are several steps.
First we define the sequence as a specific sub-type of [`Sequence`](@ref):
```julia
const DiffusionSpinEcho = Sequence{:DiffusionSpinEcho}
```
Note the duplication of the sequence name (`DiffusionSpinEcho`).
This name will have to be unique.

The next step is to define the sequence constructor:
```julia
function DiffusionSpinEcho(; scanner=DefaultScanner, parameters..., variables...)
    build_sequence(scanner) do
        (g1, g2) = dwi_gradients(...)
        return Sequence([
            :excitation => excitation_pulse(...),
            nothing,
            :gradient => g1,
            nothing,
            :refocus => refocus_pulse(...),
            nothing,
            g2,
            nothing,
            :readout => readout_event(...),
        ], name=:DiffusionSpinEcho, variables...)
    end
end
```
The crucial part here are the individual parts used to build the sequence, defined as a vector:
```julia
[
    :excitation => excitation_pulse(...),
    nothing,
    :gradient => g1,
    nothing,
    :refocus => refocus_pulse(...),
    nothing,
    :gradient2 => g2,
    nothing,
    :readout => readout_event(...),
]
```

We can see that this sequence is built in order by an excitation pulse, some unknown dead time (indicated by `nothing`), a gradient, more dead time, a refocus pulse, more dead time, another gradient, more dead time, and finally a readout.

Some of these components have been given specific names (e.g., `:excitation => ...`). This is optional, but can be useful to refer to individual components. There are [helper functions](@id helper_functions) available to create these components.

The next step is to define [summary variables](@ref variables) that the user can constrain when setting up a specific instance of this sequence:
```julia
@defvar begin
    diffusion_time(ge::DiffusionSpinEcho) = start_time(ge, :gradient2) - start_time(ge, :gradient)
    echo_time(ge::DiffusionSpinEcho) = 2 * (variables.effective_time(ge, :refocus) - variables.effective_time(ge, :excitation))
end
```
For this sequence, we can see that we define the [`diffusion_time`](@ref) as the time between the start of the first and second gradient pulse, and the [`echo_time`](@ref) as twice the time between the refocus and excitation pulses.
These variables need to be defined within a [`@defvar`](@ref) block.
They will be available as function calls within the [`variables`](@ref) object.

In addition to these sequence-specific summary variables, there are also a lot of variables already pre-defined on individual components, such as the [`slice_thickness](@ref) of the RF pulse or the [`gradient_strength`](@ref) of the gradient pulses. To access these summary variables on a sequence-level, we need to tell MRIBuilder for which RF pulses/gradients we are interested in computing these variables:
```julia
get_pulse(seq::DiffusionSpinEcho) = (excitation=seq[:excitation], refocus=seq[:refocus])
get_gradient(seq::DiffusionSpinEcho) = (gradient=seq[:gradient], gradient2=seq[:gradient2])
get_readout(seq::DiffusionSpinEcho) = seq.readout
```
Note that we can indicate we are interested in multiple RF pulses/gradients by supplying them as a named tuple (`(excitation=..., refocus=...)`).

Setting this allows us to get RF pulse or gradient-specific properties by calling the sequence, for example:
```@example
using MRIBuilder
sequence = DiffusionSpinEcho(bval=1., TE=:min, slice_thickness=2.)
flip_angle(sequence)
```
Here we can see that we get the [`flip_angle`](@ref) for each of the two RF pulses defined using [`get_pulse`](@ref) above.

The final component to defining summary variables is to define one or more default coherence pathways using [`get_pathway`]:
```julia
get_pathway(seq::DiffusionSpinEcho) = Pathway(seq, [90, 180], 1, group=:diffusion)
```
Here the coherence [`Pathway`](@ref) sets out what a specific set of spins might experience during the sequence.
In this case the sequence experiences two RF pulses and is excited by the first pulse and flipped by the second (`[90, 180]`).
It is then observed during the first readout (`1`). 
For such a [`Pathway`](@ref) we can compute:
- the time that the spin spends in each longitudinal and transverse direction, which is particularly useful in the transverse direction to compute the amount of T2-weighting ([`duration_transverse`](@ref)) and the amount of time spent dephasing to compute the amount of T2*-weighting ([`duration_dephase`](@ref)).
- the diffusion weighting experienced ([`bval`](@ref), [`bmat`](@ref), and [`qvec`](@ref))
By defining a default pathway for the sequence, the user can now put constraints on any or all of these variables.


## [Component helper functions](@id helper_functions)
```@meta
CollapsedDocStrings = true
```
```@autodocs
Modules = [
    MRIBuilder.Parts.HelperFunctions,
]
```

## Optimisation helper functions
```@autodocs
Modules = [
    MRIBuilder.BuildSequences,
]
```


## Pathways API
```@autodocs
Modules = [
    MRIBuilder.Pathways,
]
```