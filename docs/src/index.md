# MRIBuilder
[MRIBuilder](https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl) allows for the creation and optimisation of MRI sequences within [Julia](https://julialang.org).

Depending on your application, there are several levels at which you can interact with MRIBuilder.
The ones lower down require more expertise with Julia and the internals of MRIBuilder:
1. Many sequences have already been implemented and can be obtained through a simple function call (see [Using implemented sequences](@ref implemented_sequences)).
2. New sequences can be created out of pre-defined sequence components and by defining sequence-specific metrics (see [Defining new sequences](@ref defining_sequences)).
3. Finally, one can actually define new sequence components (not documented yet).

Typically, the resulting sequence will only cover a single repetition time (TR).
MRIBuilder enables the concatenation of single-TR sequences into a multi-TR sequence.
During these repeats minor adjustments can be made to the single-TR sequence.
This can be used to allow different repeats to image different lines in k-space or
excite different slices (see [Post-hoc adjustment of sequences](@ref adjust_sequences)).

The signal formation for the resulting sequence can be predicted using [MCMRSimulator](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/stable) given some representation of the imaged tissue.
MRIBuilder can be used to read/write to the [pulseq](https://pulseq.github.io) MR sequence file format.
This can be used to run the sequence on MRI scanners as described in the [pulseq homepage](https://pulseq.github.io).
Rather than just directly running the sequences from this library on the scanner, we strongly recommend to load it using the [MATLAB pulseq](https://github.com/pulseq/pulseq) or [python pypulseq](https://github.com/imr-framework/pypulseq) first as these libraries run additional checks!


## Installation
It can be run from the command line using the Julia REPL, from a Julia script, or in a [Jupyter notebook](https://jupyter.org).
Like any Julia package, Julia can be installed using the built-in [Julia package manager](https://pkgdocs.julialang.org/v1/):
```
pkg> add https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl.git
```
