# Using implemented sequences (@id implemented_sequences)

## Usage

MRIBuilder comes with several sequences pre-defined.
Each sequence can be created through a simple function call.
To get help on a specific sequence, either follow the link in the sequence list below or type `?<function_name>` in Julia.

When reading the help, you will notice that there are two type of expected inputs:
- `Parameters`: these define the type of components that will be includes, e.g., the shape of the excitation pulse or the readout strategy. These parameters have to be set to certain fixed values. If unset, they will be determined by their default value as defined in the documentation.
- `Variables`: These are a special type of parameters. In addition to being set to a fixed value, they can also be set to `:min` or `:max` to minimise or maximise the variable. If they are unset, they will be determined based on the other variables.
Internally whenever you call a sequence function it will run a non-linear optimisation given any constraints and/or objectives you have given to the `Variables`. It will give a single sequence that matches all of the constraints. If no solution exists, an error will be raised instead. If multiple solutions exist, a single (somewhat random) one is returned.

## Available sequences
- [`GradientEcho`](@ref)
- [`SpinEcho`](@ref)
- [`DWI`](@ref)/[`DW_SE`](@ref)/[`DiffusionSpinEcho`](@ref)