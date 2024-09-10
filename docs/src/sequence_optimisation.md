# [Sequence optimisation](@id sequence_optimisation)
In MRIBuilder an MR [`Sequence`](@ref) is defined as a sequence of [`BuildingBlock`](@ref) objects.
Most `BuildingBlock` objects will contain free parameters determining, for example, the duration of the block or the strength/orientation of the MR gradient.
In most MR sequence building software, the user will have to set all of these free parameters by computing the appropriate values given a desired echo time, b-value, etc.

In MRIBuilder the internal free parameters are not set directly. 
Instead, they are inferred using a non-linear, constrained optimisation.
For each sequence type, the developer defines how to compute various summary variables from the `BuildingBlock` free parameters, such as [`variables.echo_time`](@ref), [`variables.duration`](@ref), [`variables.resolution`](@ref), [`variables.gradient_strength`](@ref), [`variables.diffusion_time`](@ref), [`variables.duration_transverse`](@ref) etc.
A user can then create a specific instantiation of the sequence by fixing any of these summary variables to their desired values (or setting them to `:min`/`:max` to minimise/maximise them).
In addition to the user-defined constraints, this optimisation will also take into account any [scanner-defined constraints](@ref scanners).
Internally, MRIBuilder will then optimise the `BuildingBlock` free parameters to match any user-defined constraints and/or objectives.
This optimisation uses the [Ipopt](https://github.com/coin-or/Ipopt) optimiser accessed through the [JuMP.jl](https://jump.dev/JuMP.jl/stable/) library.

In addition to any user-defined objectives, the developer might also have defined secondary objectives (e.g., minimise the total sequence duration). 
These objective functions will only be considered if they do not affect the result of the user-defined primary objective.
More details on these developer-defined secondary objectives can be found in the section on [defining new sequences](@ref defining_sequences)

## [Summary variables](@id variables)
All variables are available as members of the [`variables`](@ref) structure.
```@meta
CollapsedDocStrings = true
```
```@docs
variables
```
```@docs; canonical=false
variables.N_left                    
variables.N_right                   
variables.TE                        
variables.TR                        
variables.all_gradient_strengths
variables.amplitude                 
variables.area_under_curve          
variables.bandwidth                 
variables.bmat                      
variables.bmat_gradient
variables.bval                      
variables.delay                     
variables.diffusion_time            
variables.duration                  
variables.duration_dephase
variables.duration_state            
variables.duration_transverse       
variables.dwell_time
variables.echo_time                 
variables.effective_time            
variables.flat_time                 
variables.flip_angle
variables.fov                       
variables.frequency                 
variables.gradient_strength         
variables.gradient_strength_norm    
variables.lobe_duration             
variables.net_dephasing             
variables.nsamples
variables.oversample                
variables.phase                     
variables.qval                      
variables.qvec
variables.ramp_overlap              
variables.readout_times             
variables.resolution                
variables.rise_time                 
variables.slew_rate
variables.slew_rate_norm            
variables.slice_thickness           
variables.spoiler
variables.time_to_center            
variables.voxel_size                
variables.Δ                         
variables.δ
```

## Variables interface
```@autodocs
Modules = [
    MRIBuilder.Variables,
]
```