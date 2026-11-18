# Resilience metrics

Composed communities in the terms of
[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl), whose metrics dispatch
on what a matrix means:

- [`linearise`](@ref linearise(::OpenCommunity, ::AbstractVector, ::Any)) gives the Jacobian of
  one step (or of a cycle) at a fixed point as a `DiscreteJacobian`, for `isstable`,
  `resilience`, `reactivity`, `max_amplification`, the median metrics and the press
  sensitivities;
- [`periodic_operator`](@ref) gives the matrices of a rare invader over a cycle of
  environments as a `PeriodicOperator`, for its `growth_rate` and `resilience` per step;
- [`MCMState`](@ref) gives a community at a fixed point the community interface, so that
  ResilienceCore's simulated experiments (`press_response`, `removal_responses`,
  `global_mortality_tolerance`, `local_mortality_tolerances`, `extinction_tolerance`) run on
  it.

A mortality press of strength `f` can mean two things, chosen with `MCMState(...; press)`:
`press = :hazard` (the default) multiplies every stage's mortality hazard by `f` (survival
`σ → σ^f`), and `press = :removal` removes a fraction `m = removal_scale × (f - 1)` of every
stage after each step (`n ← (1 - m) advance(n)`), an added loss rather than a multiple of the
baseline mortality. Both are the model itself at `f = 1`. The hazard press falls hardest on
stages that almost all die (a 10% press halves the survival of riparian cottonwood seedlings,
which survive a year with probability about 0.0009), the removal press equally on every
stage; see the resilience tutorial.

[`jacobian`](@ref) and `linearise` are ResilienceCore's functions with methods added here, so
the two packages can be loaded together. EcologicalResilience.jl, which also builds on
ResilienceCore, exports some names that differ from this package's, such as `equilibrium` and
`subcommunity`; with both loaded, qualify them (`MatrixCommunityModels.equilibrium`). (At
present the two cannot be installed in one environment: MatrixProjectionModels.jl requires
SciMLBase 2 and EcologicalResilience.jl SciMLBase 3.)

```@autodocs
Modules = [MatrixCommunityModels]
Pages = ["src/resilience.jl"]
```
