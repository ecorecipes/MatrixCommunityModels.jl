# Management and decisions

A [`ManagedCommunity`](@ref) is a closed community with what a decision analysis needs of a
model: named parameters (a `build` function maps the names to vital rates or matrix entries),
a management action, and named outcomes. [`project`](@ref project(::ManagedCommunity)) runs it
and returns a [`ManagedTrajectory`](@ref). None of this needs a decision package.

An action is plain data, a named tuple:

| field | meaning |
|:--|:--|
| `remove = target => fraction` | remove that fraction of a species or of some of its stages after each step (harvest, cull) |
| `take = target => number` | remove that number of individuals after each step, or all there are (a quota) |
| `add = target => number` | add that number after each step (supplementation, a release) |
| `parameters = (; name = value)` | replace parameter values; a function is applied to the value in force |
| `environment = (t, env) -> env′` | replace the environment of step `t` |
| `steps` | the steps after which `remove`, `take` and `add` apply |

A target is a species (`:prey`) or `(species, stages)`, for example `(:prey, :adult)`. The
empty tuple `(;)` is no intervention. Removal acts on the state after a step, as the removal
press of [`MCMState`](@ref) does.

## DecisionAnalysis.jl

[DecisionAnalysis.jl](https://github.com/ecorecipes/DecisionAnalysis.jl) analyses decisions
from the outcomes of a model: payoff tables, decision rules, and the value of information.
It reaches a model through four verbs. This package does not depend on it; when both are
loaded, the extension `MatrixCommunityModelsDecisionAnalysisExt` gives the verbs methods:

| verb | for a managed community |
|:--|:--|
| `with_parameters(m, θ)` | the draw `θ`, a named tuple, merged into `m.parameters` |
| `apply_action(m, action)` | the action set to `action` (`(;)` or `nothing`: none) |
| `simulate(m; rng)` | `project(m; rng)` |
| `outcome(run)` | `run.problem.outcomes(run)`, by default the final abundance of each species |

so that `outcomes(m, actions, draws; replicates, seed)` of DecisionAnalysis.jl runs the
community over parameter draws, actions and replicates, with the same seeds under every
action (common random numbers). `draws` is any Tables.jl table with one column per named
parameter, for example a named tuple of vectors. A deterministic problem
(`demographic = false`) ignores the seed and gives equal replicates.

DecisionAnalysis.jl needs Julia 1.12 or later; on Julia 1.11 the extension is not loaded and
the rest of the package is unaffected. It is not registered: add it, and FiniteKernels.jl
which it needs, by URL.

```julia
using Pkg
Pkg.add([PackageSpec(url = "https://github.com/ecorecipes/FiniteKernels.jl"),
         PackageSpec(url = "https://github.com/ecorecipes/DecisionAnalysis.jl")])
```

The tutorial "Which stage to harvest, and how hard" takes a harvest decision from the model
to the value of measuring a vital rate, and "Choosing a river flow regime across fish,
trees and insects under uncertainty" applies multi-objective Pareto analysis, deep-uncertainty
rules, certainty-equivalent VOI and PRIM scenario discovery to Verde River flow regimes.

```@autodocs
Modules = [MatrixCommunityModels]
Pages = ["src/management.jl"]
```
