module MatrixCommunityModelsDecisionAnalysisExt

# The adapter between this package and DecisionAnalysis.jl: methods of its adapter verbs for
# `ManagedCommunity` and `ManagedTrajectory` (src/management.jl), so that
#
#     outcomes(m, actions, draws; replicates, seed)
#
# runs a managed community over parameter draws × actions × replicates with common random numbers.
# The types and the projection are this package's own and need no decision package; the methods
# here only name them in DecisionAnalysis's vocabulary. DecisionAnalysis.jl needs Julia 1.12, so
# on Julia 1.11 the extension is never loaded.

using MatrixCommunityModels: ManagedCommunity, ManagedTrajectory, project
import DecisionAnalysis

"""
    DecisionAnalysis.with_parameters(m::ManagedCommunity, θ)

The managed community with the values of the draw `θ` (a named tuple, one row of the parameter
draws) merged into its parameters. The values reach the model through `m.build`.
"""
DecisionAnalysis.with_parameters(m::ManagedCommunity, θ) =
    ManagedCommunity(m; parameters = merge(m.parameters, NamedTuple(θ)))

"""
    DecisionAnalysis.apply_action(m::ManagedCommunity, action)

The managed community under `action`, a named tuple with the fields `remove`, `take`, `add`,
`parameters`, `environment` and `steps` documented with `ManagedCommunity`; `(;)` or `nothing`
is no intervention, the baseline.
"""
DecisionAnalysis.apply_action(m::ManagedCommunity, action) = ManagedCommunity(m; action)

"""
    DecisionAnalysis.simulate(m::ManagedCommunity; rng)

One projection of the managed community, `project(m; rng)`. A problem with `demographic = false`
is deterministic and does not touch `rng`.
"""
DecisionAnalysis.simulate(m::ManagedCommunity; rng) = project(m; rng)

"""
    DecisionAnalysis.outcome(run::ManagedTrajectory)

The named outcomes of a run: `run.problem.outcomes(run)`, by default the final abundance of each
species.
"""
DecisionAnalysis.outcome(run::ManagedTrajectory) = run.problem.outcomes(run)

end
