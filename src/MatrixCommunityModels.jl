"""
    MatrixCommunityModels

Matrix community models (Lytle & Tonkin 2023) built from the ecorecipes categorical stack.

A species is a projection net from CategoricalPopulationDynamics.jl: named stages, named
transitions, and responses saying how a transition's values depend on shared pools and
the environment. A community is a composite of species along an undirected wiring diagram
whose junctions are the pools (space, biomass, parasites per host, ...). The composite is
a `PooledProjectionSystem`, which is lowered to a coupled problem: a `CoupledIPMProblem`
solved by IntegralProjectionModels.jl when a species has a continuous stage (an integral
projection model), a `CoupledMPMProblem` solved by MatrixProjectionModels.jl otherwise.

The open communities and their composition along wiring diagrams (`species`, `community`,
`compose`, ...) are those of CategoricalPopulationDynamics.jl, re-exported here. This
package adds their projection, the analysis of the composed community, continuous stages
and kernels, and a set of published models.
"""
module MatrixCommunityModels

using LinearAlgebra
import StructuredPopulationCore
using StructuredPopulationCore: ContinuousDomain
import StructuredPopulationCore: meshpoints, step_size, quasi_extinction   # extended here, not redefined
import MatrixProjectionModels
using MatrixProjectionModels: MatrixProjectionModel, DirectIteration
import IntegralProjectionModels
import CategoricalPopulationDynamics
using CategoricalPopulationDynamics: ValuedProjectionNet, ResponsiveProjectionNet, ResponseContext,
    StageVector, Pool, PooledProjectionSystem, StateDependentMPMTarget,
    respond, inflow, after_step, births, context, expand_stage, expand_stages, stratify,
    stratum_stage, stratum_transition, transition_names, transition_matrix, to_matrix,
    evaluate, survival_fecundity_matrices, ⊕, left_kan_extension, meshes, CoupledIPMTarget
import CategoricalPopulationDynamics: stage_names, pool_names, pool_totals, advance
# open communities and their composition (CategoricalPopulationDynamics issue #4), re-exported
# below, with `species_states` (state by species), used internally
using CategoricalPopulationDynamics: OpenCommunity, species, community, compose, system, species_names,
    diagram, subcommunity, state_vector, species_states
# resilience metrics: `jacobian` and `linearise` are ResilienceCore's functions, extended here
# (resilience.jl), so that the two packages can be loaded together
import ResilienceCore
import ResilienceCore: jacobian, linearise

# the modelling vocabulary of CategoricalPopulationDynamics.jl, re-exported
export ValuedProjectionNet, ResponsiveProjectionNet, ResponseContext, StageVector, Pool,
       PooledProjectionSystem, respond, inflow, after_step, births, context, expand_stage,
       expand_stages, stratify, stratum_stage, stratum_transition, stage_names,
       transition_names, transition_matrix, to_matrix, evaluate, survival_fecundity_matrices, ⊕

export OpenCommunity, CommunityTrajectory
export species, community, compose, system, species_names, pool_names, diagram
export subcommunity, matrix_net
export ContinuousStage, ContinuousDomain, meshpoints, step_size, kernel_net, responsive, tabulated, by_trait, trait_weights
export project, advance, pool_totals, community_matrix, frozen_net, matrix_model, ipm_problem, state_vector, species_states

"""
    matrix_net(stages, name => matrix, ...)

A `ValuedProjectionNet` whose transitions are given as whole matrices (rows are
destinations, columns are sources): every non-zero entry of `matrix` becomes an entry of
the transition `name`. Use it for a discretised integral projection model, whose kernels
arrive as matrices.
"""
function matrix_net(stages::AbstractVector{Symbol}, matrices::Pair{Symbol, <:AbstractMatrix}...)
    k = length(stages)
    transitions = map(matrices) do (name, M)
        size(M) == (k, k) || throw(DimensionMismatch("matrix for :$name must be $k x $k"))
        name => [(stages[j] => stages[i]) => Float64(M[i, j]) for j in 1:k for i in 1:k if M[i, j] != 0]
    end
    return ValuedProjectionNet(collect(Symbol, stages), transitions...)
end

"""
    CommunityTrajectory

Result of [`project`](@ref). `u[t + 1, :]` is the state at time `t`; `pools[t + 1, :]`
holds the totals of the exposed pools used for the step from `t`. Index with a species
name to get its stage-by-time matrix.

A pool with a value for each of a set of classes (see [`species`](@ref)) has `NaN` in
`pools`; its values are in `pool_vectors[name][t + 1, :]`, one column per class. The last
row of both, for the final state, is `NaN`.
"""
struct CommunityTrajectory
    t::Vector{Int}
    u::Matrix{Float64}
    pools::Matrix{Float64}
    community::OpenCommunity
    pool_vectors::Dict{Symbol, Matrix{Float64}}
end
CommunityTrajectory(t, u, pools, community) = CommunityTrajectory(t, u, pools, community, Dict{Symbol, Matrix{Float64}}())

# Record the pools used for step `k`: a number in its column of `pools`, a vector in `vectors`.
function _log_pools!(pools::Matrix{Float64}, vectors::Dict{Symbol, Matrix{Float64}}, k::Int, ports, totals)
    for (j, p) in enumerate(ports)
        x = totals[p]
        if x isa Number
            pools[k, j] = x
        else
            rows = get!(() -> fill(NaN, size(pools, 1), length(x)), vectors, p)
            size(rows, 2) == length(x) || throw(DimensionMismatch("pool :$p changed its number of classes"))
            rows[k, :] = x
        end
    end
    return pools
end

Base.getindex(sol::CommunityTrajectory, sp::Symbol) = sol.u[:, sol.community.layout[sp]]

function Base.show(io::IO, sol::CommunityTrajectory)
    print(io, "CommunityTrajectory(", length(sol.t) - 1, " steps, ", size(sol.u, 2), " states)")
end

"""
    project(c, n0, envs; rng = nothing, solver = :auto)

Project a closed community from state `n0` through the environments `envs` (one per time
step). `n0` is a flat vector, or `species => stages` pairs.

The community is lowered by CategoricalPopulationDynamics.jl to a coupled problem and solved
by one of the model packages of the stack, chosen by `solver`:

- `:auto` (the default): a community with a continuous stage (a species built with
  [`kernel_net`](@ref)) is an integral projection model, lowered to a
  `IntegralProjectionModels.CoupledIPMProblem` and solved by IntegralProjectionModels.jl; a
  species without one is a component whose classes are all discrete. Any other community is
  lowered to a `MatrixProjectionModels.CoupledMPMProblem` and solved by
  MatrixProjectionModels.jl.
- `:ipm` or `:mpm`: the one or the other, whatever the stages.

Both step the same model in the same way, so they give the same trajectory.

With `rng` (a random number generator) the projection has demographic stochasticity:
individuals are whole, each survives and moves according to a multinomial draw from its
column of the survival matrix, and births and inflows are Poisson with the deterministic
values as means. This needs whole numbers in `n0`, the births of each species marked with
`births`, and survival matrices whose columns sum to at most one. A community solved by
IntegralProjectionModels.jl is then solved in the demographic mode of its coupled problem;
the others are stepped here. Both draw the numbers in the same order (for each species in
turn: survival, births, inflows), so a generator in a given state gives the same projection
whichever solver is used.
"""
function project(c::OpenCommunity, n0, envs; rng = nothing, solver::Symbol = :auto)
    solver in (:auto, :ipm, :mpm) || throw(ArgumentError("solver must be :auto, :ipm or :mpm, got :$solver"))
    start = state_vector(c, n0)
    T = length(envs)
    u = zeros(Float64, T + 1, length(c))
    pools = fill(NaN, T + 1, length(c.ports))
    vectors = Dict{Symbol, Matrix{Float64}}()
    u[1, :] = start
    T == 0 && return CommunityTrajectory([0], u, pools, c, vectors)
    by_ipm = solver === :ipm || (solver === :auto && _continuous(c))
    if by_ipm || rng === nothing
        # only the trajectory and the pools are read from the solution, so each step's matrices are not kept
        sol = if by_ipm
            stochasticity = rng === nothing ? IntegralProjectionModels.Deterministic() : IntegralProjectionModels.Demographic()
            target = CoupledIPMTarget(species_states(c, start), (0, T); p = collect(envs), stochasticity)
            IntegralProjectionModels.solve(CategoricalPopulationDynamics.lower(c.system, target),
                IntegralProjectionModels.DirectIteration(); save_matrices = false, rng)
        else
            target = StateDependentMPMTarget(species_states(c, start), (0, T); p = collect(envs))
            MatrixProjectionModels.solve(CategoricalPopulationDynamics.lower(c.system, target), DirectIteration();
                save_matrices = false)
        end
        for k in 1:(T + 1), sp in c.names
            u[k, c.layout[sp]] = sol.u[k][sp]
        end
        logged = sol.substep_log[:pools]
        for k in 1:T
            _log_pools!(pools, vectors, k, c.ports, logged[k])
        end
    else
        n = start
        for (k, env) in enumerate(envs)
            _log_pools!(pools, vectors, k, c.ports, pool_totals(c, n, env))
            n = advance(c, n, env; rng)
            u[k + 1, :] = n
        end
    end
    return CommunityTrajectory(collect(0:T), u, pools, c, vectors)
end

# Whether a community holds a species with a continuous stage, and is therefore an integral projection model.
_continuous(c::OpenCommunity) = any(sp -> !isempty(meshes(c.system[sp])), c.names)

_context(c::OpenCommunity, sp::Symbol, n::AbstractVector, env) =
    context(c.system[sp], state_vector(c, n)[c.layout[sp]], pool_totals(c.system, species_states(c, n), env), env)

"""
    frozen_net(c, species, n, env)

The `ValuedProjectionNet` of `species` with its responses evaluated in community state
`n` and environment `env`: the species frozen at one moment, every value still attached
to its named transition.
"""
frozen_net(c::OpenCommunity, sp::Symbol, n::AbstractVector, env) = evaluate(c.system[sp], _context(c, sp, n, env))

"""
    matrix_model(c, species, n, env)

The `MatrixProjectionModels.MatrixProjectionModel` of `species` frozen in community state
`n` and environment `env`, with births (see `births`) kept apart from survival, so that
the single-species analyses of MatrixProjectionModels.jl apply.
"""
function matrix_model(c::OpenCommunity, sp::Symbol, n::AbstractVector, env)
    U, F = survival_fecundity_matrices(c.system[sp], _context(c, sp, n, env))
    return MatrixProjectionModel(U, F; stage_names = stage_names(c, sp))
end

"""
    community_matrix(c, n, env)

The block-diagonal projection matrix of a closed community frozen at state `n` and
environment `env`, one block per species.
"""
function community_matrix(c::OpenCommunity, n::AbstractVector, env)
    A = zeros(Float64, length(c), length(c))
    for sp in c.names
        r = c.layout[sp]
        A[r, r] = to_matrix(c.system[sp], _context(c, sp, n, env))
    end
    return A
end

"""
    ipm_problem(c, species, n, env; tspan = (0, 1))

The `IntegralProjectionModels.IPMProblem` of `species` frozen in community state `n` and
environment `env`. Its kernel is an `IntegralProjectionModels.DiscretizedKernel`: the
species' projection matrix at that moment, with survival and births kept apart, laid out
over its stages, a continuous stage (see [`ContinuousStage`](@ref)) as one continuous state
and every other stage as a discrete state of its own. Solve it with
IntegralProjectionModels.jl, for example `solve(ipm_problem(c, :plant, n, env),
EigenAnalysis())` for the growth rate, stable state and reproductive values of the species
in that state. The starting population of the problem is the species' block of `n`.
"""
function ipm_problem(c::OpenCommunity, sp::Symbol, n::AbstractVector, env; tspan::Tuple{Int, Int} = (0, 1))
    IPM = IntegralProjectionModels
    net = c.system[sp]
    U, F = survival_fecundity_matrices(net, _context(c, sp, n, env))
    states = CategoricalPopulationDynamics.state_layout(net)
    return IPM.IPMProblem(IPM.GeneralIPM(), IPM.DiscretizedKernel(U, F, states), states, state_vector(c, n)[c.layout[sp]], tspan)
end

include("kernels.jl")

include("recruitment.jl")
export lottery, lottery_demand

include("analysis.jl")
export equilibrium, jacobian, invasion_growth_rate, invasion_sensitivity, mean_abundance, abundance_sensitivity
export equilibrium_sensitivity, stochastic_invasion_rate, periodic_point
export InteractionNetwork, links, connectance, link_density, reciprocity, out_degree, in_degree
export periodic_matrix, abundance, stage_total, CommunityEnsemble, ensemble, extinction_times, quasi_extinction

include("resilience.jl")
export linearise, periodic_operator, MCMState

include("models.jl")
export verde_fish_model, riparian_model, cypripedium_model, with_rate, VERDE_SPECIES, RIPARIAN_GUILDS
export density_response

include("bees.jl")
export day_of_year, hopomo_season, forage_index, flight_index
export honeybee, honeybee_state, varroa, bumblebee, floral_supply, floral_competition, flower_transmission
export HONEYBEE_PARAMETERS, VARROA_PARAMETERS, BUMBLEBEE_PARAMETERS

include("thistles.jl")
export platte_thistle, platte_thistle_stage, musk_thistle, flowerhead_weevil
export PLATTE_THISTLE_PARAMETERS, MUSK_THISTLE_PARAMETERS, WEEVIL_PARAMETERS

include("forests.jl")
export forest_tree, forest_stage, forest_stand, forest_vital_rates
export FOREST_PARAMETERS, FOREST_CLIMATE
include("seafans.jl")
export sea_fan, sea_fan_stage, SEA_FAN_PARAMETERS, SEA_FAN_PHASES
include("sagebrush.jl")
export sagebrush_plant, sagebrush_stage, sagebrush_community, sagebrush_crowding, sagebrush_state, SAGEBRUSH_PARAMETERS
include("lemans.jl")
export lemans_parameters, lemans_fish, lemans_community, lemans_state, lemans_environment, lemans_biomass,
       lemans_predation, LEMANS_NORTH_SEA
include("corals.jl")
export coral, coral_stage, coral_community, coral_recruits, CORAL_PARAMETERS, CORAL_RECRUITMENT_CAP
include("grasslands.jl")
export grassland_plant, grassland_community, grassland_state, synchrony, GRASSLAND_PARAMETERS
include("tribolium.jl")
export flour_beetle, flour_beetle_competition, TRIBOLIUM_EDMUNDS_2003, TRIBOLIUM_DESHARNAIS_2005
include("flowregimes.jl")
export flood_invertebrate, flow_invertebrates, invertebrate_flows, flood_pulses, tonkin_scenario, tonkin_fish_start,
       FLOW_INVERTEBRATES, TONKIN_SCENARIOS, TONKIN_FISH_INITIAL
include("predatorprey.jl")
export predator_prey, BARRAQUAND_GIMENEZ_2019, BARRAQUAND_GIMENEZ_2019_SET2
include("chytrid.jl")
export CHYTRID_SITES, CHYTRID_ZOOSPORE_SURVIVAL, chytrid_temperature, chytrid_environment, chytrid_load_stage,
       leopard_frog, chytrid_zoospores, chytrid_community
include("mizer.jl")
export MIZER_NORTH_SEA, mizer_parameters, mizer_fish, mizer_resource, mizer_community, mizer_state,
       mizer_environment, mizer_biomass, mizer_rates

end # module MatrixCommunityModels
