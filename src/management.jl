# A community under management: a model built from named parameters, projected with interventions
# on its state (removal, supplementation), its parameters or its environments, and read off as named
# outcomes. Nothing here depends on a decision package: the weak-dependency extension
# ext/MatrixCommunityModelsDecisionAnalysisExt.jl gives DecisionAnalysis.jl's adapter verbs methods
# for these types, and they are all it needs.
#
# An action is plain data, a named tuple, read by `_action` below. Removal follows the removal press
# of `MCMState` (resilience.jl): it acts on the state after a step, newborns included.

const _ACTION_FIELDS = (:remove, :take, :add, :parameters, :environment, :steps)

"""
    ManagedCommunity(build; parameters = (;), action = (;), outcomes = final abundances,
                     demographic = false)
    ManagedCommunity(community, n0, envs; kwargs...)
    ManagedCommunity(m::ManagedCommunity; parameters, action, outcomes, demographic)

A closed community with named parameters, a management action and named outcomes: what a decision
analysis needs of a model. [`project`](@ref)`(m; rng)` runs it and returns a [`ManagedTrajectory`](@ref).

**Parameters.** `build(θ)` takes a named tuple of parameter values and returns the model at those
values as a named tuple with `community`, `n0` and `envs` (further fields are ignored), which is what
the model functions of this package return or nearly so, for example

    θ -> (m = predator_prey(merge(BARRAQUAND_GIMENEZ_2019, (adult_survival = (prey = θ.prey_survival,
                                                                             predator = 0.7),)));
          (community = m.community, n0 = m.n0, envs = fill(nothing, 20)))

So the mapping from names to vital rates or matrix entries is `build` itself, written once and in the
open. `parameters` holds the values in force; names that `build` does not read are ignored by it.
`build` is called once for every projection, so keep it cheap or memoise it.

The second form is for a community whose responses read their rates from the environment (`ctx.env`,
as the Verde River fish and the ducks do): the parameters are merged into the environment of every
step, which must then be a named tuple.

**Actions.** `action` is a named tuple with any of the fields below; `(;)` or `nothing` is no
intervention.

- `remove = target => fraction`: remove that fraction of the target after each step, as harvesting
  or culling would. A target is a species (`:prey`) or `(species, stages)` with a stage or group
  name as understood by [`stage_total`](@ref), or a vector of them (`(:prey, :adult)`). Several
  removals are a vector of pairs.
- `take = target => number`: remove that number of individuals from the target after each step, or
  all of them if there are fewer, spread over the target's states in proportion to their numbers.
- `add = target => number`: add that number to the target after each step (supplementation, a
  release). A target with several states needs a vector, one number per state.
- `parameters = (; name = value, ...)`: replace parameter values; a function is applied to the value
  in force, as in `(; fecundity = x -> 0.8x)`. This is the place for an intervention that changes a
  vital rate.
- `environment = (t, env) -> env′`: replace the environment of step `t`, for an intervention that
  acts through it (a designed flow, a bag limit).
- `steps`: the steps after which `remove`, `take` and `add` apply (default: every step).

Interventions on the state happen at the census that ends a step, on the state that
[`advance`](@ref) returned, in the order `remove`, `take`, `add`; what was removed and added is
recorded by state and time in the result.

**Outcomes.** `outcomes(run)` receives the [`ManagedTrajectory`](@ref) and returns a named tuple of
numbers, for example `run -> (; prey = abundance(run.trajectory, :prey)[end], harvested =
sum(run.removed))`. The default is the final abundance of each species, by name.

**Chance.** With `demographic = true` the projection has demographic stochasticity (see
[`project`](@ref)) and needs `rng`; `remove` then takes each individual with the given probability
and `take` draws individuals without replacement. `rng` itself only seeds other generators: its
first draw seeds the one used by the interventions, and one further draw seeds the generator of each
step. These draws are made whatever the action, so two actions projected with generators in the same
state meet the same chance: every step starts from the same seed under both, and the species that are
stepped before the first one whose state differs draw the same numbers.

The third form copies a problem with some fields replaced. The type is immutable.
"""
struct ManagedCommunity{B, O}
    build::B
    parameters::NamedTuple
    action::NamedTuple
    outcomes::O
    demographic::Bool
end

function ManagedCommunity(build; parameters = (;), action = (;), outcomes = _final_abundances,
                          demographic::Bool = false)
    return ManagedCommunity(build, _named(parameters, "parameters"), _action(action), outcomes, demographic)
end

function ManagedCommunity(c::OpenCommunity, n0, envs; kwargs...)
    return ManagedCommunity(_EnvironmentParameters(c, n0, collect(Any, envs)); kwargs...)
end

function ManagedCommunity(m::ManagedCommunity; parameters = m.parameters, action = m.action, outcomes = m.outcomes,
                          demographic::Bool = m.demographic)
    return ManagedCommunity(m.build; parameters, action, outcomes, demographic)
end

Base.show(io::IO, m::ManagedCommunity) =
    print(io, "ManagedCommunity(", length(m.parameters), " parameters set; action ",
          isempty(m.action) ? "none" : join(keys(m.action), ", "), m.demographic ? "; demographic" : "", ")")

# the `build` of a community whose parameters are read from the environment
struct _EnvironmentParameters
    community::OpenCommunity
    n0::Any
    envs::Vector{Any}
end
function (b::_EnvironmentParameters)(θ::NamedTuple)
    isempty(θ) && return (community = b.community, n0 = b.n0, envs = b.envs)
    all(env -> env isa NamedTuple, b.envs) || throw(ArgumentError(
        "parameters $(keys(θ)) are merged into the environment of every step, which must then be a named tuple; " *
        "give `ManagedCommunity` a `build` function instead"))
    return (community = b.community, n0 = b.n0, envs = Any[merge(env, θ) for env in b.envs])
end

_named(x::NamedTuple, what) = x
_named(::Nothing, what) = (;)
_named(x, what) = throw(ArgumentError("$what must be a named tuple, got $(typeof(x))"))

# An action checked for its field names; the targets are checked against the community when it is projected.
function _action(a)
    action = _named(a, "an action")
    for k in keys(action)
        k in _ACTION_FIELDS ||
            throw(ArgumentError("unknown field :$k in an action; the fields are $(join(_ACTION_FIELDS, ", "))"))
    end
    haskey(action, :parameters) && _named(action.parameters, "the `parameters` of an action")
    return action
end

_final_abundances(run) = NamedTuple{Tuple(run.trajectory.community.names)}(
    Tuple(sum(run.trajectory.u[end, run.trajectory.community.layout[sp]]) for sp in run.trajectory.community.names))

"""
    ManagedTrajectory

Result of [`project`](@ref) on a [`ManagedCommunity`](@ref): `trajectory`, the
[`CommunityTrajectory`](@ref) (states after the interventions of each census); `removed` and `added`,
matrices laid out as `trajectory.u`, with in row `t + 1` what the action removed from and added to each
state at time `t` (row 1, the start, is zero); and `problem`, the `ManagedCommunity` that was run.
`run.trajectory.community.layout[species]` gives the columns of a species.
"""
struct ManagedTrajectory{M <: ManagedCommunity}
    trajectory::CommunityTrajectory
    removed::Matrix{Float64}
    added::Matrix{Float64}
    problem::M
end

Base.show(io::IO, run::ManagedTrajectory) =
    print(io, "ManagedTrajectory(", length(run.trajectory.t) - 1, " steps, ", size(run.trajectory.u, 2),
          " states; removed ", sum(run.removed), ", added ", sum(run.added), ")")

# the parameters in force under an action: its `parameters` replace values, a function acting on the value in force
function _acted_parameters(m::ManagedCommunity)
    haskey(m.action, :parameters) || return m.parameters
    θ = m.parameters
    for (name, v) in pairs(m.action.parameters)
        if v isa Function
            haskey(θ, name) || throw(ArgumentError(
                "the action changes parameter :$name by a function, and the problem has no value for it"))
            v = v(θ[name])
        end
        θ = merge(θ, NamedTuple{(name,)}((v,)))
    end
    return θ
end

# the columns of the state that a target names: a species, or (species, stage or stages)
function _target(c::OpenCommunity, target)
    sp, stages = target isa Symbol ? (target, nothing) :
                 target isa Tuple{Symbol, Any} ? target :
                 throw(ArgumentError("a target is a species or (species, stages), got $(repr(target))"))
    sp in c.names || throw(ArgumentError("the community has no species :$sp; it has $(c.names)"))
    block = c.layout[sp]
    stages === nothing && return collect(block)
    groups = c.system[sp].groups
    idx = Int[]
    for s in (stages isa Symbol ? (stages,) : stages)
        haskey(groups, s) || throw(ArgumentError("species :$sp has no stage or group :$s"))
        append!(idx, groups[s])
    end
    return collect(block)[unique(idx)]
end

_listed(x::Pair) = (x,)
_listed(x) = x
_interventions(c, action, field) =
    haskey(action, field) ? [(_target(c, target), amount) for (target, amount) in _listed(action[field])] : []

function _check_interventions(removals, takes, adds, demographic)
    for (_, f) in removals
        f isa Real && 0 <= f <= 1 || throw(ArgumentError("`remove` takes a fraction between 0 and 1, got $(repr(f))"))
    end
    for (_, k) in takes
        k isa Real && k >= 0 || throw(ArgumentError("`take` takes a non-negative number, got $(repr(k))"))
        demographic && !isinteger(k) &&
            throw(ArgumentError("with demographic stochasticity `take` needs a whole number, got $k"))
    end
    for (idx, k) in adds
        k isa Real || (k isa AbstractVector{<:Real} && length(k) == length(idx)) || throw(ArgumentError(
            "`add` takes a number, or one number for each of the $(length(idx)) states of its target"))
        k isa Real && length(idx) != 1 && throw(ArgumentError(
            "the target of `add` has $(length(idx)) states: give one number for each"))
        all(>=(0), k) || throw(ArgumentError("`add` takes non-negative numbers, got $(repr(k))"))
        demographic && !all(isinteger, k) &&
            throw(ArgumentError("with demographic stochasticity `add` needs whole numbers, got $(repr(k))"))
    end
    return nothing
end

# Apply the interventions of one census to the state `n`, recording them in `removed` and `added` (views of a row).
function _intervene!(n, removed, added, removals, takes, adds, rng)
    for (idx, f) in removals, i in idx
        x = rng === nothing ? f * n[i] : Float64(StructuredPopulationCore.rand_binomial(rng, round(Int, n[i]), f))
        n[i] -= x
        removed[i] += x
    end
    for (idx, k) in takes
        total = sum(view(n, idx))
        total > 0 || continue
        if k >= total || length(idx) == 1            # all of them, or a number from a single state: no rounding
            for i in idx
                x = k >= total ? n[i] : Float64(k)
                n[i] -= x
                removed[i] += x
            end
        elseif rng === nothing                       # in proportion to the numbers in each state
            for i in idx
                x = k / total * n[i]
                n[i] -= x
                removed[i] += x
            end
        else                                         # whole individuals from several states, without replacement
            left = total
            for _ in 1:round(Int, k)
                u, i = rand(rng) * left, first(idx)
                for j in idx
                    n[j] > 0 || continue
                    i = j
                    u -= n[j]
                    u < 0 && break
                end
                n[i] -= 1
                removed[i] += 1
                left -= 1
            end
        end
    end
    for (idx, k) in adds, (j, i) in enumerate(idx)
        x = k isa Real ? k : k[j]
        n[i] += x
        added[i] += x
    end
    return n
end

"""
    project(m::ManagedCommunity; rng = nothing)

Project a [`ManagedCommunity`](@ref): build the model at its parameters (with the action's changes to
them), then for each environment in turn [`advance`](@ref) the community one step and apply the
action's interventions to the state. Returns a [`ManagedTrajectory`](@ref).

With `m.demographic` the steps have demographic stochasticity and `rng` is required; each step
draws from a generator of its own, seeded from `rng` (see [`ManagedCommunity`](@ref)), so the run
is not the one `project(community, n0, envs; rng)` would give from the same generator. Otherwise
`rng` is not used, the projection is deterministic, and without interventions on the state it is
`project(community, n0, envs)`, step for step.
"""
function project(m::ManagedCommunity; rng = nothing)
    m.demographic && rng === nothing && throw(ArgumentError("a problem with demographic = true needs an rng"))
    model = m.build(_acted_parameters(m))
    model isa NamedTuple && haskey(model, :community) && haskey(model, :n0) && haskey(model, :envs) ||
        throw(ArgumentError("`build` must return a named tuple with `community`, `n0` and `envs`"))
    c = model.community
    c isa OpenCommunity || throw(ArgumentError("`build` must return an `OpenCommunity` as `community`"))
    action = m.action
    removals, takes, adds = (_interventions(c, action, f) for f in (:remove, :take, :add))
    _check_interventions(removals, takes, adds, m.demographic)
    T = length(model.envs)
    steps = haskey(action, :steps) ? action.steps : 1:T
    # `rng` only seeds other generators: one for the interventions, drawn first and always, then one for each
    # step, so that the demographic draws of a step do not depend on how many numbers were drawn before it
    action_rng = m.demographic ? Random.Xoshiro(rand(rng, UInt64)) : nothing
    n = state_vector(c, model.n0)
    u = zeros(Float64, T + 1, length(c))
    pools = fill(NaN, T + 1, length(c.ports))
    vectors = Dict{Symbol, Matrix{Float64}}()
    removed, added = zero(u), zero(u)
    u[1, :] = n
    for (t, env) in enumerate(model.envs)
        haskey(action, :environment) && (env = action.environment(t, env))
        _log_pools!(pools, vectors, t, c.ports, pool_totals(c, n, env))
        n = advance(c, n, env; rng = m.demographic ? Random.Xoshiro(rand(rng, UInt64)) : nothing)
        t in steps && _intervene!(n, view(removed, t + 1, :), view(added, t + 1, :), removals, takes, adds, action_rng)
        u[t + 1, :] = n
    end
    return ManagedTrajectory(CommunityTrajectory(collect(0:T), u, pools, c, vectors), removed, added, m)
end
