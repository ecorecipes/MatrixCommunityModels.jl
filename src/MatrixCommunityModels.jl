"""
    MatrixCommunityModels

Matrix community models (Lytle & Tonkin 2023) built from the ecorecipes categorical stack.

A species is a projection net from CategoricalPopulationDynamics.jl: named stages, named
transitions, and responses saying how a transition's values depend on shared pools and
the environment. A community is a composite of species along an undirected wiring diagram
whose junctions are the pools (space, biomass, parasites per host, ...). The composite is
a `PooledProjectionSystem`, which is lowered to a `CoupledMPMProblem` and solved by
MatrixProjectionModels.jl.

This package adds the composition along wiring diagrams, the analysis of the composed
community, and a set of published models.
"""
module MatrixCommunityModels

using LinearAlgebra
using Catlab: Catlab, UndirectedWiringDiagram, RelationDiagram, nparts, subpart, ports,
              junction, add_box!, add_junctions!, set_junction!
import StructuredPopulationCore
using StructuredPopulationCore: StateBlockLayout
import MatrixProjectionModels
using MatrixProjectionModels: MatrixProjectionModel, DirectIteration
import CategoricalPopulationDynamics
using CategoricalPopulationDynamics: ValuedProjectionNet, ResponsiveProjectionNet, ResponseContext,
    StageVector, Pool, PooledProjectionSystem, StateDependentMPMTarget,
    respond, inflow, after_step, births, context, expand_stage, expand_stages, stratify,
    stratum_stage, stratum_transition, transition_names, transition_matrix, to_matrix,
    evaluate, survival_fecundity_matrices, ⊕
import CategoricalPopulationDynamics: stage_names, pool_names, pool_totals, advance

# the modelling vocabulary of CategoricalPopulationDynamics.jl, re-exported
export ValuedProjectionNet, ResponsiveProjectionNet, ResponseContext, StageVector, Pool,
       PooledProjectionSystem, respond, inflow, after_step, births, context, expand_stage,
       expand_stages, stratify, stratum_stage, stratum_transition, stage_names,
       transition_names, transition_matrix, to_matrix, evaluate, survival_fecundity_matrices, ⊕

export OpenCommunity, CommunityTrajectory
export species, community, compose, system, species_names, pool_names, diagram
export subcommunity, matrix_net
export project, advance, pool_totals, community_matrix, frozen_net, matrix_model, state_vector

# A view of the pool totals under a component's own (local) pool names.
struct PoolView{P <: NamedTuple}
    pool::P
    alias::Dict{Symbol, Symbol}
end
Base.getproperty(v::PoolView, name::Symbol) =
    getfield(getfield(v, :pool), get(getfield(v, :alias), name, name))
Base.getindex(v::PoolView, name::Symbol) = getproperty(v, name)
Base.haskey(v::PoolView, name::Symbol) = haskey(getfield(v, :pool), get(getfield(v, :alias), name, name))

"""
    OpenCommunity

One or more species with exposed pools (its ports). Build one species with
[`species`](@ref), and join them with [`community`](@ref) or [`compose`](@ref).

The flattened model is a `CategoricalPopulationDynamics.PooledProjectionSystem`, available
from [`system`](@ref). The state of a community is a flat vector: one block of stages per
species, in the order of [`species_names`](@ref).
"""
struct OpenCommunity
    names::Vector{Symbol}
    nets::Dict{Symbol, ResponsiveProjectionNet}
    aliases::Dict{Symbol, Dict{Symbol, Symbol}}
    contributions::Vector{Tuple{Symbol, Symbol, Any}}
    ports::Vector{Symbol}
    diagram::Any
    system::PooledProjectionSystem
    layout::StateBlockLayout
end

function _open(names, nets, aliases, contributions, ports, diagram)
    components = map(names) do name
        net, alias = nets[name], aliases[name]
        if all(k == v for (k, v) in alias)
            name => net
        else
            # responses were written with the species' own pool names: translate them
            local_view(f) = ctx -> f(ResponseContext(ctx.n, PoolView(ctx.pool, alias), ctx.env))
            name => ResponsiveProjectionNet(net.base;
                responses = Pair{Symbol, Any}[t => local_view(f) for (t, f) in net.responses],
                inflows = Pair{Symbol, Any}[s => local_view(f) for (s, f) in net.inflows],
                after = net.after, fecundity = net.fecundity, groups = net.groups)
        end
    end
    pool_order = copy(ports)
    for (pool, _, _) in contributions
        pool in pool_order || push!(pool_order, pool)
    end
    pools = [Pool(p, Pair{Symbol, Any}[c => w for (q, c, w) in contributions if q === p]) for p in pool_order]
    sys = PooledProjectionSystem(components...; pools = pools)
    layout = StateBlockLayout(collect(Symbol, names), [length(stage_names(nets[n])) for n in names])
    return OpenCommunity(collect(Symbol, names), nets, aliases, contributions, collect(Symbol, ports), diagram, sys, layout)
end

"""
    species(name, net; occupies = [], reads = [])

A single species as an [`OpenCommunity`](@ref). `net` is a valued or responsive projection
net. `occupies` lists `pool => contribution` pairs saying how much of each pool the species
occupies; a contribution is a vector of weights (one per stage), a list of `stage =>
weight` pairs, a number (weight on every stage), or a function `(n, env) -> number`.
`reads` names pools that the species' responses read without contributing to them.

```julia
grass = species(:grass, respond(life, :seeding => ctx -> max(0, 1 - ctx.pool.space / 1000));
                occupies = [:space => [:seedling => 0.1, :adult => 1.0]])
```
"""
function species(name::Symbol, net; occupies = Pair{Symbol, Any}[], reads = Symbol[])
    r = ResponsiveProjectionNet(net)
    pools = Symbol[p.first for p in occupies]
    allunique(pools) || throw(ArgumentError("species :$name lists a pool twice"))
    ports = vcat(pools, [p for p in reads if !(p in pools)])
    contributions = Tuple{Symbol, Symbol, Any}[(p.first, name, p.second) for p in occupies]
    return _open([name], Dict(name => r), Dict(name => Dict(p => p for p in ports)), contributions, ports, nothing)
end

"""
    subcommunity(c, names...)

The named species of community `c` as a community of their own, keeping the pools they
are attached to: for example the natives of a fish community without the non-natives, or
one species on its own as an invader.
"""
function subcommunity(c::OpenCommunity, names::Symbol...)
    all(n -> n in c.names, names) || throw(ArgumentError("not all of $(names) are species of this community"))
    keep = [n for n in c.names if n in names]
    aliases = Dict(n => c.aliases[n] for n in keep)
    contributions = Tuple{Symbol, Symbol, Any}[x for x in c.contributions if x[2] in keep]
    used = unique(vcat([collect(values(aliases[n])) for n in keep]..., [x[1] for x in contributions]))
    ports = [p for p in pool_names(c.system) if p in used]
    return _open(keep, Dict(n => c.nets[n] for n in keep), aliases, contributions, ports, nothing)
end

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
    system(c)

The `PooledProjectionSystem` of a community: its species as responsive projection nets,
and its pools.
"""
system(c::OpenCommunity) = c.system

"""
    species_names(c)

Names of the species of a community, in state order.
"""
species_names(c::OpenCommunity) = c.names

"""
    stage_names(c, species)

Stage names of `species` in a community.
"""
stage_names(c::OpenCommunity, sp::Symbol) = stage_names(c.nets[sp])

"""
    pool_names(c; internal = false)

Names of the pools a community exposes; with `internal = true`, all of its pools,
including those shared only inside it.
"""
pool_names(c::OpenCommunity; internal::Bool = false) = internal ? pool_names(c.system) : c.ports

"""
    diagram(c)

The undirected wiring diagram a community was composed along, or `nothing` for a single
species.
"""
diagram(c::OpenCommunity) = c.diagram

Base.length(c::OpenCommunity) = length(c.layout)
Base.getindex(c::OpenCommunity, sp::Symbol) = c.system[sp]

function Base.show(io::IO, c::OpenCommunity)
    print(io, "OpenCommunity(", join(c.names, ", "), "; pools = ", join(c.ports, ", "), ")")
end

const _JUNCTION_COUNTER = Ref(0)

"""
    compose(diagram, components)

Compose open communities along an undirected wiring diagram. Box `b` is filled by
`components[b]`, whose pools are matched to the box's ports in order. Each junction is one
pool: its total is the sum of what every attached species contributes. Outer ports of the
diagram are the pools of the composite; junctions that are not exposed are shared only
inside it.
"""
function compose(d::UndirectedWiringDiagram, components::AbstractVector)
    nb = nparts(d, :Box)
    nb == length(components) ||
        throw(DimensionMismatch("diagram has $nb boxes, got $(length(components)) components"))
    comps = OpenCommunity[c for c in components]
    boxports = [collect(ports(d, b)) for b in 1:nb]
    for b in 1:nb
        length(boxports[b]) == length(comps[b].ports) || throw(DimensionMismatch(
            "box $b has $(length(boxports[b])) ports but its component exposes $(length(comps[b].ports))"))
    end
    nj = nparts(d, :Junction)
    portjunction = [[junction(d, p) for p in boxports[b]] for b in 1:nb]
    outer = [junction(d, q; outer = true) for q in 1:nparts(d, :OuterPort)]
    jnames = _junction_names(d, comps, portjunction, nj)
    for j in 1:nj
        # a junction that is not exposed gets a name nothing outside can clash with
        j in outer || (jnames[j] = Symbol(jnames[j], "#", _JUNCTION_COUNTER[] += 1))
    end

    names = Symbol[]
    nets = Dict{Symbol, ResponsiveProjectionNet}()
    aliases = Dict{Symbol, Dict{Symbol, Symbol}}()
    contributions = Tuple{Symbol, Symbol, Any}[]
    for b in 1:nb
        c = comps[b]
        rename = Dict(c.ports[i] => jnames[j] for (i, j) in enumerate(portjunction[b]))
        for name in c.names
            name in names && throw(ArgumentError("species :$name appears twice in the composite"))
            push!(names, name)
            nets[name] = c.nets[name]
            aliases[name] = Dict(k => get(rename, v, v) for (k, v) in c.aliases[name])
        end
        for (pool, name, w) in c.contributions
            push!(contributions, (get(rename, pool, pool), name, w))
        end
    end
    return _open(names, nets, aliases, contributions, Symbol[jnames[j] for j in outer], d)
end

function _junction_names(d, comps, portjunction, nj)
    vars = try
        subpart(d, :variable)
    catch
        nothing
    end
    vars !== nothing && eltype(vars) == Symbol && return collect(Symbol, vars)
    names = Vector{Symbol}(undef, nj)
    for (b, js) in enumerate(portjunction), (i, j) in enumerate(js)
        isassigned(names, j) || (names[j] = comps[b].ports[i])
    end
    return names
end

"""
    community(components...)

Compose open communities that share pools by name: every pool named by any component
becomes one junction, exposed as a pool of the result. The wiring diagram that was
generated is available from [`diagram`](@ref).
"""
function community(components::OpenCommunity...)
    pools = Symbol[]
    for c in components, p in c.ports
        p in pools || push!(pools, p)
    end
    d = RelationDiagram(length(pools))
    js = add_junctions!(d, length(pools); variable = pools)
    for c in components
        b = add_box!(d, length(c.ports); name = first(c.names))
        set_junction!(d, ports(d, b), [findfirst(==(p), pools) for p in c.ports])
    end
    set_junction!(d, js; outer = true)
    return compose(d, collect(components))
end

# ---------------------------------------------------------------- state

_split(c::OpenCommunity, n::AbstractVector) =
    Dict{Symbol, Vector{Float64}}(sp => collect(Float64, n[c.layout[sp]]) for sp in c.names)
_join(c::OpenCommunity, state::AbstractDict) = reduce(vcat, [state[sp] for sp in c.names])

"""
    state_vector(c, species => stages, ...)

The flat state vector of a community from one stage vector per species.
"""
function state_vector(c::OpenCommunity, blocks::Pair{Symbol}...)
    d = Dict(blocks)
    Set(keys(d)) == Set(c.names) || throw(ArgumentError("need one stage vector for each of $(c.names)"))
    n = reduce(vcat, [collect(Float64, d[sp]) for sp in c.names])
    length(n) == length(c) || throw(DimensionMismatch("stage vectors do not match the community's stages"))
    return n
end

_state(c::OpenCommunity, n::AbstractVector{<:Real}) =
    length(n) == length(c) ? collect(Float64, n) :
    throw(DimensionMismatch("state has length $(length(n)), community has $(length(c)) states"))
_state(c::OpenCommunity, n::AbstractVector{<:Pair}) = state_vector(c, n...)

"""
    pool_totals(c, n, env)

Totals of the pools a closed community exposes, in state `n`, as a `NamedTuple`.
"""
function pool_totals(c::OpenCommunity, n::AbstractVector, env; internal::Bool = false)
    all_pools = pool_totals(c.system, _split(c, _state(c, n)), env)
    internal && return all_pools
    return NamedTuple{Tuple(c.ports)}(Tuple(all_pools[p] for p in c.ports))
end

"""
    advance(c, n, env; rng = nothing)

Advance a closed community one time step. With `rng`, the step has demographic
stochasticity (see [`project`](@ref)).
"""
function advance(c::OpenCommunity, n::AbstractVector, env; rng = nothing)
    state = _split(c, _state(c, n))
    rng === nothing && return _join(c, advance(c.system, state, env))
    return _join(c, _demographic_advance(c.system, state, env, rng))
end

# One step with demographic stochasticity: individuals are whole, survival and movement are
# a multinomial draw down each column of U, births and inflows are Poisson.
function _demographic_advance(sys::PooledProjectionSystem, state, env, rng)
    pool = pool_totals(sys, state, env)
    next = Dict{Symbol, Vector{Float64}}()
    for name in CategoricalPopulationDynamics.component_names(sys)
        net, n = sys[name], state[name]
        all(isinteger, n) || throw(ArgumentError(
            "demographic stochasticity needs whole numbers of individuals in :$name"))
        ctx = context(net, n, pool, env)
        U, F = survival_fecundity_matrices(net, ctx)
        worst = maximum(sum(U; dims = 1))
        worst <= 1 + 1e-6 || throw(ArgumentError(
            "survival matrix of :$name has a column summing to $worst; demographic stochasticity needs " *
            "survival probabilities (column sums at most 1). Mark births with `births(net, ...)`"))
        counts = StructuredPopulationCore.demographic_step!(rng, zeros(Int, length(n)), round.(Int, n), U, F)
        stepped = Float64.(counts)
        for (i, r) in enumerate(CategoricalPopulationDynamics._inflows(net, ctx))
            r > 0 && (stepped[i] += StructuredPopulationCore.rand_poisson(rng, r))
        end
        next[name] = collect(Float64, CategoricalPopulationDynamics._after(net, stepped))
    end
    return next
end

"""
    CommunityTrajectory

Result of [`project`](@ref). `u[t + 1, :]` is the state at time `t`; `pools[t + 1, :]`
holds the totals of the exposed pools used for the step from `t`. Index with a species
name to get its stage-by-time matrix.
"""
struct CommunityTrajectory
    t::Vector{Int}
    u::Matrix{Float64}
    pools::Matrix{Float64}
    community::OpenCommunity
end

Base.getindex(sol::CommunityTrajectory, sp::Symbol) = sol.u[:, sol.community.layout[sp]]

function Base.show(io::IO, sol::CommunityTrajectory)
    print(io, "CommunityTrajectory(", length(sol.t) - 1, " steps, ", size(sol.u, 2), " states)")
end

"""
    project(c, n0, envs; rng = nothing)

Project a closed community from state `n0` through the environments `envs` (one per time
step). `n0` is a flat vector, or `species => stages` pairs.

The community is lowered to a `MatrixProjectionModels.CoupledMPMProblem` and solved by
MatrixProjectionModels.jl.

With `rng` (a random number generator) the projection has demographic stochasticity:
individuals are whole, each survives and moves according to a multinomial draw from its
column of the survival matrix, and births and inflows are Poisson with the deterministic
values as means. This needs whole numbers in `n0`, the births of each species marked with
`births`, and survival matrices whose columns sum to at most one.
"""
function project(c::OpenCommunity, n0, envs; rng = nothing)
    start = _state(c, n0)
    T = length(envs)
    u = zeros(Float64, T + 1, length(c))
    pools = fill(NaN, T + 1, length(c.ports))
    u[1, :] = start
    if rng === nothing
        T == 0 && return CommunityTrajectory([0], u, pools, c)
        target = StateDependentMPMTarget(_split(c, start), (0, T); p = collect(envs))
        prob = CategoricalPopulationDynamics.lower(c.system, target)
        sol = MatrixProjectionModels.solve(prob, DirectIteration())
        for k in 1:(T + 1), sp in c.names
            u[k, c.layout[sp]] = sol.u[k][sp]
        end
        logged = sol.substep_log[:pools]
        for k in 1:T, (j, p) in enumerate(c.ports)
            pools[k, j] = logged[k][p]
        end
    else
        n = start
        for (k, env) in enumerate(envs)
            totals = pool_totals(c, n, env)
            pools[k, :] = collect(totals)
            n = advance(c, n, env; rng)
            u[k + 1, :] = n
        end
    end
    return CommunityTrajectory(collect(0:T), u, pools, c)
end

_context(c::OpenCommunity, sp::Symbol, n::AbstractVector, env) =
    context(c.system[sp], _state(c, n)[c.layout[sp]], pool_totals(c.system, _split(c, _state(c, n)), env), env)

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

include("analysis.jl")
export equilibrium, jacobian, invasion_growth_rate, mean_abundance, abundance_sensitivity
export equilibrium_sensitivity, stochastic_invasion_rate
export InteractionNetwork, links, connectance, link_density, reciprocity, out_degree, in_degree
export periodic_matrix, abundance, stage_total, CommunityEnsemble, ensemble, extinction_times, quasi_extinction

include("models.jl")
export verde_fish_model, riparian_model, cypripedium_model, with_rate, VERDE_SPECIES, RIPARIAN_GUILDS
export density_response

include("bees.jl")
export day_of_year, hopomo_season, forage_index, flight_index
export honeybee, honeybee_state, varroa, bumblebee, floral_supply, floral_competition, flower_transmission
export HONEYBEE_PARAMETERS, VARROA_PARAMETERS, BUMBLEBEE_PARAMETERS

end # module MatrixCommunityModels
