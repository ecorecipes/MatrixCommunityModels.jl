"""
    MatrixCommunityModels

Matrix community models (Lytle & Tonkin 2023) built by composition: each species is an
open stage-structured model, each shared density pool (space, biomass, ...) is a junction
of an undirected wiring diagram, and a community is the composite.
"""
module MatrixCommunityModels

using LinearAlgebra
using Catlab: Catlab, UndirectedWiringDiagram, RelationDiagram, nparts, subpart, boxes, ports,
              junction, add_box!, add_junctions!, set_junction!
import StructuredPopulationCore
using StructuredPopulationCore: StateBlockLayout
import MatrixProjectionModels
using MatrixProjectionModels: MatrixProjectionModel

export Context, Transition, Species, OpenCommunity, CommunityTrajectory
export survival, fecundity, vital_matrices, projection_matrix, lower
export open_community, community, compose, pool_totals, step, project, community_matrix
export species, species_names, stage_names, pool_names, diagram

"""
    Context(n, pool, env)

What a vital rate may depend on: the species' own stage vector `n`, the totals of the pools
it is attached to (`pool`, a `NamedTuple`), and the environment of the current step (`env`).
"""
struct Context{N,P,E}
    n::N
    pool::P
    env::E
end

"""
    Transition(name, from, to, kind, rate)

One named entry of a projection matrix: individuals of stage `from` contribute to stage
`to` at `rate`. `kind` is `:survival` (goes into `U`) or `:fecundity` (goes into `F`).
`rate` is a number or a function of a [`Context`](@ref).
"""
struct Transition
    name::Symbol
    from::Symbol
    to::Symbol
    kind::Symbol
    rate::Any
    function Transition(name, from, to, kind, rate)
        kind in (:survival, :fecundity) ||
            throw(ArgumentError("kind must be :survival or :fecundity, got :$kind"))
        return new(name, from, to, kind, rate)
    end
end

"""
    survival(name, from => to, rate)

A survival/growth [`Transition`](@ref) (an entry of `U`).
"""
survival(name::Symbol, ft::Pair{Symbol,Symbol}, rate) = Transition(name, ft.first, ft.second, :survival, rate)

"""
    fecundity(name, from => to, rate)

A reproduction [`Transition`](@ref) (an entry of `F`).
"""
fecundity(name::Symbol, ft::Pair{Symbol,Symbol}, rate) = Transition(name, ft.first, ft.second, :fecundity, rate)

_value(x::Number, ctx) = x
_value(f, ctx) = f(ctx)

"""
    Species(name, stages, transitions; occupancy = [], recruitment = [], threshold = nothing)

A stage-structured species whose vital rates may depend on shared pools and the environment.

- `occupancy`: `pool => weights` pairs saying how much of each pool the species occupies.
  `weights` is a vector (one weight per stage), or a function `(n, env) -> Real`.
- `recruitment`: `stage => f` pairs; `f(ctx)` individuals are added to `stage` each step,
  independently of the species' own abundance.
- `threshold`: an optional map applied to the new stage vector (for example
  quasi-extinction).

One step is `threshold(A(ctx) * n + r(ctx))` with `A = U + F`.
"""
struct Species
    name::Symbol
    stages::Vector{Symbol}
    transitions::Vector{Transition}
    occupancy::Vector{Pair{Symbol,Any}}
    recruitment::Vector{Pair{Symbol,Any}}
    threshold::Any
end

function Species(name::Symbol, stages::AbstractVector{Symbol}, transitions::AbstractVector{Transition};
                 occupancy = Pair{Symbol,Any}[], recruitment = Pair{Symbol,Any}[], threshold = nothing)
    allunique(stages) || throw(ArgumentError("duplicate stage names in species :$name"))
    for t in transitions
        (t.from in stages && t.to in stages) ||
            throw(ArgumentError("transition :$(t.name) of :$name uses an unknown stage"))
    end
    for (pool, w) in occupancy
        w isa AbstractVector && length(w) != length(stages) &&
            throw(DimensionMismatch("occupancy of pool :$pool needs $(length(stages)) weights"))
    end
    for (stage, _) in recruitment
        stage in stages || throw(ArgumentError("recruitment into unknown stage :$stage of :$name"))
    end
    return Species(name, collect(Symbol, stages), collect(Transition, transitions),
                   Pair{Symbol,Any}[p.first => p.second for p in occupancy],
                   Pair{Symbol,Any}[p.first => p.second for p in recruitment], threshold)
end

_stage_index(sp::Species, s::Symbol) = findfirst(==(s), sp.stages)

"""
    vital_matrices(species, ctx) -> (U, F)

The survival and fecundity matrices of `species` evaluated at a [`Context`](@ref).
Transitions sharing a matrix entry are added in the order they were declared.
"""
function vital_matrices(sp::Species, ctx::Context)
    k = length(sp.stages)
    U = zeros(Float64, k, k)
    F = zeros(Float64, k, k)
    for t in sp.transitions
        i, j = _stage_index(sp, t.to), _stage_index(sp, t.from)
        M = t.kind === :survival ? U : F
        M[i, j] += _value(t.rate, ctx)
    end
    return U, F
end

"""
    projection_matrix(species, ctx)

`U + F` for `species` at a [`Context`](@ref).
"""
function projection_matrix(sp::Species, ctx::Context)
    U, F = vital_matrices(sp, ctx)
    return U + F
end

"""
    lower(species, ctx)

Freeze `species` at a [`Context`](@ref) and return an ordinary
`MatrixProjectionModels.MatrixProjectionModel` with the `U`/`F` split preserved.
"""
function lower(sp::Species, ctx::Context)
    U, F = vital_matrices(sp, ctx)
    return MatrixProjectionModel(U, F; stage_names = sp.stages)
end

function _species_step(sp::Species, n, pool, env)
    ctx = Context(n, pool, env)
    n1 = projection_matrix(sp, ctx) * n
    for (stage, f) in sp.recruitment
        n1[_stage_index(sp, stage)] += _value(f, ctx)
    end
    return sp.threshold === nothing ? n1 : sp.threshold(n1)
end

_occupancy(w::AbstractVector, n, env) = dot(w, n)
_occupancy(f, n, env) = f(n, env)

"""
    OpenCommunity

A community with exposed pools. It has a state layout (one block of stages per species),
`ports` (pool names), a readout giving its contribution to each pool, and an update taking
the pool totals. Build one from a [`Species`](@ref) with [`open_community`](@ref), or from
several with [`community`](@ref) or [`compose`](@ref).
"""
struct OpenCommunity
    layout::StateBlockLayout
    stages::Dict{Symbol,Vector{Symbol}}
    ports::Vector{Symbol}
    readout::Any
    update::Any
    parts::Vector{Any}
    diagram::Any
end

"""
    species_names(c)

Names of the species (state blocks) of an [`OpenCommunity`](@ref), in state order.
"""
species_names(c::OpenCommunity) = collect(keys(c.layout))

"""
    stage_names(c, species)

Stage names of `species` in an [`OpenCommunity`](@ref).
"""
stage_names(c::OpenCommunity, sp::Symbol) = c.stages[sp]

"""
    pool_names(c)

Names of the pools an [`OpenCommunity`](@ref) exposes.
"""
pool_names(c::OpenCommunity) = c.ports

"""
    diagram(c)

The undirected wiring diagram a composite [`OpenCommunity`](@ref) was built from, or
`nothing` for a single species.
"""
diagram(c::OpenCommunity) = c.diagram

Base.length(c::OpenCommunity) = length(c.layout)

function Base.show(io::IO, c::OpenCommunity)
    print(io, "OpenCommunity(", join(species_names(c), ", "), "; pools = ", join(c.ports, ", "), ")")
end

"""
    open_community(species)

View a single [`Species`](@ref) as an [`OpenCommunity`](@ref) with one port per pool it
occupies.
"""
function open_community(sp::Species)
    pools = Symbol[p.first for p in sp.occupancy]
    allunique(pools) || throw(ArgumentError("species :$(sp.name) lists a pool twice"))
    occ = Any[p.second for p in sp.occupancy]
    readout = (n, env) -> Float64[_occupancy(w, n, env) for w in occ]
    names = Tuple(pools)
    update = (n, totals, env) -> _species_step(sp, n, NamedTuple{names}(Tuple(totals)), env)
    layout = StateBlockLayout([sp.name], [length(sp.stages)])
    return OpenCommunity(layout, Dict(sp.name => sp.stages), pools, readout, update, Any[sp], nothing)
end

open_community(c::OpenCommunity) = c

"""
    compose(diagram, components)

Compose open communities (or species) along an undirected wiring diagram. Box `b` is
filled by `components[b]`, whose ports are matched to the box's ports in order. The total
of each junction is the sum, in port order, of the readouts attached to it; every component
is then updated with the totals of its junctions, all from the same time-`t` state. Outer
ports of the diagram become the pools of the composite.
"""
function compose(d::UndirectedWiringDiagram, components::AbstractVector)
    comps = OpenCommunity[open_community(c) for c in components]
    nb = nparts(d, :Box)
    nb == length(comps) || throw(DimensionMismatch("diagram has $nb boxes, got $(length(comps)) components"))
    boxports = [collect(ports(d, b)) for b in 1:nb]
    for b in 1:nb
        length(boxports[b]) == length(comps[b].ports) || throw(DimensionMismatch(
            "box $b has $(length(boxports[b])) ports but its component exposes $(length(comps[b].ports))"))
    end
    names = Symbol[]
    sizes = Int[]
    stages = Dict{Symbol,Vector{Symbol}}()
    for c in comps, s in species_names(c)
        s in names && throw(ArgumentError("species :$s appears twice in the composite"))
        push!(names, s)
        push!(sizes, length(c.stages[s]))
        stages[s] = c.stages[s]
    end
    layout = StateBlockLayout(names, sizes)
    offsets = cumsum([0; [length(c) for c in comps]])
    ranges = [(offsets[b] + 1):offsets[b + 1] for b in 1:nb]
    nj = nparts(d, :Junction)
    portjunction = [[junction(d, p) for p in boxports[b]] for b in 1:nb]
    outer = [junction(d, q; outer = true) for q in 1:nparts(d, :OuterPort)]
    jnames = _junction_names(d, comps, portjunction, nj)
    outernames = Symbol[jnames[j] for j in outer]

    function internal(n, env)
        T = zeros(Float64, nj)
        for b in 1:nb
            y = comps[b].readout(view(n, ranges[b]), env)
            for (i, j) in enumerate(portjunction[b])
                T[j] += y[i]
            end
        end
        return T
    end
    readout = (n, env) -> internal(n, env)[outer]
    function update(n, totals, env)
        T = internal(n, env)
        for (q, j) in enumerate(outer)
            T[j] = totals[q]
        end
        n1 = similar(n, Float64)
        for b in 1:nb
            n1[ranges[b]] = comps[b].update(n[ranges[b]], T[portjunction[b]], env)
        end
        return n1
    end
    return OpenCommunity(layout, stages, outernames, readout, update, Any[comps...], d)
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
    community(species...)

Compose species that share pools by name: every pool named in any species' `occupancy`
becomes one junction, exposed as an outer port. Returns the composite
[`OpenCommunity`](@ref); the generated wiring diagram is available from [`diagram`](@ref).
"""
function community(species...)
    comps = OpenCommunity[open_community(s) for s in species]
    pools = Symbol[]
    for c in comps, p in c.ports
        p in pools || push!(pools, p)
    end
    d = RelationDiagram(length(pools))
    js = add_junctions!(d, length(pools); variable = pools)
    for c in comps
        b = add_box!(d, length(c.ports); name = first(species_names(c)))
        set_junction!(d, ports(d, b), [findfirst(==(p), pools) for p in c.ports])
    end
    set_junction!(d, js; outer = true)
    return compose(d, comps)
end

"""
    pool_totals(c, n, env)

Totals of the exposed pools of a closed community `c` in state `n`, as a `NamedTuple`.
"""
function pool_totals(c::OpenCommunity, n::AbstractVector, env)
    return NamedTuple{Tuple(c.ports)}(Tuple(c.readout(n, env)))
end

"""
    step(c, n, env)

Advance a closed community one time step: pool totals are the community's own readouts.
"""
step(c::OpenCommunity, n::AbstractVector, env) = c.update(n, c.readout(n, env), env)

"""
    CommunityTrajectory

Result of [`project`](@ref). `u[t + 1, :]` is the state at time `t` and `pools[t + 1, :]`
the pool totals used for the step from `t`. Index with a species name to get its
stage-by-time matrix.
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
    project(c, n0, envs)

Project a closed community from state `n0` through the environments `envs` (one per time
step). `n0` is a flat vector in layout order, or `species => vector` pairs.
"""
function project(c::OpenCommunity, n0::AbstractVector{<:Real}, envs)
    length(n0) == length(c) || throw(DimensionMismatch("n0 has length $(length(n0)), community has $(length(c)) states"))
    T = length(envs)
    u = zeros(Float64, T + 1, length(c))
    pools = fill(NaN, T + 1, length(c.ports))
    n = collect(Float64, n0)
    u[1, :] = n
    for (k, env) in enumerate(envs)
        totals = c.readout(n, env)
        pools[k, :] = totals
        n = c.update(n, totals, env)
        u[k + 1, :] = n
    end
    return CommunityTrajectory(collect(0:T), u, pools, c)
end

function project(c::OpenCommunity, n0::AbstractVector{<:Pair}, envs)
    d = Dict(n0)
    return project(c, reduce(vcat, [collect(Float64, d[s]) for s in species_names(c)]), envs)
end

_species(c::OpenCommunity) = reduce(vcat, [p isa Species ? Species[p] : _species(p) for p in c.parts])

"""
    species(c)
    species(c, name)

The [`Species`](@ref) definitions inside an [`OpenCommunity`](@ref), in state order, or
the one called `name`.
"""
species(c::OpenCommunity) = _species(c)
species(c::OpenCommunity, name::Symbol) = only(s for s in _species(c) if s.name === name)

"""
    community_matrix(c, n, env)

The block-diagonal projection matrix of a closed community frozen at state `n` and
environment `env`, one block per species in layout order. Pool totals are those of `n`.
Only communities whose species all share the top-level pools are supported.
"""
function community_matrix(c::OpenCommunity, n::AbstractVector, env)
    pool = pool_totals(c, n, env)
    A = zeros(Float64, length(c), length(c))
    for sp in _species(c)
        r = c.layout[sp.name]
        names = Tuple(p.first for p in sp.occupancy)
        ctx = Context(n[r], NamedTuple{names}(Tuple(pool[q] for q in names)), env)
        A[r, r] = projection_matrix(sp, ctx)
    end
    return A
end

include("models.jl")
export verde_fish_model, riparian_model, VERDE_SPECIES, RIPARIAN_GUILDS

end # module MatrixCommunityModels
