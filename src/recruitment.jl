# Recruitment into a limited resource: a lottery for free space.
#
# With synchronous updates, recruitment limited only by the free space at the start of a
# step is not limited within the step: many recruits can arrive in a pulse and overshoot
# the space. The lottery fixes this. The recruits trying to establish in a step (the
# demand) share the free space: each establishes with a probability that depends on the
# free space and the demand, so that no more establish than there is room for.

# A response added by `lottery`, marked so that the demand on a transition can be computed
# from its other responses.
struct _LotteryResponse{F} <: Function
    f::F
end
(r::_LotteryResponse)(ctx) = r.f(ctx)

# The pools seen when computing a demand as a pool contribution, which cannot read pools.
struct _NoPools end
_no_pool(p) = throw(ArgumentError(
    "a response on a lottery transition reads the pool :$p, so its demand cannot be a pool contribution " *
    "(contributions see only the stages and the environment); use `demand = :own` in `lottery` instead"))
Base.getproperty(::_NoPools, p::Symbol) = _no_pool(p)
Base.getindex(::_NoPools, p::Symbol) = _no_pool(p)
Base.haskey(::_NoPools, ::Symbol) = false

# The compiled entries of the named transitions of a responsive net.
function _lottery_entries(net::ResponsiveProjectionNet, transitions)
    isempty(transitions) && throw(ArgumentError("name at least one transition"))
    for t in transitions
        t in transition_names(net) ||
            throw(ArgumentError("unknown transition :$t; transitions are $(transition_names(net))"))
    end
    return [entry for entry in net.compiled if entry[1] in transitions]
end

# Individuals flowing through `entries` in context `ctx`: every value with the responses of
# its transition evaluated, except lotteries, times the abundance of its source stage.
function _flow(net::ResponsiveProjectionNet, entries, ctx)
    k = length(stage_names(net))
    total = 0.0
    for entry in entries
        t, rows, cols, vals, froms, tos = entry
        v = copy(vals)
        for f in get(net.by_transition, t, ())
            f isa _LotteryResponse && continue
            CategoricalPopulationDynamics._scale!(v, f(ctx), t, k, rows, cols, froms, tos)
        end
        for e in eachindex(v)
            total += v[e] * ctx.n[cols[e]]
        end
    end
    return total
end

# The free space: `capacity` less what is held, read from a pool or from the species' own stages.
function _free_space(space::Symbol, K::Float64)
    return ctx -> K - Float64(getproperty(ctx.pool, space))
end
function _free_space(space::AbstractVector, K::Float64)
    weights = Pair{Symbol, Float64}[s isa Pair ? (Symbol(s.first) => Float64(s.second)) : (Symbol(s) => 1.0)
                                   for s in space]
    return function (ctx)
        free = K
        for (s, w) in weights
            free -= w * ctx.n[s]
        end
        return free
    end
end
_free_space(space, K) = throw(ArgumentError(
    "space must be a pool name, or the stages that hold space (names or `stage => weight` pairs); got $space"))

_lottery_rule(rule::Function) = rule
function _lottery_rule(rule::Symbol)
    rule === :lottery && return (free, demand) -> min(1.0, free / demand)
    rule === :saturating && return (free, demand) -> free / (free + demand)
    throw(ArgumentError("rule must be :lottery, :saturating or a function (free, demand) -> probability; got :$rule"))
end

# Establishment probability from the free space and the demand. With no demand, a single
# recruit establishes if there is any free space (the limit of both rules as demand → 0).
function _establishment(rule, free::Float64, demand::Float64)
    free = max(free, 0.0)
    demand > 0 || return free > 0 ? 1.0 : 0.0
    return Float64(rule(free, demand))
end

"""
    lottery(net, transitions...; capacity, space = :space, demand = :demand, rule = :lottery)

Make the named recruitment `transitions` of `net` compete for a limited resource, such as
space, within each time step. The individuals passing through these transitions in a step
(the demand) share the free space, `capacity` less the space held: each establishes with
probability `min(1, free / demand)` (`rule = :lottery`: every recruit has the same chance,
and no more establish than there is room for) or `free / (free + demand)` (`rule =
:saturating`: recruits and free sites meet at random, and the number established, `demand ×
free / (free + demand)`, rises smoothly towards the free space). `rule` can also be a
function `(free, demand) -> probability`. With no demand the probability is one if there is
free space and zero otherwise; a negative free space counts as none. Returns a responsive
net with this response attached to each of `transitions`.

This is the device for pulse recruitment into a finite resource. With synchronous updates,
recruitment limited only by the free space at the start of a step is not limited within
the step, and a pulse of recruits can overshoot the space; the lottery cannot.

- `space`: where the space held is read. A pool name (the space held by every species
  attached to that pool), or the species' own stages that hold space, as stage names or
  `stage => weight` pairs (the species competes only with itself).
- `demand`: the recruits competing in a step. A pool name: the recruits of every species
  attached to that pool, each of which contributes its own with
  [`lottery_demand`](@ref) in [`species`](@ref). Or `:own`: the species' own recruits,
  computed from the transitions' values in the step.

In both cases the demand counts, for each entry of `transitions`, its value with the
transition's other responses (those attached before `lottery`) evaluated, times the
abundance of its source stage. Attach `lottery` after `expand_stage` or `stratify`, since
it refers to the net's stages as they are.

Two annual plants sharing space by lottery, as in the seasons vignette (germinating seeds
leave the seed bank and establish if they win a site):

```julia
plant = respond(life, :germination => ctx -> ctx.env == germinate ? 1.0 : 0.0)
plant = lottery(plant, :germination; capacity = 1000.0, space = :space, demand = :demand)
species(name, plant; occupies = [:space => [:plant => 1.0], :demand => lottery_demand(plant, :germination)])
```

A weed whose seedlings compete for the sites its own rosettes leave open:

```julia
lottery(net, :emerging, :establishing; capacity = sites, space = [:small, :medium, :large],
        demand = :own, rule = :saturating)
```
"""
function lottery(net, transitions::Symbol...; capacity::Real, space = :space, demand = :demand,
                 rule = :lottery)
    r = ResponsiveProjectionNet(net)
    entries = _lottery_entries(r, transitions)
    free = _free_space(space, Float64(capacity))
    chance = _lottery_rule(rule)
    demand isa Symbol || throw(ArgumentError("demand must be a pool name or :own; got $demand"))
    f = if demand === :own
        _LotteryResponse(ctx -> _establishment(chance, free(ctx), _flow(r, entries, ctx)))
    else
        _LotteryResponse(ctx -> _establishment(chance, free(ctx), Float64(getproperty(ctx.pool, demand))))
    end
    return respond(r, (t => f for t in transitions)...)
end

"""
    lottery_demand(net, transitions...)

A species' contribution to the demand pool of a [`lottery`](@ref): a function `(n, env) ->
number` giving the individuals that pass through the named `transitions` of `net` in a
step, for use in `occupies` of [`species`](@ref). Each entry of the transitions counts its
value, with the transition's responses evaluated (except lotteries), times the abundance of
its source stage. A contribution sees only the stages and the environment, so a response on
these transitions that reads a pool is an error; use `demand = :own` in `lottery` for a
species that competes only with itself.
"""
function lottery_demand(net, transitions::Symbol...)
    r = ResponsiveProjectionNet(net)
    entries = _lottery_entries(r, transitions)
    return function (n, env)
        stages = n isa StageVector ? n : StageVector(n, r.groups)
        return _flow(r, entries, ResponseContext(stages, _NoPools(), env))
    end
end
