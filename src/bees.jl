# Honey bees, Varroa mites, deformed wing virus (DWV) and bumblebees.
#
# Each species is a short life-cycle net. Stage durations come from `expand_stages`,
# infection status from `stratify`, and every dependence on the rest of the community is a
# named response to a pool. The time step is one day; the environment of a step is the day
# number (an integer, day 1 = 1 January of the first year).

"""
    day_of_year(env)

Day of the year (1 to 365) for an environment given as a day number.
"""
day_of_year(env::Integer) = mod1(env, 365)

"""
    hopomo_season(day)

The seasonal factor of the HoPoMo honey bee model (Schmickl & Crailsheim 2007, eq. 1 with
their parameters 385, 30, 36, 155, 30): close to 0 in midwinter and close to 1 in early
June. The queen's laying rate and the availability of forage follow it.
"""
function hopomo_season(day::Real)
    x1, x2, x3, x4, x5 = 385.0, 30.0, 36.0, 155.0, 30.0
    return 1 - max(1 - 1 / (1 + x1 * exp(-2 * day / x2)), 1 / (1 + x3 * exp(-2 * (day - x4) / x5)))
end

"""
    forage_index(day)

Availability of forage on a day of the year, between 0 and 1 (`min(1, 1.5 ×
hopomo_season)`, Schmickl & Crailsheim 2007).
"""
forage_index(day::Real) = min(1.0, 1.5 * hopomo_season(day))

_temperature(day) = 14 + 11 * sin(2pi * (day - 105) / 365)

"""
    flight_index(day; threshold = 10, full = 18)

Fraction of a day of the year on which bees can fly, between 0 and 1: flight rises
linearly as the daily maximum temperature goes from `threshold` to `full` °C. The
temperature is an illustrative temperate climate (14 ± 11 °C, warmest in late July), not
fitted to a site. Bumblebees fly in cooler weather (MacQueen et al. 2025); use
`threshold = 6, full = 12` for them.
"""
flight_index(day::Real; threshold = 10, full = 18) = clamp((_temperature(day) - threshold) / (full - threshold), 0.0, 1.0)

# ---------------------------------------------------------------- honey bees

"""
    HONEYBEE_PARAMETERS

Default parameters of [`honeybee`](@ref); sources are in the
comments of `src/bees.jl`. Mite and virus values were checked against BEEHAVE's parameter list and code.
"""
const HONEYBEE_PARAMETERS = (
    max_laying = 1600.0,          # eggs per day (HoPoMo, as quoted by BEEHAVE)
    brood_per_nurse = 3.0,        # brood items one in-hive bee can tend (BEEHAVE)
    forager_nursing = 0.2,        # nursing by a forager relative to a hive bee (BEEHAVE)
    egg_survival = 0.97,          # per day (HoPoMo)
    larva_survival = 0.99,        # per day, before cannibalism (HoPoMo)
    cannibalism = (0.23, 0.30, 0.58, 0.06, 0.0),   # of larvae aged 1..5 d when unfed (HoPoMo)
    pupa_survival = 0.999,        # per day (HoPoMo)
    hive_mortality = 0.004,       # in-hive adults per day (BEEHAVE; Russell et al. 2013 use 0.007)
    forager_mortality = 0.10,     # extra per day on a full flying day (Russell et al. 2013: 0.10 to 0.15)
    brood_signal = 1000.0,        # larvae at which recruitment to foraging is half its maximum (proposed)
    recruit_min = 0.25,           # hive bee to forager, per day (Khoury et al. 2013)
    recruit_max = 0.25,           # additional recruitment when food is short (Khoury et al. 2013)
    inhibition = 0.75,            # social inhibition of recruitment by foragers (Khoury et al. 2011)
    food_half = 500.0,            # g of honey at which food shortage is half felt (Khoury et al. 2013)
    nectar_per_forager = 0.14,    # g honey per nectar forager per day at full forage; 0.1 per forager of any kind in
                                  # Khoury et al. (2013) and Russell et al. (2013), raised so a healthy colony gains stores
    pollen_per_forager = 0.12,    # g pollen per pollen forager per day: 15 mg x 10 trips x 0.8 (HoPoMo)
    honey_per_adult = 0.007,      # g per adult per day (Khoury et al. 2013)
    honey_per_larva = 0.0326,     # g per larva per day (Russell et al. 2013)
    pollen_per_adult = 0.0015,    # g per adult per day (BEEHAVE, after Pernal & Currie 2000)
    pollen_per_larva = 0.0284,    # g per larva per day: 142 mg over 5 d (HoPoMo)
    pollen_reserve_days = 6.0,    # days of pollen need the colony tries to hold (HoPoMo)
    max_honey = 60_000.0,         # g the hive can hold (proposed; BEEHAVE uses 50 kg)
    max_pollen = 3_000.0,         # g the hive can hold (proposed)
    starvation_honey = 100.0,     # g of honey below which brood rearing stops (Russell et al. 2013)
    starvation_mortality = 0.2,   # adult mortality per day when honey cannot meet the day's need (proposed)
    critical_adults = 1500.0,     # adults below which the colony fails (Russell et al. 2013)
    failure_mortality = 0.7,      # adult mortality per day below the critical size (Russell et al. 2013: 0.7 for hive bees, 0.8 for foragers, both arbitrary)
    infection_per_cell = 0.89,    # probability a pupa is infected by an infectious mite in its cell (Martin 2001)
    deformed = 0.2,               # fraction of infected pupae that die or emerge deformed (Martin 2001)
    infected_mortality = 0.012,   # in-hive mortality per day of adults infected as pupae (BEEHAVE)
    # drones (all from BEEHAVE's parameter list, which gives the sources named)
    drone_eggs = 0.04,            # drone eggs per worker egg (Wilkinson & Smith 2002, from Allen 1963, 1965)
    drone_season = (115, 240),    # first and last day of drone egg laying (Allen 1963)
    drone_brood = (3, 7, 14),     # days as egg, larva, pupa (Jay 1963; Winston 1987)
    drone_survival = (0.936, 0.956, 0.995),   # per day as egg, larva, pupa (Fukuda & Ohtani 1977)
    drone_mortality = 0.05,       # adult drones per day (Fukuda & Ohtani 1977)
    honey_per_drone = 0.010,      # g per adult drone per day
    honey_per_drone_larva = 0.0192,   # g per drone larva per day
    pollen_per_drone = 0.002,     # g per adult drone per day (an estimate in BEEHAVE)
    pollen_per_drone_larva = 0.050)   # g per drone larva per day (an estimate in BEEHAVE; never measured)

# Fractions of phoretic mites entering worker and drone brood cells in a day (Boot et al.
# 1995 and Calis et al. 1999, as used by Martin 1998 and BEEHAVE): the rates are 0.56 and
# 6.49 x cells being capped / grams of adult worker bees (0.1 g each); a mite enters a cell
# with probability 1 - exp(-(sum of rates)), and the cell is of either kind in proportion to its rate
function _invasion(worker_cells, drone_cells, adult_bees)
    adult_bees >= 1 || return (0.0, 0.0)
    rw, rd = 0.56 * worker_cells / (0.1 * adult_bees), 6.49 * drone_cells / (0.1 * adult_bees)
    rw + rd > 0 || return (0.0, 0.0)
    entering = 1 - exp(-(rw + rd))
    return (entering * rw / (rw + rd), entering * rd / (rw + rd))
end

# A cell takes at most `most` mites (four in Martin 2001 and BEEHAVE); the rest stay on the bees. With
# mites arriving at cells at random (Poisson, `mean` per cell), this is the fraction that gets in.
function _admitted(mean, most = 4)
    mean > 1e-9 || return 1.0
    q, inside, below = exp(-mean), 0.0, 0.0
    for k in 0:(most - 1)
        inside += k * q
        below += q
        q *= mean / (k + 1)
    end
    return (inside + most * (1 - below)) / mean
end

# Fractions of phoretic mites that enter and stay in worker and drone cells, given the cells being
# capped, the adult worker bees and the phoretic mites
function _entering(worker_cells, drone_cells, adult_bees, phoretic)
    w, d = _invasion(worker_cells, drone_cells, adult_bees)
    return (worker_cells > 0 ? w * _admitted(w * phoretic / worker_cells) : 0.0,
            drone_cells > 0 ? d * _admitted(d * phoretic / drone_cells) : 0.0)
end

"""
    honeybee(; virus = false, drones = true, parameters = HONEYBEE_PARAMETERS, forage = nothing,
               flowers = nothing, oral_susceptibility = 0.0, name = :honeybee)

The bees of a honey bee (*Apis mellifera*) colony on a daily time step.

The workers' life cycle is `egg → larva → pupa → hive → forager`, with eggs, larvae and
pupae lasting 3, 5 and 12 days. The queen's laying follows [`hopomo_season`](@ref) and is
capped by the brood the adults can tend; larvae are lost when nurses or pollen are short;
hive bees become foragers at the rate of Khoury et al. (2013) on flying days while there is
brood, so bees emerging after brood rearing stops become long-lived winter bees.

With `drones = true` a second life cycle, `drone_egg → drone_larva → drone_pupa → drone`
(3, 7 and 14 days), is merged in: from late April to late August the queen lays `drone_eggs`
drone eggs for each worker egg. Drones eat and need nursing but do no work. They matter
because mites prefer, and breed better in, drone cells.

The colony's honey and pollen (grams) are two further entries of its state, added as a
module: each day the stores lose what the bees eat and gain what the foragers bring in,
and the share of foragers collecting pollen rises when the pollen store is low.

With `virus = true` the colony has three strata: `healthy`, `carrier` and `infected`. A
larva becomes `infected` when an infectious mite enters its cell at capping (pool
`viral_phoretic`, from [`varroa`](@ref)); a fraction of infected pupae die, and the adults
they become die faster. A healthy in-hive bee becomes a `carrier` when an infectious mite
takes it as a host (pools `viral_phoretic` and `viral_emerging`); as in BEEHAVE, a bee
infected as an adult lives as long as a healthy one but passes the virus to mites.
`oral_susceptibility` above zero lets adults become carriers from virus on flowers (see
[`flower_transmission`](@ref), passed as `flowers`).

`forage` is the result of [`floral_competition`](@ref), if the colony shares flowers.
"""
function honeybee(; virus::Bool = false, drones::Bool = true, parameters = HONEYBEE_PARAMETERS, forage = nothing,
                  flowers = nothing, oral_susceptibility::Real = 0.0, name::Symbol = :honeybee)
    p = parameters
    stages = [:egg, :larva, :pupa, :hive, :forager]
    drones && append!(stages, [:drone_egg, :drone_larva, :drone_pupa, :drone])
    life = ValuedProjectionNet(stages,
        :hatching => [(:egg => :larva) => p.egg_survival],
        :capping => [(:larva => :pupa) => p.larva_survival],
        :emerging => [(:pupa => :hive) => p.pupa_survival],
        :hive_life => [(:hive => :hive) => 1 - p.hive_mortality],
        :recruiting => [(:hive => :forager) => 1 - p.hive_mortality],
        :foraging => [(:forager => :forager) => 1.0])
    durations = Pair{Symbol, Int}[:egg => 3, :larva => 5, :pupa => 12]
    ageing = Dict(:egg => p.egg_survival, :larva => p.larva_survival, :pupa => p.pupa_survival)
    if drones
        # the drones' life cycle is a module on its own stages, merged with the workers'
        life = life ⊕ ValuedProjectionNet(stages,
            :drone_hatching => [(:drone_egg => :drone_larva) => p.drone_survival[1]],
            :drone_capping => [(:drone_larva => :drone_pupa) => p.drone_survival[2]],
            :drone_emerging => [(:drone_pupa => :drone) => p.drone_survival[3]],
            :drone_life => [(:drone => :drone) => 1 - p.drone_mortality])
        append!(durations, [:drone_egg => p.drone_brood[1], :drone_larva => p.drone_brood[2], :drone_pupa => p.drone_brood[3]])
        merge!(ageing, Dict(:drone_egg => p.drone_survival[1], :drone_larva => p.drone_survival[2],
                            :drone_pupa => p.drone_survival[3]))
    end
    bees = expand_stages(ResponsiveProjectionNet(life), durations...; ageing)
    met(kind) = forage === nothing ? (ctx -> 1.0) : getfield(forage.access, kind)
    last_larva, last_pupa = :larva5, :pupa12
    # mites enter worker cells on the last day as a larva and drone cells a day earlier (Boot et al. 1992)
    drone_entry, last_drone_pupa = Symbol(:drone_larva, p.drone_brood[2] - 1), Symbol(:drone_pupa, p.drone_brood[3])

    day(ctx) = day_of_year(ctx.env)
    adults(n) = n[:hive] + n[:forager]
    brood(n) = n[:egg] + n[:larva] + n[:pupa] + (drones ? n[:drone_egg] + n[:drone_larva] + n[:drone_pupa] : 0.0)
    nursing(n) = p.brood_per_nurse * (n[:hive] + p.forager_nursing * n[:forager])
    honey_need(n) = p.honey_per_adult * adults(n) + p.honey_per_larva * n[:larva] +
                    (drones ? p.honey_per_drone * n[:drone] + p.honey_per_drone_larva * n[:drone_larva] : 0.0)
    pollen_need(n) = p.pollen_per_adult * adults(n) + p.pollen_per_larva * n[:larva] +
                     (drones ? p.pollen_per_drone * n[:drone] + p.pollen_per_drone_larva * n[:drone_larva] : 0.0)
    starving(ctx) = ctx.n[:honey] < p.starvation_honey
    failing(ctx) = adults(ctx.n) < p.critical_adults
    pollen_reserve(n) = min(1.0, n[:pollen] / (p.pollen_reserve_days * pollen_need(n) + 1e-9))
    fed(ctx) = pollen_reserve(ctx.n) * min(1.0, nursing(ctx.n) / max(brood(ctx.n), 1.0))
    # more foragers collect pollen when the pollen store is low (proposed)
    pollen_share(n) = min(0.8, 0.05 + 3 * (1 - pollen_reserve(n)))
    nectar_demand(n, env) = p.nectar_per_forager * forage_index(day_of_year(env)) * (1 - pollen_share(n)) * n[:forager]
    pollen_demand(n, env) = p.pollen_per_forager * forage_index(day_of_year(env)) * pollen_share(n) * n[:forager]
    # adults die when the colony is too small to function or the stores cannot feed them
    hardship(ctx) = max(0.0, 1 - (failing(ctx) ? p.failure_mortality : 0.0) -
                             p.starvation_mortality * (1 - min(1.0, ctx.n[:honey] / max(honey_need(ctx.n), 1e-9))))
    # hive bee to forager (Khoury et al. 2013), on flying days while there is brood to rear
    function recruit(ctx)
        n = ctx.n
        R = p.recruit_min + p.recruit_max * p.food_half^2 / (p.food_half^2 + ctx.n[:honey]^2) -
            p.inhibition * n[:forager] / max(adults(n), 1.0)
        return flight_index(day(ctx)) * n[:larva] / (n[:larva] + p.brood_signal) * (1 - exp(-max(R, 0.0)))
    end
    age(stage) = parse(Int, match(r"larva(\d)", String(stage))[1])     # works for larva3 and larva3_infected
    # the queen lays whole eggs (as in BEEHAVE), so there are days in midwinter with no laying and no brood
    function eggs(ctx)
        (starving(ctx) || failing(ctx)) && return 0.0
        laid = min(p.max_laying * hopomo_season(day(ctx)), max(0.0, nursing(ctx.n) - brood(ctx.n)))
        return laid < 0.5 ? 0.0 : laid
    end

    bees = respond(bees,
        :larva_ageing => ctx -> starving(ctx) ? 0.0 : ((from, to) -> 1 - p.cannibalism[age(from)] * (1 - fed(ctx))),
        :capping => ctx -> starving(ctx) ? 0.0 : 1.0,
        :hive_life => ctx -> hardship(ctx) * (1 - recruit(ctx)),
        :recruiting => ctx -> hardship(ctx) * recruit(ctx),
        :foraging => ctx -> hardship(ctx) * max(0.0, 1 - p.hive_mortality - p.forager_mortality * flight_index(day(ctx))))
    bees = inflow(bees, :egg1 => eggs)

    reads = Symbol[]
    brood_stages = drones ? [:egg, :larva, :pupa, :drone_egg, :drone_larva, :drone_pupa] : [:egg, :larva, :pupa]
    occupies = Pair{Symbol, Any}[
        :adult_bees => [:hive => 1.0, :forager => 1.0],
        :brood => [s => 1.0 for s in brood_stages],
        :capping => [last_larva => 1.0],
        :emerging => [last_pupa => 1.0]]
    if drones
        # unfed drone larvae are given up as readily as the most vulnerable worker larvae (proposed)
        bees = respond(bees,
            :drone_larva_ageing => ctx -> starving(ctx) ? 0.0 : 1 - maximum(p.cannibalism) * (1 - fed(ctx)),
            :drone_capping => ctx -> starving(ctx) ? 0.0 : 1.0,
            :drone_life => hardship)
        bees = inflow(bees, :drone_egg1 => ctx -> p.drone_season[1] <= day(ctx) <= p.drone_season[2] ?
                                                  p.drone_eggs * eggs(ctx) : 0.0)
        push!(occupies, :drone_capping => [drone_entry => 1.0], :drone_emerging => [last_drone_pupa => 1.0])
    end

    if virus
        # three strata. Larvae cross to `infected` at capping; adults cross to `carrier` (matrices are [to, from])
        to_carrier = [1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 0.0 1.0]
        to_infected = [1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 0.0 1.0]
        coupling = Dict{Symbol, Any}(:capping => to_infected, :hive_life => to_carrier, :recruiting => to_carrier,
                                     :foraging => to_carrier)
        drones && (coupling[:drone_larva_ageing] = to_infected)
        bees = stratify(bees, [:healthy, :carrier, :infected]; coupling)
        entering(ctx) = _entering(ctx.n[last_larva], drones ? ctx.n[drone_entry] : 0.0, adults(ctx.n), ctx.pool.phoretic)
        # a pupa is infected if any of the viral mites in its cell passes the virus on; mites enter cells at random
        function at_capping(kind, entry)
            return function (ctx)
                cells = ctx.n[entry]
                cells > 0 || return 0.0
                return 1 - exp(-p.infection_per_cell * entering(ctx)[kind] * ctx.pool.viral_phoretic / cells)
            end
        end
        # an in-hive bee becomes a carrier if a viral mite takes it as a host: mites that have just left a
        # cell, and mites whose host died or became a forager (BEEHAVE)
        turnover(ctx) = 1 - hardship(ctx) * (1 - p.hive_mortality) * (1 - recruit(ctx))
        function by_mite(ctx)
            hosts = ctx.n[:hive]
            hosts > 0 || return 0.0
            return 1 - exp(-(turnover(ctx) * ctx.pool.viral_phoretic + ctx.pool.viral_emerging) / hosts)
        end
        by_mouth(ctx) = flowers === nothing ? 0.0 : 1 - exp(-oral_susceptibility * flowers.contamination(ctx))
        in_hive(ctx) = 1 - (1 - by_mite(ctx)) * (1 - by_mouth(ctx))
        faster = (1 - p.infected_mortality) / (1 - p.hive_mortality)
        flying(ctx) = p.forager_mortality * flight_index(day(ctx))
        bees = respond(bees,
            :capping_healthy => ctx -> 1 - at_capping(1, last_larva)(ctx),
            :capping_healthy_to_infected => at_capping(1, last_larva),
            :emerging_infected => ctx -> 1 - p.deformed,
            (Symbol(t, :_healthy) => (ctx -> 1 - in_hive(ctx)) for t in (:hive_life, :recruiting))...,
            (Symbol(t, :_healthy_to_carrier) => in_hive for t in (:hive_life, :recruiting))...,
            :foraging_healthy => ctx -> 1 - by_mouth(ctx),
            :foraging_healthy_to_carrier => by_mouth,
            :hive_life_infected => ctx -> faster,
            :recruiting_infected => ctx -> faster,
            :foraging_infected => ctx -> max(0.0, 1 - p.infected_mortality - flying(ctx)) / max(1 - p.hive_mortality - flying(ctx), 1e-9))
        if drones
            # drones infected as pupae die faster in the same proportion as workers (BEEHAVE); mites do not ride on drones
            sick_drone = min(1.0, p.infected_mortality * p.drone_mortality / p.hive_mortality)
            entered(stage) = startswith(String(stage), String(drone_entry))
            bees = respond(bees,
                # only the larvae whose cells the mites are entering today can be infected
                :drone_larva_ageing_healthy => ctx -> ((from, to) -> entered(from) ? 1 - at_capping(2, drone_entry)(ctx) : 1.0),
                :drone_larva_ageing_healthy_to_infected => ctx -> ((from, to) -> entered(from) ? at_capping(2, drone_entry)(ctx) : 0.0),
                :drone_emerging_infected => ctx -> 1 - p.deformed,
                :drone_life_infected => ctx -> (1 - sick_drone) / (1 - p.drone_mortality))
        end
        push!(occupies, :infected_adults => [:hive_infected => 1.0, :forager_infected => 1.0])
        # what a mite taking a new host needs to know: in-hive bees, how many carry the virus, and
        # how many leave that class each day (by dying or becoming foragers)
        push!(occupies, :hive_bees => [:hive => 1.0])
        push!(occupies, :infected_hive_bees => [:hive_infected => 1.0, :hive_carrier => 1.0])
        push!(occupies, :host_turnover => (n, env) -> n[:hive] * turnover((n = n, env = env)))
        append!(reads, [:phoretic, :viral_phoretic, :viral_emerging])
    end
    # the stores: a module on two more "stages", merged into the colony
    stores = ValuedProjectionNet(vcat(stage_names(bees), [:honey, :pollen]),
        :honey_kept => [(:honey => :honey) => 1.0], :pollen_kept => [(:pollen => :pollen) => 1.0])
    stores = respond(stores,
        :honey_kept => ctx -> max(0.0, 1 - honey_need(ctx.n) / max(ctx.n[:honey], 1e-9)),
        :pollen_kept => ctx -> max(0.0, 1 - pollen_need(ctx.n) / max(ctx.n[:pollen], 1e-9)))
    stores = inflow(stores,
        :honey => ctx -> min(nectar_demand(ctx.n, ctx.env) * met(:hb_nectar)(ctx), max(0.0, p.max_honey - ctx.n[:honey])),
        :pollen => ctx -> min(pollen_demand(ctx.n, ctx.env) * met(:hb_pollen)(ctx), max(0.0, p.max_pollen - ctx.n[:pollen])))
    bees = _add_stages(bees, :honey, :pollen) ⊕ stores
    if forage !== nothing
        push!(occupies, :hb_nectar => (n, env) -> forage.hives * forage.removal * nectar_demand(n, env))
        push!(occupies, :hb_pollen => (n, env) -> forage.hives * pollen_demand(n, env))
        append!(reads, forage.pools)
    end
    if flowers !== nothing
        infected = virus ? :forager_infected : nothing
        push!(occupies, :infectious_visits => (n, env) -> infected === nothing ? 0.0 :
                                                          flowers.hives * flowers.honeybee_visits * n[infected])
    end
    # a colony with less than one adult is dead: the bees go, the stores stay
    adult_idx = bees.groups[:hive] ∪ bees.groups[:forager]
    store_idx = bees.groups[:honey] ∪ bees.groups[:pollen]
    bees = after_step(bees, n -> sum(n[adult_idx]) < 1 ? [i in store_idx ? n[i] : 0.0 for i in eachindex(n)] : n)
    return species(name, bees; occupies, reads = unique(reads))
end

"""
    honeybee_state(c; adults = 10_000, honey = 25_000, pollen = 500, mites = 0, viral = 0.0, queens = 0)

A starting state on 1 January for a community containing a honey bee colony (and
optionally [`varroa`](@ref) and a [`bumblebee`](@ref) population): `adults` in-hive bees,
no brood, stores in grams (the BEEHAVE default is 10,000 workers and 25 kg of honey),
`mites` phoretic mites of which a fraction `viral` carry virus, and `queens` hibernating
bumblebee queens.
"""
function honeybee_state(c::OpenCommunity; adults::Real = 10_000, honey::Real = 25_000, pollen::Real = 500,
                        mites::Real = 0, viral::Real = 0.0, queens::Real = 0)
    n = zeros(length(c))
    function put!(sp, stage, value)
        sp in c.names || return
        groups = c.system[sp].groups
        haskey(groups, stage) || return
        n[c.layout[sp][first(groups[stage])]] = value
    end
    put!(:honeybee, :hive, adults)                   # the first (healthy) hive class
    put!(:honeybee, :honey, honey)
    put!(:honeybee, :pollen, pollen)
    if :varroa in c.names
        clean = haskey(c.system[:varroa].groups, :phoretic_clean) ? :phoretic_clean : :phoretic
        put!(:varroa, clean, mites * (1 - viral))
        viral > 0 && put!(:varroa, :phoretic_viral, mites * viral)
    end
    put!(:bumblebee, :hibernating, queens)
    return n
end

# ---------------------------------------------------------------- Varroa

"""
    VARROA_PARAMETERS

Default parameters of [`varroa`](@ref). The values are those of Martin (1998), checked
against the parameter list and code of BEEHAVE (Becher et al. 2014), which cite them.
"""
const VARROA_PARAMETERS = (
    offspring = 1.01,             # viable daughters per foundress per cycle in a worker cell (Martin 1994, 1998)
    drone_offspring = 2.91,       # the same in a drone cell (Martin 1998)
    emergence_loss = 0.30,        # mites that fall and die when the bee emerges from a worker cell (Martin 1998)
    drone_emergence_loss = 0.20,  # the same for a drone cell (Martin 1998)
    # daughters per foundress, relative to a foundress alone, with one to four foundresses in a cell (Martin 1998, Table 4)
    sharing = (1.0, 0.91, 0.86, 0.60),
    drone_sharing = (1.0, 0.84, 0.65, 0.66),
    mortality_brood = 0.006,      # phoretic mites per day while the colony has brood (Martin 1998; Fries et al. 1994)
    mortality_broodless = 0.002)  # phoretic mites per day while it has none (Martin 1998)

"""
    varroa(; virus = false, drones = true, parameters = VARROA_PARAMETERS, treatment = day -> 0.0, name = :varroa)

The mite *Varroa destructor* in a honey bee colony. A mite is `phoretic` (on an adult bee)
or in a brood cell, a worker `cell` for 12 days or a `drone_cell` for 15, until the bee
emerges; it then leaves with its daughters. Each day a fraction of the phoretic mites
enters the cells about to be capped; the fractions depend on the pools `capping`,
`drone_capping` and `adult_bees` reported by [`honeybee`](@ref), and a cell takes at most
four mites. Mites enter drone cells almost twelve times as readily as worker cells and
leave nearly three times as many daughters there.
Phoretic mites die faster while the colony has brood (pool `brood`). `drones` must be the
same as in [`honeybee`](@ref).

With `virus = true` mites are stratified into `clean` and `viral`. A clean mite becomes
viral by entering a cell whose pupa is infected, or by taking a bee that carries the virus
as its host (pools `host_turnover`, `hive_bees`, `infected_hive_bees`). Taking a host
follows BEEHAVE: a clean mite that takes a new host, because it has just left a cell or
because its host died or became a forager, becomes viral with a probability equal to the
fraction of in-hive bees carrying the virus. A viral mite stays viral. Its daughters are
viral if the pupa in their cell was infected and lived; if it died there are no daughters,
and if it escaped infection they are clean.

`treatment(day_of_year)` is an extra daily mortality of phoretic mites. `infection_per_cell`,
`pupal_death` and `brood_survival` (the chance that worker and drone brood survive the days
a mite spends in the cell) repeat values of the honey bee model and default to them.
"""
function varroa(; virus::Bool = false, drones::Bool = true, parameters = VARROA_PARAMETERS, treatment = day -> 0.0,
                infection_per_cell = HONEYBEE_PARAMETERS.infection_per_cell,
                pupal_death = HONEYBEE_PARAMETERS.deformed,
                brood_survival = (HONEYBEE_PARAMETERS.pupa_survival^12,
                                  HONEYBEE_PARAMETERS.drone_survival[2] * HONEYBEE_PARAMETERS.drone_survival[3]^14),
                name::Symbol = :varroa)
    p = parameters
    # a foundress whose brood dies in the cell leaves it without daughters (BEEHAVE)
    worker_daughters, drone_daughters = p.offspring * brood_survival[1], p.drone_offspring * brood_survival[2]
    stages = drones ? [:phoretic, :cell, :drone_cell] : [:phoretic, :cell]
    life = ValuedProjectionNet(stages,
        :waiting => [(:phoretic => :phoretic) => 1.0],
        :invading => [(:phoretic => :cell) => 1.0],
        :leaving => [(:cell => :phoretic) => 1 - p.emergence_loss],
        :daughters => [(:cell => :phoretic) => worker_daughters * (1 - p.emergence_loss)])
    durations = Pair{Symbol, Int}[:cell => 12]
    # one entry per kind of cell: (invading, leaving, daughters, last day in the cell, cells capped, cells emerging, table)
    kinds = [(:invading, :leaving, :daughters, :cell12, :capping, :emerging, p.sharing)]
    if drones
        life = life ⊕ ValuedProjectionNet(stages,
            :invading_drone => [(:phoretic => :drone_cell) => 1.0],
            :leaving_drone => [(:drone_cell => :phoretic) => 1 - p.drone_emergence_loss],
            :daughters_drone => [(:drone_cell => :phoretic) => drone_daughters * (1 - p.drone_emergence_loss)])
        push!(durations, :drone_cell => 15)
        push!(kinds, (:invading_drone, :leaving_drone, :daughters_drone, :drone_cell15, :drone_capping, :drone_emerging,
                      p.drone_sharing))
    end
    mites = births(expand_stages(ResponsiveProjectionNet(life), durations...), (k[3] for k in kinds)...)

    invade(ctx) = _entering(ctx.pool.capping, drones ? ctx.pool.drone_capping : 0.0, ctx.pool.adult_bees, ctx.n[:phoretic])
    function survive(ctx)
        ctx.pool.adult_bees >= 1 || return 0.0            # no hosts left
        mortality = ctx.pool.brood >= 1 ? p.mortality_brood : p.mortality_broodless
        return (1 - mortality) * (1 - treatment(day_of_year(ctx.env)))
    end
    # Foundresses sharing a cell leave fewer daughters each. With mites spread over cells at random, a
    # foundress has j companions with the Poisson probability for the mean number of mites per cell;
    # more than three companions count as three.
    function shared(last, cells, table)
        return function (ctx)
            emerging = getproperty(ctx.pool, cells)
            emerging > 0 || return 1.0
            mean = ctx.n[last] / emerging
            q, total, rest = exp(-mean), 0.0, 1.0
            for j in 1:(length(table) - 1)
                total += q * table[j]
                rest -= q
                q *= mean / j
            end
            return total + max(rest, 0.0) * table[end]
        end
    end
    mites = respond(mites,
        :waiting => ctx -> survive(ctx) * (1 - sum(invade(ctx))),
        (k[1] => (ctx -> survive(ctx) * invade(ctx)[i]) for (i, k) in enumerate(kinds))...,
        (k[3] => shared(k[4], k[6], k[7]) for k in kinds)...)

    occupies = Pair{Symbol, Any}[:phoretic => [:phoretic => 1.0]]
    reads = [:capping, :adult_bees, :emerging, :brood]
    drones && append!(reads, [:drone_capping, :drone_emerging])
    if virus
        one_way, both_ways = [1.0 0.0; 1.0 1.0], [1.0 1.0; 1.0 1.0]      # [to, from]
        coupling = Dict{Symbol, Any}(:waiting => one_way)
        for k in kinds
            coupling[k[1]] = one_way
            coupling[k[2]] = one_way
            coupling[k[3]] = both_ways
        end
        mites = stratify(mites, [:clean, :viral]; coupling)
        # a clean mite meets the virus in a cell if any viral mite sharing it infects the pupa
        function in_cell(i, cells)
            return function (ctx)
                capped = getproperty(ctx.pool, cells)
                capped > 0 || return 0.0
                return 1 - exp(-infection_per_cell * invade(ctx)[i] * ctx.n[:phoretic_viral] / capped)
            end
        end
        # a mite taking a new host meets the virus with the prevalence among in-hive bees (as in BEEHAVE)
        prevalence(ctx) = ctx.pool.hive_bees > 0 ? clamp(ctx.pool.infected_hive_bees / ctx.pool.hive_bees, 0.0, 1.0) : 0.0
        switching(ctx) = ctx.pool.hive_bees > 0 ? clamp(ctx.pool.host_turnover / ctx.pool.hive_bees, 0.0, 1.0) : 0.0
        on_adults(ctx) = switching(ctx) * prevalence(ctx)
        # the pupa of a viral foundress is infected and lives (viral daughters), is infected and dies (none),
        # or escapes (clean daughters, which may then take an infected host)
        viral_daughters = infection_per_cell * (1 - pupal_death)
        mites = respond(mites, :waiting_clean => ctx -> 1 - on_adults(ctx), :waiting_clean_to_viral => on_adults)
        for (i, k) in enumerate(kinds)
            enter, leave, born = k[1], k[2], k[3]
            mites = respond(mites,
                Symbol(enter, :_clean) => ctx -> 1 - in_cell(i, k[5])(ctx), Symbol(enter, :_clean_to_viral) => in_cell(i, k[5]),
                Symbol(leave, :_clean) => ctx -> 1 - prevalence(ctx), Symbol(leave, :_clean_to_viral) => prevalence,
                Symbol(born, :_clean) => ctx -> 1 - prevalence(ctx), Symbol(born, :_clean_to_viral) => prevalence,
                Symbol(born, :_viral) => ctx -> viral_daughters + (1 - infection_per_cell) * prevalence(ctx),
                Symbol(born, :_viral_to_clean) => ctx -> (1 - infection_per_cell) * (1 - prevalence(ctx)))
        end
        push!(occupies, :viral_phoretic => [:phoretic_viral => 1.0])
        # viral mites that will leave their cells today and look for a host (crowding in the cell is ignored here)
        leaving = [(Symbol(:cell12, :_viral), (1 - p.emergence_loss) * (1 + worker_daughters * viral_daughters))]
        drones && push!(leaving, (Symbol(:drone_cell15, :_viral),
                                  (1 - p.drone_emergence_loss) * (1 + drone_daughters * viral_daughters)))
        push!(occupies, :viral_emerging => [s => w for (s, w) in leaving])
        append!(reads, [:hive_bees, :infected_hive_bees, :host_turnover])
    end
    mites = after_step(mites, n -> [x < 1e-6 ? 0.0 : x for x in n])
    return species(name, mites; occupies, reads)
end

# ---------------------------------------------------------------- bumblebees

"""
    BUMBLEBEE_PARAMETERS

Default parameters of [`bumblebee`](@ref); sources are in the comments of `src/bees.jl`.
"""
const BUMBLEBEE_PARAMETERS = (
    emergence_mean = 91.0,        # day of year queens leave hibernation (Bumble-BEEHAVE, B. terrestris)
    emergence_sd = 28.0,          # (Bumble-BEEHAVE)
    last_emergence = 200,         # no emergence after this day: later queens are this year's (proposed)
    season_end = 305,             # day on which all bees not in hibernation die (Bumble-BEEHAVE)
    winter_survival = 0.64,       # of hibernating queens (Bumble-BEEHAVE maximum, after Beekman et al. 1998)
    nest_finding = 0.2,           # per day (Bumble-BEEHAVE)
    nest_sites = 30.0,            # nest sites per square kilometre (Knight et al. 2005: 29 nests per km2; Bumble-BEEHAVE simulates 34)
    search_mortality = 0.155,     # per day while searching (derived from Bumble-BEEHAVE)
    founding_days = 28,           # nest founding to first workers (measured: 21 and 33 d, as cited by Becher et al. 2018)
    founding_mortality = 0.04,    # per day for a queen foraging alone (calibrated to establishment of about 0.12)
    first_workers = 8.0,          # workers in the first brood (Cresswell 2017: 8.25)
    queen_mortality = 0.0154,     # per day in the nest (Banks et al. 2017)
    eggs = 1.6,                   # eggs per colony per day = eggs x (workers per colony)^husbandry (after Cresswell 2017)
    husbandry = 0.54,             # (Cresswell 2017)
    worker_brood = (4, 9, 9),     # days as egg, larva, pupa (Banks et al. 2020)
    male_brood = (4, 11, 11),
    gyne_brood = (4, 13, 13),
    brood_survival = 0.995,       # per day when fully fed (proposed)
    starvation = 0.3,             # larvae lost per day when unfed (proposed)
    gyne_sensitivity = 2.0,       # how much more a gyne larva suffers from a shortfall (proposed)
    switch = 0.08,                # per day, colony starts rearing males and gynes (Bumble-BEEHAVE: 0.13 once conditions are met)
    ending = 1 / 9,               # per day, reproductive colony stops rearing brood (about 9 d, derived from Becher et al. 2018, Table 3)
    gyne_fraction = 0.45,         # of eggs in the reproductive phase; set to the share of gynes among reproductives reared (Becher et al. 2018, Table 3)
    gyne_survival = 0.9,          # new gyne mates and enters hibernation (proposed)
    worker_mortality = 0.04,      # per day (Banks et al. 2020: 0.05; Cresswell 2017 field 0.048, not checked; lower in the nest)
    male_mortality = 0.05,        # per day (proposed)
    shortage_mortality = 1.0,     # proportional rise in forager mortality when nectar is unobtainable (proposed)
    forager_fraction = 0.5,       # of workers (Haeussler et al. 2017)
    queen_foraging = 3.0,         # nectar and pollen a founding queen collects, relative to a worker (crop volume, Bumble-BEEHAVE)
    nectar_per_forager = 0.31,    # g sugar per day: 0.6 ml at 1.5 mol/l sucrose, the top of the range for flowers in Bumble-BEEHAVE.
                                  # Banks et al. (2020) give 0.6 ml per worker, not per forager: see nectar_per_worker
    pollen_per_forager = 0.11,    # g per day: 0.018 g load (Bumble-BEEHAVE) x 5.9 trips (Minahan & Brunet 2018)
    nectar_per_colony = 0.6,      # g sugar per day to incubate the brood (about 10 kJ per day, Bumble-BEEHAVE)
    nectar_per_worker = 0.077,    # g sugar per day (0.15 ml, Banks et al. 2017, not checked). Banks et al. (2020) use 0.35 ml with
                                  # 0.6 ml collected per worker; that pair (0.18 here, and 0.62 per forager) gives nearly the same
                                  # results, but 0.18 with the collection rate left at 0.31 starves every colony
    pollen_per_worker = 0.15,     # g to rear a worker or male (Bumble-BEEHAVE mass balance)
    pollen_per_gyne = 0.8,        # g to rear a gyne (Bumble-BEEHAVE mass balance)
    infected_hazard = 1.0,        # ratio of mortality of DWV-infected to healthy workers (evidence conflicting)
    visits = 1000.0)              # flower visits per forager per day (proposed)

"""
    bumblebee(; virus = false, parameters = BUMBLEBEE_PARAMETERS, forage = nothing,
                flowers = nothing, name = :bumblebee)

A population of bumblebee (*Bombus terrestris*) colonies on a daily time step; numbers can
be read as per square kilometre.

Queens hibernate, emerge in spring, search for a nest site and found a nest alone; after
`founding_days` the nest has its first workers. A colony first rears workers (`growing`),
then males and new queens (`reproducing`), then stops (`ending`). The number of nests is
limited by `nest_sites`. Worker, male and gyne
brood take 22, 26 and 30 days. New queens mate and stay dormant, joining the hibernating
queens when the season ends; everything else dies then. Colonies are pooled: workers per colony is the total number
of workers over the total number of colonies.

Bumblebees keep almost no stores, so each day's intake is compared with each day's need:
larvae are lost when pollen or nectar is short (gyne larvae more readily), and foragers
die faster when nectar is hard to find.

With `virus = true` there are infected workers, infected at flowers with the daily
probability given by `flowers` (from [`flower_transmission`](@ref)); they die
`infected_hazard` times as fast. There is no transmission within the nest, and no mites.

`forage` is the result of [`floral_competition`](@ref), if the population shares flowers.
"""
function bumblebee(; virus::Bool = false, parameters = BUMBLEBEE_PARAMETERS, forage = nothing,
                   flowers = nothing, name::Symbol = :bumblebee)
    p = parameters
    stay = 1 - p.queen_mortality
    stages = [:hibernating, :searching, :founding, :growing, :reproducing, :ending,
              :worker_brood, :worker, :male_brood, :male, :gyne_brood, :new_queen]
    queens = ValuedProjectionNet(stages,
        :dormancy => [(:hibernating => :hibernating) => 1.0],
        :waiting => [(:new_queen => :new_queen) => 1.0],
        :overwintering => [(:new_queen => :hibernating) => 1.0],
        :emergence => [(:hibernating => :searching) => p.winter_survival],
        :searching => [(:searching => :searching) => 1 - p.nest_finding - p.search_mortality],
        :nesting => [(:searching => :founding) => p.nest_finding],
        :first_workers => [(:founding => :growing) => 1.0],
        :colony_cycle => [(:growing => :growing) => stay * (1 - p.switch), (:growing => :reproducing) => stay * p.switch,
                          (:reproducing => :reproducing) => stay * (1 - p.ending), (:reproducing => :ending) => stay * p.ending,
                          (:ending => :ending) => stay])
    rearing = ValuedProjectionNet(stages,
        :first_brood => [(:founding => :worker) => p.first_workers],
        :worker_eggs => [(:growing => :worker_brood) => p.eggs],
        :sexual_eggs => [(:reproducing => :gyne_brood) => p.gyne_fraction * p.eggs,
                         (:reproducing => :male_brood) => (1 - p.gyne_fraction) * p.eggs],
        :eclosion => [(:worker_brood => :worker) => p.brood_survival, (:male_brood => :male) => p.brood_survival,
                      (:gyne_brood => :new_queen) => p.brood_survival * p.gyne_survival])
    adults = ValuedProjectionNet(stages,
        :worker_life => [(:worker => :worker) => 1.0], :male_life => [(:male => :male) => 1 - p.male_mortality])
    colony = births(ResponsiveProjectionNet(queens ⊕ rearing ⊕ adults), :first_brood, :worker_eggs, :sexual_eggs)
    if virus
        illness = ValuedProjectionNet(vcat(stages, :infected_worker),
            :infection => [(:worker => :infected_worker) => 1.0],
            :infected_life => [(:infected_worker => :infected_worker) => 1.0])
        colony = _add_stages(colony, :infected_worker) ⊕ illness
    end
    brood = (:worker_brood => p.worker_brood, :male_brood => p.male_brood, :gyne_brood => p.gyne_brood)
    colony = expand_stages(colony, :founding => p.founding_days, (s => sum(d) for (s, d) in brood)...;
                           ageing = Dict(s => p.brood_survival for (s, _) in brood))

    day(ctx) = day_of_year(ctx.env)
    alive(ctx) = day(ctx) < p.season_end ? 1.0 : 0.0
    workers(n) = n[:worker] + (virus ? n[:infected_worker] : 0.0)
    colonies(n) = n[:growing] + n[:reproducing] + n[:ending]
    foragers(n) = p.forager_fraction * workers(n) + p.queen_foraging * n[:founding]
    larvae(n, s, d) = sum(n[Symbol(s, i)] for i in (d[1] + 1):(d[1] + d[2]))
    is_larva = Dict(Symbol(s, i) => (d[1] < i <= d[1] + d[2]) for (s, d) in brood for i in 1:sum(d))
    flying(env) = flight_index(day_of_year(env); threshold = 6, full = 12)
    nectar_demand(n, env) = p.nectar_per_forager * flying(env) * foragers(n)
    pollen_demand(n, env) = p.pollen_per_forager * flying(env) * foragers(n)
    nectar_need(n) = p.nectar_per_colony * (colonies(n) + n[:founding]) + p.nectar_per_worker * workers(n)
    pollen_need(n) = p.pollen_per_worker / p.worker_brood[2] * larvae(n, :worker_brood, p.worker_brood) +
                     p.pollen_per_worker / p.male_brood[2] * larvae(n, :male_brood, p.male_brood) +
                     p.pollen_per_gyne / p.gyne_brood[2] * larvae(n, :gyne_brood, p.gyne_brood)
    met(kind) = forage === nothing ? (ctx -> 1.0) : getfield(forage.access, kind)
    ratio(intake, need) = need > 0 ? min(1.0, intake / need) : 1.0
    nectar_fed(ctx) = ratio(nectar_demand(ctx.n, ctx.env) * met(:bb_nectar)(ctx), nectar_need(ctx.n))
    pollen_fed(ctx) = ratio(pollen_demand(ctx.n, ctx.env) * met(:bb_pollen)(ctx), pollen_need(ctx.n))
    hunger(ctx) = 1 - min(nectar_fed(ctx), pollen_fed(ctx))
    # daily chance that a hibernating queen emerges (logistic approximation to a normal distribution of dates)
    emerge(ctx) = day(ctx) > p.last_emergence ? 0.0 :
                  min(1.0, 1.702 / p.emergence_sd / (1 + exp(-1.702 * (day(ctx) - p.emergence_mean) / p.emergence_sd)))
    hard_times(ctx) = 1 + p.shortage_mortality * (1 - nectar_fed(ctx))
    # a queen can only settle where there is a free nest site
    free_sites(ctx) = max(0.0, 1 - (colonies(ctx.n) + ctx.n[:founding]) / p.nest_sites)
    founding(ctx) = alive(ctx) * max(0.0, 1 - p.founding_mortality * hard_times(ctx))
    working(hazard) = ctx -> alive(ctx) * max(0.0, 1 - hazard * p.worker_mortality * hard_times(ctx))
    husbandry(ctx) = max(workers(ctx.n) / max(colonies(ctx.n), 1e-9), 1.0)^p.husbandry
    caught(ctx) = (virus && flowers !== nothing) ? flowers.infection(ctx) : 0.0
    # unfed larvae die; gyne larvae, which need most, die first
    rationing(sensitivity) = ctx -> ((from, to) -> alive(ctx) *
        (is_larva[from] ? max(0.0, 1 - min(1.0, sensitivity * p.starvation * hunger(ctx))) : 1.0))

    colony = respond(colony,
        :dormancy => ctx -> 1 - emerge(ctx), :emergence => emerge,
        # this year's new queens stay dormant until the season ends, then join the hibernating queens
        :waiting => alive, :overwintering => ctx -> 1 - alive(ctx),
        # a queen that finds no free site keeps searching
        :searching => ctx -> alive(ctx) * (1 - p.search_mortality - p.nest_finding * free_sites(ctx)) /
                             (1 - p.search_mortality - p.nest_finding),
        :nesting => ctx -> alive(ctx) * free_sites(ctx), :colony_cycle => alive, :male_life => alive, :eclosion => alive,
        :founding_ageing => founding, :first_workers => founding, :first_brood => founding,
        :worker_eggs => husbandry, :sexual_eggs => husbandry,
        :worker_brood_ageing => rationing(1.0), :male_brood_ageing => rationing(1.0),
        :gyne_brood_ageing => rationing(p.gyne_sensitivity),
        :worker_life => ctx -> working(1.0)(ctx) * (1 - caught(ctx)))
    virus && (colony = respond(colony, :infection => ctx -> working(1.0)(ctx) * caught(ctx),
                               :infected_life => working(p.infected_hazard)))

    occupies = Pair{Symbol, Any}[]
    reads = Symbol[]
    if forage !== nothing
        push!(occupies, :bb_nectar => nectar_demand, :bb_pollen => pollen_demand)
        append!(reads, forage.pools)
    end
    if flowers !== nothing
        push!(occupies, :infectious_visits => (n, env) -> virus ?
              flowers.bumblebee_deposit * p.visits * p.forager_fraction * n[:infected_worker] : 0.0)
    end
    return species(name, colony; occupies, reads = unique(reads))
end

# A responsive net with more (empty) stages at the end, so that a module using them can be merged in.
function _add_stages(net::ResponsiveProjectionNet, stages::Symbol...)
    base = ValuedProjectionNet(vcat(stage_names(net), collect(stages)),
        (t => net.base.transition_values[t] for t in transition_names(net))...)
    return ResponsiveProjectionNet(base; responses = net.responses, inflows = net.inflows, after = net.after,
                                   fecundity = net.fecundity, groups = copy(net.groups))
end

# ---------------------------------------------------------------- shared flowers

"""
    floral_supply(peak; start = 80, finish = 280, ramp = 20)

A seasonal supply function `day -> amount per day`: zero before `start` and after
`finish`, rising to `peak` over `ramp` days and falling over the last `ramp` days.
"""
floral_supply(peak::Real; start::Real = 80, finish::Real = 280, ramp::Real = 20) =
    day -> peak * clamp(min((day - start) / ramp, (finish - day) / ramp), 0.0, 1.0)

"""
    floral_competition(mode = :scramble; nectar, pollen, hives = 1, removal = 1,
                       overlap = 0.5, honeybee_share = 0.5, refuge = 0.2, carryover = 0.5)

Link honey bees and bumblebees through floral nectar and pollen. `nectar` and `pollen` are
supply functions of the day of the year (grams of sugar and of pollen per day, see
[`floral_supply`](@ref)). `hives` is the number of honey bee colonies drawing on the
supply (per square kilometre, if bumblebees are per square kilometre), each identical to
the one modelled; `removal` converts a gram of honey gained into
grams of floral sugar removed (a colony collects about 120 kg of nectar a year, some 48 kg
of sugar, which is the sugar in about 60 kg of honey: roughly 0.8 to 1).

The default supplies peak at 3 kg of sugar and 3 kg of pollen per square kilometre per day.
That is of the order measured for intensively farmed land; flower-rich habitat offers 3 to
60 times as much (Hicks et al. 2016; Baude et al. 2016), with pollen and sugar in roughly
equal mass.

Each species reports its demand to the pools `hb_nectar`, `hb_pollen`, `bb_nectar` and
`bb_pollen`, and receives the fraction of its demand that is met. `mode` decides how:

- `:scramble`: everyone gets the same fraction, `min(1, supply / total demand)`.
- `:standing_crop`: the supply is a stock (a `flowers` species, pools `nectar_stock` and
  `pollen_stock`) that is secreted daily and of which a fraction `carryover` of what is
  not removed persists, so heavy use one day lowers what is there the next.
- `:partitioned`: a fraction `overlap` of the supply is open to both species; the rest is
  exclusive, a share `honeybee_share` of it to honey bees and the remainder to bumblebees.
- `:preemptive`: honey bees take what they want first, except from a fraction `refuge` of
  the supply that only bumblebees can use (deep flowers); bumblebees get what is left.

Pass the result as `forage` to [`honeybee`](@ref) and [`bumblebee`](@ref), and add its `species` (the flowers, if any) to the community.
"""
function floral_competition(mode::Symbol = :scramble; nectar = floral_supply(3000.0), pollen = floral_supply(3000.0),
                            hives::Real = 1, removal::Real = 1, overlap::Real = 0.5, honeybee_share::Real = 0.5,
                            refuge::Real = 0.2, carryover::Real = 0.5)
    mode in (:scramble, :standing_crop, :partitioned, :preemptive) || throw(ArgumentError("unknown mode :$mode"))
    frac(have, want) = want > 0 ? min(1.0, max(have, 0.0) / want) : 1.0
    # fractions of honey bee and bumblebee demand met, given both demands and what is on offer
    function met(dh, db, s)
        if mode === :scramble || mode === :standing_crop
            f = frac(s, dh + db)
            return f, f
        elseif mode === :partitioned
            eh, eb = (1 - overlap) * honeybee_share * s, (1 - overlap) * (1 - honeybee_share) * s
            rh, rb = max(dh - eh, 0.0), max(db - eb, 0.0)
            f = frac(overlap * s, rh + rb)
            return frac(min(dh, eh) + f * rh, dh), frac(min(db, eb) + f * rb, db)
        else
            open = (1 - refuge) * s
            return frac(open, dh), frac(refuge * s + max(open - dh, 0.0), db)
        end
    end
    offer(ctx, supply, stock) = mode === :standing_crop ? getproperty(ctx.pool, stock) : supply(day_of_year(ctx.env))
    nectar_met(ctx) = met(ctx.pool.hb_nectar, ctx.pool.bb_nectar, offer(ctx, nectar, :nectar_stock))
    pollen_met(ctx) = met(ctx.pool.hb_pollen, ctx.pool.bb_pollen, offer(ctx, pollen, :pollen_stock))
    access = (hb_nectar = ctx -> nectar_met(ctx)[1], bb_nectar = ctx -> nectar_met(ctx)[2],
              hb_pollen = ctx -> pollen_met(ctx)[1], bb_pollen = ctx -> pollen_met(ctx)[2])
    pools = [:hb_nectar, :hb_pollen, :bb_nectar, :bb_pollen]
    extra = OpenCommunity[]
    if mode === :standing_crop
        append!(pools, [:nectar_stock, :pollen_stock])
        left(stock, dh, db) = ctx -> carryover * max(0.0, 1 - (getproperty(ctx.pool, dh) + getproperty(ctx.pool, db)) /
                                                           max(getproperty(ctx.pool, stock), 1e-9))
        flowers = ValuedProjectionNet([:nectar, :pollen],
            :nectar_left => [(:nectar => :nectar) => 1.0], :pollen_left => [(:pollen => :pollen) => 1.0])
        flowers = respond(flowers, :nectar_left => left(:nectar_stock, :hb_nectar, :bb_nectar),
                          :pollen_left => left(:pollen_stock, :hb_pollen, :bb_pollen))
        flowers = inflow(flowers, :nectar => ctx -> (1 - carryover) * nectar(day_of_year(ctx.env)),
                         :pollen => ctx -> (1 - carryover) * pollen(day_of_year(ctx.env)))
        push!(extra, species(:flowers, flowers; occupies = [:nectar_stock => [:nectar => 1.0], :pollen_stock => [:pollen => 1.0]],
                             reads = [:hb_nectar, :hb_pollen, :bb_nectar, :bb_pollen]))
    end
    return (mode = mode, hives = hives, removal = removal, access = access, pools = pools, species = extra)
end

"""
    flower_transmission(; flowers = floral_supply(5.0e6), hives = 1, honeybee_visits = 500,
                        deposit = 0.02, bumblebee_deposit = 0.0, acquisition = 1.0)

Link honey bees and bumblebees through virus left on shared flowers. Infected foragers
report their visits to the pool `infectious_visits` (bumblebee visits weighted by
`bumblebee_deposit`, zero by default: bumblebees pick the virus up but do not pass it on).
The fraction of flowers contaminated on a day is `1 - exp(-deposit × infectious visits /
flowers(day))`, so more flowers dilute the virus. A bumblebee worker is infected with
daily probability `1 - exp(-acquisition × contamination)`; a honey bee colony's exposure
is the contamination itself, which matters only if its `oral_susceptibility` is above zero.

None of these rates has been measured. `deposit` is set so that an infected forager leaves
virus on about ten flowers a day (the model fitted by Burnham et al. 2021 has four a day,
lasting ten days), and `acquisition` so that roughly a tenth of bumblebee workers carry the
virus beside an untreated colony, as in field surveys (Fürst et al. 2014; Alger et al.
2019; Manley et al. 2019). Only their product matters to the bumblebees.

Pass the result as `flowers` to [`honeybee`](@ref) and [`bumblebee`](@ref).
"""
function flower_transmission(; flowers = floral_supply(5.0e6), hives::Real = 1, honeybee_visits::Real = 500,
                             deposit::Real = 0.02, bumblebee_deposit::Real = 0.0, acquisition::Real = 1.0)
    function contamination(ctx)
        f = flowers(day_of_year(ctx.env))
        return f > 0 ? 1 - exp(-deposit * ctx.pool.infectious_visits / f) : 0.0
    end
    return (hives = hives, honeybee_visits = honeybee_visits, bumblebee_deposit = bumblebee_deposit,
            contamination = contamination, infection = ctx -> 1 - exp(-acquisition * contamination(ctx)))
end
