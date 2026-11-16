# The length-based multispecies fish model LeMans (Hall et al. 2006), North Sea version (Thorpe et al. 2015).
#
# Fish of each species are counted in length classes of equal width, shared by all species. A time step is a
# tenth of a year. In each step, in this order:
#
#   1. at the first step of a year, recruits (hockey-stick in the spawning stock biomass) join the smallest
#      length class;
#   2. fish die of background mortality, of fishing and of predation, all at once: a share exp(-Z) survives,
#      with Z the sum of the three per-step rates;
#   3. a fixed share of the survivors of each class grows into the next, as von Bertalanffy growth requires.
#
# Predation is the coupling between species. A predator of a given species and length class needs a ration
# of food each step (what it takes to grow as von Bertalanffy says), and takes it from the prey it finds
# suitable: suitability is a log-normal preference for prey weighing about 10^-2.25 of the predator's
# weight, allowed only for the species pairs of a diet matrix, and scaled to sum to one over all prey. It
# also eats "other food", a constant amount of prey outside the model. Prey of each class die at the rate
#
#   M2(prey) = Σ_predator classes  ration × N(predator) × suitability(predator, prey) / (Σ_prey suitability × weight × N + other)
#
# (a Holling type II response). Both sums run over the length classes of every species, so the model needs a
# pool with a value for each predator class: the food the predators of that class need (`ration`) and the
# prey biomass they find suitable (`prey`). Each species adds its own fish to both, as a predator and as prey,
# and reads both to know how fast it is eaten. No species refers to another.
#
# The equations, defaults, units and the order of operations are those of the R package LeMaRns 0.1.2
# (Spence et al. 2020), whose data sets NS_par, NS_tau, NS_other, NS_mixed_fish, NS_eta and NS_L50 are in
# data/lemans (GPL-3; see data/lemans/LICENSE.txt). Lengths are in cm, weights in g, numbers are of fish in
# the North Sea.

const _LEMANS_DIR = joinpath(DATA_DIR, "lemans")

_lemans_symbol(name::AbstractString) = Symbol(lowercase(replace(strip(name), " " => "_")))

function _lemans_north_sea(dir)
    table(file) = (lines = [split(l, ",") for l in readlines(joinpath(dir, file)) if !isempty(strip(l))];
        (String.(lines[1]), lines[2:end]))
    num(x) = parse(Float64, x)
    sh, rows = table("north_sea_species.csv")
    col(name) = (j = findfirst(==(name), sh); j === nothing ? throw(ArgumentError("column $name missing")) : [r[j] for r in rows])
    names = String.(col("species"))
    dh, drows = table("north_sea_diet.csv")
    dh[2:end] == names || throw(ArgumentError("the diet matrix does not list the species in order"))
    [String(r[1]) for r in drows] == names || throw(ArgumentError("the diet matrix does not list the predators in order"))
    diet = [parse(Int, drows[i][j + 1]) == 1 for i in eachindex(names), j in eachindex(names)]
    _, orows = table("north_sea_other_food.csv")
    return (species = _lemans_symbol.(names), names = names,
            Linf = num.(col("Linf")), W_a = num.(col("W_a")), W_b = num.(col("W_b")), k = num.(col("k")),
            Lmat = num.(col("Lmat")), a = num.(col("a")), b = num.(col("b")),
            gear = _lemans_symbol.(col("gear")), max_catchability = num.(col("max_catchability")),
            eta = num.(col("eta")), L50 = num.(col("L50")), diet = diet, other = num(only(only(orows))))
end

"""
    LEMANS_NORTH_SEA

The North Sea data of the LeMans model for [`lemans_parameters`](@ref): the 21 species of
Thorpe et al. (2015), as shipped with the R package LeMaRns 0.1.2 (Spence et al. 2020;
GPL-3) and read from `data/lemans`. A `NamedTuple` of vectors, one entry per species:

- `species` (symbols such as `:cod`, `:poor_cod`) and `names` ("Cod", "Poor cod");
- `Linf` (asymptotic length, cm) and `k` (per year), the von Bertalanffy growth curve;
- `W_a`, `W_b`: weight in g is `W_a L^W_b` for length `L` in cm;
- `Lmat`: the length (cm) at which half the fish are mature;
- `a`, `b`: hockey-stick recruitment, `min(a S, b)` million recruits a year for a spawning
  stock biomass of `S` thousand tonnes (LeMaRns's `NS_par`);
- `gear`, `max_catchability`, `eta`, `L50`: the fleet (`:industrial`, `:otter`, `:beam` or
  `:pelagic`) that catches each species, with a logistic selectivity of slope `eta` per cm
  and half-selection length `L50` cm, as in LeMaRns's mixed-fishery data (`NS_mixed_fish`,
  `NS_eta`, `NS_L50`);

and two entries for the community: `diet`, a matrix whose entry `[i, j]` is `true` when
species `i` eats species `j` (`NS_tau`), and `other`, the biomass in g of the food outside
the model (`NS_other`).
"""
const LEMANS_NORTH_SEA = _lemans_north_sea(_LEMANS_DIR)

"""
    lemans_parameters(data = LEMANS_NORTH_SEA; classes = 32, dt = 0.1, pred_mu = -2.25,
                      pred_sigma = 0.5, other = data.other, kappa = 10.0, Nmort = 0.8,
                      prop = 0.75, growth_eff = 0.5, growth_eff_decay = 0.11, eps = 1e-5,
                      selectivity_at = :upper)

The quantities of the LeMans model derived from the species' data, as LeMaRns's
`LeMansParam` derives them with its defaults (hockey-stick recruitment, growth efficiency
falling with length, background mortality `"std_RNM"`). A `NamedTuple` with

- `species`, `names`, `gears` and `dt` (the time step in years);
- `lower`, `upper`, `mid`: the bounds and midpoints (cm) of the `classes` length classes,
  of equal width from 0 to 1.01 times the largest `Linf`;
- `fish`: one `NamedTuple` for each species (`p.fish.cod`), with its
  [`ContinuousStage`](@ref) `stage`, whose cells are its length classes up to the one that
  holds its `Linf`, and for each of those classes:
  - `weight` (g): at the midpoint; in the last class, at the midpoint between its lower
    bound and `Linf`;
  - `leaving`: the share of the fish that grow into the next class in a step, `dt` over the
    time von Bertalanffy growth takes to cross the class (zero in the last class);
  - `ration` (g): the food a fish needs in a step, its growth in weight over the step
    divided by a growth efficiency `growth_eff (1 - (w/W∞)^growth_eff_decay)`;
  - `mature`: the share mature, logistic in length with slope `kappa` (one in the last
    class, zero below `eps`); `spawners` is `mature × weight`;
  - `natural`: background mortality per step, `Nmort dt` in the classes from
    `floor(prop × last)` on, zero below;
  - `catchability`: fishing mortality per year at a fishing effort of one, logistic in
    length (see [`LEMANS_NORTH_SEA`](@ref)) evaluated at the upper bound of the class,
    divided by its value at the largest length class of the model, zero below 10⁻⁵. With
    `selectivity_at = :mid` it is evaluated at the midpoint instead. LeMaRns 0.1.2 uses the
    upper bound: `LeMansParam` passes the midpoints and bounds to `calc_Q` by position, in
    an order that gives its argument `mid` the upper bounds;
  - `suitability`: a matrix with one row for each predator class of the community (the
    rows of `predators`) and one column for each class of this species, the share of the
    diet of that predator this class would make up if every prey class were equally
    abundant: `exp(-(log10(w_prey / w_predator) - pred_mu)² / (2 pred_sigma))` where the
    diet matrix allows it, scaled to sum to one over all prey of the predator. Predators
    in the smallest length class eat nothing;
  and its `gear`, recruitment parameters `a` and `b`, and `Linf`;
- `predators`: the predator classes of the community, `(species, class)` pairs, species by
  species and class by class; the pools of [`lemans_fish`](@ref) have one value for each;
- `other` (g), the food outside the model, and `eps`, added to every mortality rate as in
  LeMaRns.

`pred_sigma` divides the squared log-ratio directly, as in LeMaRns, so it is a variance
rather than a standard deviation. LeMaRns computes the growth efficiency of the last length
class with a fixed 0.5 in place of `growth_eff`; the two agree at the default.
"""
function lemans_parameters(data = LEMANS_NORTH_SEA; classes::Integer = 32, dt::Real = 0.1, pred_mu::Real = -2.25,
                           pred_sigma::Real = 0.5, other::Real = data.other, kappa::Real = 10.0, Nmort::Real = 0.8,
                           prop::Real = 0.75, growth_eff::Real = 0.5, growth_eff_decay::Real = 0.11, eps::Real = 1e-5,
                           selectivity_at::Symbol = :upper)
    dt > 0 || throw(ArgumentError("the time step dt must be positive"))
    selectivity_at in (:upper, :mid) || throw(ArgumentError("selectivity_at must be :upper or :mid, got :$selectivity_at"))
    pred_sigma > 0 || throw(ArgumentError("pred_sigma must be positive"))
    nsp = length(data.species)
    # the length classes, computed as LeMaRns computes them (R's seq)
    largest = maximum(data.Linf) * 1.01
    width = largest / classes
    lower = [(j - 1) * width for j in 1:classes]
    upper = [width + (j - 1) * width for j in 1:classes]
    mid = lower .+ (upper .- lower) ./ 2
    W(i, L) = data.W_a[i] * L^data.W_b[i]

    function one(i)
        Linf, k = data.Linf[i], data.k[i]
        m = findlast(<(Linf), lower)
        crossing = [j < m ? log((Linf - lower[j]) / (Linf - upper[j])) / k : 0.0 for j in 1:m]
        leaving = [t > 0 ? dt / t : 0.0 for t in crossing]
        # the last class holds fish between its lower bound and Linf: its weight is taken half way
        last_mid = lower[m] + (Linf - lower[m]) / 2
        length_of = [j < m ? mid[j] : last_mid for j in 1:m]
        weight = [W(i, L) for L in length_of]
        grown = [L + (Linf - L) * (1 - exp(-k * dt)) for L in length_of]
        efficiency = [(1 - (weight[j] / W(i, Linf))^growth_eff_decay) * growth_eff for j in 1:m]
        ration = [(W(i, grown[j]) - weight[j]) / efficiency[j] for j in 1:m]
        mature = [j == m ? 1.0 : 1 / (1 + exp(-kappa * (mid[j] - data.Lmat[i]))) for j in 1:m]
        mature[mature .< eps] .= 0.0
        first_dying = floor(Int, m * prop)
        natural = [j >= first_dying ? Nmort * dt : 0.0 for j in 1:m]
        selected = [1 / (1 + exp(-data.eta[i] * (L - data.L50[i]))) for L in (selectivity_at === :upper ? upper : mid)]
        selected ./= maximum(selected)
        selected[selected .< 1e-5] .= 0.0
        catchability = data.max_catchability[i] .* selected[1:m]
        return (name = data.names[i], stage = ContinuousStage(:length, 0.0, m * width, m), classes = m,
                Linf = Linf, weight = weight, leaving = leaving, ration = ration, mature = mature,
                spawners = mature .* weight, natural = natural, catchability = catchability,
                gear = data.gear[i], a = data.a[i], b = data.b[i])
    end
    fish = [one(i) for i in 1:nsp]
    predators = [(data.species[i], j) for i in 1:nsp for j in 1:fish[i].classes]

    # suitability of each prey class for each predator class
    preference(w_predator, w_prey) = exp(-(log10(w_prey / w_predator) - pred_mu)^2 / (2 * pred_sigma))
    columns = [sum(f.classes for f in fish[1:(s - 1)]; init = 0) .+ (1:fish[s].classes) for s in 1:nsp]
    S = zeros(length(predators), sum(f.classes for f in fish))
    for (r, (sp, j)) in enumerate(predators)
        i = findfirst(==(sp), data.species)
        for s in 1:nsp
            data.diet[i, s] || continue
            S[r, columns[s]] .= [preference(fish[i].weight[j], w) for w in fish[s].weight]
        end
        total = sum(S[r, :])
        total > 0 && (S[r, :] ./= total)
        j == 1 && (S[r, :] .= 0.0)     # the smallest predators eat nothing
    end
    fish = [merge(fish[s], (suitability = S[:, columns[s]],)) for s in 1:nsp]
    return (species = data.species, names = data.names, gears = unique(data.gear), dt = Float64(dt),
            lower = lower, upper = upper, mid = mid,
            fish = NamedTuple{Tuple(data.species)}(Tuple(fish)), predators = predators,
            other = Float64(other), eps = Float64(eps))
end

_lemans_fish(p, sp::Symbol) = haskey(p.fish, sp) ? p.fish[sp] :
    throw(ArgumentError("unknown species :$sp; the species are $(join(p.species, ", "))"))

# hockey-stick recruitment: recruits (numbers) for a spawning stock biomass in g, as LeMaRns's rec_hockey
_lemans_recruits(f, ssb) = min(f.a * ssb / 1e9, f.b) * 1e6

"""
    lemans_fish(p, species)

One species of the LeMans model with parameters `p` (see [`lemans_parameters`](@ref)), as an
[`OpenCommunity`](@ref) that meets the others through two pools with a value for each
predator class of the community (`p.predators`):

- `ration`: the food in g that the predators of each class need in a step, `ration × N`;
- `prey`: the prey biomass in g that the predators of each class find suitable,
  `Σ suitability × weight × N` over the length classes of every species.

The species adds its own fish to both and reads both. Its stages are its length classes,
the cells of `p.fish[species].stage` (`length1`, `length2`, ...), and it has three
transitions:

- `staying`: the fish of a class that survive the step and stay in their class;
- `growing`: those that survive and grow into the next class, a share `leaving` of each;
- `spawning`: the recruits, which join the smallest class at the first step of a year,
  then die and grow as the fish of that class do in the same step. Marked as births.

The share that survives a step is `exp(-Z)`, with `Z` the sum of the background mortality
`natural`, the fishing mortality `effort × catchability × dt`, the predation mortality

    M2 = suitabilityᵀ (pool.ration ./ (pool.prey .+ p.other))

and `p.eps`. Recruitment is `min(a S / 10⁹, b) × 10⁶` fish for a spawning stock biomass of
`S` g (`Σ spawners × N`), shared among the spawners in proportion to their weight. The
recruits are prey in the step they join, so they count in the pool `prey` of that step.

The environment of a step is a `NamedTuple` with `recruiting` (whether this is the first
step of a year) and `effort` (a `NamedTuple` of fishing effort by fleet, missing fleets
not fishing); [`lemans_environment`](@ref) makes the sequence for a run.
"""
function lemans_fish(p, sp::Symbol)
    f = _lemans_fish(p, sp)
    fish = f.stage
    cells = stage_names(fish)
    m = length(cells)
    life = ValuedProjectionNet(cells,
        :staying => [(c => c) => 1 - g for (c, g) in zip(cells, f.leaving)],
        :growing => [(cells[j] => cells[j + 1]) => f.leaving[j] for j in 1:(m - 1)],
        # each g of spawner gives a recruit at the base value; the recruits then stay or grow as the smallest fish do
        :spawning => vcat([(c => cells[1]) => w * (1 - f.leaving[1]) for (c, w) in zip(cells, f.spawners) if w > 0],
                          [(c => cells[2]) => w * f.leaving[1] for (c, w) in zip(cells, f.spawners) if w > 0]))

    # predation mortality per step of each length class, from the two pools
    predation(pool) = f.suitability' * (pool.ration ./ (pool.prey .+ p.other))
    fishing(env) = get(env.effort, f.gear, 0.0) * p.dt .* f.catchability
    surviving(ctx) = exp.(-(f.natural .+ fishing(ctx.env) .+ predation(ctx.pool) .+ p.eps))
    spawning_stock(n) = dot(f.spawners, n[fish])
    recruits(n, env) = env.recruiting ? _lemans_recruits(f, spawning_stock(n)) : 0.0
    # recruits per g of spawner, which has a limit as the spawning stock falls to zero
    function per_gram(ctx)
        ctx.env.recruiting || return 0.0
        return min(f.a * 1e-3, f.b * 1e6 / spawning_stock(ctx.n))
    end

    net = births(respond(life,
        :staying => ctx -> by_trait(fish, surviving(ctx)),
        :growing => ctx -> by_trait(fish, surviving(ctx)),
        :spawning => ctx -> per_gram(ctx) * first(surviving(ctx))), :spawning)

    # the species' rows among the predator classes of the community
    rows = [findfirst(==((sp, j)), p.predators) for j in 1:m]
    needs = zeros(length(p.predators), m)
    for j in 1:m
        needs[rows[j], j] = f.ration[j]
    end
    # prey on offer: the biomass of each class weighted by each predator class's suitability; this step's
    # recruits join the smallest class before predation
    offer = f.suitability .* transpose(f.weight)
    prey(n, env) = offer * n[fish] .+ recruits(n, env) .* offer[:, 1]
    return species(sp, net; occupies = [:ration => needs, :prey => prey])
end

"""
    lemans_community(p = lemans_parameters(); species = p.species)

The species of the LeMans model with parameters `p` composed through their predation
pools (see [`lemans_fish`](@ref)), or those named in `species`. A species left out is
neither predator nor prey, but the predators keep the suitabilities of the full
community, so they eat more of the food outside the model.
"""
lemans_community(p = lemans_parameters(); species = p.species) =
    community((lemans_fish(p, sp) for sp in species)...)

"""
    lemans_state(p = lemans_parameters(); intercept = 1e10, slope = -5, species = p.species)

The starting state of LeMaRns (its `get_N0`), as `species => numbers` pairs, one number
for each length class, for [`project`](@ref) or `state_vector`: the whole community holds
`intercept × L^slope` fish in the class with midpoint `L`, shared among the species in
proportion to their share of that class when each species, spread over its own classes as
the community is, has a biomass of one.
"""
function lemans_state(p = lemans_parameters(); intercept::Real = 1e10, slope::Real = -5, species = p.species)
    total = intercept .* p.mid .^ slope
    # each species spread over its own classes in proportion to `total`, then scaled to a biomass of one
    shape = Dict(sp => total[1:p.fish[sp].classes] ./ sum(total[1:p.fish[sp].classes]) for sp in p.species)
    per_biomass = Dict(sp => shape[sp] ./ dot(shape[sp], p.fish[sp].weight) for sp in p.species)
    in_class(j) = sum(j <= length(v) ? v[j] : 0.0 for v in values(per_biomass))
    return [sp => [total[j] * per_biomass[sp][j] / in_class(j) for j in eachindex(per_biomass[sp])] for sp in species]
end

"""
    lemans_environment(years; effort = (;), dt = 0.1)

The environments of a run of `years` years of the LeMans model in steps of `dt` years
(`floor(years / dt)` steps), as [`lemans_fish`](@ref) reads them: for step `k`,
`(recruiting, effort, year)`. Step `k` is in year `ceil(k dt)`, and recruits when that
differs from `ceil((k - 1) dt)`, as in LeMaRns's `run_LeMans`. `effort` is the fishing
effort by fleet, `(otter = 1.0, beam = 0.5)` say, the same in every year; or a vector with
one such `NamedTuple` for each year; or a function from the year to one.
"""
function lemans_environment(years::Real; effort = NamedTuple(), dt::Real = 0.1)
    steps = floor(Int, years / dt)
    by_year(y) = effort isa NamedTuple ? effort : effort isa AbstractVector ? effort[y] : effort(y)
    return [(recruiting = ceil(k * dt) != ceil((k - 1) * dt), effort = by_year(ceil(Int, k * dt)), year = ceil(Int, k * dt))
            for k in 1:steps]
end

"""
    lemans_biomass(p, sol; spawning = false)

The biomass in g of each species through a run of the LeMans model (a
[`CommunityTrajectory`](@ref)): a matrix with one row per time (as `sol.u`) and one column
per species of the community, in the order of `species_names(sol.community)`. With
`spawning = true`, the spawning stock biomass.
"""
function lemans_biomass(p, sol::CommunityTrajectory; spawning::Bool = false)
    kinds = species_names(sol.community)
    w(sp) = spawning ? p.fish[sp].spawners : p.fish[sp].weight
    return reduce(hcat, [sol[sp] * w(sp) for sp in kinds])
end

"""
    lemans_predation(p, c, n, env)

The predation mortality per step of each length class of each species of the LeMans
community `c` in state `n` and environment `env`, as a `Dict` from species to a vector:
the `M2` of [`lemans_fish`](@ref), read from the pools of that state. In the first step of a
year the pools include the year's recruits.
"""
function lemans_predation(p, c::OpenCommunity, n::AbstractVector, env)
    pools = pool_totals(c, n, env)
    pressure = pools.ration ./ (pools.prey .+ p.other)
    return Dict(sp => p.fish[sp].suitability' * pressure for sp in species_names(c))
end
