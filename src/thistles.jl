# A native thistle (an integral projection model), an invasive thistle (a matrix model) and the
# flower-head weevil they share. Each species is a net; they meet only through three pools:
#
#   heads    flower heads on offer to egg-laying weevils, in Platte thistle heads
#   rearing  the number of weevils those heads can rear
#   eggs     the eggs the weevils lay
#
# Numbers are individuals on 720 m², the area of the plots in which the native thistle was counted.

_logistic(u) = 1 / (1 + exp(-u))
_normal_density(y, mean, sd) = exp(-0.5 * ((y - mean) / sd)^2) / (sd * sqrt(2pi))

# ---------------------------------------------------------------- Platte thistle

"""
    PLATTE_THISTLE_PARAMETERS

Default parameters of [`platte_thistle`](@ref): the integral projection model of Platte
thistle (*Cirsium canescens*) of Rose, Louda & Rees (2005, their Tables 1 and 2). Size is
the natural logarithm of root crown diameter in mm, and the time step is a year. Sources
of the other values are in the comments of `src/thistles.jl`.
"""
const PLATTE_THISTLE_PARAMETERS = (
    size_range = (-1.0, 4.5),     # wide enough to hold the population (the paper used 0.9 x the smallest to 1.1 x the largest size seen, not given)
    growth = (0.83, 0.69),        # size next year = a + b x (Table 2)
    growth_variance = 0.19,
    survival = (-0.62, 0.85),     # logit of survival
    flowering = (-10.22, 4.25),   # logit of flowering; a plant flowers once and dies
    seeds = (0.37, 2.02),         # seeds of a flowering plant without insects = exp(a + b x) (Eq. 4)
    native_insects = 0.65,        # share of seeds lost to native insects (Eq. 4 at mean damage, "a reduction of 65%")
    seedling_size = (0.75, 0.17), # mean and variance
    recruitment = 0.67,           # recruits = (seeds produced)^0.67, on the plots (Eq. 5)
    establishment = 0.067,        # mean chance that a seed establishes; used only without seed limitation
    greatest_establishment = 0.47, # the chance of establishment cannot exceed this: seeds sown at low density at the same
                                  # site gave 0.47 seedlings per seed, and 0.29 at a second site (Rand et al. 2020, Table 2)
    weevil_damage = 0.02,         # seeds x exp(-0.02 x weevil eggs on the plant) (Eq. 4)
    egg_aggregation = 0.32,       # k of the negative binomial distribution of eggs among plants of one size (Table 1)
    # flower heads initiated by a flowering plant = a x diameter^b. b is the size dependence of the egg load (Table 1);
    # a makes the average flowering plant initiate 10.4 heads (West & Louda 2018, plants under ambient herbivory)
    heads = (0.119, 1.71),
    # weevils that one initiated head can rear: 10 to 13 in a head that flowers (Louda et al. 2011), and 3.8 of the
    # 10.4 heads flower (West & Louda 2018)
    weevils_per_head = 4.2)

"""
    platte_thistle_stage(; parameters = PLATTE_THISTLE_PARAMETERS, cells = 100)

The [`ContinuousStage`](@ref) of [`platte_thistle`](@ref): plants by size, the natural
logarithm of root crown diameter in mm. `meshpoints` of it are the sizes of the cells.
"""
platte_thistle_stage(; parameters = PLATTE_THISTLE_PARAMETERS, cells::Integer = 100) =
    ContinuousStage(:plant, parameters.size_range[1], parameters.size_range[2], cells)

"""
    platte_thistle(; parameters = PLATTE_THISTLE_PARAMETERS, cells = 100, native_insects = true,
                     seed_limited = true, aggregated = true, eggs = nothing, name = :platte_thistle)

Platte thistle (*Cirsium canescens*), a monocarpic thistle of Nebraska sand prairie, as the
integral projection model of Rose, Louda & Rees (2005). Plants are structured by size on
`cells` mesh cells (see [`platte_thistle_stage`](@ref)); numbers are plants on 720 m².
Transitions are `growing` (survive, do not flower, change size) and `seeding` (survive,
flower, set seed, and the seedlings that result).

- `native_insects`: native flower-head insects destroy a fixed share of the seeds.
- `seed_limited`: recruits are a power of the seeds produced, as fitted, so a seed is less
  likely to establish in a large crop; the chance is capped at the value measured for
  seeds sown at low density. With `false` every seed has the same chance, which gives a
  growth rate.
- The flower-head weevil reduces seed set according to the eggs on the plant. With
  `eggs = nothing` these follow from the pools: the plant offers its flower heads to the
  pool `heads` and receives a share of the pool `eggs` in proportion, so that larger
  plants carry more eggs. It also adds the weevils its heads can rear to the pool
  `rearing`. With `eggs` a function of size, that many eggs are imposed instead.
- `aggregated`: eggs are spread unevenly among plants of the same size (a negative
  binomial distribution), as observed; with `false` every plant of a size carries the mean.

Alone, or with nothing laying eggs, the thistle has no weevil damage. Compose it with
[`flowerhead_weevil`](@ref) and [`musk_thistle`](@ref).
"""
function platte_thistle(; parameters = PLATTE_THISTLE_PARAMETERS, cells::Integer = 100, native_insects::Bool = true,
                        seed_limited::Bool = true, aggregated::Bool = true, eggs = nothing, name::Symbol = :platte_thistle)
    p = parameters
    plant = platte_thistle_stage(; parameters, cells)
    z = meshpoints(plant)
    survival(x) = _logistic(p.survival[1] + p.survival[2] * x)
    flowering(x) = _logistic(p.flowering[1] + p.flowering[2] * x)
    growth(y, x) = _normal_density(y, p.growth[1] + p.growth[2] * x, sqrt(p.growth_variance))
    seeds(x) = exp(p.seeds[1] + p.seeds[2] * x)
    seedling(y) = _normal_density(y, p.seedling_size[1], sqrt(p.seedling_size[2]))
    heads(x) = p.heads[1] * exp(p.heads[2] * x)

    # the life cycle without insects, every seed establishing
    life = kernel_net(plant,
        :growing => (y, x) -> survival(x) * (1 - flowering(x)) * growth(y, x),
        :seeding => (y, x) -> survival(x) * flowering(x) * seeds(x) * seedling(y))

    spared = native_insects ? 1 - p.native_insects : 1.0
    k, per_egg = p.egg_aggregation, 1 - exp(-p.weevil_damage)
    # share of seeds escaping the weevil on plants that carry `mean` eggs on average
    escape(mean) = aggregated ? (k / (k + mean * per_egg))^k : exp(-p.weevil_damage * mean)
    function egg_load(ctx)
        eggs === nothing || return eggs
        per_head = ctx.pool.heads > 0 ? ctx.pool.eggs / ctx.pool.heads : 0.0
        return x -> per_head * heads(x)
    end
    potential = [survival(x) * flowering(x) * seeds(x) for x in z]       # seeds per plant, before insects
    function seeding(ctx)
        load = egg_load(ctx)
        kept = [spared * escape(load(x)) for x in z]
        establishing = if seed_limited
            crop = sum(potential[i] * kept[i] * ctx.n[i] for i in eachindex(z))
            crop > 0 ? min(crop^(p.recruitment - 1), p.greatest_establishment) : p.greatest_establishment
        else
            p.establishment
        end
        return by_trait(plant, kept .* establishing)
    end
    net = births(respond(life, :seeding => seeding), :seeding)

    offered(x) = survival(x) * flowering(x) * heads(x)
    occupies = Pair{Symbol, Any}[:heads => trait_weights(plant, offered),
                                 :rearing => trait_weights(plant, x -> p.weevils_per_head * offered(x))]
    return species(name, net; occupies, reads = eggs === nothing ? [:eggs] : Symbol[])
end

# ---------------------------------------------------------------- musk thistle

"""
    MUSK_THISTLE_PARAMETERS

Default parameters of [`musk_thistle`](@ref). The vital rates are those of the matrix model
of Jongejans et al. (2008) for a natural population in Kansas (Lee & Hamrick 1983), as
tabulated by Zhang, Jongejans & Shea (2011, Table S2); its growth rate is 1.748. Sources of
the other values are in the comments of `src/thistles.jl`.
"""
const MUSK_THISTLE_PARAMETERS = (
    seed_survival = 0.2597,                # a seed stays alive in the seed bank
    emergence = 0.03,                      # a seed in the seed bank becomes a rosette
    survival = (0.2182, 0.5633, 0.7474),   # small, medium and large rosettes
    bolting = (0.0180, 0.3596, 1.0),       # a surviving rosette bolts, flowers and dies
    seeds = (100.0, 195.0, 1473.0),        # potential seeds of a bolting plant
    to_seed_bank = 0.2333,                 # a new seed enters the seed bank
    establishment = 0.03,                  # a new seed becomes a rosette within the year
    seedling_size = (0.8373, 0.1627, 0.0), # a new rosette is small, medium or large
    growth_small = (0.8028, 0.1268),       # a surviving small rosette that does not bolt becomes medium, large
    growth_medium = 0.3824,                # a surviving medium rosette that does not bolt becomes large
    seeds_per_head = 448.0,                # a head of 25 mm (Popay et al. 1984); 447 in Surles et al. (1974, abstract)
    seeds_per_larva = 29.6,                # seeds destroyed by one weevil larva (Popay et al. 1984)
    weevils_per_head = 25.0,               # larvae a head can hold: means of 18 to 28 at 27 to 134 eggs per head (Rees 1978)
    attraction = 2.3,                      # eggs laid on a musk thistle head per egg on a Platte thistle head of the same
                                           # diameter (2.8 against about 1.2, read from Fig. 1 of Arnett & Louda 2002)
    exposed = 0.5,                         # share of heads open while weevils lay (proposed: 86% of first heads were infested,
                                           # 24% of second and 3% of third; McNeill & Fletcher 2005)
    sites = 6000.0)                        # rosettes the open ground can hold. Proposed; it sets how abundant the weed is, and
                                           # 6000 gives 0.4 bolting plants per m2 with the weevil present, near the median of
                                           # 0.46 in Nebraska patches (Russell et al. 2007, Appendix A)

"""
    musk_thistle(; parameters = MUSK_THISTLE_PARAMETERS, weevil_survival = WEEVIL_PARAMETERS.survival,
                   name = :musk_thistle)

Musk thistle (*Carduus nutans*), an invasive monocarpic thistle, as the matrix model of
Jongejans et al. (2008) written from its vital rates: a seed bank and small, medium and
large rosettes, counted at the end of autumn. Transitions are `dormancy`, `emerging` (from
the seed bank), `growing` (rosettes that survive without bolting), and `seed_rain` and
`establishing` (the seeds of plants that bolt, entering the seed bank or becoming rosettes).

Two responses are added to the published model:

- New rosettes need open ground: the seedlings of a year compete for the `sites` that hold
  no rosette, and of `seedlings` contending for `free` sites the share
  `free / (free + seedlings)` establishes (a response on `emerging` and `establishing`).
  This bounds a population that would otherwise grow 1.75-fold a year, and `sites` sets
  how abundant the weed is.
- The flower-head weevil destroys seed. The published model has a multiplier for the share
  of seeds escaping floral herbivory; here it follows from the pools. Bolting plants offer
  their early heads to the pool `heads` (each worth `attraction` Platte thistle heads) and
  the weevils these can rear to the pool `rearing`; each weevil larva that develops
  destroys `seeds_per_larva` seeds. `weevil_survival` is the chance that an egg yields a
  developing larva, repeated from [`flowerhead_weevil`](@ref).
"""
function musk_thistle(; parameters = MUSK_THISTLE_PARAMETERS, weevil_survival::Real = WEEVIL_PARAMETERS.survival,
                      name::Symbol = :musk_thistle)
    p = parameters
    rosettes = [:small, :medium, :large]
    bolters = [p.survival[j] * p.bolting[j] for j in 1:3]          # rosettes that survive the winter and bolt
    staying = [p.survival[j] * (1 - p.bolting[j]) for j in 1:3]    # rosettes that survive and do not
    to_medium, to_large = p.growth_small
    sized(value) = [(rosettes[k], value * p.seedling_size[k]) for k in 1:3 if p.seedling_size[k] > 0]
    life = ValuedProjectionNet([:seed_bank, :small, :medium, :large],
        :dormancy => [(:seed_bank => :seed_bank) => p.seed_survival],
        :emerging => [(:seed_bank => to) => v for (to, v) in sized(p.emergence)],
        :growing => [(:small => :small) => staying[1] * (1 - to_medium - to_large),
                     (:small => :medium) => staying[1] * to_medium, (:small => :large) => staying[1] * to_large,
                     (:medium => :medium) => staying[2] * (1 - p.growth_medium),
                     (:medium => :large) => staying[2] * p.growth_medium,
                     (:large => :large) => staying[3]],
        :seed_rain => [(rosettes[j] => :seed_bank) => bolters[j] * p.seeds[j] * p.to_seed_bank for j in 1:3],
        :establishing => [(rosettes[j] => to) => v for j in 1:3 for (to, v) in sized(bolters[j] * p.seeds[j] * p.establishment)])

    # share of seeds escaping the weevil: larvae develop in the exposed heads only, and die young when crowded
    function escaping(ctx)
        (ctx.pool.heads > 0 && ctx.pool.rearing > 0) || return 1.0
        eggs_per_head = p.attraction * ctx.pool.eggs / ctx.pool.heads
        larvae = weevil_survival * eggs_per_head / (1 + weevil_survival * ctx.pool.eggs / ctx.pool.rearing)
        return 1 - p.exposed * min(1.0, larvae * p.seeds_per_larva / p.seeds_per_head)
    end
    # seedlings compete for the sites that hold no rosette: of those that would establish, the share
    # free / (free + seedlings) gets a site, so the rosettes never outnumber the sites
    function open_ground(ctx)
        free = max(0.0, p.sites - ctx.n[:small] - ctx.n[:medium] - ctx.n[:large])
        seedlings = p.emergence * ctx.n[:seed_bank] +
                    escaping(ctx) * p.establishment * sum(bolters[j] * p.seeds[j] * ctx.n[rosettes[j]] for j in 1:3)
        return seedlings > 0 ? free / (free + seedlings) : 1.0
    end
    net = respond(life, :emerging => open_ground, :establishing => open_ground,
                  :seed_rain => escaping, :establishing => escaping)
    net = births(net, :seed_rain, :establishing)

    # early heads of the plants that bolt, in heads of `seeds_per_head` seeds
    early = [bolters[j] * p.exposed * p.seeds[j] / p.seeds_per_head for j in 1:3]
    occupies = Pair{Symbol, Any}[:heads => [rosettes[j] => p.attraction * early[j] for j in 1:3],
                                 :rearing => [rosettes[j] => p.weevils_per_head * early[j] for j in 1:3]]
    return species(name, net; occupies, reads = [:eggs])
end

# ---------------------------------------------------------------- the weevil

"""
    WEEVIL_PARAMETERS

Default parameters of [`flowerhead_weevil`](@ref); sources are in the comments of
`src/thistles.jl`.
"""
const WEEVIL_PARAMETERS = (
    eggs = 50.0,          # per adult: 100.5 per female in the laboratory (Surles 1972), taking half the adults to be females
    survival = 0.32,      # egg to new adult in a head that is not crowded (field, musk thistle: 1128 eggs, 360 adults; Surles 1972)
    # new adult to the next spring. Cages gave 0.03 to 0.29 (Smith 1982) and 0.48 (Kok 1976, abstract). Calibrated within
    # that range so that the weevils reared in a Platte thistle head that flowers settle at 5, the number found at Arapaho
    # Prairie in 1995 and 1996 (Louda 1998); 3.3 were found over 108 sites (Rand & Louda 2006)
    overwintering = 0.11)

"""
    flowerhead_weevil(; parameters = WEEVIL_PARAMETERS, name = :weevil)

The flower-head weevil *Rhinocyllus conicus*, released against musk thistle and now also
on native thistles. It has one generation a year: adults that have overwintered lay eggs
on flower heads in spring, larvae feed in the head, and new adults overwinter. The state
is the number of adults in spring.

The weevil lays its eggs into the pool `eggs`; thistles share them out among the heads they
offer. Larvae compete within heads: the more eggs there are for the weevils that the heads
can rear (the pool `rearing`, supplied by the thistles), the smaller the share that
develops. No thistle is named: any species that offers `heads` and `rearing` is a host.
"""
function flowerhead_weevil(; parameters = WEEVIL_PARAMETERS, name::Symbol = :weevil)
    p = parameters
    life = ValuedProjectionNet([:adult], :breeding => [(:adult => :adult) => p.eggs * p.survival * p.overwintering])
    function crowding(ctx)
        ctx.pool.rearing > 0 || return 0.0
        return 1 / (1 + p.survival * p.eggs * ctx.n[:adult] / ctx.pool.rearing)
    end
    net = births(respond(life, :breeding => crowding), :breeding)
    return species(name, net; occupies = [:eggs => [p.eggs]], reads = [:rearing])
end
