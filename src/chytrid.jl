# Leopard frogs and the chytrid fungus Batrachochytrium dendrobatidis (Bd): the full host-parasite integral
# projection model of Wilber, Ohmer, Altman et al. (2022, Ecology 103:e3759), after Wilber et al. (2016, Methods
# in Ecology and Evolution 7:1182-1194).
#
# One host species at one location. Tadpoles develop for the larval period and metamorphose into susceptible
# post-metamorphic frogs; a susceptible frog in the water is infected by zoospores, and an infected frog is
# described by its Bd load, the natural logarithm of its zoospore equivalents. Each week the load changes, the
# infection may be lost, and an infected frog in the water sheds zoospores into a pool that decays at a rate
# set by the water temperature. Temperature follows a sinusoid through a year of 364 days, and breeding is one
# pulse in the middle of the breeding season.
#
# The equations, the order of operations and every value are those of the authors' Python code (full_model.py
# and ccs_simulation.py, https://github.com/mqwilber/bd_seasonality, Zenodo doi:10.5281/zenodo.6407088), which
# carries no licence and is not copied here; its values are facts, cited. The fitted values are the medians of
# the MCMC chains of the reduced IPM that ccs_simulation.py loads into the full IPM (density assumption K = 4,
# intensity-dependent loss of infection, an external zoospore reservoir during fitting). Time is in days, a
# step is 7 days, numbers are totals in an arena of `area` m³.

"""
    CHYTRID_SITES

The four locations of Wilber et al. (2022), for [`leopard_frog`](@ref): `PA` (Pennsylvania)
and `VT` (Vermont), with the northern leopard frog *Rana pipiens*, and `LA`
(Louisiana) and `TN` (Tennessee), with the southern leopard frog *R. sphenocephala*. Each
is a `NamedTuple` with

- `location`, `host`, `model` (the fitted model whose medians are used: `:a_temp`, Bd
  growth depending on temperature, or `:both_model`, Bd growth and transmission depending
  on temperature) and `temperature`, the yearly minimum and maximum water temperature
  (°C) of the sinusoid;
- `s0`, the weekly survival of an uninfected post-metamorphic frog, and `sI`, the
  survival of an infected frog relative to it;
- `a`, `a_temp`, `b`, `sigma`: the log load next week is normal with mean
  `a + a_temp (T - 15) + b x` for log load `x` now and standard deviation `sigma`; a new
  infection starts with the load `a + a_temp (T - 15)` and the same deviation;
- `mu_l`, `sigma_l`: an infected frog of log load `x` loses its infection within the week
  with chance `1 - Φ((x - mu_l) / sigma_l)`;
- `beta`, `beta_temp`: a susceptible frog in the water is infected within the week with
  chance `1 - exp(-beta exp(beta_temp (T - 15)) Z / area)` for `Z` zoospores;
- `shedding`: zoospores shed in a week by a frog of log load `x` in the water, `shedding × e^x`;
- `fecundity`: tadpoles per post-metamorphic frog at the breeding pulse;
- `K`: recruitment of metamorphs is multiplied by `exp(-e^K N / area)` for `N`
  post-metamorphic frogs;
- `larval_period` (days), `breeding` (first and last day of the breeding season),
  `aquatic` (whether post-metamorphic frogs live in the water all year), `hibernation`,
  `hibernation_aquatic` and `hibernation_period` (first and last day, which may wrap round
  the year's end).

`a`, `a_temp`, `sigma`, `mu_l` and `beta` (and `beta_temp` for LA) are medians of the
posterior fitted by the authors; the rest are fixed values of their parameter files.
Sources of the values are in the comments of `src/chytrid.jl`.
"""
const CHYTRID_SITES = (
    # the medians, computed as ccs_simulation.load_reduced_parameters computes them from the authors' pickled MCMC
    # results (results/pickled_results/*_allsites_ext_beta_prior_K=4_loss_load_omega=1.pkl) by
    # test/reference/python/run_bd_reference.py, and written to full precision; the fixed values are those of the
    # same files (base_params) and of code/model_params/*_params/spp_params_*.yml. The temperature ranges are those of
    # ccs_simulation.py (PA, LA, TN) and of dod_dict in model_analysis.ipynb (VT).
    PA = (location = "Pennsylvania", host = "northern leopard frog, Rana pipiens", model = :a_temp,
          temperature = (4.0, 25.0), s0 = 0.995, sI = 1.0, a = 3.841640984857781, a_temp = -0.1991081815756661,
          b = 0.5, sigma = 3.2578534407782054, mu_l = 3.7563388668845272, sigma_l = 1.15,
          beta = 4.484275304986494e-6, beta_temp = 0.0, shedding = 1.0, fecundity = 600.0, K = 4.0,
          larval_period = 90.0, breeding = (75, 150), aquatic = false, hibernation = true,
          hibernation_aquatic = true, hibernation_period = (315, 75)),
    LA = (location = "Louisiana", host = "southern leopard frog, Rana sphenocephala", model = :both_model,
          temperature = (4.0, 30.0), s0 = 0.995, sI = 1.0, a = 4.735916942098988, a_temp = -0.18545828052298247,
          b = 0.5, sigma = 3.1932051958809797, mu_l = 1.1366182325776504, sigma_l = 1.15,
          beta = 4.599782063576056e-5, beta_temp = -1.1619552789892564, shedding = 1.0, fecundity = 600.0, K = 4.0,
          larval_period = 57.0, breeding = (30, 120), aquatic = false, hibernation = false,
          hibernation_aquatic = true, hibernation_period = (315, 105)),
    TN = (location = "Tennessee", host = "southern leopard frog, Rana sphenocephala", model = :a_temp,
          temperature = (4.0, 27.0), s0 = 0.995, sI = 1.0, a = 4.338111459822141, a_temp = -0.12211014135261514,
          b = 0.5, sigma = 2.2745905992529565, mu_l = 3.276549945988139, sigma_l = 1.15,
          beta = 2.5608012692579595e-5, beta_temp = 0.0, shedding = 1.0, fecundity = 600.0, K = 4.0,
          larval_period = 90.0, breeding = (30, 150), aquatic = false, hibernation = false,
          hibernation_aquatic = true, hibernation_period = (315, 105)),
    VT = (location = "Vermont", host = "northern leopard frog, Rana pipiens", model = :a_temp,
          temperature = (4.0, 25.0), s0 = 0.995, sI = 1.0, a = 4.524309246640505, a_temp = -0.1716974641899735,
          b = 0.5, sigma = 2.8620117581975477, mu_l = 1.232608088353262, sigma_l = 1.15,
          beta = 3.3442725913792966e-8, beta_temp = 0.0, shedding = 1.0, fecundity = 600.0, K = 4.0,
          larval_period = 90.0, breeding = (105, 150), aquatic = false, hibernation = true,
          hibernation_aquatic = true, hibernation_period = (315, 105)))

"""
    CHYTRID_ZOOSPORE_SURVIVAL

The weekly survival of Bd zoospores in water at temperature `T` (°C) used by
[`chytrid_zoospores`](@ref): `maxsurv / (1 + exp((T - ld50) / s))`, a `NamedTuple` with
`ld50`, `s` and `maxsurv`. The authors' code fits this curve by least squares (`scipy`'s
`fmin`, from `(20, 1)`) to the weekly survival `exp(-7 × 24 × rate)` at the four
temperatures 4, 14.5, 23 and 28 °C, with the death rates per hour 0.0027, 0.0035, 0.0160 and
0.0410 that it attributes to Woodhams et al. (2008); `maxsurv` is the largest of the four.
The values are those of that fit, read back from the authors' function.
"""
const CHYTRID_ZOOSPORE_SURVIVAL = (ld50 = 18.547852060025306, s = 2.082483488695515, maxsurv = 0.63533681715262)

const _CHYTRID_YEAR = 364        # days in the authors' year (the default year_length of full_model.Community)
const _CHYTRID_STEP = 7          # days in a step
const _CHYTRID_BASE_TEMPERATURE = 15.0   # temperature (°C) at which the temperature effects are zero

_chytrid_site(site::Symbol) = haskey(CHYTRID_SITES, site) ? CHYTRID_SITES[site] :
    throw(ArgumentError("unknown site :$site; the sites are $(join(keys(CHYTRID_SITES), ", "))"))
_chytrid_site(site::NamedTuple) = site

# Whether `day` of the year lies in the window (first, last), which wraps round the end of the year when first > last.
_chytrid_within(day, (first, last)) = first <= last ? first <= day <= last : !(last <= day <= first)

# Whether post-metamorphic frogs are in the water on `day`, where they meet zoospores and shed them: all year if
# aquatic, otherwise while breeding and, if they hibernate in the water, while hibernating.
function _chytrid_in_water(s, day)
    d = mod(day, _CHYTRID_YEAR)
    s.aquatic && return true
    _chytrid_within(d, s.breeding) && return true
    return s.hibernation && s.hibernation_aquatic && _chytrid_within(d, s.hibernation_period)
end

# Whether the breeding pulse falls in the step that starts on `day`: the middle of the breeding season, as the
# authors' get_repro_time computes it (with its 365-day year for a season that wraps), in (day - 7, day] of the year.
function _chytrid_breeding_now(s, day)
    first, last = s.breeding
    length_of = last >= first ? last - first : 365 - (first - last)
    pulse = first + length_of / 2
    pulse > 365 && (pulse -= 365)
    d = mod(day, _CHYTRID_YEAR)
    return d - _CHYTRID_STEP < pulse <= d
end

"""
    chytrid_temperature(site, day)

The water temperature (°C) on `day` (counted from the start of a run, which is the
coldest day) at `site` (a key of [`CHYTRID_SITES`](@ref) or a site `NamedTuple`): the
authors' sinusoid, `min + (max - min) (1 - cos(2π day / 364)) / 2`.
"""
function chytrid_temperature(site, day::Real)
    s = _chytrid_site(site)
    low, high = s.temperature
    return (high - low) / 2 * (1 - cos(2pi * day / _CHYTRID_YEAR)) + low
end

"""
    chytrid_environment(site, steps; start = 0)

The environments of `steps` weekly steps at `site`, as [`leopard_frog`](@ref) and
[`chytrid_zoospores`](@ref) read them: for step `k`, `(day, temperature)`, with
`day = start + 7 (k - 1)` the day the step starts and `temperature` that of
[`chytrid_temperature`](@ref) on that day. Day 0 is the coldest day of the year.
"""
chytrid_environment(site, steps::Integer; start::Integer = 0) =
    [(day = start + _CHYTRID_STEP * (k - 1), temperature = chytrid_temperature(site, start + _CHYTRID_STEP * (k - 1)))
     for k in 1:steps]

"""
    chytrid_load_stage(; cells = 300, range = (-20.0, 25.0))

The [`ContinuousStage`](@ref) of infected frogs in [`leopard_frog`](@ref): the natural
logarithm of the Bd load (zoospore equivalents) from -20 to 25 in 300 cells, as in the
authors' simulations. Its cells are `load1`, `load2`, ...
"""
chytrid_load_stage(; cells::Integer = 300, range = (-20.0, 25.0)) = ContinuousStage(:load, range[1], range[2], cells)

# The chance that a normal variable of mean `m` and deviation `sd` falls in [lo, hi], taken on the side away from
# the upper tail so that it keeps its precision there.
function _chytrid_cell_chance(lo, hi, m, sd)
    u, l = (hi - m) / sd, (lo - m) / sd
    return l > 0 ? _normal_cdf(-l) - _normal_cdf(-u) : _normal_cdf(u) - _normal_cdf(l)
end

"""
    leopard_frog(site; area = 500, cells = 300, name = :frog, crowding = :frogs)

The leopard frog at `site` (a key of [`CHYTRID_SITES`](@ref) such as `:PA`, or a site
`NamedTuple`) in the full host-parasite integral projection model of Wilber et al. (2022),
as an [`OpenCommunity`](@ref) of one species. Numbers are totals in an arena of `area` m³
and the time step is 7 days.

The stages are `tadpole`, `susceptible` (uninfected post-metamorphic frogs) and the cells of
[`chytrid_load_stage`](@ref) (`load1`, ...), infected frogs by the logarithm of their Bd
load. The transitions are

- `developing` and `metamorphosing`: a tadpole survives the week with chance
  `0.2^(7 / larval_period)` (a fifth survive the larval period) and metamorphoses with
  chance `7 / larval_period`, the metamorphs damped by `exp(-e^K N / area)`, where `N` is
  the pool named by `crowding` (`frogs`);
- `escaping` and `infecting`: a susceptible frog survives (`s0`) and is infected with the
  chance given in [`CHYTRID_SITES`](@ref), which reads the pool `zoospores`, while in the
  water; an infection starts at a load drawn from the distribution of a new infection;
- `growing`: an infected frog survives (`sI s0`), keeps its infection and its load changes;
- `clearing`: an infected frog survives and loses its infection, becoming susceptible;
- `breeding`: at the breeding pulse every post-metamorphic frog, infected or not, produces
  `fecundity` tadpoles. Marked as births.

As in the authors' code the load kernels are integrated over each cell of the next load (a
difference of normal distribution functions), not evaluated at its middle. The environment
of a step is `(day, temperature)` (see [`chytrid_environment`](@ref)).

The species adds its post-metamorphic frogs to the pool named by `crowding`, and the
zoospores its infected frogs shed in the week, `shedding × Σ e^x`, to the pool `shedding`
while they are in the water. It reads the pool `zoospores`, which
[`chytrid_zoospores`](@ref) provides; [`chytrid_community`](@ref) composes the two. In the
authors' model the recruitment of a species is damped by its own frogs only, so two hosts
composed into one community should be given pools of their own (`crowding =
:pipiens_frogs`, say); with the same pool, each host's recruitment would be damped by
both.
"""
function leopard_frog(site; area::Real = 500.0, cells::Integer = 300, name::Symbol = :frog, crowding::Symbol = :frogs)
    area > 0 || throw(ArgumentError("area must be positive"))
    crowding in (:zoospores, :shedding) && throw(ArgumentError("the crowding pool cannot be :$crowding"))
    s = _chytrid_site(site)
    load = chytrid_load_stage(; cells)
    h = step_size(load)
    weeks = s.larval_period / _CHYTRID_STEP
    metamorphosis = 1 / weeks
    tadpole_survival = 0.2^(1 / weeks)
    centred(ctx) = ctx.env.temperature - _CHYTRID_BASE_TEMPERATURE
    start_load(ctx) = s.a + s.a_temp * centred(ctx)              # mean log load of a new infection, and of growth at x = 0
    # chance (per unit of log load) of a log load in the cell around y next week
    in_cell(y, m) = _chytrid_cell_chance(y - h / 2, y + h / 2, m, s.sigma) / h
    losing(x) = _normal_cdf(-(x - s.mu_l) / s.sigma_l)          # chance of losing the infection at log load x
    infected_survival = s.sI * s.s0
    water(ctx) = _chytrid_in_water(s, ctx.env.day)
    infection(ctx) = water(ctx) ? 1 - exp(-(s.beta / area) * exp(s.beta_temp * centred(ctx)) * ctx.pool.zoospores) : 0.0
    damping(ctx) = exp(-(exp(s.K) / area) * ctx.pool[crowding])

    life = kernel_net([:tadpole, :susceptible, load],
        :developing => [(:tadpole => :tadpole) => tadpole_survival * (1 - metamorphosis)],
        :metamorphosing => [(:tadpole => :susceptible) => tadpole_survival * metamorphosis],
        :escaping => [(:susceptible => :susceptible) => s.s0],
        :infecting => [(:susceptible => :load) => responsive(ctx -> (y -> s.s0 * in_cell(y, start_load(ctx))))],
        :growing => [(:load => :load) =>
            responsive(ctx -> ((y, x) -> infected_survival * (1 - losing(x)) * in_cell(y, start_load(ctx) + s.b * x)))],
        :clearing => [(:load => :susceptible) => x -> infected_survival * losing(x)],
        :breeding => [(:susceptible => :tadpole) => s.fecundity, (:load => :tadpole) => x -> s.fecundity])
    net = births(respond(life,
        :metamorphosing => damping,
        :escaping => ctx -> 1 - infection(ctx),
        :infecting => infection,
        :breeding => ctx -> _chytrid_breeding_now(s, ctx.env.day) ? 1.0 : 0.0), :breeding)

    frogs = vcat([:susceptible => 1.0], trait_weights(load, x -> 1.0))
    shed = s.shedding .* exp.(meshpoints(load))
    shedding(n, env) = _chytrid_in_water(s, env.day) ? dot(shed, n[load]) : 0.0
    return species(name, net; occupies = [crowding => frogs, :shedding => shedding], reads = [:zoospores])
end

"""
    chytrid_zoospores(; survival = CHYTRID_ZOOSPORE_SURVIVAL, name = :zoospores)

The zoospores of Bd in the water, as an [`OpenCommunity`](@ref) of one stage,
`zoospores`, for [`chytrid_community`](@ref). Each week a share
`survival.maxsurv / (1 + exp((T - survival.ld50) / survival.s))` of them survives, at the
water temperature `T` of the step's environment (see [`CHYTRID_ZOOSPORE_SURVIVAL`](@ref)),
and the zoospores shed by the hosts, the pool `shedding`, are added (an inflow). Its
zoospores are the pool `zoospores`, which hosts read to know how many they meet.
"""
function chytrid_zoospores(; survival = CHYTRID_ZOOSPORE_SURVIVAL, name::Symbol = :zoospores)
    alive(T) = survival.maxsurv / (1 + exp((T - survival.ld50) / survival.s))
    pool = respond(ValuedProjectionNet([:zoospores], :surviving => [(:zoospores => :zoospores) => 1.0]),
                   :surviving => ctx -> alive(ctx.env.temperature))
    pool = inflow(pool, :zoospores => ctx -> ctx.pool.shedding)
    return species(name, pool; occupies = [:zoospores => [:zoospores => 1.0]], reads = [:shedding])
end

"""
    chytrid_community(site; area = 500, cells = 300)

The leopard frog of `site` ([`leopard_frog`](@ref)) and the Bd zoospores in its water
([`chytrid_zoospores`](@ref)), composed through the pools `zoospores` and `shedding`: the
full IPM of Wilber et al. (2022) for one host species. The pool `frogs` stays exposed, so
that another host can be added.
"""
chytrid_community(site; area::Real = 500.0, cells::Integer = 300) =
    community(leopard_frog(site; area, cells), chytrid_zoospores())
