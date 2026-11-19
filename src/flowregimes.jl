# Designed flow regimes for a whole river (Tonkin et al. 2021): riparian plants (riparian_model), fish
# (verde_fish_model) and three aquatic invertebrates, run under the flow regimes designed for each group and under
# the natural regime of the upper Verde River. The invertebrates are the flood-disturbance logistic models of
# McMullen et al. (2017) as coded in GitHub jdtonkin/flow-tradeoffs (MIT): one population each, a time step of a
# day, a carrying capacity that a flood resets and that then relaxes back, and a flood mortality. They share no
# pool; each has its own carrying capacity, scaled to 100.

"""
    FLOW_INVERTEBRATES

Parameters of the three invertebrates of Tonkin et al. (2021), from their code (`dod-3taxa-inverts-*.R`): `r`
intrinsic rate of increase per day, `h` strength of flood mortality (a flood of `Q` m³/s leaves `exp(-h Q)` of
the population), `g` rate per day at which the carrying capacity relaxes after a flood, `Kd` the carrying capacity
just after a strong flood and `Kb` the baseline without floods. The mayfly *Fallceon* and the dragonfly
*Progomphus* gain habitat from floods (`Kd > Kb`); the ostracod loses it. `Qmin = 5` m³/s is the smallest flow
that counts as a flood and `a = 100` m³/s the half-saturation of the flood effect on carrying capacity.
"""
const FLOW_INVERTEBRATES = (
    mayfly = (r = 0.23, h = 0.02, g = 0.01, Kd = 100.0, Kb = 40.0),      # Fallceon, fast life cycle
    dragonfly = (r = 0.08, h = 0.01, g = 0.01, Kd = 100.0, Kb = 40.0),   # Progomphus, slow life cycle
    ostracod = (r = 0.16, h = 0.05, g = 0.01, Kd = 40.0, Kb = 100.0),    # flood-averse
    Qmin = 5.0, a = 100.0)

_flood_effect(Q, Qmin, a) = Q < Qmin ? 0.0 : (Q - Qmin) / (a + (Q - Qmin))

# Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch)
const _GL = let m = 10
    β = [k / sqrt(4k^2 - 1) for k in 1:(m - 1)]
    F = eigen(SymTridiagonal(zeros(m), β))
    ((F.values .+ 1) ./ 2, F.vectors[1, :] .^ 2)
end

"""
    flood_invertebrate(name; parameters = FLOW_INVERTEBRATES[name], Qmin = 5, a = 100)

One aquatic invertebrate of Tonkin et al. (2021) as a net with a single stage, `N`, its abundance on a scale on
which the carrying capacity is at most 100, and a time step of a day. It has one transition, `living`, whose value
depends on the day's environment (see [`invertebrate_flows`](@ref)):

- on the first day of a flood (`env.t == 1`), `exp(-h Q)`: the flood kills;
- on any other day, the exact one-day solution of logistic growth `dN/dt = r N (1 - N / K(t))`, whose carrying
  capacity relaxes from its post-flood value `K0` to `Kb` as `K(t) = Kb + (K0 - Kb) exp(-g t)`, `t` being the day
  since the flood. `K0 = Kpre + (Kd - Kpre) Qf(Q)` moves the carrying capacity before the flood, `env.Kpre[name]`,
  towards `Kd` by the flood effect `Qf(Q) = (Q - Qmin) / (a + Q - Qmin)`.

Stepping a day at a time gives the closed form of McMullen et al. (2017) used in the authors' code, which solves
each inter-flood interval at once.
"""
function flood_invertebrate(name::Symbol; parameters = FLOW_INVERTEBRATES[name], Qmin::Real = FLOW_INVERTEBRATES.Qmin,
                            a::Real = FLOW_INVERTEBRATES.a)
    p = parameters
    life = ValuedProjectionNet([:N], :living => [(:N => :N) => 1.0])
    function today(ctx)
        env = ctx.env
        env.t == 1 && return exp(-p.h * env.Q)
        K0 = env.Kpre[name] + (p.Kd - env.Kpre[name]) * _flood_effect(env.Q, Qmin, a)
        K(y) = p.Kb + (K0 - p.Kb) * exp(-p.g * y)
        # 1/N(t) = exp(-r) / N(t - 1) + J, J = ∫ r / K(y) exp(-r (t - y)) dy over the day
        J = sum(w * p.r / K(env.t - 1 + x) * exp(-p.r * (1 - x)) for (x, w) in zip(_GL...))
        return 1 / (exp(-p.r) + ctx.n[:N] * J)
    end
    return species(name, respond(life, :living => today))
end

"""
    invertebrate_flows(t, Q; invertebrates = FLOW_INVERTEBRATES, corrected = false)

One environment per day for [`flood_invertebrate`](@ref)s: `t` is the day since the last flood (1 on the day of a
flood) and `Q` the flood's peak flow in m³/s on Bill Williams River scale, both one value per day. Each environment
holds `day`, `t`, `Q` and `Kpre`, the carrying capacity of each invertebrate just before the current flood; the
first flood finds every invertebrate at its baseline `Kb`.

The authors' code carries the carrying capacity from one flood to the next as `Kpre + (K0 - Kpre) exp(-g t)`, which
relaxes towards the value before the flood rather than towards `Kb`, although within an interval the population
grows under a carrying capacity that relaxes towards `Kb`. `corrected = false` (the default) reproduces the code;
`corrected = true` carries forward the carrying capacity that the population grew under.
"""
function invertebrate_flows(t::AbstractVector{<:Integer}, Q::AbstractVector{<:Real}; invertebrates = FLOW_INVERTEBRATES,
                            corrected::Bool = false)
    length(t) == length(Q) || throw(DimensionMismatch("need one flow per day"))
    taxa = Tuple(k for k in keys(invertebrates) if invertebrates[k] isa NamedTuple)
    Qmin, a = invertebrates.Qmin, invertebrates.a
    Kpre = Dict(k => invertebrates[k].Kb for k in taxa)
    last = Dict(k => invertebrates[k].Kb for k in taxa)       # carrying capacity carried to the next flood
    envs = map(eachindex(t)) do k
        if t[k] == 1 && k > 1
            for x in taxa
                Kpre[x] = last[x]
            end
        end
        for x in taxa
            p = invertebrates[x]
            K0 = Kpre[x] + (p.Kd - Kpre[x]) * _flood_effect(Q[k], Qmin, a)
            last[x] = corrected ? p.Kb + (K0 - p.Kb) * exp(-p.g * t[k]) : Kpre[x] + (K0 - Kpre[x]) * exp(-p.g * t[k])
        end
        (day = k, t = Int(t[k]), Q = Float64(Q[k]), Kpre = NamedTuple{taxa}(Tuple(Kpre[x] for x in taxa)))
    end
    return envs
end

"""
    flood_pulses(every, Q, days)

Days since the last flood and peak flow for `days` days with a flood of `Q` m³/s every `every` days, the first on
day 1: the arguments `t` and `Q` of [`invertebrate_flows`](@ref). `every` is truncated to whole days, as in the
authors' code (four pulses a year come every 91 days).
"""
function flood_pulses(every::Real, Q::Real, days::Integer)
    n = floor(Int, every)
    return [mod(k - 1, n) + 1 for k in 1:days], fill(Float64(Q), days)
end

"""
    TONKIN_SCENARIOS

The four flow regimes of Tonkin et al. (2021): `:plants` (the riparian prescription: a spring flood every sixth
year, preceded by a drought year), `:fish` (a spring flood every year), `:invertebrates` (four small pulses a year,
below the flood threshold of the plant and fish models, so that every year is a year without events for them) and
`:natural` (the upper Verde River near Paulden, 1964-2008, repeated).
"""
const TONKIN_SCENARIOS = (:plants, :fish, :invertebrates, :natural)

"""
    TONKIN_FISH_INITIAL

Negative binomial distributions, `(size, mu)`, from which the authors' code draws the initial adult biomass (g) of
each fish species in each of its 100 iterations (juveniles and subadults start at zero).
"""
const TONKIN_FISH_INITIAL = (CACL = (1.52, 5284.0), GIRO = (0.44, 2376.0), LECY = (0.34, 164.0),
                             CAIN = (1.33, 34068.0), MIDO = (0.66, 4202.0), CYLU = (1.78, 238.0), AMNA = (0.36, 1306.0))

"""
    tonkin_scenario(scenario; years = 200, invertebrate_years = 30, corrected = false)

The environments of one flow regime of Tonkin et al. (2021) (see [`TONKIN_SCENARIOS`](@ref)), for each group:

- `riparian`: one per year for the community of [`riparian_model`](@ref), with a flood, a drought or neither; seed
  release and drawdown are held at the values of the authors' code (`cprop = 0.9`, `tprop = 0.4`, no decline);
- `fish`: one per year for [`verde_fish_model`](@ref), with the year types of Rogosch et al. (2019);
- `invertebrates`: one per day for [`flood_invertebrate`](@ref)s (see [`invertebrate_flows`](@ref)).

Returns `(riparian, fish, invertebrates, riparian_n0, invertebrate_n0)`: the riparian start is 3978 m² of every
stage of every guild, the invertebrates start at 50. Peak flows of the natural regime are the mean spring and
autumn peaks at Paulden (1134 and 459 cfs), scaled by 0.35 to m³/s on the Bill Williams River, on which the
invertebrate models were parameterised. The factor is not a unit conversion (1 m³/s is about 35 cfs) but the
authors' match between corresponding floods on the two rivers: a bankfull flood is 220 cfs on the Verde and 87 m³/s
on the Bill Williams, a large flood 700 cfs and 215 m³/s, and 0.35 lies between the two ratios. The designed floods
use it too (`700 * 0.35`).
"""
function tonkin_scenario(scenario::Symbol; years::Integer = 200, invertebrate_years::Integer = 30, corrected::Bool = false,
                         path = joinpath(DATA_DIR, "verde_fish", "naturalflow_long.csv"))
    scenario in TONKIN_SCENARIOS || throw(ArgumentError("scenario must be one of $(TONKIN_SCENARIOS)"))
    cycle(v) = [v[mod(k - 1, length(v)) + 1] for k in 1:years]
    if scenario === :plants
        flood, drought = cycle([0, 0, 0, 0, 0, 1]), cycle([0, 0, 0, 0, 1, 0])
        fishyears = [(SP_highflood = f, SU_highflood = 0, medflood = 0, drought = d, nonevent = 1 - f - d)
                     for (f, d) in zip(flood, drought)]
    elseif scenario === :fish
        flood, drought = fill(1, years), fill(0, years)
        fishyears = fill((SP_highflood = 1, SU_highflood = 0, medflood = 0, drought = 0, nonevent = 0), years)
    elseif scenario === :invertebrates
        flood, drought = fill(0, years), fill(0, years)
        fishyears = fill((SP_highflood = 0, SU_highflood = 0, medflood = 0, drought = 0, nonevent = 1), years)
    else
        rows = [parse.(Int, split(replace(l, "\"" => ""), ",")) for l in readlines(path)[2:end]]
        rows = cycle([r for r in rows if r[1] <= 2008])
        fishyears = [(SP_highflood = r[3], SU_highflood = r[4], medflood = r[5], drought = r[6], nonevent = r[7])
                     for r in rows]
        flood, drought = [r[3] for r in rows], [r[6] for r in rows]
    end
    riparian = [(year = k, flood = Float64(flood[k]), drought = Float64(drought[k]),
                 nonflood = flood[k] + drought[k] == 1 ? 0.0 : 1.0,
                 cprop = 0.9, tprop = 0.4, decline = 1.0, tdecline = 1.0) for k in 1:years]
    days = 365 * invertebrate_years
    t, Q = scenario === :invertebrates ? flood_pulses(365 / 4, 50, days) :
           scenario === :fish ? flood_pulses(365, 700 * 0.35, days) :
           scenario === :plants ? flood_pulses(365 * 6, 700 * 0.35, days) :
           ([mod(k - 1, 365) < 200 ? mod(k - 1, 365) + 1 : mod(k - 1, 365) - 199 for k in 1:days],
            [mod(k - 1, 365) < 200 ? 1134 * 0.35 : 459 * 0.35 for k in 1:days])
    riparian_n0 = [g => 3978 .* RIPARIAN_GUILDS[g].den for g in (:C, :T, :W, :S, :M)]
    return (riparian = riparian, fish = fishyears, invertebrates = invertebrate_flows(t, Q; corrected),
            riparian_n0 = riparian_n0, invertebrate_n0 = [:mayfly => [50.0], :dragonfly => [50.0], :ostracod => [50.0]])
end

"""
    flow_invertebrates(; invertebrates = FLOW_INVERTEBRATES)

The mayfly, dragonfly and ostracod of Tonkin et al. (2021), each a [`flood_invertebrate`](@ref), side by side in
one community. They share no pool.
"""
flow_invertebrates(; invertebrates = FLOW_INVERTEBRATES) =
    community((flood_invertebrate(x; parameters = invertebrates[x], Qmin = invertebrates.Qmin, a = invertebrates.a)
               for x in (:mayfly, :dragonfly, :ostracod))...)

# a gamma variate of shape k and scale 1 (Marsaglia & Tsang 2000, with the boost for k < 1)
function _rand_gamma(rng, k::Real)
    k < 1 && return _rand_gamma(rng, k + 1) * rand(rng)^(1 / k)
    d = k - 1 / 3
    c = 1 / sqrt(9d)
    while true
        x = randn(rng)
        v = (1 + c * x)^3
        v <= 0 && continue
        u = rand(rng)
        log(u) < x^2 / 2 + d - d * v + d * log(v) && return d * v
    end
end

"""
    tonkin_fish_start(rng)

One draw of the initial fish biomass of the authors' code: for each species of [`VERDE_SPECIES`](@ref), adults
from the negative binomial distribution of [`TONKIN_FISH_INITIAL`](@ref) (a Poisson draw with a gamma-distributed
mean), juveniles and subadults zero. Returns `species => [0, 0, adults]` pairs for [`verde_fish_model`](@ref).
"""
function tonkin_fish_start(rng)
    return [Symbol(sp) => begin
                size, mu = TONKIN_FISH_INITIAL[Symbol(sp)]
                [0.0, 0.0, Float64(StructuredPopulationCore.rand_poisson(rng, mu / size * _rand_gamma(rng, size)))]
            end for sp in VERDE_SPECIES]
end
