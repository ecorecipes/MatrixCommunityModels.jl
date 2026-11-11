# Published matrix community models, rebuilt by composition.

const DATA_DIR = joinpath(dirname(@__DIR__), "data")

_csv(path) = [split(replace(l, "\"" => ""), ",") for l in readlines(path)[2:end]]

# ---------------------------------------------------------------- Verde River fish

"""
    VERDE_SPECIES

Species codes of the Verde River fish model, in the order of Rogosch et al. (2019).
"""
const VERDE_SPECIES = ["CACL", "GIRO", "LECY", "CAIN", "MIDO", "CYLU", "AMNA"]

function _verde_species(sp, rate, modifier, K)
    r(code) = rate[code]
    yeartypes = (:SUHF => :SU_highflood, :SPHF => :SP_highflood, :MF => :medflood,
                 :NE => :nonevent, :DR => :drought)
    # year-type mortality: one factor per year type, multiplied in the original order
    function mort(base, stage, env)
        x = base
        for (yt, flag) in yeartypes
            x *= 1 - getfield(env, flag) * r("STMort$sp") * modifier["$(sp)_$(stage)_$yt"]
        end
        return x
    end
    G1 = ctx -> mort(r("a$(sp)1") * r("den$(sp)J") * (1 / r("den$(sp)2")), "J", ctx.env)
    G2 = ctx -> mort(r("a$(sp)2") * r("den$(sp)2") * (1 / r("den$(sp)3")), "A", ctx.env)
    P3 = ctx -> mort(1 - r("a$(sp)3"), "A", ctx.env)
    # the dependency assumption: egg survival falls to zero as total biomass approaches K
    F = mult -> (ctx -> ((0.5 * mult * r("GSI.$sp") * (1 - r("S0Mort$sp"))) *
                         max((K - ctx.pool.biomass) / K, 0)) * r("den$(sp)1") * (1 / r("den$(sp)J")))
    reproducing = sp in ("CACL", "GIRO", "LECY") ? [:s2 => 1, :s3 => 1] :
                  sp == "CYLU" ? [:s1 => 1, :s2 => 2, :s3 => 2] : [:s3 => 1]
    transitions = [survival(:growth1, :s1 => :s2, G1), survival(:growth2, :s2 => :s3, G2),
                   survival(:stasis3, :s3 => :s3, P3)]
    for (stage, mult) in reproducing
        push!(transitions, fecundity(Symbol(:eggs_, stage), stage => :s1, F(mult)))
    end
    return Species(Symbol(sp), [:s1, :s2, :s3], transitions; occupancy = [:biomass => ones(3)])
end

"""
    verde_fish_model(; K = 47660.0, dir = joinpath(DATA_DIR, "verde_fish"))

The seven-species Verde River fish community of Rogosch et al. (2019): three-stage
biomass models linked by one pool, total community biomass, which limits every species'
fecundity through `max(0, 1 - biomass / K)`.

Returns `(community, n0, envs, years, density)`: the composed community, the initial
biomass (the stage-3 values of the vital-rate table), one environment per year of the
1964-2017 flow record (year-type flags), the years, and the individuals-per-gram
conversion for each species' stages.
"""
function verde_fish_model(; K = 47660.0, dir = joinpath(DATA_DIR, "verde_fish"))
    rate = Dict(String(r[2]) => parse(Float64, r[3]) for r in _csv(joinpath(dir, "vital-rates.csv")))
    modifier = Dict(String(r[4]) => parse(Float64, r[5]) for r in _csv(joinpath(dir, "modifiers-all-spp.csv")))
    flow = _csv(joinpath(dir, "naturalflow_long.csv"))
    envs = [(SP_highflood = parse(Int, r[3]), SU_highflood = parse(Int, r[4]), medflood = parse(Int, r[5]),
             drought = parse(Int, r[6]), nonevent = parse(Int, r[7])) for r in flow]
    com = community((_verde_species(sp, rate, modifier, K) for sp in VERDE_SPECIES)...)
    n0 = [Symbol(sp) => [0.0, 0.0, rate["biom$(sp)3"]] for sp in VERDE_SPECIES]
    density = Dict(Symbol(sp) => [rate["den$(sp)J"], rate["den$(sp)2"], rate["den$(sp)3"]] for sp in VERDE_SPECIES)
    return (community = com, n0 = n0, envs = envs, years = [parse(Int, r[1]) for r in flow], density = density)
end

# ---------------------------------------------------------------- Riparian guilds

"""
    RIPARIAN_GUILDS

Parameters of the five riparian guilds of Lytle et al. (2017): cottonwood `C`, tamarisk
`T`, willow `W`, sagebrush `S`, meadow `M`. `den` is stem density per stage (per m²), `a`
the maturation probabilities, `S` and `D` the flood and drought mortalities, `stasis` the
stages with a stasis term in the original projection matrix, and `area` the initial area
per stage (m²).
"""
const RIPARIAN_GUILDS = (
    C = (den = [350, 10, 1, 0.91, 0.6, 0.12], a = [1, 1, 1, 1, 0.167, 0.03],
         S = [0.97, 0.33, 0.224, 0.19, 0.073, 0.02], D = [0.49, 0.16, 0.083, 0.05, 0.05, 0.05],
         stasis = 5:6, area = [11816, 11816, 11144, 11144, 13819, 59605]),
    T = (den = [400, 29, 4.5, 1.4, 1.3, 1.3], a = [1, 1, 0.25, 0.11, 0.07, 0.05],
         S = [0.9, 0.55, 0.25, 0.05, 0.01, 0.01], D = [0.5, 0.15, 0.05, 0.025, 0.025, 0.025],
         stasis = 3:6, area = [1000, 0, 0, 0, 0, 0]),
    W = (den = [350, 35, 1, 1, 1, 1], a = [1, 1, 1, 1, 0.167, 0.01],
         S = [0.49, 0.17, 0.11, 0.10, 0.04, 0.01], D = [0.75, 0.58, 0.54, 0.53, 0.51, 0.51],
         stasis = 5:6, area = fill(100, 6)),
    S = (den = [350, 35, 1, 1, 1, 1], a = [1, 1, 1, 1, 0.167, 0.03],
         S = [0.99, 0.67, 0.61, 0.60, 0.54, 0.51], D = [0.24, 0.08, 0.042, 0.025, 0.005, 0.005],
         stasis = 6:6, area = fill(100, 6)),
    M = (den = [350, 35, 1, 1, 1, 1], a = [1, 1, 1, 1, 0.167, 0.03],
         S = [0.99, 0.67, 0.61, 0.19, 0.073, 0.02], D = [0.49, 0.16, 0.083, 0.025, 0.005, 0.005],
         stasis = 6:6, area = fill(100, 6)))

const _RIP_STAGES = [:s1, :s2, :s3, :s4, :s5, :s6]
const _RIP_K = 119344 * 350.0
const _DEN_C1 = 350.0

_adult(x) = x > 0.99999 ? 1.0 : 0.0
_quasi(n) = [x < 1 ? 0.0 : x for x in n]
_quasiten(n) = [x < 10 ? 10.0 : x for x in n]

function _riparian_environments(path)
    rows = [parse.(Float64, split(l, ",")) for l in readlines(path)[2:end]]
    return map(rows) do r
        mag, fdate, slopecot, slopetam, bdate = r[2], r[3], r[4], r[5], r[6]
        flood = mag > 9888 ? 1.0 : 0.0
        (year = Int(r[1]), flood = flood, drought = mag < 7416 ? 1.0 : 0.0, nonflood = 1 - flood,
         cprop = max(0.0, min(bdate, 298) - max(fdate, 260)) / 38,
         tprop = max(0.0, min(bdate, 336) - max(fdate, 280)) / 56,
         decline = 0.944 * exp(-0.5 * (log(-slopecot / 1.279) / 0.987)^2),
         tdecline = 0.917 * exp(-0.5 * ((-slopetam - 1.8) / 3.4)^2))
    end
end

function _riparian_guild(name::Symbol, corrected::Bool, guilds)
    g = guilds[name]
    b = [g.den[i + 1] / g.den[i] for i in 1:5]                 # self-thinning
    G(i, env) = g.a[i] * b[i] * (1 - env.flood * g.S[i]) * (1 - env.drought * g.D[i])
    P(i, env) = (1 - g.a[i]) * (1 - env.flood * g.S[i]) * (1 - env.drought * g.D[i])
    cum(k) = prod(b[j] for j in k:-1:1)
    # corrected: sagebrush and meadow keep their non-maturing stage-5 individuals
    stasis = corrected && name in (:S, :M) ? (5:6) : g.stasis
    # the original counts the space of a stage-5 stasis term for S and M regardless
    spaceP = name in (:S, :M) ? (5:6) : g.stasis

    # space occupied after this year's mortality, converted to cottonwood-seedling units
    function space(n, env)
        x = n[1] * G(1, env) / b[1] + n[2] * G(2, env) / cum(2)
        for i in 3:5
            gterm = n[i] * G(i, env) / cum(i)
            pterm = i in spaceP ? n[i] * P(i, env) / cum(i - 1) : nothing
            if name === :T
                x += gterm
                pterm === nothing || (x += pterm)
            else
                pterm === nothing || (x += pterm)
                x += gterm
            end
        end
        x += n[6] * P(6, env) / cum(5)
        return x * (_DEN_C1 / g.den[1])
    end

    transitions = [survival(Symbol(:grow, i), _RIP_STAGES[i] => _RIP_STAGES[i + 1], ctx -> G(i, ctx.env))
                   for i in 1:5]
    for i in stasis
        push!(transitions, survival(Symbol(:stay, i), _RIP_STAGES[i] => _RIP_STAGES[i], ctx -> P(i, ctx.env)))
    end
    free(ctx) = _RIP_K - ctx.pool.space
    seedlings(ctx, adults, gate) = max(adults * gate * (g.den[1] / _DEN_C1) * free(ctx), 0.0)
    occupancy = Pair{Symbol,Any}[:space => space]
    recruitment = Pair{Symbol,Any}[]
    if name === :C
        # one tree can seed all free space, so the per-capita rate carries 1/N6
        nonind(x) = x == 0 ? 666.0 : x
        push!(transitions, fecundity(:seedlings, :s6 => :s1,
              ctx -> max(_adult(ctx.n[6]) * (1 / nonind(ctx.n[6])) * ctx.env.flood * ctx.env.cprop *
                         ctx.env.decline * free(ctx), 0.0)))
    elseif name === :T
        push!(recruitment, :s1 => ctx -> max(_adult(ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]) * ctx.env.flood *
                                             ctx.env.tprop * ctx.env.tdecline * (g.den[1] / _DEN_C1) * free(ctx), 0.0))
    elseif name === :W
        if corrected
            push!(recruitment, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]), ctx.env.flood))
        else
            # willow announces whether it has adults, because the original meadow guild reads it
            push!(occupancy, :willow_adults => (n, env) -> _adult(n[3] + n[4] + n[5] + n[6]))
            push!(recruitment, :s1 => ctx -> seedlings(ctx, ctx.pool.willow_adults, ctx.env.flood))
        end
    elseif name === :S
        push!(recruitment, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[2] + ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]),
                                                   ctx.env.nonflood))
    elseif name === :M
        if corrected
            push!(recruitment, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[2] + ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]),
                                                       ctx.env.flood))
        else
            # the original script gives the meadow guild the willow recruitment
            push!(occupancy, :willow_adults => (n, env) -> 0.0)
            push!(recruitment, :s1 => ctx -> seedlings(ctx, ctx.pool.willow_adults, ctx.env.flood))
        end
    end
    return Species(name, _RIP_STAGES, transitions; occupancy, recruitment,
                   threshold = name === :S ? _quasiten : _quasi)
end

"""
    riparian_model(; corrected = false, guilds = RIPARIAN_GUILDS, path = ...)

The five-guild riparian vegetation community of Lytle et al. (2017). Guilds are linked by
one pool, `space`: the area occupied after the year's flood and drought mortality, in
cottonwood-seedling units. Recruitment of every guild is proportional to the free space
`K - space`.

With `corrected = false` the model reproduces the published R script exactly, including
two behaviours that look unintended: the meadow guild receives the willow guild's
recruitment (expressed here as a second pool, `willow_adults`), and sagebrush and meadow
lose the stage-5 individuals that do not mature. `corrected = true` gives meadow its own
recruitment, gated on its own adults, and keeps those stage-5 individuals. Both variants
keep the script's floor of 10 individuals per sagebrush stage.

`guilds` supplies the parameters; pass a copy changed with [`with_rate`](@ref) to
perturb a vital rate.

Returns `(community, n0, envs, K)`; `envs` holds one environment per year of the
1916-1998 Yampa River record at Maybell.
"""
function riparian_model(; corrected::Bool = false, guilds = RIPARIAN_GUILDS,
                        path = joinpath(DATA_DIR, "riparian", "maybell.csv"))
    names = (:C, :T, :W, :S, :M)
    com = community((_riparian_guild(g, corrected, guilds) for g in names)...)
    n0 = [g => guilds[g].area .* guilds[g].den for g in names]
    return (community = com, n0 = n0, envs = _riparian_environments(path), K = _RIP_K)
end

"""
    with_rate(guilds, guild, field, stage, value)

A copy of a guild parameter table such as [`RIPARIAN_GUILDS`](@ref) in which
`guilds[guild][field][stage]` is `value`; for example
`with_rate(RIPARIAN_GUILDS, :C, :S, 6, 0.05)` sets the flood mortality of stage-6
cottonwood.
"""
function with_rate(guilds::NamedTuple, guild::Symbol, field::Symbol, stage::Integer, value::Real)
    g = guilds[guild]
    v = collect(Float64, g[field])
    v[stage] = value
    return merge(guilds, NamedTuple{(guild,)}((merge(g, NamedTuple{(field,)}((v,))),)))
end

# ---------------------------------------------------------------- Cypripedium (adapt3)

"""
    cypripedium_model(; equivalence = true, integeronly = false, dir = joinpath(DATA_DIR, "cypripedium"))

The two-orchid community of the adapt3 community-projection vignette (Shefferson 2026):
*Cypripedium candidum* (5 annual matrices) and *C. parviflorum* (17), nine stages each,
cycled through in order. The species are linked by one pool, `density`, the number of
individuals of both species, not counting dormant seed when `equivalence` is true.
Reproduction into dormant seed (`SD`) and protocorms (`P1`) is multiplied by Ricker
terms `0.5 exp(-density)` and `1.2 exp(-2 density)`. With `integeronly` the stage vectors
are floored each step, as in adapt3.

The annual matrices were estimated in R with lefko3 and are shipped as data. Returns
`(community, n0, envs)`; the environment of step `t` is the integer `t - 1`, which selects
each species' matrix for that year.
"""
function cypripedium_model(; equivalence::Bool = true, integeronly::Bool = false,
                           dir = joinpath(DATA_DIR, "cypripedium"))
    stages = Symbol.(readlines(joinpath(dir, "stages.txt")))
    k = length(stages)
    mats = Dict{Int,Dict{Int,Matrix{Float64}}}()
    for row in _csv(joinpath(dir, "matrices.csv"))
        sp, yr, i, j = parse(Int, row[1]), parse(Int, row[2]), parse(Int, row[4]), parse(Int, row[5])
        M = get!(get!(mats, sp, Dict{Int,Matrix{Float64}}()), yr, zeros(k, k))
        M[i, j] = parse(Float64, row[6])
    end
    weights = equivalence ? [0.0; ones(k - 1)] : ones(k)
    reproductive = (:XSm, :Sm, :Md, :Lg, :XLg)
    function orchid(name, annual)
        year(ctx) = mod(ctx.env, length(annual)) + 1
        transitions = Transition[]
        for (j, from) in enumerate(stages), (i, to) in enumerate(stages)
            any(M[i, j] != 0 for M in annual) || continue
            label = Symbol(from, :_to_, to)
            if to === :SD && from in reproductive
                push!(transitions, fecundity(label, from => to,
                      ctx -> annual[year(ctx)][i, j] * 0.5 * exp(-1.0 * ctx.pool.density)))
            elseif to === :P1 && from in reproductive
                # adapt3 types this element as survival, but it is production of protocorms
                push!(transitions, fecundity(label, from => to,
                      ctx -> annual[year(ctx)][i, j] * 1.2 * exp(-2.0 * ctx.pool.density)))
            else
                push!(transitions, survival(label, from => to, ctx -> annual[year(ctx)][i, j]))
            end
        end
        return Species(name, stages, transitions; occupancy = [:density => weights],
                       threshold = integeronly ? (n -> floor.(n)) : nothing)
    end
    com = community(orchid(:candidum, [mats[1][y] for y in 1:length(mats[1])]),
                    orchid(:parviflorum, [mats[2][y] for y in 1:length(mats[2])]))
    n0 = [:candidum => [500.0; 100.0; 200.0; zeros(k - 3)], :parviflorum => [5000.0; 1000.0; 2000.0; zeros(k - 3)]]
    return (community = com, n0 = n0, envs = 0:99)
end
