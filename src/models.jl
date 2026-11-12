# Published matrix community models. Each species is a projection net: a life cycle of
# named transitions, with responses saying how each transition depends on the year and on
# the shared pool.

const DATA_DIR = joinpath(dirname(@__DIR__), "data")

_csv(path) = [split(replace(l, "\"" => ""), ",") for l in readlines(path)[2:end]]

"""
    density_response(pool, f)

A response that applies a density-response function from MatrixProjectionModels.jl (for
example `RickerDensity(0.0, 1.0)` or `BevertonHoltDensity(1.0, 0.01)`) to the total of
`pool`: use it as `respond(net, :transition => density_response(:density, f))`.
"""
density_response(pool::Symbol, f) = ctx -> f(getproperty(ctx.pool, pool))

# ---------------------------------------------------------------- Verde River fish

"""
    VERDE_SPECIES

Species codes of the Verde River fish model, in the order of Rogosch et al. (2019).
"""
const VERDE_SPECIES = ["CACL", "GIRO", "LECY", "CAIN", "MIDO", "CYLU", "AMNA"]

const _VERDE_YEAR_TYPES = (:SUHF => :SU_highflood, :SPHF => :SP_highflood, :MF => :medflood,
                           :NE => :nonevent, :DR => :drought)

function _verde_species(sp, rate, modifier, K)
    r(code) = rate[code]
    spawners = sp in ("CACL", "GIRO", "LECY") ? [:s2 => 1, :s3 => 1] :
               sp == "CYLU" ? [:s1 => 1, :s2 => 2, :s3 => 2] : [:s3 => 1]
    # the life cycle in a benign year with an empty river; growth is in biomass
    life = ValuedProjectionNet([:s1, :s2, :s3],
        :juvenile_growth => [(:s1 => :s2) => r("a$(sp)1") * r("den$(sp)J") * (1 / r("den$(sp)2"))],
        :adult_growth => [(:s2 => :s3) => r("a$(sp)2") * r("den$(sp)2") * (1 / r("den$(sp)3"))],
        :adult_stasis => [(:s3 => :s3) => 1 - r("a$(sp)3")],
        :spawning => [(stage => :s1) => 0.5 * mult * r("GSI.$sp") * (1 - r("S0Mort$sp")) for (stage, mult) in spawners])
    # each kind of year takes its toll on juveniles and adults
    toll(stage) = [ctx -> 1 - getfield(ctx.env, flag) * r("STMort$sp") * modifier["$(sp)_$(stage)_$yt"]
                   for (yt, flag) in _VERDE_YEAR_TYPES]
    # the dependency assumption: egg survival falls to zero as total biomass approaches K
    room(ctx) = max((K - ctx.pool.biomass) / K, 0)
    eggs_to_biomass = (ctx -> r("den$(sp)1"), ctx -> 1 / r("den$(sp)J"))
    fish = respond(life,
        (:juvenile_growth .=> toll("J"))..., (:adult_growth .=> toll("A"))..., (:adult_stasis .=> toll("A"))...,
        :spawning => room, (:spawning .=> eggs_to_biomass)...)
    return species(Symbol(sp), births(fish, :spawning); occupies = [:biomass => 1.0])
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
    stage(i) = _RIP_STAGES[i]
    index(s) = findfirst(==(s), _RIP_STAGES)
    # corrected: sagebrush and meadow keep their non-maturing stage-5 individuals
    stasis = corrected && name in (:S, :M) ? (5:6) : g.stasis

    # the life cycle in a year with neither flood nor drought
    life = ValuedProjectionNet(_RIP_STAGES,
        :growing => [(stage(i) => stage(i + 1)) => g.a[i] * b[i] for i in 1:5],
        :staying => [(stage(i) => stage(i)) => 1 - g.a[i] for i in stasis])
    # floods and droughts kill a stage-specific fraction
    flood_survival(ctx) = (from, to) -> 1 - ctx.env.flood * g.S[index(from)]
    drought_survival(ctx) = (from, to) -> 1 - ctx.env.drought * g.D[index(from)]
    guild = respond(life, :growing => flood_survival, :growing => drought_survival,
                    :staying => flood_survival, :staying => drought_survival)

    # the dependency assumption: space occupied after this year's mortality, in
    # cottonwood-seedling units; seedlings can only establish on what is left
    G(i, env) = g.a[i] * b[i] * (1 - env.flood * g.S[i]) * (1 - env.drought * g.D[i])
    P(i, env) = (1 - g.a[i]) * (1 - env.flood * g.S[i]) * (1 - env.drought * g.D[i])
    cum(k) = prod(b[j] for j in k:-1:1)
    spaceP = name in (:S, :M) ? (5:6) : g.stasis     # the original counts a stage-5 stasis for S and M regardless
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
    free(ctx) = _RIP_K - ctx.pool.space
    seedlings(ctx, adults, gate) = max(adults * gate * (g.den[1] / _DEN_C1) * free(ctx), 0.0)
    occupies = Pair{Symbol, Any}[:space => space]
    reads = Symbol[]

    if name === :C
        # one tree can seed all free space, so the per-capita rate carries 1/N6
        nonind(x) = x == 0 ? 666.0 : x
        seeding = ValuedProjectionNet(_RIP_STAGES, :seeding => [(:s6 => :s1) => 1.0])
        guild = births(guild ⊕ respond(seeding, :seeding => ctx ->
            max(_adult(ctx.n[6]) * (1 / nonind(ctx.n[6])) * ctx.env.flood * ctx.env.cprop * ctx.env.decline * free(ctx), 0.0)),
            :seeding)
    elseif name === :T
        guild = inflow(guild, :s1 => ctx -> max(_adult(ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]) * ctx.env.flood *
                                                ctx.env.tprop * ctx.env.tdecline * (g.den[1] / _DEN_C1) * free(ctx), 0.0))
    elseif name === :W
        if corrected
            guild = inflow(guild, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]), ctx.env.flood))
        else
            # willow announces whether it has adults, because the original meadow guild reads it
            push!(occupies, :willow_adults => (n, env) -> _adult(n[3] + n[4] + n[5] + n[6]))
            guild = inflow(guild, :s1 => ctx -> seedlings(ctx, ctx.pool.willow_adults, ctx.env.flood))
        end
    elseif name === :S
        guild = inflow(guild, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[2] + ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]),
                                                      ctx.env.nonflood))
    elseif name === :M
        if corrected
            guild = inflow(guild, :s1 => ctx -> seedlings(ctx, _adult(ctx.n[2] + ctx.n[3] + ctx.n[4] + ctx.n[5] + ctx.n[6]),
                                                          ctx.env.flood))
        else
            # the original script gives the meadow guild the willow recruitment
            push!(reads, :willow_adults)
            guild = inflow(guild, :s1 => ctx -> seedlings(ctx, ctx.pool.willow_adults, ctx.env.flood))
        end
    end
    guild = after_step(guild, name === :S ? _quasiten : _quasi)
    return species(name, guild; occupies, reads)
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
Production of dormant seed (`SD`) and of protocorms (`P1`) by flowering plants is
multiplied by the Ricker terms `0.5 exp(-density)` and `1.2 exp(-2 density)`, using
`RickerDensity` from MatrixProjectionModels.jl. With `integeronly` the stage vectors are
floored each step, as in adapt3.

The annual matrices were estimated in R with lefko3 and are shipped as data. Returns
`(community, n0, envs)`; the environment of step `t` is the integer `t - 1`, which selects
each species' matrix for that year.
"""
function cypripedium_model(; equivalence::Bool = true, integeronly::Bool = false,
                           dir = joinpath(DATA_DIR, "cypripedium"))
    stages = Symbol.(readlines(joinpath(dir, "stages.txt")))
    k = length(stages)
    position = Dict(s => i for (i, s) in enumerate(stages))
    mats = Dict{Int, Dict{Int, Matrix{Float64}}}()
    for row in _csv(joinpath(dir, "matrices.csv"))
        sp, yr, i, j = parse(Int, row[1]), parse(Int, row[2]), parse(Int, row[4]), parse(Int, row[5])
        M = get!(get!(mats, sp, Dict{Int, Matrix{Float64}}()), yr, zeros(k, k))
        M[i, j] = parse(Float64, row[6])
    end
    weights = equivalence ? [0.0; ones(k - 1)] : ones(k)
    flowering = (:XSm, :Sm, :Md, :Lg, :XLg)
    function orchid(name, annual)
        used(i, j) = any(M[i, j] != 0 for M in annual)
        entries(keep) = [(from => to) => 1.0 for (j, from) in enumerate(stages) for (i, to) in enumerate(stages)
                         if used(i, j) && keep(from, to)]
        seed(from, to) = to === :SD && from in flowering
        protocorm(from, to) = to === :P1 && from in flowering
        # every element estimated from the field data is a transition; three kinds are told apart
        life = ValuedProjectionNet(stages,
            :life_cycle => entries((from, to) -> !seed(from, to) && !protocorm(from, to)),
            :seed_production => [e.first => 0.5 for e in entries(seed)],
            :protocorm_production => [e.first => 1.2 for e in entries(protocorm)])
        # this year's estimate of each element
        this_year(ctx) = (from, to) -> annual[mod(ctx.env, length(annual)) + 1][position[to], position[from]]
        plant = respond(life,
            :life_cycle => this_year, :seed_production => this_year, :protocorm_production => this_year,
            :seed_production => density_response(:density, MatrixProjectionModels.RickerDensity(0.0, 1.0)),
            :protocorm_production => density_response(:density, MatrixProjectionModels.RickerDensity(0.0, 2.0)))
        plant = births(plant, :seed_production, :protocorm_production)
        integeronly && (plant = after_step(plant, n -> floor.(n)))
        return species(name, plant; occupies = [:density => weights])
    end
    com = community(orchid(:candidum, [mats[1][y] for y in 1:length(mats[1])]),
                    orchid(:parviflorum, [mats[2][y] for y in 1:length(mats[2])]))
    n0 = [:candidum => [500.0; 100.0; 200.0; zeros(k - 3)], :parviflorum => [5000.0; 1000.0; 2000.0; zeros(k - 3)]]
    return (community = com, n0 = n0, envs = 0:99)
end
