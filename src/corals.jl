# Three coral genera on the fore reefs of Mo'orea as an open-population integral projection model.
#
# The model is that of Kayal, Lenihan, Brooks, Holbrook, Schmitt & Kendall (2018, Ecology Letters 21:1790-1799)
# for Acropora, Pocillopora and Porites recolonising four fore reefs (LTER 1, 2, 4 and 5) after a predator
# outbreak and a cyclone. Size is z = log10(S + 1), S the surface area of a colony's living tissue in cm², on
# [0, 7] in 100 cells; the time step is a year; numbers are colonies per 10 m² of reef.
#
# Survival is logistic and growth normal in z, with the parameters of Table S1. Recruits come from outside
# (an open population): their number depends on local crowding, and they arrive with z normal (mean 0.6, sd
# 0.2). Recruitment of Acropora rises with the density of Acropora, that of Porites (barely) with the density
# of Porites, and that of Pocillopora falls with the summed log sizes of all three genera. The species meet
# only through pools:
#
#   acropora_colonies, pocillopora_colonies, porites_colonies      colonies per 10 m²
#   acropora_log_size, pocillopora_log_size, porites_log_size      Σ z over colonies, per 10 m²
#
# Two readings of the recruitment function, log r = R0 + Rγ γ (Appendix S3, eq. iii), are offered: the one written
# in the authors' script (`:as_written`, the default) and one that reproduces the published figures (`:calibrated`).
# See `coral`.
#
# The parameters are read from data/corals, typed from Table S1; the colonies observed in 2013 and 2015, from which the
# paper's updated runs start (`coral_observed_state`), are binned from MCR LTER data (CC BY 4.0), and also totalled by
# transect (`coral_observed_transects`). See
# data/corals/LICENSE.txt.

const _CORAL_TAXA = (acropora = "Acr", pocillopora = "Poc", porites = "Por")
const _CORAL_NAMES = (acropora = "Acropora", pocillopora = "Pocillopora", porites = "Porites")
const _CORAL_REEFS = (lter1 = 1, lter2 = 2, lter4 = 4, lter5 = 5)
const _CORAL_READINGS = (:as_written, :calibrated)

function _coral_parameters(file)
    lines = [split(l, ",") for l in readlines(file)]
    header, rows = lines[1], lines[2:end]
    col(name) = findfirst(==(name), header)
    function one(reef, taxon)
        r = only(row for row in rows if parse(Int, row[col("reef")]) == _CORAL_REEFS[reef] && row[col("taxon")] == _CORAL_TAXA[taxon])
        x(name) = parse(Float64, r[col(name)])
        return (survival = (intercept = x("survival_intercept"), slope = x("survival_slope")),
                growth = (intercept = x("growth_intercept"), slope = x("growth_slope"), sd = x("growth_sd")),
                recruitment = (intercept = x("recruitment_intercept"), slope = x("recruitment_slope"),
                               crowding = taxon === :pocillopora ? :assemblage_log_size : :own_density))
    end
    reefs = keys(_CORAL_REEFS)
    return NamedTuple{reefs}(Tuple(NamedTuple{keys(_CORAL_TAXA)}(Tuple(one(r, t) for t in keys(_CORAL_TAXA))) for r in reefs))
end

"""
    CORAL_PARAMETERS

Demographic parameters of three coral genera on four fore reefs of Mo'orea, French
Polynesia, for [`coral`](@ref): Table S1 of Kayal et al. (2018), estimated from colonies
followed in 2011-2013 (survival and growth, by generalised linear mixed models) and from
recruits counted in 2011-2012 (the intercept of recruitment); the slopes of recruitment come
from an earlier study (Kayal et al. 2015).

`CORAL_PARAMETERS[reef][taxon]`, with `reef` one of `:lter1`, `:lter2` (north shore),
`:lter4` (east), `:lter5` (west) and `taxon` one of `:acropora`, `:pocillopora`,
`:porites`, holds:

- `survival`: `intercept` and `slope` of the logit of survival in size `z`;
- `growth`: `intercept` and `slope` of the mean size next year, and its standard deviation
  `sd`;
- `recruitment`: `intercept` R0 and `slope` Rγ of log recruitment, and `crowding`, the
  quantity γ it responds to (`:own_density`, the colonies of the taxon itself, for
  *Acropora* and *Porites*; `:assemblage_log_size`, the summed sizes of all colonies of the
  three taxa, for *Pocillopora*).

Size `z` is log10(S + 1), with S the surface area of a colony's living tissue in cm².
"""
const CORAL_PARAMETERS = _coral_parameters(joinpath(DATA_DIR, "corals", "kayal2018_table_s1.csv"))

"""
    CORAL_RECRUITMENT_CAP

The largest number of recruits per 10 m² in a year, 348, the highest recruitment measured
for a population by Kayal et al. (2018, Appendix S3), who capped recruitment there.
"""
const CORAL_RECRUITMENT_CAP = 348.0

"""
    coral_stage(; cells = 100)

The [`ContinuousStage`](@ref) of [`coral`](@ref): colonies by size z = log10(S + 1), S the
surface area of living tissue in cm², from 0 to 7 in `cells` cells (100 in Kayal et al.
2018, Appendix S3).
"""
coral_stage(; cells::Integer = 100) = ContinuousStage(:colony, 0.0, 7.0, cells)

_coral_pool(taxon::Symbol, what::Symbol) = Symbol(taxon, :_, what)

_coral_check(name, value, allowed) = value in allowed ? value :
    throw(ArgumentError("unknown $name :$value; expected one of $(join(allowed, ", "))"))

"""
    coral_recruits(taxon, γ; reef = :lter1, parameters = CORAL_PARAMETERS,
                   reading = :as_written, cap = CORAL_RECRUITMENT_CAP)

Recruits per 10 m² in a year of `taxon` on `reef`, when its crowding is `γ` per 10 m²
(colonies of the taxon itself for *Acropora* and *Porites*, summed sizes z of the colonies of
all three taxa for *Pocillopora*). The recruitment function of Kayal et al. (2018, Appendix
S3) is log r = R0 + Rγ γ; `reading` says how its units are read:

- `:as_written` (the default): the natural logarithm, as in the authors' R script, with r and γ
  per 10 m². Read this way, *Acropora* and *Porites* hold far more of the recovering
  assemblage than the paper's figures show.
- `:calibrated`: the logarithm is to base 10, R0 is for the 25 m² plots in
  which recruits were counted (so r = 0.4 × 10^(R0 + Rγ γ/10) per 10 m²), and γ is per m².
  This reading reproduces the compositions and densities published for the simulations
  (Fig. 3, Fig. S2 and Fig. S4 of the paper); it is a reconstruction (see the vignette).

Either is capped at `cap`.
"""
function coral_recruits(taxon::Symbol, γ::Real; reef::Symbol = :lter1, parameters = CORAL_PARAMETERS,
                        reading::Symbol = :as_written, cap::Real = CORAL_RECRUITMENT_CAP)
    _coral_check("reef", reef, keys(parameters))
    _coral_check("taxon", taxon, keys(_CORAL_TAXA))
    _coral_check("reading", reading, _CORAL_READINGS)
    r = parameters[reef][taxon].recruitment
    value = reading === :calibrated ? 0.4 * 10^(r.intercept + r.slope * γ / 10) : exp(r.intercept + r.slope * γ)
    return min(Float64(cap), value)
end

"""
    coral(taxon; reef = :lter1, parameters = CORAL_PARAMETERS, reading = :as_written,
          recruit_size = (0.6, 0.2), cap = CORAL_RECRUITMENT_CAP, cells = 100, name = taxon)

A coral genus (`:acropora`, `:pocillopora` or `:porites`) on one fore reef of Mo'orea
(`:lter1`, `:lter2`, `:lter4`, `:lter5`) as the open-population integral projection model
of Kayal et al. (2018), with the parameters of [`CORAL_PARAMETERS`](@ref). Colonies are
structured by size z = log10(S + 1) ([`coral_stage`](@ref)); numbers are colonies per
10 m² and the time step is a year.

The net has one transition, `growing`: a colony of size z survives with probability
logit⁻¹(s₀ + s_z z) and grows to a size that is normal with mean G₀ + G_z z and standard
deviation G_SD. Recruits are not born of the local colonies but settle from the plankton:
each year [`coral_recruits`](@ref) of them arrive, a number set by crowding, with sizes
normal with mean and standard deviation `recruit_size` (0.6 and 0.2 in the initial runs of
the paper; 1.2 and 0.4 for *Acropora* in its later runs). They are added as inflows to the
cells of the stage.

The species adds its colonies and their summed sizes z to the pools `<taxon>_colonies` and
`<taxon>_log_size`. *Acropora* and *Porites* read their own colonies; *Pocillopora* reads
the summed sizes of all three taxa (`acropora_log_size`, `pocillopora_log_size`,
`porites_log_size`), which a community of the three supplies (see
[`coral_community`](@ref)). The environment of a step is not used.
"""
function coral(taxon::Symbol; reef::Symbol = :lter1, parameters = CORAL_PARAMETERS, reading::Symbol = :as_written,
               recruit_size = (0.6, 0.2), cap::Real = CORAL_RECRUITMENT_CAP, cells::Integer = 100, name::Symbol = taxon)
    _coral_check("reef", reef, keys(parameters))
    _coral_check("taxon", taxon, keys(_CORAL_TAXA))
    _coral_check("reading", reading, _CORAL_READINGS)
    p = parameters[reef][taxon]
    s, g = p.survival, p.growth
    colony = coral_stage(; cells)
    life = kernel_net(colony,
        :growing => (y, x) -> _logistic(s.intercept + s.slope * x) * _normal_density(y, g.intercept + g.slope * x, g.sd))
    # recruits settle in every cell: the cell's share of the recruit-size density, times the cell width
    μ, σ = recruit_size
    share = [step_size(colony) * _normal_density(z, μ, σ) for z in meshpoints(colony)]
    pools = keys(_CORAL_TAXA)
    crowding = p.recruitment.crowding === :own_density ?
        (ctx -> getproperty(ctx.pool, _coral_pool(taxon, :colonies))) :
        (ctx -> sum(getproperty(ctx.pool, _coral_pool(k, :log_size)) for k in pools))
    recruits(ctx) = coral_recruits(taxon, crowding(ctx); reef, parameters, reading, cap)
    net = inflow(life, (cell => (ctx -> share[i] * recruits(ctx)) for (i, cell) in enumerate(stage_names(colony)))...)
    occupies = [_coral_pool(taxon, :colonies) => trait_weights(colony, z -> 1.0),
                _coral_pool(taxon, :log_size) => trait_weights(colony, z -> z)]
    reads = p.recruitment.crowding === :own_density ? Symbol[] : [_coral_pool(k, :log_size) for k in pools if k != taxon]
    return species(name, net; occupies, reads)
end

const _CORAL_YEARS = (2013, 2015)

function _coral_observed(file)
    lines = [split(l, ",") for l in readlines(file)]
    header, rows = lines[1], lines[2:end]
    col(name) = findfirst(==(name), header)
    code = Dict(v => k for (k, v) in pairs(_CORAL_TAXA))
    reef_of = Dict(v => k for (k, v) in pairs(_CORAL_REEFS))
    states = Dict((r, y) => NamedTuple{keys(_CORAL_TAXA)}(Tuple(zeros(100) for _ in _CORAL_TAXA))
                  for r in keys(_CORAL_REEFS), y in _CORAL_YEARS)
    for row in rows
        reef, year = reef_of[parse(Int, row[col("reef")])], parse(Int, row[col("year")])
        taxon = code[row[col("taxon")]]
        # colonies per 10 m²: each genus was mapped over its own area (20 m² for Acropora, 4 m² for the others)
        getproperty(states[(reef, year)], taxon)[parse(Int, row[col("cell")])] +=
            parse(Int, row[col("colonies")]) * 10 / parse(Float64, row[col("sampled_area_m2")])
    end
    return states
end

const _CORAL_OBSERVED = _coral_observed(joinpath(DATA_DIR, "corals", "mcr_lter_observed_colonies.csv"))

"""
    coral_observed_state(reef, year)

The coral colonies observed on the fore reef `reef` (`:lter1`, `:lter2`, `:lter4`,
`:lter5`) of Mo'orea in `year` (2013 or 2015), binned onto the mesh of
[`coral_stage`](@ref) (100 cells), in colonies per 10 m² in each cell: a vector of
`taxon => stages` pairs, one per genus (`:acropora`, `:pocillopora`, `:porites`), which
[`project`](@ref) and `state_vector` accept as the state of a [`coral_community`](@ref)
of the default 100 cells. Kayal et al. (2018) started their updated runs from these
states (Fig. 4: from 2013 to 2015, and from 2015 to 2020):

```julia
c = coral_community(:lter1; recruit_size = (acropora = (1.2, 0.4), pocillopora = (0.6, 0.2), porites = (0.6, 0.2)))
sol = project(c, coral_observed_state(:lter1, 2015), fill(nothing, 5))   # 2020 in the last row
```

The colonies are those mapped in four permanent 5 m × 1 m transects per reef by the
Mo'orea Coral Reef LTER (Lenihan et al., EDI data package knb-lter-mcr.4009, revision 5,
CC BY 4.0; see `data/corals/LICENSE.txt`). *Acropora* was mapped over the whole of each
transect (20 m² per reef), *Pocillopora* and *Porites* over its first metre (4 m²), and
the counts are scaled to 10 m² accordingly. A colony's size is z = log10(S + 1), with its
surface area S estimated from its length, width and height as in Appendix S2 of the paper
(a sphere for *Pocillopora*, a cylinder for the other two).
"""
function coral_observed_state(reef::Symbol, year::Integer)
    _coral_check("reef", reef, keys(_CORAL_REEFS))
    year in _CORAL_YEARS || throw(ArgumentError("no observed state for $year; expected one of $(join(_CORAL_YEARS, ", "))"))
    observed = _CORAL_OBSERVED[(reef, year)]
    return [taxon => copy(getproperty(observed, taxon)) for taxon in keys(_CORAL_TAXA)]
end

function _coral_transects(file)
    lines = [split(l, ",") for l in readlines(file)]
    header, rows = lines[1], lines[2:end]
    col(name) = findfirst(==(name), header)
    code = Dict(v => k for (k, v) in pairs(_CORAL_TAXA))
    reef_of = Dict(v => k for (k, v) in pairs(_CORAL_REEFS))
    out = Dict((r, y) => NamedTuple{keys(_CORAL_TAXA)}(Tuple((colonies = Float64[], log_size = Float64[]) for _ in _CORAL_TAXA))
               for r in keys(_CORAL_REEFS), y in _CORAL_YEARS)
    for row in rows
        reef, year = reef_of[parse(Int, row[col("reef")])], parse(Int, row[col("year")])
        t = getproperty(out[(reef, year)], code[row[col("taxon")]])
        scale = 10 / parse(Float64, row[col("sampled_area_m2")])
        push!(t.colonies, parse(Int, row[col("colonies")]) * scale)
        push!(t.log_size, parse(Float64, row[col("sum_log_size")]) * scale)
    end
    return out
end

const _CORAL_TRANSECTS = _coral_transects(joinpath(DATA_DIR, "corals", "mcr_lter_observed_transects.csv"))

"""
    coral_observed_transects(reef, year)

The colonies of [`coral_observed_state`](@ref) totalled by transect: for each genus
(`:acropora`, `:pocillopora`, `:porites`), a `NamedTuple` of `colonies` (colonies per
10 m²) and `log_size` (summed sizes z per 10 m²), each a vector with one entry for each of the
four transects mapped on `reef` in `year` (2013 or 2015). Their means are the observed
abundances of Figs. S2, S5 and S7 of Kayal et al. (2018), and their standard errors those
behind the figures' confidence intervals and z-scores.
"""
function coral_observed_transects(reef::Symbol, year::Integer)
    _coral_check("reef", reef, keys(_CORAL_REEFS))
    year in _CORAL_YEARS || throw(ArgumentError("no observed transects for $year; expected one of $(join(_CORAL_YEARS, ", "))"))
    t = _CORAL_TRANSECTS[(reef, year)]
    return map(x -> (colonies = copy(x.colonies), log_size = copy(x.log_size)), t)
end

"""
    coral_community(reef = :lter1; taxa = (:acropora, :pocillopora, :porites), kwargs...)

The three coral genera of [`coral`](@ref) on one reef, composed through their pools. The
keyword arguments are passed to `coral` (`reading`, `recruit_size`, `cap`, `cells`,
`parameters`); `recruit_size` may also be a `NamedTuple` with one `(mean, sd)` per taxon.
The reef of the paper's first simulations starts empty: `zeros(length(c))`; its updated runs
start from the assemblages observed in 2013 or 2015, [`coral_observed_state`](@ref).
"""
function coral_community(reef::Symbol = :lter1; taxa = keys(_CORAL_TAXA), recruit_size = (0.6, 0.2), kwargs...)
    size_of(t) = recruit_size isa NamedTuple ? getproperty(recruit_size, t) : recruit_size
    return community((coral(t; reef, recruit_size = size_of(t), kwargs...) for t in taxa)...)
end
