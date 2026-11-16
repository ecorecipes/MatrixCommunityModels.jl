# Four perennial plants of the sagebrush steppe as integral projection models coupled by crowding.
#
# The model is the multispecies IPM of Adler, Ellner & Levine (2010) for the shrub Artemisia tripartita and
# the grasses Hesperostipa comata, Poa secunda and Pseudoroegneria spicata, mapped in permanent quadrats at
# the U.S. Sheep Experiment Station near Dubois, Idaho (the fits have effects for 22 years). The fitted parameters and the
# functions that build the kernels are those of the R package IPMdoit (Tredennick & Adler,
# https://github.com/atredennick/IPMdoit, commit 94124a4, MIT licence), whose Idaho example runs this
# model. Its recruitment coefficients are those published by Chu & Adler (2015, Appendix H); its growth and
# survival coefficients are a fit that keeps the positive interactions Chu & Adler removed.
#
# Size is the logarithm of the basal area of a genet in cm², the time step is a year and numbers are genets
# per m² (a quadrat). Each species is a net; the species meet only through four pools:
#
#   artr_cover, heco_cover, pose_cover, pssp_cover   the share of the ground covered by each species
#
# A plant is crowded by its neighbours. Crowding by another species depends only on that species' cover;
# crowding by its own species depends on the sizes of its conspecifics too, because genets of one species
# do not overlap (IPMdoit's crowd_overlap). That part is read from the species' own state.
#
# The parameters are read from data/sagebrush, copied from or extracted from IPMdoit; the notice of its
# licence is in data/sagebrush/LICENSE.txt.

const _SAGEBRUSH_CODES = (artr = "ARTR", heco = "HECO", pose = "POSE", pssp = "PSSP")
const _SAGEBRUSH_NAMES = (artr = "Artemisia tripartita", heco = "Hesperostipa comata",
                          pose = "Poa secunda", pssp = "Pseudoroegneria spicata")
# the largest genet in the data (cm²) and the number of mesh cells, from IPMdoit's Idaho example
const _SAGEBRUSH_MESH = (artr = (3000.0, 75), heco = (202.0, 75), pose = (260.0, 50), pssp = (225.0, 50))
const _SAGEBRUSH_QUADRAT = 10_000.0     # cm²: crowding and cover are computed for a quadrat of 1 m²
const _SAGEBRUSH_SMALLEST = 0.2         # cm²: the mesh starts at log(0.2)
const _SAGEBRUSH_TRUNCATION = -1.61     # recruit sizes are a normal distribution truncated below at this log area

function _sagebrush_parameters(dir)
    kinds = keys(_SAGEBRUSH_CODES)
    column(rows, header, name) = (j = findfirst(==(name), header);
        j === nothing ? throw(ArgumentError("column $name missing")) : [r[j] for r in rows])
    num(x) = parse(Float64, x)
    table(file) = (lines = [split(replace(l, "\"" => ""), ",") for l in readlines(joinpath(dir, file))];
        (lines[1], lines[2:end]))
    rh, rrows = table("recruitment.csv")
    yh, yrows = table("recruitment_years.csv")
    function one(kind)
        code = _SAGEBRUSH_CODES[kind]
        gh, g = table("Growth_params_$code.csv")
        sh, s = table("Surv_params_$code.csv")
        r = only(row for row in rrows if row[1] == code)
        rcol(name) = num(r[findfirst(==(name), rh)])
        crowding(h, rows) = [num(column(rows, h, "crowd$i")[1]) for i in 1:4]
        growth = (intercept = num(column(g, gh, "Intercept")[1]), slope = num(column(g, gh, "logarea.t0")[1]),
                  crowding = crowding(gh, g), alpha = num.(column(g, gh, "alpha")[1:4]),
                  sigma_a = num(column(g, gh, "sigma.a")[1]), sigma_b = num(column(g, gh, "sigma.b")[1]),
                  intercept_year = num.(column(g, gh, "Intercept.yr")), slope_year = num.(column(g, gh, "logarea.t0.yr")))
        survival = (intercept = num(column(s, sh, "Intercept")[1]), slope = num(column(s, sh, "logarea")[1]),
                    crowding = crowding(sh, s), alpha = num.(column(s, sh, "alpha")[1:4]),
                    intercept_year = num.(column(s, sh, "Intercept.yr")), slope_year = num.(column(s, sh, "logarea.yr")))
        recruitment = (intercept = rcol("intercept"),
                       intercept_year = num.(column(yrows, yh, code)),
                       crowding = [rcol("crowding_" * c) for c in values(_SAGEBRUSH_CODES)],
                       size_mean = rcol("size_mean"), size_variance = rcol("size_variance"))
        return (name = _SAGEBRUSH_NAMES[kind], code = code, largest = _SAGEBRUSH_MESH[kind][1],
                cells = _SAGEBRUSH_MESH[kind][2], growth = growth, survival = survival, recruitment = recruitment)
    end
    return NamedTuple{kinds}(Tuple(one(k) for k in kinds))
end

"""
    SAGEBRUSH_PARAMETERS

Fitted vital rates of four perennial plants of the sagebrush steppe in Idaho, for
[`sagebrush_plant`](@ref): the shrub *Artemisia tripartita* (`artr`) and the bunchgrasses
*Hesperostipa comata* (`heco`), *Poa secunda* (`pose`) and *Pseudoroegneria spicata*
(`pssp`). They are the parameters of the Idaho example of the R package IPMdoit (Tredennick
& Adler; MIT licence), read from `data/sagebrush`.

Size `u` is the logarithm of a genet's basal area in cm². For each species:

- `growth`: the log area next year is normal with mean `intercept + slope u + Σⱼ
  crowding[j] wⱼ` and variance `sigma_a exp(sigma_b mean)`;
- `survival`: the logit of survival is `intercept + slope u + Σⱼ crowding[j] wⱼ`;
- `recruitment`: recruits per m² are `C exp(intercept + Σⱼ crowding[j] √Cⱼ)`, where `Cⱼ` is
  the percent cover of species `j` and `C` that of the species itself; recruits have a log
  area that is normal with `size_mean` and `size_variance`, truncated below at -1.61;
- `alpha`: the spatial scale of crowding, by neighbour species, in cm⁻² (`wⱼ` sums the
  areas of neighbours of species `j`, each weighted by `exp(-alpha r²)` at distance `r`);
- `intercept_year`, `slope_year`: the effects of the 22 years of the data, used when the
  environment of a step is a year (see [`sagebrush_plant`](@ref)).

`crowding[j]` is the effect of species `j`, in the order `artr`, `heco`, `pose`, `pssp`.
Also `largest`, the largest genet in the data (cm²), and `cells`, the number of mesh cells.
"""
const SAGEBRUSH_PARAMETERS = _sagebrush_parameters(joinpath(DATA_DIR, "sagebrush"))

_sagebrush_species(parameters, kind::Symbol) = haskey(parameters, kind) ? parameters[kind] :
    throw(ArgumentError("unknown plant :$kind; the parameters hold $(join(keys(parameters), ", "))"))

_cover_pool(kind::Symbol) = Symbol(kind, :_cover)

"""
    sagebrush_stage(kind; parameters = SAGEBRUSH_PARAMETERS, cells = parameters[kind].cells)

The [`ContinuousStage`](@ref) of [`sagebrush_plant`](@ref): genets of species `kind` by the
logarithm of their basal area in cm², from log 0.2 to 1.1 times the logarithm of the
largest genet in the data, the mesh of IPMdoit (75 cells for `artr` and `heco`, 50 for
`pose` and `pssp`).
"""
function sagebrush_stage(kind::Symbol; parameters = SAGEBRUSH_PARAMETERS, cells::Integer = _sagebrush_species(parameters, kind).cells)
    p = _sagebrush_species(parameters, kind)
    return ContinuousStage(:plant, log(_SAGEBRUSH_SMALLEST), 1.1 * log(p.largest), cells)
end

# The natural cubic spline through (x, y), as R's splinefun(method = "natural") makes it: second derivative
# zero at both ends, linear beyond them. Returns its coefficients on each interval.
function _natural_spline(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    d = diff(x)
    slope = diff(y) ./ d
    # second derivatives M at the knots, M[1] = M[n] = 0, from the tridiagonal system
    M = zeros(n)
    if n > 2
        T = SymTridiagonal(2 .* (d[1:(end - 1)] .+ d[2:end]), d[2:(end - 1)])
        M[2:(n - 1)] = T \ (6 .* diff(slope))
    end
    return (x = x, y = y, d = d, slope = slope, M = M)
end

# value of the spline at t (t at or above x[1] - anything below is extrapolated linearly)
function _spline_value(s, t::Float64, i::Int)
    if i == 0       # left of the first knot: linear, with the slope of the spline there
        b = s.slope[1] - s.d[1] * (2 * s.M[1] + s.M[2]) / 6
        return s.y[1] + b * (t - s.x[1])
    end
    a, h = t - s.x[i], s.d[i]
    return s.y[i] + a * (s.slope[i] - h * (2 * s.M[i] + s.M[i + 1]) / 6) + a^2 * s.M[i] / 2 +
           a^3 * (s.M[i + 1] - s.M[i]) / (6h)
end

function _gauss_legendre(n::Int)
    F = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]))
    return F.values, 2 .* F.vectors[1, :] .^ 2
end
const _GL_NODES, _GL_WEIGHTS = _gauss_legendre(8)

# Crowding of a genet of radius r by the other genets of its own species, which cannot overlap it: the
# function wrijG/wrijS of IPMdoit. `C` is the spline of the cover (cm² per
# quadrat) of the conspecifics with radius below ρ, `reach` the radius of the largest genet of the mesh,
# `total` their whole cover. Integrated by Gauss-Legendre on each piece of the spline, where IPMdoit uses
# R's integrate.
function _own_crowding(r::Float64, alpha::Float64, C, reach::Float64, total::Float64)
    f(ρ, i) = 2pi * (r + ρ) * exp(-alpha * (r + ρ)^2) * _spline_value(C, ρ, i)
    edges = [0.0; C.x]
    integral = 0.0
    for i in 0:(length(C.x) - 1)
        lo, hi = edges[i + 1], edges[i + 2]
        hi > lo || continue
        mid, half = (lo + hi) / 2, (hi - lo) / 2
        integral += half * sum(w * f(mid + half * t, i) for (t, w) in zip(_GL_NODES, _GL_WEIGHTS))
    end
    return (integral + pi * total * exp(-alpha * (r + reach)^2) / alpha) / _SAGEBRUSH_QUADRAT
end

# The coefficients of the neighbours of species i, with or without niche differences. Without them a species is
# affected by the cover of every species as by its own: in recruitment, which reads the covers, every coefficient
# becomes its own; in growth and survival, where crowding by species j is π Cⱼ / alpha[j] on the scale of j, the
# coefficient of j becomes the species' own times alpha[j] / alpha[i].
function _sagebrush_effects(effects, i, niche::Bool; alpha = nothing, interactions::Bool = true)
    β = if niche
        collect(Float64, effects)
    else
        scale = alpha === nothing ? ones(length(effects)) : alpha ./ alpha[i]
        Float64(effects[i]) .* scale
    end
    # without interspecific interactions a species responds to its own species alone
    interactions || (β = [j == i ? b : 0.0 for (j, b) in enumerate(β)])
    return β
end

"""
    sagebrush_crowding(kind, plants, cover; parameters = SAGEBRUSH_PARAMETERS, cells = ...)

The crowding experienced by genets of species `kind` of each size: a `NamedTuple` with
`growth` and `survival`, each a matrix with one row for each cell of
[`sagebrush_stage`](@ref) and one column for each neighbour species (`artr`, `heco`,
`pose`, `pssp`). `plants` are the numbers of the species' own genets per m² in the cells
of its stage, and `cover` the cover of every species (a `NamedTuple` with fields
`artr_cover`, ..., as the pools of [`sagebrush_plant`](@ref)).

Crowding by another species `j` is `π Cⱼ / alpha`, with `Cⱼ` its cover as a share of the
ground: its genets are spread at random and may overlap the focal genet. Crowding by the
species' own genets, which do not overlap the focal genet, depends on their sizes and on
the focal genet's radius. Both follow IPMdoit (`crowd_overlap`).
"""
function sagebrush_crowding(kind::Symbol, plants::AbstractVector{<:Real}, cover; parameters = SAGEBRUSH_PARAMETERS,
                            cells::Integer = _sagebrush_species(parameters, kind).cells)
    stage = sagebrush_stage(kind; parameters, cells)
    return _sagebrush_crowding(kind, stage, collect(Float64, plants), cover, parameters)
end

function _sagebrush_crowding(kind::Symbol, stage::ContinuousStage, plants::Vector{Float64}, cover, parameters)
    kinds = keys(parameters)
    i = findfirst(==(kind), kinds)
    p = parameters[kind]
    z, h = meshpoints(stage), step_size(stage)
    length(plants) == length(z) || throw(DimensionMismatch("need one number of genets for each of the $(length(z)) cells"))
    lower, upper = stage.domain.lower, stage.domain.upper
    edges = sqrt.(exp.(range(lower, upper; length = length(z) + 1)) ./ pi)       # radii at the cell edges
    held = exp.(z) .* plants                                                      # cm² in each cell, per m²
    C = _natural_spline(collect(edges), [0.0; cumsum(held)])
    reach, total = edges[end], sum(held)
    radius = sqrt.(exp.(z) ./ pi)
    function matrix(alpha)
        W = zeros(length(z), length(kinds))
        for (j, other) in enumerate(kinds)
            W[:, j] .= j == i ? [_own_crowding(r, alpha[j], C, reach, total) for r in radius] :
                                pi * getproperty(cover, _cover_pool(other)) / alpha[j]
        end
        return W
    end
    return (growth = matrix(p.growth.alpha), survival = matrix(p.survival.alpha))
end

const _SAGEBRUSH_RATES = (:growth, :survival, :recruitment)

# the vital rates that keep their niche differences
function _sagebrush_kept(niche)
    niche === true && return _SAGEBRUSH_RATES
    niche === false && return ()
    (niche isa Union{Tuple, AbstractVector} && all(in(_SAGEBRUSH_RATES), niche)) && return Tuple(niche)
    throw(ArgumentError("niche_differences is true, false or a list of vital rates among $(_SAGEBRUSH_RATES), got $niche"))
end

_sagebrush_year(env, n::Int) = env === nothing ? 0 :
    (env isa Integer && 1 <= env <= n) ? Int(env) :
    throw(ArgumentError("the environment of a sagebrush step is nothing (the mean year) or a year from 1 to $n, got $env"))

"""
    sagebrush_plant(kind; parameters = SAGEBRUSH_PARAMETERS, cells = parameters[kind].cells,
                      niche_differences = true, name = kind)

A perennial plant of the sagebrush steppe (`:artr`, `:heco`, `:pose` or `:pssp`; see
[`SAGEBRUSH_PARAMETERS`](@ref)) as the integral projection model of Adler, Ellner & Levine
(2010), with the parameters and numerical scheme of the R package IPMdoit. Genets are
structured by the logarithm of their area (see [`sagebrush_stage`](@ref)); numbers are
genets per m² and the time step is a year.

The species adds its basal area, as a share of the ground, to the pool `<kind>_cover`
(`artr_cover`, ...) and reads the cover of the other three. It has two transitions:

- `growing`: a genet survives and grows. Both depend on the crowding of the genet by
  neighbours of each species ([`sagebrush_crowding`](@ref)), which depends on the genet's
  size, so the kernel is rebuilt each year ([`responsive`](@ref)).
- `recruiting`: the species' recruits, `C exp(intercept + Σⱼ crowding[j] √Cⱼ)` per m² for
  percent covers `C`, are shared among its genets in proportion to their area. They are
  marked as births.

The environment of a time step is `nothing`, for the mean year (the random year effects
set to zero, and the recruitment intercept at its posterior mean), or a year `1` to `22`
of the data, for that year's effects on growth, survival and recruitment.

With `niche_differences = false` every species is affected by the cover of each species as
by its own, which is an assumption of this package (see the vignette). In recruitment,
which reads the covers, the coefficient of each neighbour species becomes the species'
coefficient on itself. In growth and survival, where crowding by species `j` is measured
on the spatial scale `alpha[j]`, it becomes the species' own coefficient times
`alpha[j] / alpha[i]`, so that a unit of cover of any species, placed at random, crowds a
plant as a unit of its own species' cover would at its own scale. The rule that
conspecific genets do not overlap is kept. `niche_differences` can also list the vital
rates that keep their niche differences, for example `(:growth, :survival)` to remove them
from recruitment alone.

With `interactions = false` every coefficient of another species is set to zero, in all
three vital rates, so that the species responds to its own crowding and cover alone (the
"no competition" simulations of Tredennick et al. 2017).

`parameters` can be any set in the form of `SAGEBRUSH_PARAMETERS`, such as a site of
[`GRASSLAND_PARAMETERS`](@ref); `kind` is then one of its species.
"""
function sagebrush_plant(kind::Symbol; parameters = SAGEBRUSH_PARAMETERS,
                         cells::Integer = _sagebrush_species(parameters, kind).cells,
                         niche_differences = true, interactions::Bool = true, name::Symbol = kind)
    p = _sagebrush_species(parameters, kind)
    kinds = collect(keys(parameters))
    i = findfirst(==(kind), kinds)
    plant = sagebrush_stage(kind; parameters, cells)
    g, s, r = p.growth, p.survival, p.recruitment
    kept = _sagebrush_kept(niche_differences)
    βg = _sagebrush_effects(g.crowding, i, :growth in kept; alpha = g.alpha, interactions)
    βs = _sagebrush_effects(s.crowding, i, :survival in kept; alpha = s.alpha, interactions)
    βr = _sagebrush_effects(r.crowding, i, :recruitment in kept; interactions)
    years = length(r.intercept_year)
    z = meshpoints(plant)
    sd = sqrt(r.size_variance)
    # recruit sizes: a normal density truncated below at log area -1.61, as in IPMdoit
    recruit_size(y) = _normal_density(y, r.size_mean, sd) / (1 - _normal_cdf((_SAGEBRUSH_TRUNCATION - r.size_mean) / sd))

    function growing(ctx)
        k = _sagebrush_year(ctx.env, years)
        W = _sagebrush_crowding(kind, plant, ctx.n[plant], ctx.pool, parameters)
        a = g.intercept + (k == 0 ? 0.0 : g.intercept_year[k])
        b = g.slope + (k == 0 ? 0.0 : g.slope_year[k])
        c = s.intercept + (k == 0 ? 0.0 : s.intercept_year[k])
        d = s.slope + (k == 0 ? 0.0 : s.slope_year[k])
        cell = Dict(x => j for (j, x) in enumerate(z))
        mean_size = a .+ b .* z .+ W.growth * βg
        alive = _logistic.(c .+ d .* z .+ W.survival * βs)
        return function (y, x)
            j = cell[x]
            μ = mean_size[j]
            return alive[j] * _normal_density(y, μ, sqrt(g.sigma_a * exp(g.sigma_b * μ)))
        end
    end
    life = kernel_net(plant,
        :growing => responsive(growing),
        # recruits for each cm² of the species' basal area; the response below supplies their number
        :recruiting => (y, x) -> exp(x) * recruit_size(y))

    # recruits per m² per cm² of the species' basal area: IPMdoit's get_rpa, written without dividing by the
    # species' cover, so that it also holds for a species without genets
    function recruiting(ctx)
        k = _sagebrush_year(ctx.env, years)
        percent = [100 * getproperty(ctx.pool, _cover_pool(other)) for other in kinds]
        return 100 * exp((k == 0 ? r.intercept : r.intercept_year[k]) + sum(βr .* sqrt.(max.(percent, 0.0)))) / _SAGEBRUSH_QUADRAT
    end
    net = births(respond(life, :recruiting => recruiting), :recruiting)
    own = _cover_pool(kind)
    return species(name, net; occupies = [own => trait_weights(plant, u -> exp(u) / _SAGEBRUSH_QUADRAT)],
                   reads = [_cover_pool(k) for k in kinds if k != kind])
end

"""
    sagebrush_community(; parameters = SAGEBRUSH_PARAMETERS, niche_differences = true,
                          interactions = true, kinds = keys(parameters))

The four plants of [`sagebrush_plant`](@ref) composed through their cover pools, or those
named in `kinds`.
"""
sagebrush_community(; parameters = SAGEBRUSH_PARAMETERS, niche_differences = true, interactions::Bool = true,
                    kinds = keys(parameters)) =
    community((sagebrush_plant(k; parameters, niche_differences, interactions) for k in kinds)...)

"""
    sagebrush_state(kind; density = 0.001, cover = nothing, parameters = SAGEBRUSH_PARAMETERS,
                      cells = parameters[kind].cells)

A stage vector for `sagebrush_plant(kind; cells)`: genets spread evenly over the log areas
of the mesh, `density` genets per m² for each unit of log area, which is how the Idaho
example of IPMdoit starts its runs. With `cover` (a share of the ground) the same shape is
scaled to that cover.
"""
function sagebrush_state(kind::Symbol; density::Real = 0.001, cover = nothing, parameters = SAGEBRUSH_PARAMETERS,
                         cells::Integer = _sagebrush_species(parameters, kind).cells)
    stage = sagebrush_stage(kind; parameters, cells)
    n = fill(Float64(density) * step_size(stage), length(stage))
    cover === nothing && return n
    return n .* (cover / sum(exp.(meshpoints(stage)) .* n ./ _SAGEBRUSH_QUADRAT))
end
