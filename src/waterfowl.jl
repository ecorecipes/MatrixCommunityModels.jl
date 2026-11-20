# Canvasbacks and redheads: the process model of the integrated population model of Gibson, Arnold, Buderman &
# Koons (2025, Ecological Monographs 95:e1642), two diving ducks that share the breeding ponds of the Prairie
# Pothole Region; redheads lay eggs in canvasback nests. Each species is counted at the May breeding survey in four
# classes, yearling (second-year, SY) and older (ASY) females and males. One step is a year:
#
#   ducklings      F_t (N_ASY,F + ω N_SY,F), half of each sex
#   yearlings      ducklings × natural winter survival × (1 - harvest) of hatch-years × √(summer survival, t+1)
#   older birds    √(summer survival, t) × natural winter survival × (1 - harvest)^r × √(summer survival, t+1)
#                  × (yearlings + older birds of that sex), where r = exp(β_indirect) turns the direct harvest
#                  hazard of adults into the indirect one (IPM_code.R)
#
# The published model fitted the yearly rates with covariates; the densities and the canvasback-to-redhead
# ratio that drive fecundity entered as data, not as the projected numbers. Here (a) the yearly posterior means
# are the environment of each step and the two species do not meet in the projection, as in the fitted model;
# (b) an extension, not in the paper, closes the loop: fecundity responds to pools of the females of both species,
# through their departure from the posterior-mean female numbers of the fit.
#
#   canvasback_females, redhead_females   breeding females (SY + ASY) of each species, millions
#
# The equations are those of the authors' code, IPM_code.R (Zenodo doi:10.5281/zenodo.8076163, MIT); the transient
# LTRE transcribes LTRE_Code.R. The yearly rates are posterior means of the draws on Dryad
# (doi:10.5061/dryad.hqbzkh1n9, CC0), made by test/reference/R/run_waterfowl_reference.R (data/waterfowl).

const _WATERFOWL_DIR = joinpath(DATA_DIR, "waterfowl")
const _WATERFOWL_RATES = (:fecundity, :sw_hy_f, :sw_ahy_f, :sw_hy_m, :sw_ahy_m, :kappa_hy_f, :kappa_ahy_f,
                          :kappa_hy_m, :kappa_ahy_m, :ss_f, :ss_m)

function _waterfowl_table(file)
    lines = readlines(file)
    header = Symbol.(split(lines[1], ","))
    rows = [split(l, ",") for l in lines[2:end] if !isempty(l)]
    return header, rows
end
_waterfowl_number(s) = s == "NA" ? NaN : parse(Float64, s)

# the yearly rates of each species from a table with columns species, year and _WATERFOWL_RATES (rows of `pick`)
function _waterfowl_rates(header, rows; pick = r -> true)
    col(name) = findfirst(==(name), header)
    out = map((:canvasback, :redhead)) do sp
        rs = [r for r in rows if r[col(:species)] == String(sp) && pick(r)]
        sp => NamedTuple{_WATERFOWL_RATES}(Tuple([_waterfowl_number(r[col(v)]) for r in rs] for v in _WATERFOWL_RATES))
    end
    return (; out...)
end

function _waterfowl_parameters(dir)
    vh, vrows = _waterfowl_table(joinpath(dir, "vital_rates.csv"))
    ph, prows = _waterfowl_table(joinpath(dir, "parameters.csv"))
    ah, arows = _waterfowl_table(joinpath(dir, "abundance.csv"))
    ch, crows = _waterfowl_table(joinpath(dir, "covariates.csv"))
    pcol(sp, name) = _waterfowl_number(prows[findfirst(r -> r[1] == String(sp), prows)][findfirst(==(name), ph)])
    bysp(name) = (canvasback = pcol(:canvasback, name), redhead = pcol(:redhead, name))
    acol(name) = findfirst(==(name), ah)
    abundance(sp) = (rs = [r for r in arows if r[1] == String(sp)];
                     (year = [parse(Int, r[2]) for r in rs],
                      n = [_waterfowl_number(r[acol(c)]) for r in rs, c in (:sy_f, :asy_f, :sy_m, :asy_m)],
                      total_q05 = [_waterfowl_number(r[acol(:total_q05)]) for r in rs],
                      total_q95 = [_waterfowl_number(r[acol(:total_q95)]) for r in rs]))
    ccol(name) = [_waterfowl_number(r[findfirst(==(name), ch)]) for r in crows]
    ab = (canvasback = abundance(:canvasback), redhead = abundance(:redhead))
    return (species = (:canvasback, :redhead),
            years = 1961:2020,
            rates = _waterfowl_rates(vh, vrows),
            omega = bysp(:omega),
            indirect_harvest = bysp(:indirect_harvest),
            bp_ltre = bysp(:bp_ltre),
            n0 = (canvasback = ab.canvasback.n[1, :], redhead = ab.redhead.n[1, :]),
            abundance = ab,
            covariates = (year = Int.(ccol(:year)),
                          females = (canvasback = ccol(:density_canvasback), redhead = ccol(:density_redhead)),
                          ratio = ccol(:ratio),
                          bpop = (canvasback = ccol(:bpop_canvasback), redhead = ccol(:bpop_redhead))),
            density_mean = bysp(:density_mean), density_sd = bysp(:density_sd),
            ratio_mean = pcol(:canvasback, :ratio_mean), ratio_sd = pcol(:canvasback, :ratio_sd),
            beta_density = (canvasback = -0.12, redhead = -0.09),
            beta_ratio = (canvasback = 0.15, redhead = -0.16))
end

"""
    WATERFOWL_GIBSON_2025

The canvasback and redhead model of Gibson, Arnold, Buderman & Koons (2025, *Ecol. Monogr.* 95:e1642) at the
posterior means of the authors' 8000 draws (Dryad doi:10.5061/dryad.hqbzkh1n9, CC0; see `data/waterfowl`):

- `rates.<species>`: for 1961 to 2020, the yearly per-female fecundity (ducklings of both sexes per older
  female) and, for females (`_f`) and males (`_m`), natural winter survival of hatch-years (`sw_hy`) and adults
  (`sw_ahy`), direct harvest probability (`kappa_hy`, `kappa_ahy`) and natural summer survival (`ss`).
- `omega`: breeding success of yearling relative to older females, and `indirect_harvest`, the ratio
  `exp(β_indirect)` of the indirect to the direct harvest hazard of adults. The draws do not store these; each draw's
  value is recovered from its abundances and the posterior mean taken.
- `bp_ltre`: the mean of the `bp` column used by the authors' LTRE script (it is not aligned with the draws).
- `n0`, `abundance`: posterior means of the breeding numbers (millions) of yearling and older females and males,
  1961 to 2021, with 5% and 95% quantiles of the total.
- `covariates`: the breeding female numbers and the canvasback-to-redhead female ratio that the fitted model took
  as data, with the means and standard deviations it standardised them by (`density_mean`, `density_sd`,
  `ratio_mean`, `ratio_sd`), and the breeding survey estimates (`bpop`), 1961 to 2021 (no survey in 2020 and 2021).
- `beta_density`, `beta_ratio`: the effects of conspecific female density and of female canvasbacks per female
  redhead on log fecundity, as printed in Appendix S1, Table S6 (used only by the feedback extension).
"""
const WATERFOWL_GIBSON_2025 = _waterfowl_parameters(_WATERFOWL_DIR)

"""
    waterfowl_environment(parameters = WATERFOWL_GIBSON_2025; years = 1961:2019)

The environment of each yearly step from May of `year` to May of `year + 1`: a named tuple with `year` and, for
each species, its rates of that year (see [`WATERFOWL_GIBSON_2025`](@ref)) with `ss_f_next` and `ss_m_next`, the
summer survival of the next year, and `reference`, the posterior-mean numbers of breeding females of both species
that year (millions), from which the feedback extension measures its departures. The rates end in 2020, so the
last step starts in 2019.
"""
function waterfowl_environment(parameters = WATERFOWL_GIBSON_2025; years = 1961:2019)
    p = parameters
    first(years) >= 1961 && last(years) <= 2019 ||
        throw(ArgumentError("the rates cover 1961 to 2020, so steps can start from 1961 to 2019"))
    return map(years) do y
        i = y - 1960
        sp(r) = merge(NamedTuple{_WATERFOWL_RATES}(Tuple(getproperty(r, v)[i] for v in _WATERFOWL_RATES)),
                      (ss_f_next = r.ss_f[i + 1], ss_m_next = r.ss_m[i + 1]))
        females(a) = a.n[i, 1] + a.n[i, 2]
        (year = y, canvasback = sp(p.rates.canvasback), redhead = sp(p.rates.redhead),
         reference = (canvasback_females = females(p.abundance.canvasback),
                      redhead_females = females(p.abundance.redhead)))
    end
end

"""
    waterfowl_duck(name; omega, indirect_harvest, feedback = nothing)

One duck species of the model of Gibson et al. (2025), with stages `sy_female`, `asy_female`, `sy_male` and
`asy_male` (breeding yearlings and older birds, counted in May) and transitions `female_recruitment` and
`male_recruitment` (females to yearlings of each sex: half the ducklings, `omega` per yearling female and one per
older female, marked as births) and `female_survival` and `male_survival` (yearlings and older birds to older
birds). The rates of each year come from the environment of the step, the field `name` of
[`waterfowl_environment`](@ref): recruitment is fecundity times natural winter survival times `1 - harvest` of
hatch-years times the square root of next year's summer survival; survival is the square root of this summer's
survival times natural winter survival times `(1 - harvest)^indirect_harvest` times the square root of next
year's summer survival.

With `feedback = (beta_density, beta_ratio, density_sd, ratio_sd)`, an extension that is not the published model,
fecundity is multiplied by `exp(β_N (N - N_ref) / sd_N + β_C/R (C/R - C_ref/R_ref) / sd_C/R)`: `N` is the number
of breeding females of this species, read from the pool `<name>_females`, `C/R` the ratio of the pools
`canvasback_females` and `redhead_females`, and `N_ref`, `C_ref`, `R_ref` the posterior-mean numbers of the fit that
year (`reference` in the environment). The slopes are those of the fitted fecundity model, per standard deviation
of the covariate, so the published fecundity is recovered when the projection follows the fit. The species then
occupies `<name>_females` and reads the other species' pool.
"""
function waterfowl_duck(name::Symbol; omega::Real, indirect_harvest::Real, feedback = nothing)
    life = ValuedProjectionNet([:sy_female, :asy_female, :sy_male, :asy_male],
        :female_recruitment => [(:sy_female => :sy_female) => 0.5 * omega, (:asy_female => :sy_female) => 0.5],
        :male_recruitment => [(:sy_female => :sy_male) => 0.5 * omega, (:asy_female => :sy_male) => 0.5],
        :female_survival => [(:sy_female => :asy_female) => 1.0, (:asy_female => :asy_female) => 1.0],
        :male_survival => [(:sy_male => :asy_male) => 1.0, (:asy_male => :asy_male) => 1.0])
    r = indirect_harvest
    rates(ctx) = getproperty(ctx.env, name)
    fecundity(ctx) = rates(ctx).fecundity * _waterfowl_feedback(ctx, name, feedback)
    net = respond(life,
        :female_recruitment => ctx -> (e = rates(ctx);
            fecundity(ctx) * e.sw_hy_f * (1 - e.kappa_hy_f) * sqrt(e.ss_f_next)),
        :male_recruitment => ctx -> (e = rates(ctx);
            fecundity(ctx) * e.sw_hy_m * (1 - e.kappa_hy_m) * sqrt(e.ss_m_next)),
        :female_survival => ctx -> (e = rates(ctx);
            sqrt(e.ss_f) * e.sw_ahy_f * (1 - e.kappa_ahy_f)^r * sqrt(e.ss_f_next)),
        :male_survival => ctx -> (e = rates(ctx);
            sqrt(e.ss_m) * e.sw_ahy_m * (1 - e.kappa_ahy_m)^r * sqrt(e.ss_m_next)))
    net = births(net, :female_recruitment, :male_recruitment)
    feedback === nothing && return species(name, net)
    own = Symbol(name, :_females)
    other = name === :canvasback ? :redhead_females : :canvasback_females
    return species(name, net; occupies = [own => [:sy_female => 1.0, :asy_female => 1.0]], reads = [other])
end

_waterfowl_feedback(ctx, name, ::Nothing) = 1.0
function _waterfowl_feedback(ctx, name, f)
    ref = ctx.env.reference
    zN = (getproperty(ctx.pool, Symbol(name, :_females)) - getproperty(ref, Symbol(name, :_females))) / f.density_sd
    zR = (ctx.pool.canvasback_females / ctx.pool.redhead_females - ref.canvasback_females / ref.redhead_females) /
         f.ratio_sd
    return exp(f.beta_density * zN + f.beta_ratio * zR)
end

"""
    waterfowl_community(parameters = WATERFOWL_GIBSON_2025; feedback = false)

Canvasbacks and redheads ([`waterfowl_duck`](@ref)) at the parameters of [`WATERFOWL_GIBSON_2025`](@ref). Project
through [`waterfowl_environment`](@ref). With `feedback = false` (the published model) the species do not meet in
the projection: their yearly rates, which the paper fitted with the densities and the canvasback-to-redhead ratio
as data, are the environment. With `feedback = true`, an extension that is not the published model, fecundity
responds to the projected numbers of breeding females of both species through the pools `canvasback_females` and
`redhead_females`, with the slopes of Appendix S1, Table S6 (see [`waterfowl_duck`](@ref)). Change the slopes with
`merge`, e.g. `merge(WATERFOWL_GIBSON_2025, (beta_ratio = (canvasback = 0.0, redhead = 0.0),))` for conspecific
density dependence alone.

Returns `(community, n0)`, `n0` the posterior-mean breeding numbers of 1961 (millions).
"""
function waterfowl_community(parameters = WATERFOWL_GIBSON_2025; feedback::Bool = false)
    p = parameters
    duck(sp) = waterfowl_duck(sp; omega = getproperty(p.omega, sp), indirect_harvest = getproperty(p.indirect_harvest, sp),
        feedback = feedback ? (beta_density = getproperty(p.beta_density, sp), beta_ratio = getproperty(p.beta_ratio, sp),
                               density_sd = getproperty(p.density_sd, sp), ratio_sd = p.ratio_sd) : nothing)
    com = community(duck(:canvasback), duck(:redhead))
    return (community = com, n0 = [:canvasback => p.n0.canvasback, :redhead => p.n0.redhead])
end

"""
    waterfowl_ltre(parameters = WATERFOWL_GIBSON_2025)

The transient life table response experiment of Gibson et al. (2025, Equation 35), transcribed from the authors'
LTRE_Code.R and evaluated on the posterior means instead of each draw. For each species, the 18 yearly series of
the authors' script (fecundity of yearling and older females `F1`, `F2`; female and male hatch-year and adult
winter (`_w`), summer (`_s`) and harvest (`_h`) survival; the proportions of the breeding population in each
class, `nf1`, `nf2`, `nm1`, `nm2`) are formed as the script forms them, sensitivities of the realised growth rate
are evaluated at their temporal means with the script's formulas, and the contribution of each series to the
variance of the growth rate is the sum over the others of covariance times both sensitivities.

Returns a named tuple per species with `names` (the 18 series), `contribution`, `sensitivity`, `total` and
`lambda`, the geometric mean growth rate of the whole breeding population. On the posterior means the
contributions differ from the paper's, which average this calculation over the draws (vignette 28).
"""
function waterfowl_ltre(parameters = WATERFOWL_GIBSON_2025)
    p = parameters
    out = map(p.species) do sp
        r = getproperty(p.rates, sp); n = getproperty(p.abundance, sp).n
        T = length(r.fecundity)                                   # 60 years, 1961 to 2020
        tot = vec(sum(n; dims = 2))
        prop = n[1:T, :] ./ tot[1:T]
        ms, msm = sum(r.ss_f) / T, sum(r.ss_m) / T
        adult_summer(ss, m) = [t == 1 ? m^(1 / 2) * ss[1]^(1 / 2) : ss[t - 1]^(1 / 2) * ss[t]^(1 / 2) for t in 1:T]
        series = hcat(r.fecundity .* getproperty(p.bp_ltre, sp), r.fecundity,
                      r.sw_hy_f, r.sw_ahy_f, r.ss_f .^ (1 / 2), adult_summer(r.ss_f, ms),
                      1 .- r.kappa_hy_f, 1 .- r.kappa_ahy_f, prop[:, 1], prop[:, 2],
                      r.sw_hy_m, r.sw_ahy_m, r.ss_m .^ (1 / 2), adult_summer(r.ss_m, msm),
                      1 .- r.kappa_hy_m, 1 .- r.kappa_ahy_m, prop[:, 3], prop[:, 4])
        μ = vec(sum(series; dims = 1)) ./ T
        F1, F2, Jw, Aw, Js, As, Jh, Ah, nf1, nf2, Mjw, Maw, Mjs, Mas, Mjh, Mah, nm1, nm2 = μ
        psr = 0.5
        JF, JM, AF, AM = Jw * Jh * Js, Mjw * Mjh * Mjs, Aw * Ah * As, Maw * Mah * Mas
        Fn = F1 * nf1 + F2 * nf2
        Λ = psr * JF * Fn + (1 - psr) * JM * Fn + AF + AM
        sens = [psr * JF * nf1 + (1 - psr) * JM * nf1, psr * JF * nf2 + (1 - psr) * JM * nf2,
                psr * Fn * Jh * Js, As * Ah, psr * Fn * Jh * Jw, Aw * Ah, psr * Fn * Jw * Js, As * Aw,
                (psr * F1 * JF + (1 - psr) * F1 * JM + AF) - Λ, (psr * F2 * JF + (1 - psr) * F2 * JM + AF) - Λ,
                (1 - psr) * Fn * Mjh * Mjs, Mas * Mah, (1 - psr) * Fn * Mjh * Mjw, Maw * Mah,
                (1 - psr) * Fn * Mjw * Mjs, Mas * Maw, AM - Λ, AM - Λ]
        dev = series .- μ'
        C = (dev' * dev) ./ (T - 1)
        contribution = [sum(C[l, m] * sens[l] * sens[m] for m in 1:18) for l in 1:18]
        sp => (names = [:F1, :F2, :Sfj_w, :Sfa_w, :Sfj_s, :Sfa_s, :Sfj_h, :Sfa_h, :nf1, :nf2,
                        :Smj_w, :Sma_w, :Smj_s, :Sma_s, :Smj_h, :Sma_h, :nm1, :nm2],
               contribution = contribution, sensitivity = sens, total = sum(contribution),
               lambda = (tot[T + 1] / tot[1])^(1 / T))
    end
    return (; out...)
end
