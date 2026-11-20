# Annual plants with a seed bank competing for one season: the model of Godoy, Kraft & Levine (2014, Ecology
# Letters 17:836-844) for 18 annual species of a California serpentine grassland, with the vital rates and
# interaction coefficients the authors deposited on figshare (data/annual_plants; see LICENSE.txt there).
#
# The paper's model counts the seeds N in the soil before germination:
#
#   N[i, t+1] = (1 - g[i]) s[i] N[i, t] + g[i] lambda[i] N[i, t] / (1 + sum_j alpha[i, j] g[j] N[j, t])
#
# Here each species has two stages counted just after germination: `seeds`, the (1 - g) N seeds that stay
# dormant in the soil, and `plants`, the g N germinants. A germinant sets lambda / (1 + crowding) seeds, and
# every seed in the soil next autumn, old or new, germinates with probability g; a dormant seed first survives
# the year with probability s. The two stages sum to the paper's N, step for step.
#
# The species meet through one pool per species, `<code>_crowding`: the germinants of every species j weighted
# by alpha[i, j], the crowding felt by a germinant of species i.

const _ANNUAL_PLANT_NAMES = Dict(
    "AGHE" => "Agoseris heterophylla", "AGRE" => "Agoseris retrorsa", "AMME" => "Amsinckia menziesii",
    "ANAR" => "Anagallis arvensis", "CEME" => "Centaurea melitensis", "CLPU" => "Clarkia purpurea",
    "ERBO" => "Erodium botrys", "ERCI" => "Erodium cicutarium", "EUPE" => "Euphorbia peplus",
    "GECA" => "Geranium carolinianum", "HECO" => "Hemizonia congesta", "LACA" => "Lasthenia californica",
    "LOPU" => "Lotus purshianus", "LOWR" => "Lotus wrangelianus", "MEPO" => "Medicago polymorpha",
    "NAAT" => "Navarretia atractyloides", "PLER" => "Plantago erecta", "SACA" => "Salvia columbariae")

function _annual_plant_parameters(dir)
    table(file) = [split(l, ",") for l in readlines(joinpath(dir, file)) if !isempty(l)]
    vital = table("species_vital_rates.csv")
    head = vital[1]
    col(name) = findfirst(==(name), head)
    codes = [String(r[1]) for r in vital[2:end]]
    num(x) = x == "NA" ? NaN : parse(Float64, x)
    alpha_rows = table("alpha_estimates_row_is_target.csv")
    [String(x) for x in alpha_rows[1][2:end]] == codes || error("alpha columns do not match the species")
    [String(r[1]) for r in alpha_rows[2:end]] == codes || error("alpha rows do not match the species")
    alpha = permutedims(reduce(hcat, [num.(r[2:end]) for r in alpha_rows[2:end]]))
    common = table("godoy2014_appendix_s2_table2.csv")
    [String(r[1]) for r in common[2:end]] == codes || error("Table 2 rows do not match the species")
    return (codes = Tuple(Symbol(lowercase(c)) for c in codes),
            names = Tuple(_ANNUAL_PLANT_NAMES[c] for c in codes),
            lambda = [num(r[col("lambda")]) for r in vital[2:end]],
            g = [num(r[col("g")]) for r in vital[2:end]],
            s = [num(r[col("s")]) for r in vital[2:end]],
            alpha = alpha,
            effect = [num(r[2]) for r in common[2:end]],
            response = [num(r[4]) for r in common[2:end]])
end

"""
    ANNUAL_PLANTS_GODOY_2014

Parameters of the annual plant model of Godoy, Kraft & Levine (2014) for 18 annual species
of the serpentine grassland of the Sedgwick Reserve, California, read from
`data/annual_plants` (figshare doi:10.6084/m9.figshare.4793488, CC0):

- `codes`: the species as lower-case four-letter codes, `:aghe` (*Agoseris heterophylla*) to
  `:saca` (*Salvia columbariae*); `names` their Latin names;
- `lambda`: viable seeds set per germinant without competitors;
- `g`: the fraction of seeds in the soil that germinate;
- `s`: the annual survival of ungerminated seed in the soil;
- `alpha`: the interaction coefficients, `alpha[i, j]` the per capita effect of a germinant
  of species `j` on the seed set of species `i`, in the order of `codes`; `NaN` where the
  authors had little confidence in the estimate (66 of 324, none on the diagonal);
- `effect`, `response`: each species' common competitive effect `e_j` and response `r_i`
  (Appendix S2, Table 2 of the paper), with which the paper fills a missing coefficient as
  `alpha[i, j] = r_i e_j` (see [`annual_plant_alphas`](@ref)).

The `g` and `s` agree with the paper's Appendix S2, Table 3 at its printed precision; `lambda`
differs from that table by up to a quarter (vignette 29).
"""
const ANNUAL_PLANTS_GODOY_2014 = _annual_plant_parameters(joinpath(DATA_DIR, "annual_plants"))

function _annual_plant_index(p, code::Symbol)
    i = findfirst(==(code), p.codes)
    i === nothing && throw(ArgumentError("unknown species :$code; the species are $(join(p.codes, ", "))"))
    return i
end

"""
    annual_plant_alphas(parameters = ANNUAL_PLANTS_GODOY_2014; fill = true)

The matrix of interaction coefficients, `alpha[i, j]` the effect of a germinant of species `j`
on species `i`. With `fill = true` a missing coefficient is replaced, as in Godoy et al.
(2014, Appendix S4), by the product `r_i e_j` of the responding species' common competitive
response and the other species' common competitive effect; with `fill = false` it stays
`NaN`.
"""
function annual_plant_alphas(p = ANNUAL_PLANTS_GODOY_2014; fill::Bool = true)
    A = copy(p.alpha)
    fill || return A
    for i in axes(A, 1), j in axes(A, 2)
        isnan(A[i, j]) && (A[i, j] = p.response[i] * p.effect[j])
    end
    return A
end

"""
    annual_plant(name; g, s, lambda, alpha = [])

One annual plant with a seed bank, counted just after germination in two stages: `seeds`,
the ungerminated seeds in the soil, and `plants`, the germinants. Its transitions are

- `seed_survival`: a dormant seed survives the year (`s`) and next year germinates (`g`,
  to `plants`) or stays dormant (`1 - g`, to `seeds`);
- `seed_production`: a germinant sets `lambda / (1 + crowding)` seeds, which next year
  germinate or stay dormant in the same shares (marked as births).

`crowding` is the pool `<name>_crowding`. The plant adds its germinants to the crowding of
each species it competes with: `alpha` lists `target => coefficient`, the per capita effect
of this species' germinants on the seed set of `target`. Include `name` itself for the
competition within the species.

The seeds in the soil before germination, the `N` of Godoy et al. (2014), are the sum of the
two stages, and a state with `N` seeds is `[(1 - g) N, g N]` (see [`annual_plant_state`](@ref)).
"""
function annual_plant(name::Symbol; g::Real, s::Real, lambda::Real, alpha = Pair{Symbol, Float64}[])
    (0 <= g <= 1 && 0 <= s <= 1) || throw(ArgumentError("g and s must lie in [0, 1]"))
    life = ValuedProjectionNet([:seeds, :plants],
        :seed_survival => [(:seeds => :seeds) => s * (1 - g), (:seeds => :plants) => s * g],
        :seed_production => [(:plants => :seeds) => lambda * (1 - g), (:plants => :plants) => lambda * g])
    own = Symbol(name, :_crowding)
    net = births(respond(life, :seed_production => ctx -> 1 / (1 + getproperty(ctx.pool, own))), :seed_production)
    occupies = Pair{Symbol, Any}[Symbol(target, :_crowding) => [:plants => Float64(a)] for (target, a) in alpha]
    return species(name, net; occupies, reads = any(o -> o.first === own, occupies) ? Symbol[] : [own])
end

"""
    annual_plant_community(parameters = ANNUAL_PLANTS_GODOY_2014; codes = parameters.codes,
                           fill = true)

The species `codes` of `parameters` (by default all 18 of Godoy et al. 2014) as
[`annual_plant`](@ref)s composed through their crowding pools, each species crowding every
other with the coefficients of [`annual_plant_alphas`](@ref) (`fill` as there). With
`fill = false` every coefficient among `codes` must be known. The state of each species is
`[seeds, plants]`.
"""
function annual_plant_community(p = ANNUAL_PLANTS_GODOY_2014; codes = p.codes, fill::Bool = true)
    A = annual_plant_alphas(p; fill)
    idx = [_annual_plant_index(p, c) for c in codes]
    plants = map(idx) do j
        alpha = [p.codes[i] => A[i, j] for i in idx]
        any(isnan(a.second) for a in alpha) &&
            throw(ArgumentError("a coefficient of :$(p.codes[j]) is missing; use fill = true"))
        annual_plant(p.codes[j]; g = p.g[j], s = p.s[j], lambda = p.lambda[j], alpha)
    end
    return length(plants) == 1 ? only(plants) : community(plants...)
end

"""
    annual_plant_state(code, seeds; parameters = ANNUAL_PLANTS_GODOY_2014)

The stage vector `[(1 - g) N, g N]` of species `code` with `N = seeds` seeds in the soil
before germination.
"""
function annual_plant_state(code::Symbol, seeds::Real; parameters = ANNUAL_PLANTS_GODOY_2014)
    g = parameters.g[_annual_plant_index(parameters, code)]
    return [(1 - g) * seeds, g * seeds]
end

"""
    annual_plant_eta(parameters = ANNUAL_PLANTS_GODOY_2014)

The demographic ratio `eta_i = lambda_i g_i / (1 - (1 - g_i) s_i)` of each species, the seeds
produced per seed lost from the soil (to death or germination) without competitors. A species
persists alone if `eta_i > 1`; it then settles with `(eta_i - 1) / alpha_ii` germinants.
"""
annual_plant_eta(p = ANNUAL_PLANTS_GODOY_2014) = @. p.lambda * p.g / (1 - (1 - p.g) * p.s)

"""
    annual_plant_pairs(parameters = ANNUAL_PLANTS_GODOY_2014; fill = true)

The niche and fitness differences of Godoy et al. (2014) for every pair of species `i < j`,
in the order of `parameters.codes`, with the coefficients of [`annual_plant_alphas`](@ref):

- `rho`: the niche overlap `sqrt(alpha_ij alpha_ji / (alpha_ii alpha_jj))` (eqn. 3), and
  `niche_difference = 1 - rho`;
- `fitness_ratio`: `kappa_j / kappa_i = (eta_j - 1) / (eta_i - 1) sqrt(alpha_ij alpha_ii /
  (alpha_jj alpha_ji))` (eqn. 4), above one when `j` is the fitter;
- `invasion_i`, `invasion_j`: the growth rate of each species when rare, with the other at
  its equilibrium alone, `(1 - g_i) s_i + g_i lambda_i / (1 + alpha_ij (eta_j - 1) / alpha_jj)`;
- `outcome`: `:coexistence` if both invade (the condition `rho < kappa_j / kappa_i < 1 / rho`
  of eqn. 5), `:exclusion` if one does (`winner` is its code), `:priority_effect` if neither;
- `filled`: how many of the four coefficients were missing and filled.

Each is a `NamedTuple`; the result is a vector of 153 of them for the 18 species.
"""
function annual_plant_pairs(p = ANNUAL_PLANTS_GODOY_2014; fill::Bool = true)
    A = annual_plant_alphas(p; fill)
    eta = annual_plant_eta(p)
    k = length(p.codes)
    invade(i, j) = (1 - p.g[i]) * p.s[i] + p.g[i] * p.lambda[i] / (1 + A[i, j] * (eta[j] - 1) / A[j, j])
    pairs = map([(i, j) for i in 1:k for j in (i + 1):k]) do (i, j)
        rho = sqrt(A[i, j] * A[j, i] / (A[i, i] * A[j, j]))
        ratio = (eta[j] - 1) / (eta[i] - 1) * sqrt(A[i, j] * A[i, i] / (A[j, j] * A[j, i]))
        ri, rj = invade(i, j), invade(j, i)
        outcome = ri > 1 && rj > 1 ? :coexistence : ri > 1 || rj > 1 ? :exclusion : :priority_effect
        winner = outcome === :exclusion ? (ri > 1 ? p.codes[i] : p.codes[j]) : nothing
        (i = p.codes[i], j = p.codes[j], rho = rho, niche_difference = 1 - rho, fitness_ratio = ratio,
         invasion_i = ri, invasion_j = rj, outcome = outcome, winner = winner,
         filled = count(isnan, (p.alpha[i, i], p.alpha[j, j], p.alpha[i, j], p.alpha[j, i])))
    end
    return pairs
end
