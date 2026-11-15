# Forest trees as integral projection models, coupled through the basal area of the stand: the model of
# the R package matreex 0.4.0 (Jaunatre et al., https://github.com/gowachin/matreex, MIT licence), which
# packages the models of Kunstler et al. (2021) and Guyennon et al. (2023). Each species is a net; the
# species meet only through one pool:
#
#   basal_area   the cross-sectional area of the trunks of all trees of the stand, in m² per hectare
#
# Size is trunk diameter at breast height in mm, the time step is a year, and numbers are trees per hectare.
#
# Left out: the harvesting rules (uneven- and even-aged stands), disturbances and a climate that changes
# through time, which matreex has, and its regional pool of seed, which is in development there. A constant
# harvest rate is kept.
#
# The fitted coefficients in FOREST_PARAMETERS and the climate in FOREST_CLIMATE are copied from matreex,
# which is distributed under the MIT licence with this notice:
#
#   Copyright (c) 2022 maxime jaunatre
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
#   associated documentation files (the "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
#   following conditions:
#
#   The above copyright notice and this permission notice shall be included in all copies or substantial
#   portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
#   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
#   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#   IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
#   THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
    FOREST_CLIMATE

Default climate of [`forest_tree`](@ref): the yearly sum of degree days above 5.5 °C
(`sgdd`) and the water availability index (`wai`: precipitation minus potential
evapotranspiration, divided by potential evapotranspiration) at the middle of the climatic
range of Norway spruce. The values are those of the R package matreex (its table
`climate_species`, *Picea abies*, `N == 2`), the climate of its introductory vignette.
"""
const FOREST_CLIMATE = (sgdd = 1444.6666281571565, wai = 0.45193869236972689)

"""
    FOREST_PARAMETERS

Fitted vital rates of three European trees, for [`forest_tree`](@ref): Norway spruce
(`spruce`, *Picea abies*), silver fir (`fir`, *Abies alba*) and European beech (`beech`,
*Fagus sylvatica*). The numbers are those shipped with the R package matreex 0.4.0 (objects
`fit_Picea_abies`, `fit_Abies_alba` and `fit_Fagus_sylvatica`; MIT licence), copied with all
their digits.

For each species: `largest`, the largest trunk diameter in the data (mm); `delay`, the
years a recruit takes to reach the smallest diameter recorded; and three regressions,
each a list of `term => coefficient`:

- `growth`: the mean of the logarithm of a tree's diameter growth in a year (mm), which is
  normally distributed with standard deviation `sigma`;
- `death`: the complementary log-log of the chance that a tree dies in a year;
- `recruitment`: the logarithm of the number of recruits in a year on a plot of 300 m².

A term is a product of variables, written as in matreex: `size` (diameter, mm), `logsize`,
`BATOTcomp` (basal area of the stand, m² per hectare), `BATOTSP` and `logBATOTSP` (basal
area of the species), `BATOTNonSP` (basal area of the other species), and the climate
variables `sgdd`, `sgdd2` (its square), `sgddb` (its inverse), `wai`, `wai2` (its square)
and `waib` (the inverse of `1 + wai`); see [`FOREST_CLIMATE`](@ref).

Growth and survival were fitted to European forest inventories by Kunstler et al. (2021)
and are the means of their 100 fits to resampled data. The recruitment coefficients are
those of matreex 0.4.0, which differ from the fits published with Guyennon et al. (2023).
"""
const FOREST_PARAMETERS = (
    spruce = (
        name = "Picea abies", largest = 1095.0, delay = 14, sigma = 0.58130213846066026,
        growth = [
            "BATOTcomp" => -0.0083643119931635784,
            "BATOTcomp:sgdd" => -8.4021405783187631e-06,
            "BATOTcomp:wai" => -0.00086724876056916683,
            "intercept" => 1.4799924646074496,
            "logsize" => -0.51591322405654239,
            "logsize:sgdd" => 0.0010326789994132208,
            "logsize:wai" => 0.28370187720113732,
            "sgdd" => -0.0022893707204971313,
            "sgdd2" => -6.4282670253797942e-07,
            "size" => 0.0014770792781632625,
            "size:sgdd" => -2.1829923219315674e-06,
            "size:wai" => -0.00030836441467356821,
            "wai" => -1.5432696702962785,
            "wai2" => 0.066758268372484941],
        death = [
            "BATOTcomp" => 0.036146046770128663,
            "BATOTcomp:sgdd" => -3.3335140555622886e-06,
            "BATOTcomp:wai" => -0.010027596993363584,
            "intercept" => -21.040711505115855,
            "logsize" => 3.298094780865489,
            "logsize:sgdd" => -0.0029264038797028511,
            "logsize:wai" => -1.6002025419697066,
            "sgdd" => 0.012414559174795189,
            "sgdd2" => 5.4920807017687575e-07,
            "size" => -0.0049277324140054278,
            "size:sgdd" => 5.1289082755233611e-06,
            "size:wai" => 0.0029559877671803745,
            "wai" => 7.7586097097048849,
            "wai2" => 0.26128370231316472],
        recruitment = [
            "BATOTNonSP" => -0.058326133795597999,
            "BATOTSP" => -0.08261184324459675,
            "intercept" => 0.15680121076377346,
            "logBATOTSP" => 0.52155184190267923,
            "sgddb" => -1507.8944482000584,
            "waib" => -0.54710947711694757]),
    fir = (
        name = "Abies alba", largest = 1420.0, delay = 28, sigma = 0.59672620418878219,
        growth = [
            "BATOTcomp" => -0.018168922665520226,
            "intercept" => -5.3320914534028176,
            "logsize" => 1.3423510131707401,
            "logsize:sgdd" => 0.00011356180690487391,
            "logsize:sgddb" => -316.94307740370215,
            "logsize:wai" => -0.13377359558180868,
            "logsize:waib" => 0.14696823830803502,
            "sgdd" => -0.00015080745973269842,
            "sgdd2" => -7.6035933890937168e-08,
            "sgddb" => 1381.3817302653613,
            "size" => -0.0024181315776910608,
            "size:sgdd" => -2.0703057312844438e-07,
            "size:sgddb" => 0.48707262204664836,
            "size:wai" => 0.00082702922595691872,
            "size:waib" => -0.00098844193175090154,
            "wai" => 0.50521936225376241,
            "wai2" => -0.042814517923865085,
            "waib" => -0.48473073052513843],
        death = [
            "BATOTcomp" => 0.02141555902528116,
            "BATOTcomp:sgddb" => -2.8022087554720714,
            "BATOTcomp:waib" => -0.00023447190404364855,
            "intercept" => 4.1835866309762144,
            "logsize" => -2.4311280464907719,
            "logsize:sgddb" => -0.86473685006000212,
            "logsize:waib" => -0.013810039377160698,
            "sgddb" => 306.77931698289029,
            "size" => 0.0048515858782366928,
            "size:sgddb" => 0.077714351636184481,
            "size:waib" => 4.9914950556838527e-05,
            "waib" => 1.483629155648202],
        recruitment = [
            "BATOTNonSP" => -0.037110469246740593,
            "BATOTSP" => -0.036782425223471429,
            "intercept" => -0.95748084517838428,
            "logBATOTSP" => 0.062390333902184118,
            "sgdd" => 2.9381629584017514e-05,
            "sgddb" => -730.48824161497873,
            "wai" => -0.26928366450279329,
            "waib" => 0.21046865549752861]),
    beech = (
        name = "Fagus sylvatica", largest = 1649.0, delay = 48, sigma = 0.62292439355385043,
        growth = [
            "BATOTcomp" => -0.015356653309507052,
            "intercept" => -4.9649034554548868,
            "logsize" => 0.96987942989529441,
            "logsize:sgdd" => 0.00010155063157377876,
            "logsize:wai" => 0.11746194772024279,
            "sgdd" => 0.00051204099026076701,
            "sgdd2" => -1.2621957152595005e-07,
            "size" => -0.0014675320155576636,
            "size:sgdd" => -3.4381004499189121e-07,
            "size:wai" => 6.9563548570258539e-05,
            "wai" => -0.32831176589709954,
            "wai2" => -0.18407630101019054],
        death = [
            "BATOTcomp" => 0.0078866408988526628,
            "intercept" => 11.901376489252545,
            "logsize" => -3.9369322626453283,
            "logsize:sgdd" => 0.00070656660236307592,
            "logsize:sgddb" => -300.24090225085342,
            "logsize:wai" => -0.25471510935200625,
            "logsize:waib" => 0.58900275676363978,
            "sgdd" => -0.0027036923160317491,
            "sgdd2" => -1.4951797689625871e-07,
            "sgddb" => 1280.8317883852974,
            "size" => 0.0088337008614907502,
            "size:sgdd" => -1.3761466273677924e-06,
            "size:sgddb" => 0.65499413049487398,
            "size:wai" => 0.0005258198498984555,
            "size:waib" => -0.0012707726653347507,
            "wai" => 0.63437271950329055,
            "wai2" => 0.39393880289560101,
            "waib" => -2.6136924557714867],
        recruitment = [
            "BATOTNonSP" => -0.042610558754950852,
            "BATOTSP" => -0.0098180103718472616,
            "intercept" => -2.264983660004583,
            "logBATOTSP" => -0.10356062648905924,
            "sgddb" => 989.57270868132275,
            "wai" => -0.37025130022946368,
            "wai2" => 0.52905057955783175]))

const _FOREST_SMALLEST = 90.0   # mm: the lower end of the mesh (matreex: L = 90); the inventories record trees from 100 mm
const _FOREST_RANGE = 1.1       # the mesh ends at 1.1 times the largest diameter in the data (matreex)
const _FOREST_PLOT = 0.03       # ha: recruitment was fitted as recruits on a plot of 300 m²
const _FOREST_TAIL = 6.0        # growth beyond this many standard deviations of its logarithm is left out (chance 1e-9)

_trunk_area(d) = pi * (d / 2000)^2      # m², for a diameter in mm

# The climate variables of the fitted models: the two that are measured, their squares, and the transforms
# used for responses that level off (as tabulated in matreex's `climate_species`).
_climate_variables(c) = Dict("sgdd" => c.sgdd, "sgdd2" => c.sgdd^2, "sgddb" => 1 / c.sgdd,
                             "wai" => c.wai, "wai2" => c.wai^2, "waib" => 1 / (1 + c.wai))

# A fitted regression at one climate. Each term is a coefficient times a product of variables
# ("logsize:sgdd"). The climate variables are replaced by their values and the coefficients are added up
# by the variable that is left: one of `variables`, or none ("").
function _at_climate(terms, climate, variables)
    known = _climate_variables(climate)
    total = Dict{String, Float64}(v => 0.0 for v in variables)
    total[""] = 0.0
    for (term, coefficient) in terms
        value, left = Float64(coefficient), ""
        for factor in split(term, ":")
            if factor == "intercept"
                continue
            elseif haskey(known, factor)
                value *= known[factor]
            elseif left == "" && factor in variables
                left = String(factor)
            else
                throw(ArgumentError("the term \"$term\" of a fitted model is not understood"))
            end
        end
        total[left] += value
    end
    return total
end

_forest_species(parameters, kind::Symbol) = haskey(parameters, kind) ? parameters[kind] :
    throw(ArgumentError("unknown tree :$kind; the parameters hold $(join(keys(parameters), ", "))"))

"""
    forest_vital_rates(kind; parameters = FOREST_PARAMETERS, climate = FOREST_CLIMATE)

The fitted vital rates of the tree `kind` (`:spruce`, `:fir` or `:beech`) at one climate, as
functions of trunk diameter `d` (mm) and basal area (m² per hectare):

- `log_increment(d, basal_area)`: the mean of the logarithm of the diameter growth of a
  tree in a year (mm). The logarithm is normally distributed with standard deviation
  `sigma`, so the median growth is `exp(log_increment(d, basal_area))`.
- `death(d, basal_area)`: the chance that a tree dies in a year.
- `recruits(own, others)`: recruits per hectare and year, for a species that holds the
  basal area `own` among other species that hold `others`.

`climate` is a `NamedTuple` with `sgdd` and `wai`; see [`FOREST_CLIMATE`](@ref).
"""
function forest_vital_rates(kind::Symbol; parameters = FOREST_PARAMETERS, climate = FOREST_CLIMATE)
    p = _forest_species(parameters, kind)
    g = _at_climate(p.growth, climate, ("size", "logsize", "BATOTcomp"))
    m = _at_climate(p.death, climate, ("size", "logsize", "BATOTcomp"))
    r = _at_climate(p.recruitment, climate, ("BATOTSP", "BATOTNonSP", "logBATOTSP"))
    g0, g1, g2, g3 = g[""], g["size"], g["logsize"], g["BATOTcomp"]
    m0, m1, m2, m3 = m[""], m["size"], m["logsize"], m["BATOTcomp"]
    r0, r1, r2, r3 = r[""], r["BATOTSP"], r["BATOTNonSP"], r["logBATOTSP"]
    log_increment(d, basal_area) = g0 + g1 * d + g2 * log(d) + g3 * basal_area
    # complementary log-log link, kept within (0, 1) as R's binomial family does
    death(d, basal_area) = clamp(-expm1(-exp(m0 + m1 * d + m2 * log(d) + m3 * basal_area)), eps(), 1 - eps())
    recruits(own, others) = exp(r0 + r1 * own + r2 * others + r3 * log(own)) / _FOREST_PLOT
    return (log_increment = log_increment, sigma = Float64(p.sigma), death = death, recruits = recruits)
end

"""
    forest_stage(kind; parameters = FOREST_PARAMETERS, cells = 700)

The [`ContinuousStage`](@ref) of [`forest_tree`](@ref): the trees of species `kind` by
trunk diameter in mm, from 90 mm to 1.1 times the largest diameter in the data, in `cells`
mesh cells. `meshpoints` of it are the diameters at the middle of the cells.
"""
function forest_stage(kind::Symbol; parameters = FOREST_PARAMETERS, cells::Integer = 700)
    p = _forest_species(parameters, kind)
    return ContinuousStage(:tree, _FOREST_SMALLEST, _FOREST_RANGE * p.largest, cells)
end

# The standard normal distribution function: the rational approximations of Hart (1968, Computer
# Approximations) as arranged by West (2005, Wilmott Magazine, "Better approximations to cumulative normal
# functions"). Against R's pnorm its absolute error is below 2e-16 (test/test_forests.jl); its relative error
# reaches 1e-8 far in the tails, where the values are too small to matter here.
function _normal_cdf(x::Real)
    a = abs(Float64(x))
    a > 37 && return x > 0 ? 1.0 : 0.0
    e = exp(-a^2 / 2)
    tail = if a < 7.07106781186547
        top = ((((((3.52624965998911e-2 * a + 0.700383064443688) * a + 6.37396220353165) * a + 33.912866078383) * a +
                 112.079291497871) * a + 221.213596169931) * a + 220.206867912376)
        bottom = (((((((8.83883476483184e-2 * a + 1.75566716318264) * a + 16.064177579207) * a + 86.7807322029461) * a +
                     296.564248779674) * a + 637.333633378831) * a + 793.826512519948) * a + 440.413735824752)
        e * top / bottom
    else
        e / (a + 1 / (a + 2 / (a + 3 / (a + 4 / (a + 0.65))))) / 2.506628274631
    end
    return x > 0 ? 1 - tail : tail
end

# The chance that a year's growth lies between `lo` and `hi` (mm), its logarithm being normal.
function _increment_chance(lo, hi, mean_log, sigma)
    hi > 0 || return 0.0
    upper = (log(hi) - mean_log) / sigma
    lo > 0 || return _normal_cdf(upper)
    lower = (log(lo) - mean_log) / sigma
    # far in the upper tail both values are close to one: take the difference on the other side instead
    return lower > 0 ? _normal_cdf(-lower) - _normal_cdf(-upper) : _normal_cdf(upper) - _normal_cdf(lower)
end

# The mean of `f` over the cell of width `h` centred on `x`, by the three-point Gauss-Legendre rule: the
# rule matreex uses for the place of a tree within its cell.
const _GAUSS3 = sqrt(3 / 5) / 2
_cell_mean(f, x, h) = (5 * (f(x - _GAUSS3 * h) + f(x + _GAUSS3 * h)) + 8 * f(x)) / 18

"""
    forest_tree(kind; parameters = FOREST_PARAMETERS, climate = FOREST_CLIMATE, cells = 700,
                  delay = parameters[kind].delay, harvest = 0.006, conspecific = nothing, name = kind,
                  tabulate = nothing)

A tree species (`:spruce`, `:fir` or `:beech`) as the integral projection model of the R
package matreex (Kunstler et al. 2021; Guyennon et al. 2023). Trees are structured by trunk
diameter on `cells` mesh cells (see [`forest_stage`](@ref)); numbers are trees per hectare
and the time step is a year. The species contributes the cross-sectional area of its trunks
to the pool `basal_area` (m² per hectare) and reads the pool's total: a species is affected
by the others only through it.

Stages are `seedling`, `sapling` and `tree`. Transitions are

- `growing`: a tree survives, is not felled, and grows. Growth is log-normal; its mean
  falls and the chance of death changes as the basal area of the stand rises, so the kernel
  changes shape with the pool and is rebuilt each year ([`responsive`](@ref)).
- `recruiting`: the species gains recruits according to its own basal area and that of the
  other species. They are shared among its trees in proportion to their basal area.
- `establishing` and `ingrowth`: a recruit joins the smallest trees `delay` years later,
  without loss on the way. As in matreex, half of a year's recruits are a year ahead: they
  start as saplings, the other half as seedlings.

Keywords:

- `climate`: a `NamedTuple` with `sgdd` and `wai`; see [`FOREST_CLIMATE`](@ref).
- `delay`: years from recruitment to the smallest diameter, at least 2.
- `harvest`: the share of the trees felled each year, after growth and before
  recruitment. The default is that of matreex.
- `conspecific`: the basal area that recruitment takes the species to hold, as part of the
  pool's total. With `nothing` it is the basal area of the species' own trees, and a
  species without trees has no recruits. Give a number to hold it fixed whatever the
  numbers of trees: the species is then linear, and has a growth rate
  ([`invasion_growth_rate`](@ref)). That describes a rare invader, whose recruits for each
  unit of basal area depend on how rare it is (Guyennon et al. 2023 used 0.1); it also
  finds the basal area at which a species alone replaces itself, where the growth rate
  with `conspecific` equal to the pool's total is one.
- `tabulate`: `nothing` (the default) rebuilds the growth kernel every year at the basal
  area of the stand. Give a grid of basal areas, such as `0:1:80`, to build it once at each
  of them and interpolate linearly in between, as matreex does ([`tabulated`](@ref)); a
  basal area outside the grid gets the kernel built for it. On 700 cells a year of
  `advance` is then about five times faster. A kernel is kept for each grid value the
  stand reaches, 0.2 to 0.9 MB each (at most about 40 MB per species for `0:1:80`). With
  a grid of whole basal areas the basal area of the stand stays within 1e-4 of itself of
  that of the exact run.

Growth is integrated exactly over the cell a tree arrives in, and by a three-point rule
over its place in the cell it leaves, which is the scheme of matreex. Trees that outgrow
the mesh are lost.
"""
function forest_tree(kind::Symbol; parameters = FOREST_PARAMETERS, climate = FOREST_CLIMATE, cells::Integer = 700,
                     delay::Integer = _forest_species(parameters, kind).delay, harvest::Real = 0.006,
                     conspecific = nothing, name::Symbol = kind, tabulate = nothing)
    delay >= 2 || throw(ArgumentError("the delay must be at least 2 years, got $delay"))
    0 <= harvest < 1 || throw(ArgumentError("harvest is a share of the trees, got $harvest"))
    (conspecific === nothing || conspecific > 0) || throw(ArgumentError("conspecific must be a positive basal area"))
    rates = forest_vital_rates(kind; parameters, climate)
    tree = forest_stage(kind; parameters, cells)
    sizes, h, sigma = meshpoints(tree), step_size(tree), rates.sigma
    kept = 1 - harvest

    # the kernel of the trees that survive, are not felled and grow, in a stand of basal area `stand`
    function growing(stand)
        mean_log(d) = rates.log_increment(d, stand)
        reach = exp(maximum(mean_log, sizes) + _FOREST_TAIL * sigma)      # more growth than this is left out
        # What depends on the size `x` a tree leaves alone (its chance of surviving, and the mean of its growth
        # at the three points of the rule) is kept from one entry to the next, since the kernel is evaluated
        # column by column. The values are those of `_cell_mean`, computed in the same way.
        # (x, alive, d1, d2, mean_log(x), mean_log(d1), mean_log(d2)) for the last x seen
        column = fill(NaN, 7)
        return function (y, x)
            (y < x - h / 2 || y - x - h > reach) && return 0.0            # trees do not shrink
            if x != column[1]
                d1, d2 = x - _GAUSS3 * h, x + _GAUSS3 * h
                column .= (x, 1 - _cell_mean(d -> rates.death(d, stand), x, h), d1, d2, mean_log(x), mean_log(d1), mean_log(d2))
            end
            alive, d1, d2, m0, m1, m2 = column[2], column[3], column[4], column[5], column[6], column[7]
            chance(d, m) = _increment_chance(y - h / 2 - d, y + h / 2 - d, m, sigma)
            arriving = (5 * (chance(d1, m1) + chance(d2, m2)) + 8 * chance(x, m0)) / 18     # _cell_mean
            return kept * alive * arriving / h
        end
    end
    growth = if tabulate === nothing
        responsive(ctx -> growing(ctx.pool.basal_area))
    else
        tabulated(growing; pool = :basal_area, grid = tabulate, outside = :exact)
    end

    life = kernel_net([:seedling, :sapling, tree],
        :growing => [(:tree => :tree) => growth],
        # recruits for each m² per hectare of the species; the response below supplies their number
        :recruiting => [(:tree => :seedling) => d -> _trunk_area(d) / 2, (:tree => :sapling) => d -> _trunk_area(d) / 2],
        :establishing => [(:seedling => :sapling) => 1.0],
        # saplings join the smallest cell, where they can be felled like other trees
        :ingrowth => [(:sapling => :tree) => d -> d < first(sizes) + h / 2 ? kept / h : 0.0])
    delay > 2 && (life = expand_stages(life, :sapling => delay - 1))

    areas = _trunk_area.(sizes)
    function recruiting(ctx)
        # the species' own basal area, from its trees unless it is imposed
        own = conspecific === nothing ? sum(areas .* ctx.n[tree]) : Float64(conspecific)
        own > 0 || return 0.0
        return rates.recruits(own, max(0.0, ctx.pool.basal_area - own)) / own
    end
    net = births(respond(life, :recruiting => recruiting), :recruiting)
    return species(name, net; occupies = [:basal_area => trait_weights(tree, _trunk_area)])
end

"""
    forest_stand(kind; basal_area = 1.0, diameters = (90, 150), parameters = FOREST_PARAMETERS,
                   cells = 700, delay = parameters[kind].delay)

A stage vector for `forest_tree(kind; cells, delay)`: a young stand without seedlings or
saplings, whose trees are spread evenly over the mesh cells with a diameter (mm) within
`diameters` and add up to `basal_area` (m² per hectare). Its entries are in the order of
the species' stages: seedlings and saplings by year, then trees by cell, which is also the
order of the classes of matreex.
"""
function forest_stand(kind::Symbol; basal_area::Real = 1.0, diameters = (90, 150), parameters = FOREST_PARAMETERS,
                      cells::Integer = 700, delay::Integer = _forest_species(parameters, kind).delay)
    sizes = meshpoints(forest_stage(kind; parameters, cells))
    present = [diameters[1] <= d < diameters[2] for d in sizes]
    any(present) || throw(ArgumentError("no mesh cell has a diameter within $diameters"))
    trees = present .* (basal_area / sum(_trunk_area.(sizes[present])))
    return [zeros(delay); trees]              # seedlings and saplings, then the trees by cell
end
