# The dominant perennials of five western North American grasslands as integral projection models coupled by
# crowding: the models of Chu & Adler (2015, Ecological Monographs 85:373-392) as run by Tredennick, de
# Mazancourt, Loreau & Adler (2017, Ecology 98:971-981).
#
# The models have the form of the sagebrush steppe model of Adler, Ellner & Levine (2010) (src/sagebrush.jl):
# genets by log basal area, growth and survival crowded by neighbours of each species, recruitment set by the
# cover of each species, year effects for each year of the data. Only the parameters differ, so a site's plants
# are `sagebrush_plant`s built with that site's parameters, and they meet through the same cover pools
# (`<species>_cover`).
#
# The parameters are read from data/grasslands, written out from the authors' repository; see
# data/grasslands/LICENSE.txt.

const _GRASSLAND_SITES = (arizona = "Arizona", idaho = "Idaho", kansas = "Kansas", montana = "Montana",
                          newmexico = "NewMexico")
const _GRASSLAND_SPECIES_NAMES = Dict(
    "ARTR" => "Artemisia tripartita", "BOCU" => "Bouteloua curtipendula", "BOER" => "Bouteloua eriopoda",
    "BOGR" => "Bouteloua gracilis", "BOHI" => "Bouteloua hirsuta", "BORO" => "Bouteloua rothrockii",
    "HECO" => "Hesperostipa comata", "PASM" => "Pascopyrum smithii", "POSE" => "Poa secunda",
    "PSSP" => "Pseudoroegneria spicata", "SCSC" => "Schizachyrium scoparium", "SPFL" => "Sporobolus flexuosus")

function _grassland_parameters(dir)
    read_table(file) = (lines = [split(l, ","; keepempty = true) for l in readlines(joinpath(dir, file))];
                        (lines[1], lines[2:end]))
    sh, srows = read_table("species.csv")
    yh, yrows = read_table("years.csv")
    function site_parameters(site)
        rows = [r for r in srows if r[1] == site]
        n = length(rows)
        function one(r)
            x(name) = parse(Float64, r[findfirst(==(name), sh)])
            v(prefix) = [x(prefix * string(j)) for j in 1:n]
            ys = [y for y in yrows if y[1] == site && y[2] == r[2]]
            sort!(ys; by = y -> parse(Int, y[3]))
            yearly(name) = [parse(Float64, y[findfirst(==(name), yh)]) for y in ys]
            growth = (intercept = x("growth_intercept"), slope = x("growth_slope"), crowding = v("growth_crowding_"),
                      alpha = v("growth_alpha_"), sigma_a = x("growth_sigma_a"), sigma_b = x("growth_sigma_b"),
                      intercept_year = yearly("growth_intercept"), slope_year = yearly("growth_slope"))
            survival = (intercept = x("survival_intercept"), slope = x("survival_slope"), crowding = v("survival_crowding_"),
                        alpha = v("survival_alpha_"), intercept_year = yearly("survival_intercept"),
                        slope_year = yearly("survival_slope"))
            recruitment = (intercept = x("recruitment_intercept"), intercept_year = yearly("recruitment_intercept"),
                           crowding = v("recruitment_crowding_"), size_mean = x("recruit_size_mean"),
                           size_variance = x("recruit_size_variance"))
            code = String(r[2])
            return (name = _GRASSLAND_SPECIES_NAMES[code], code = code, largest = x("largest"),
                    cells = parse(Int, r[findfirst(==("cells"), sh)]), growth = growth, survival = survival,
                    recruitment = recruitment)
        end
        kinds = Tuple(Symbol(lowercase(r[2])) for r in rows)
        return NamedTuple{kinds}(Tuple(one(r) for r in rows))
    end
    return NamedTuple{keys(_GRASSLAND_SITES)}(Tuple(site_parameters(s) for s in values(_GRASSLAND_SITES)))
end

"""
    GRASSLAND_PARAMETERS

Fitted vital rates of the dominant perennials of five semiarid grasslands of the western
United States, the models of Chu & Adler (2015) as run by Tredennick et al. (2017), read
from `data/grasslands`:

- `arizona` (Sonoran desert): `boer` *Bouteloua eriopoda*, `boro` *B. rothrockii*;
- `idaho` (sagebrush steppe): `artr` *Artemisia tripartita*, `heco` *Hesperostipa comata*,
  `pose` *Poa secunda*, `pssp` *Pseudoroegneria spicata*;
- `kansas` (southern mixed prairie): `bocu` *Bouteloua curtipendula*, `bohi` *B. hirsuta*,
  `scsc` *Schizachyrium scoparium*;
- `montana` (northern mixed prairie): `bogr` *Bouteloua gracilis*, `heco`, `pasm`
  *Pascopyrum smithii*, `pose`;
- `newmexico` (Chihuahuan desert): `boer`, `spfl` *Sporobolus flexuosus*.

Each site is a set of parameters in the form of [`SAGEBRUSH_PARAMETERS`](@ref) (see there
for the meaning of each field), so [`sagebrush_plant`](@ref) builds its species; the
crowding coefficients and spatial scales are by neighbour, in the order of the site's
species. [`grassland_plant`](@ref) and [`grassland_community`](@ref) do this for a site.

The Idaho set is the one of Tredennick et al. (2017). It differs from `SAGEBRUSH_PARAMETERS`
(the R package IPMdoit), notably in its recruitment coefficients and in the mesh of the
shrub (50 cells rather than 75).
"""
const GRASSLAND_PARAMETERS = _grassland_parameters(joinpath(DATA_DIR, "grasslands"))

_grassland_site(site::Symbol) = haskey(GRASSLAND_PARAMETERS, site) ? GRASSLAND_PARAMETERS[site] :
    throw(ArgumentError("unknown site :$site; the sites are $(join(keys(GRASSLAND_PARAMETERS), ", "))"))

"""
    grassland_plant(site, kind; niche_differences = true, interactions = true, cells, name = kind)

Species `kind` of a grassland `site` of [`GRASSLAND_PARAMETERS`](@ref), for example
`grassland_plant(:kansas, :scsc)`: [`sagebrush_plant`](@ref) with the site's parameters. It
adds its basal area to the pool `<kind>_cover` and reads the covers of the site's other
species. The environment of a step is `nothing` (the mean year) or a year of the site's
data, `1` to the number of years of its fits (18 in Arizona, 22 in Idaho, 34 in Kansas, 13
in Montana, 32 in New Mexico).
"""
function grassland_plant(site::Symbol, kind::Symbol; kwargs...)
    parameters = _grassland_site(site)
    return sagebrush_plant(kind; parameters, kwargs...)
end

"""
    grassland_community(site; kinds = all of the site's species, niche_differences = true,
                        interactions = true)

The species of a grassland `site` composed through their cover pools. With
`interactions = false` each species responds to its own species alone, as in the "no
competition" simulations of Tredennick et al. (2017).
"""
function grassland_community(site::Symbol; kinds = keys(_grassland_site(site)), niche_differences = true,
                             interactions::Bool = true)
    return sagebrush_community(; parameters = _grassland_site(site), kinds, niche_differences, interactions)
end

"""
    grassland_state(site, kind; density = 0.001, cover = nothing)

A stage vector for `grassland_plant(site, kind)`: `density` genets per m² for each unit of
log area in every cell, the starting state of the authors' simulations, or that shape
scaled to `cover` (a share of the ground).
"""
grassland_state(site::Symbol, kind::Symbol; density::Real = 0.001, cover = nothing) =
    sagebrush_state(kind; density, cover, parameters = _grassland_site(site))

"""
    synchrony(X)

The community-wide synchrony of Loreau & de Mazancourt (2008) of the columns of `X`, one
column per species and one row per time: the variance of the row sums divided by the
square of the sum of the columns' standard deviations. It is 1 when the species move in
perfect step and near 0 when they compensate one another. Applied to yearly per capita
growth rates, `log.(C[2:end, :] ./ C[1:(end - 1), :])` for covers `C`, it is the synchrony
of growth rates studied by Tredennick et al. (2017).
"""
function synchrony(X::AbstractMatrix{<:Real})
    size(X, 1) >= 2 || throw(ArgumentError("synchrony needs at least two rows"))
    sd(x) = sqrt(sum(abs2, x .- sum(x) / length(x)) / (length(x) - 1))
    return sd(vec(sum(X; dims = 2)))^2 / sum(sd(c) for c in eachcol(X))^2
end
