# Writes data/corals/mcr_lter_observed_colonies.csv: the fore-reef coral colonies observed in 2013 and 2015 on the
# four reefs of Kayal et al. (2018), binned onto the mesh of `coral_stage()` (z = log10(S + 1) on [0, 7] in 100 cells
# of width 0.07); and data/corals/mcr_lter_observed_transects.csv: the same colonies totalled by transect (their
# number and summed z), from which the means and standard errors of Figs. S5 and S7 of the paper follow.
#
# Usage, from the package root:
#
#   julia data/corals/bin_observed_colonies.jl path/to/MCR_LTER_Lenihan_mapped_corals_20220413.csv
#
# The input is the table "Mapped Corals long format" of the EDI data package knb-lter-mcr.4009, revision 5
# (doi:10.6073/pasta/322a76c6f197bfac60c91d66a3e96fb5), SHA-1 cb91e282e52f84ea47085431bed6aa6be4274c76. It is not
# included in this package; see data/corals/LICENSE.txt. Only Base Julia is used.
#
# Rows kept: Habitat "OR" (fore reef), Year 2013 or 2015, Taxa Acr, Poc or Por, with Length, Width and Height all
# measured (in the long format a colony absent in a year, not yet settled or dead, has NA dimensions). Acropora was
# mapped along the whole of each 5 m x 1 m transect (X from 0 to 5 m), Pocillopora and Porites only in its first metre
# (X up to 1 m): four transects per reef, so 20 m² for Acropora and 4 m² for the other two. The few Pocillopora and
# Porites recorded beyond the first metre are dropped, so that every colony counted lies in the area sampled for its
# genus. Surface area S (cm²) is that of Appendix S2 of the paper: a sphere for Pocillopora,
# S = 4π((l + w + h)/6)², and a cylinder for Acropora and Porites, S = 2π((l + w)/4)² + 2π((l + w)/4) h.
#
# A transect counts in a year if at least one colony of the three genera was measured on it that year: four per reef
# and year. A transect of LTER 2 is labelled T4(18)Old in 2013 and T4(18) in 2015, and one of LTER 5 T2 Old in 2013 and
# T2 New in 2015; the labels are kept as they are. A genus absent from a counted transect has a row of zeros.

const SHA1 = "cb91e282e52f84ea47085431bed6aa6be4274c76"
const WIDTH = 0.07
const CELLS = 100
const SAMPLED = Dict("Acr" => (x_max = 5.0, area = 20.0), "Poc" => (x_max = 1.0, area = 4.0),
                     "Por" => (x_max = 1.0, area = 4.0))

function surface_area(taxon, l, w, h)
    taxon == "Poc" && return 4π * ((l + w + h) / 6)^2
    return 2π * ((l + w) / 4)^2 + 2π * ((l + w) / 4) * h
end

function bin_colonies(input, output, transects_output)
    lines = readlines(input)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    counts = Dict{Tuple{Int, Int, String, Int}, Int}()
    by_transect = Dict{Tuple{Int, Int, String, String}, Tuple{Int, Float64}}()
    dropped = 0
    for line in lines[2:end]
        r = split(line, ",")
        r[col("Habitat")] == "OR" || continue
        year = parse(Int, r[col("Year")])
        year in (2013, 2015) || continue
        taxon = r[col("Taxa")]
        haskey(SAMPLED, taxon) || continue
        dims = (r[col("Length")], r[col("Width")], r[col("Height")])
        any(d -> d in ("NA", "", " "), dims) && continue
        if parse(Float64, r[col("X")]) > SAMPLED[taxon].x_max
            dropped += 1
            continue
        end
        l, w, h = parse.(Float64, dims)
        z = log10(surface_area(taxon, l, w, h) + 1)
        cell = min(CELLS, floor(Int, z / WIDTH) + 1)
        reef = parse(Int, replace(r[col("Site")], "LTER" => ""))
        key = (reef, year, taxon, cell)
        counts[key] = get(counts, key, 0) + 1
        tkey = (reef, year, taxon, r[col("Transect")])
        n, s = get(by_transect, tkey, (0, 0.0))
        by_transect[tkey] = (n + 1, s + z)
    end
    open(output, "w") do io
        println(io, "reef,year,taxon,cell,colonies,sampled_area_m2")
        for k in sort(collect(keys(counts)))
            println(io, join((k..., counts[k], SAMPLED[k[3]].area), ","))
        end
    end
    transects = sort(unique((k[1], k[2], k[4]) for k in keys(by_transect)))
    open(transects_output, "w") do io
        println(io, "reef,year,taxon,transect,colonies,sum_log_size,sampled_area_m2")
        for (reef, year, transect) in transects, taxon in sort(collect(keys(SAMPLED)))
            n, s = get(by_transect, (reef, year, taxon, transect), (0, 0.0))
            println(io, join((reef, year, taxon, transect, n, round(s; digits = 6), SAMPLED[taxon].area / 4), ","))
        end
    end
    return dropped
end

if abspath(PROGRAM_FILE) == @__FILE__
    input = only(ARGS)
    digest = strip(read(`shasum -a 1 $input`, String))
    startswith(digest, SHA1) || @warn "SHA-1 of $input is not that of knb-lter-mcr.4009.5" digest
    dropped = bin_colonies(input, joinpath(@__DIR__, "mcr_lter_observed_colonies.csv"),
                           joinpath(@__DIR__, "mcr_lter_observed_transects.csv"))
    println("dropped $dropped Pocillopora and Porites outside the first metre of their transect")
end
