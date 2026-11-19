# The seeded ensembles of whole plants and weevils of vignette 14 (a native thistle, an invasive thistle and their
# weevil), as the vignette prints them.
using Random

@testset "thistles and weevil, ensembles (vignette 14)" begin
    P, W = PLATTE_THISTLE_PARAMETERS, WEEVIL_PARAMETERS
    plant = platte_thistle_stage()
    seedling(y) = exp(-0.5 * (y - 0.75)^2 / 0.17) / sqrt(2pi * 0.17)
    start = 600 .* step_size(plant) .* seedling.(meshpoints(plant))         # 600 seedlings
    native, weevil, weed = platte_thistle(), flowerhead_weevil(), musk_thistle()
    pair, biocontrol = community(native, weevil), community(weed, weevil)
    trio = community(native, weed, weevil)
    # the settled communities of the vignette
    alone = equilibrium(native, start, nothing)
    arrivals = exp(-3.74) / P.heads[1] * pool_totals(native, alone.n, nothing).heads / W.eggs
    sol = project(pair, [alone.n; arrivals], fill(nothing, 60))
    settled = equilibrium(pair, sol.u[end, :], nothing)
    stand = equilibrium(weed, [1000.0, 10.0, 5.0, 1.0], nothing)
    held = equilibrium(biocontrol, [stand.n; 5.0], nothing)
    all3 = equilibrium(trio, [settled.n[1:100]; held.n], nothing)
    @test alone.converged && settled.converged && stand.converged && held.converged && all3.converged

    function whole(n, rng)                # whole plants drawn from the size distribution n
        cells, cumulative = zeros(length(n)), cumsum(n) ./ sum(n)
        for _ in 1:round(Int, sum(n))
            cells[searchsortedfirst(cumulative, rand(rng))] += 1
        end
        return cells
    end
    function outcome(ens)
        plants, weevils = abundance(ens, :platte_thistle), abundance(ens, :weevil)
        return (median_plants = long_median(plants[:, end]), nine_in_ten_between = long_quantile.(Ref(plants[:, end]), [0.05, 0.95]),
                thistle_lost = count(==(0), plants[:, end]), weevil_lost = count(==(0), weevils[:, end]))
    end
    # 200 replicates over a century, drawn with MersenneTwister(1); the starting plants with MersenneTwister(14).
    # Seeded, so the vignette's values are exact on any platform that draws the same numbers.
    rng = MersenneTwister(14)
    century = fill(nothing, 100)
    two = ensemble(pair, [whole(settled.n[1:100], rng); round(settled.n[101])], century; replicates = 200, rng = MersenneTwister(1))
    three = ensemble(trio, [whole(all3.n[1:100], rng); round.(all3.n[101:105])], century; replicates = 200, rng = MersenneTwister(1))
    even = ensemble(community(platte_thistle(aggregated = false), weevil), [whole(alone.n, rng); 2.0], century;
                    replicates = 200, rng = MersenneTwister(1))
    for (label, ens, printed) in
        (("thistle and weevil", two, (144.5, [96.0, 219.4], 0, 10)),
         ("with the weed", three, (74.0, [38.95, 112.1], 1, 0)),
         ("eggs spread evenly", even, (572.5, [0.0, 703.25], 24, 200)))
        o = outcome(ens)
        println("thistles, ensemble, $label: ", o)
        @test o.median_plants == printed[1]
        @test all(abs.(o.nine_in_ten_between .- printed[2]) .<= 1e-9)
        @test (o.thistle_lost, o.weevil_lost) == printed[3:4]
    end
end
