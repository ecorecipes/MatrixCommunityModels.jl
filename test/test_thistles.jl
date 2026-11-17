# A native thistle (an integral projection model), an invasive thistle (a matrix model) and their weevil
using LinearAlgebra
using Random

@testset "thistles and weevil" begin
    dominant(sp) = maximum(abs, eigvals(community_matrix(sp, zeros(length(sp)), nothing)))
    # a stand of seedlings, as a starting point
    function seedlings(total; cells = 100)
        stage = platte_thistle_stage(; cells)
        return total .* step_size(stage) .* [exp(-0.5 * (x - 0.75)^2 / 0.17) / sqrt(2pi * 0.17) for x in meshpoints(stage)]
    end
    P, M, W = PLATTE_THISTLE_PARAMETERS, MUSK_THISTLE_PARAMETERS, WEEVIL_PARAMETERS

    @testset "Platte thistle (Rose et al. 2005)" begin
        native = platte_thistle()
        @test stage_names(native, :platte_thistle) == stage_names(platte_thistle_stage())
        @test Set(pool_names(native)) == Set([:heads, :rearing, :eggs])
        # Growth rates without density dependence. Published: 1.38 without insects and 1.05 with native insects. The
        # tabulated chance of establishment, 0.067, gives 1.33 and 1.02 (vignette 14 prints them in full); 0.077 gives
        # both published values, which are printed to two decimals.
        @test dominant(platte_thistle(native_insects = false, seed_limited = false)) ≈ 1.3293045416613871 rtol = 1e-12  # achieved 0 (2026-11-17)
        @test dominant(platte_thistle(seed_limited = false)) ≈ 1.019338009700462 rtol = 1e-12                           # achieved 0 (2026-11-17)
        more = merge(P, (establishment = 0.077,))
        @test dominant(platte_thistle(parameters = more, native_insects = false, seed_limited = false)) ≈ 1.38 atol = 5e-3
        @test dominant(platte_thistle(parameters = more, seed_limited = false)) ≈ 1.05 atol = 5e-3
        # Equilibria with seed-limited recruitment. Published: about 5000 plants without insects, about 600 with;
        # vignette 14 prints 5079 and 603.
        bare = equilibrium(platte_thistle(native_insects = false), seedlings(600), nothing)
        alone = equilibrium(native, seedlings(600), nothing)
        @test bare.converged && alone.converged
        @test sum(bare.n) ≈ 5079 atol = 0.5
        @test sum(alone.n) ≈ 603 atol = 0.5
        # Flowering plants are about 13 mm across, as observed (vignette: 12.8 mm, 0.061 flowering plants per m²)
        z = meshpoints(platte_thistle_stage())
        flowering = [1 / (1 + exp(-(P.survival[1] + P.survival[2] * x))) / (1 + exp(-(P.flowering[1] + P.flowering[2] * x))) for x in z] .* alone.n
        @test sum(flowering .* exp.(z)) / sum(flowering) ≈ 12.8 atol = 0.05
        @test sum(flowering) / 720 ≈ 0.061 atol = 5e-4
        # and the average flowering plant offers 10.4 heads, which can rear 4.2 weevils each
        pools = pool_totals(native, alone.n, nothing)
        @test pools.heads / sum(flowering) ≈ 10.4 atol = 0.05      # achieved 10.42 (2026-11-17)
        @test pools.rearing / pools.heads ≈ P.weevils_per_head
        @test pools.eggs == 0

        # The weevil as observed in 2000 (Table 1 of the paper): eggs on a plant = exp(-0.52 + 1.71 x). Published: a
        # growth rate below one, 69% of plants lost by 2002, up to 95% if eggs were spread evenly; vignette 14 prints
        # 0.9, 67% and 93%.
        eggs = x -> exp(-0.52 + 1.71x)
        @test dominant(platte_thistle(seed_limited = false, eggs = eggs)) ≈ 0.90 atol = 5e-3
        uneven = equilibrium(platte_thistle(eggs = eggs), seedlings(600), nothing)
        even = equilibrium(platte_thistle(eggs = eggs, aggregated = false), seedlings(600), nothing)
        @test 1 - sum(uneven.n) / sum(alone.n) ≈ 0.67 atol = 5e-3
        @test 1 - sum(even.n) / sum(alone.n) ≈ 0.93 atol = 5e-3
        @test isempty(pool_names(platte_thistle(eggs = eggs)) ∩ [:eggs])
        # the chance of establishment is capped for small seed crops
        few = zeros(100)
        few[70] = 1e-6
        net = system(native)[:platte_thistle]
        _, births_few = survival_fecundity_matrices(net, context(net, few, (heads = 0.0, rearing = 0.0, eggs = 0.0), nothing))
        _, births_all = survival_fecundity_matrices(system(platte_thistle(seed_limited = false))[:platte_thistle],
                                                    context(net, few, (heads = 0.0, rearing = 0.0, eggs = 0.0), nothing))
        @test sum(births_few) / sum(births_all) ≈ P.greatest_establishment / P.establishment
        # the mesh does not matter
        coarse = equilibrium(platte_thistle(cells = 50), seedlings(600; cells = 50), nothing)
        @test sum(coarse.n) ≈ sum(alone.n) rtol = 1e-4             # achieved 8.3e-6 (2026-11-17), the error of the mesh
    end

    @testset "musk thistle (Jongejans et al. 2008, Kansas rates)" begin
        weed = musk_thistle()
        @test stage_names(weed, :musk_thistle) == [:seed_bank, :small, :medium, :large]
        @test Set(transition_names(system(weed)[:musk_thistle])) == Set([:dormancy, :emerging, :growing, :seed_rain, :establishing])
        # the published matrix and its growth rate, 1.748 (Zhang, Jongejans & Shea 2011)
        A = community_matrix(weed, zeros(4), nothing)
        # (printed to four decimals; every entry rounds to the table's, the nearest to the edge being 0.0249505 for 0.0250)
        @test maximum(abs, A .- [0.2597 0.0916 9.2153 256.8447; 0.0251 0.0250 0.9922 27.6540; 0.0049 0.1739 0.4156 5.3736; 0 0.0272 0.1379 0]) <= 5e-5
        @test dominant(weed) ≈ 1.748 atol = 5e-4
        @test dominant(weed) ≈ 1.7482683163797041 rtol = 1e-12     # as vignette 14 prints it; achieved 0 (2026-11-17)
        # open ground bounds it: rosettes never outnumber the sites
        sol = project(weed, [1e6, 0.0, 0.0, 0.0], fill(nothing, 200))
        rosettes = stage_total(sol, :musk_thistle, [:small, :medium, :large])
        @test maximum(rosettes) <= M.sites
        @test rosettes[end] ≈ rosettes[end - 1] rtol = 1e-6
        @test 0.4 < rosettes[end] / M.sites < 0.6
        # twice the sites, twice the weed
        double = equilibrium(musk_thistle(parameters = merge(M, (sites = 2 * M.sites,))), [1000.0, 10.0, 5.0, 1.0], nothing)
        @test sum(double.n[2:4]) ≈ 2 * rosettes[end] rtol = 1e-9   # achieved 4.1e-12 (2026-11-17)
    end

    @testset "the weevil and its hosts" begin
        native, weed, weevil = platte_thistle(cells = 50), musk_thistle(), flowerhead_weevil()
        alone = equilibrium(native, seedlings(600; cells = 50), nothing)
        stand = equilibrium(weed, [1000.0, 10.0, 5.0, 1.0], nothing)
        rare = W.eggs * W.survival * W.overwintering
        # a rare weevil multiplies at the same rate on either host, and dies out without one
        @test invasion_growth_rate(weevil, native, alone.n, nothing) ≈ rare
        @test invasion_growth_rate(weevil, weed, stand.n, nothing) ≈ rare
        @test invasion_growth_rate(weevil, native, zeros(50), nothing) == 0

        # the weed alone, as vignette 14 prints it: 72469 seeds in the bank, 2993 rosettes, 0.45 bolting plants per m²
        @test round(Int, stand.n[1]) == 72469 && round(Int, sum(stand.n[2:4])) == 2993
        @test sum(M.survival[j] * M.bolting[j] * stand.n[j + 1] for j in 1:3) / 720 ≈ 0.45 atol = 5e-3

        # native thistle and weevil: the thistle falls to about a quarter. Vignette 14 prints 136.93 thistles and 20.68
        # weevils on 50 cells (and on 100 and 200).
        pair = community(native, weevil)
        @test Set(pool_names(pair)) == Set([:heads, :rearing, :eggs])
        settled = equilibrium(pair, [alone.n; 2.0], nothing)
        @test settled.converged
        @test sum(settled.n[1:50]) ≈ 136.93 atol = 0.005
        @test settled.n[51] ≈ 20.68 atol = 0.005
        pools = pool_totals(pair, settled.n, nothing)
        reared = W.survival * pools.eggs / (1 + W.survival * pools.eggs / pools.rearing)
        @test reared / pools.rearing ≈ 1 - 1 / rare                           # weevils reared, as a share of what the heads can hold
        @test reared / pools.heads * 10.4 / 3.8 ≈ 5.0 atol = 0.05             # about five in each head that flowers (achieved 4.96)
        # with eggs spread evenly the weevil destroys the thistle, and then itself
        doomed = project(community(platte_thistle(cells = 50, aggregated = false), weevil), [alone.n; 2.0], fill(nothing, 60))
        @test abundance(doomed, :platte_thistle)[end] < 1e-3 && abundance(doomed, :weevil)[end] < 1e-3

        # weed and weevil: a third of the seed is lost, and the weed barely notices (vignette 14: 2690 rosettes, 274
        # weevils, 59.4 eggs in a weed head, 10.8 larvae, 36% of the seed lost)
        biocontrol = community(weed, weevil)
        held = equilibrium(biocontrol, [stand.n; 5.0], nothing)
        @test held.converged
        @test round(Int, sum(held.n[2:4])) == 2690 && round(Int, held.n[5]) == 274
        on_weed = pool_totals(biocontrol, held.n, nothing)
        larvae = W.survival * M.attraction * on_weed.eggs / on_weed.heads / (1 + W.survival * on_weed.eggs / on_weed.rearing)
        @test M.attraction * on_weed.eggs / on_weed.heads ≈ 59.4 atol = 0.05
        @test larvae ≈ 10.8 atol = 0.05
        @test M.exposed * min(1, larvae * M.seeds_per_larva / M.seeds_per_head) ≈ 0.36 atol = 5e-3

        # all three: the weed's heads rear more weevils for the eggs they draw, so the native carries more eggs
        # (vignette 14: 10.0 eggs in a native head with the weevil alone, 24.2 beside the weed, 81 native plants
        # rather than 137)
        @test M.weevils_per_head / M.attraction > P.weevils_per_head
        @test on_weed.eggs / on_weed.heads > 2 * pools.eggs / pools.heads
        trio = community(native, weed, weevil)
        all3 = equilibrium(trio, [settled.n[1:50]; held.n[1:4]; held.n[5]], nothing)
        @test all3.converged
        together = pool_totals(trio, all3.n, nothing)
        @test pools.eggs / pools.heads ≈ 10.0 atol = 0.05
        @test together.eggs / together.heads ≈ 24.2 atol = 0.05
        @test pools.eggs / pools.heads < together.eggs / together.heads < on_weed.eggs / on_weed.heads
        @test round(Int, sum(all3.n[1:50])) == 81                            # apparent competition
        # but in the year a weed stand comes within reach, the same weevils spread their eggs over more heads (1.6)
        first_year = pool_totals(trio, [settled.n[1:50]; held.n[1:4]; settled.n[51]], nothing)
        @test first_year.eggs / first_year.heads ≈ 1.6 atol = 0.05

        # can the native thistle return when rare? Yes if eggs are spread unevenly, as observed; no if evenly. Vignette
        # 14 prints 1.2992491766060492 and 0.6196951337217237 on 100 cells; 50 cells differ by 2e-8.
        @test invasion_growth_rate(native, biocontrol, held.n, nothing) ≈ 1.2992491766060492 rtol = 1e-6
        @test invasion_growth_rate(platte_thistle(cells = 50, aggregated = false), biocontrol, held.n, nothing) ≈ 0.6196951337217237 rtol = 1e-6
        @test invasion_growth_rate(weed, pair, settled.n, nothing) > 1.5

        # whole plants and weevils
        n0 = [round.(alone.n); round.(stand.n); 20.0]
        random = project(trio, n0, fill(nothing, 10); rng = MersenneTwister(2))
        @test all(isinteger, random.u) && all(>=(0), random.u)
        @test abundance(random, :weevil)[end] > 0
    end
end
