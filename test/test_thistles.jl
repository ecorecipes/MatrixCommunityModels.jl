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
        # tabulated chance of establishment, 0.067, gives 1.33 and 1.02; 0.077 gives both published values.
        @test dominant(platte_thistle(native_insects = false, seed_limited = false)) ≈ 1.329 atol = 2e-3
        @test dominant(platte_thistle(seed_limited = false)) ≈ 1.019 atol = 2e-3
        more = merge(P, (establishment = 0.077,))
        @test dominant(platte_thistle(parameters = more, native_insects = false, seed_limited = false)) ≈ 1.38 atol = 5e-3
        @test dominant(platte_thistle(parameters = more, seed_limited = false)) ≈ 1.05 atol = 5e-3
        # Equilibria with seed-limited recruitment. Published: about 5000 plants without insects, about 600 with.
        bare = equilibrium(platte_thistle(native_insects = false), seedlings(600), nothing)
        alone = equilibrium(native, seedlings(600), nothing)
        @test bare.converged && alone.converged
        @test 4900 < sum(bare.n) < 5200
        @test 590 < sum(alone.n) < 615
        # Flowering plants are about 13 mm across, as observed
        z = meshpoints(platte_thistle_stage())
        flowering = [1 / (1 + exp(-(P.survival[1] + P.survival[2] * x))) / (1 + exp(-(P.flowering[1] + P.flowering[2] * x))) for x in z] .* alone.n
        @test 12 < sum(flowering .* exp.(z)) / sum(flowering) < 13.5
        # and the average flowering plant offers 10.4 heads, which can rear 4.2 weevils each
        pools = pool_totals(native, alone.n, nothing)
        @test pools.heads / sum(flowering) ≈ 10.4 rtol = 0.01
        @test pools.rearing / pools.heads ≈ P.weevils_per_head
        @test pools.eggs == 0

        # The weevil as observed in 2000 (Table 1 of the paper): eggs on a plant = exp(-0.52 + 1.71 x)
        eggs = x -> exp(-0.52 + 1.71x)
        @test dominant(platte_thistle(seed_limited = false, eggs = eggs)) < 0.95
        uneven = equilibrium(platte_thistle(eggs = eggs), seedlings(600), nothing)
        even = equilibrium(platte_thistle(eggs = eggs, aggregated = false), seedlings(600), nothing)
        @test 0.6 < 1 - sum(uneven.n) / sum(alone.n) < 0.72        # published: 69% by 2002
        @test 0.9 < 1 - sum(even.n) / sum(alone.n) < 0.96          # published: up to 95%
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
        @test sum(coarse.n) ≈ sum(alone.n) rtol = 1e-3
    end

    @testset "musk thistle (Jongejans et al. 2008, Kansas rates)" begin
        weed = musk_thistle()
        @test stage_names(weed, :musk_thistle) == [:seed_bank, :small, :medium, :large]
        @test Set(transition_names(system(weed)[:musk_thistle])) == Set([:dormancy, :emerging, :growing, :seed_rain, :establishing])
        # the published matrix and its growth rate, 1.748 (Zhang, Jongejans & Shea 2011)
        A = community_matrix(weed, zeros(4), nothing)
        @test A ≈ [0.2597 0.0916 9.2153 256.8447; 0.0251 0.0250 0.9922 27.6540; 0.0049 0.1739 0.4156 5.3736; 0 0.0272 0.1379 0] atol = 1e-4
        @test dominant(weed) ≈ 1.748 atol = 5e-4
        # open ground bounds it: rosettes never outnumber the sites
        sol = project(weed, [1e6, 0.0, 0.0, 0.0], fill(nothing, 200))
        rosettes = stage_total(sol, :musk_thistle, [:small, :medium, :large])
        @test maximum(rosettes) <= M.sites
        @test rosettes[end] ≈ rosettes[end - 1] rtol = 1e-6
        @test 0.4 < rosettes[end] / M.sites < 0.6
        # twice the sites, twice the weed
        double = equilibrium(musk_thistle(parameters = merge(M, (sites = 2 * M.sites,))), [1000.0, 10.0, 5.0, 1.0], nothing)
        @test sum(double.n[2:4]) ≈ 2 * rosettes[end] rtol = 1e-3
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

        # native thistle and weevil
        pair = community(native, weevil)
        @test Set(pool_names(pair)) == Set([:heads, :rearing, :eggs])
        settled = equilibrium(pair, [alone.n; 2.0], nothing)
        @test settled.converged
        @test 0.7 < 1 - sum(settled.n[1:50]) / sum(alone.n) < 0.85            # the thistle falls to about a quarter
        pools = pool_totals(pair, settled.n, nothing)
        reared = W.survival * pools.eggs / (1 + W.survival * pools.eggs / pools.rearing)
        @test reared / pools.rearing ≈ 1 - 1 / rare                           # weevils reared, as a share of what the heads can hold
        @test 4.5 < reared / pools.heads * 10.4 / 3.8 < 5.5                   # about five in each head that flowers
        # with eggs spread evenly the weevil destroys the thistle, and then itself
        doomed = project(community(platte_thistle(cells = 50, aggregated = false), weevil), [alone.n; 2.0], fill(nothing, 60))
        @test abundance(doomed, :platte_thistle)[end] < 1e-3 && abundance(doomed, :weevil)[end] < 1e-3

        # weed and weevil: a third of the seed is lost, and the weed barely notices
        biocontrol = community(weed, weevil)
        held = equilibrium(biocontrol, [stand.n; 5.0], nothing)
        @test held.converged
        @test 0.85 < sum(held.n[2:4]) / sum(stand.n[2:4]) < 0.95
        on_weed = pool_totals(biocontrol, held.n, nothing)

        # all three: the weed's heads rear more weevils for the eggs they draw, so the native carries more eggs
        @test M.weevils_per_head / M.attraction > P.weevils_per_head
        @test on_weed.eggs / on_weed.heads > 2 * pools.eggs / pools.heads
        trio = community(native, weed, weevil)
        all3 = equilibrium(trio, [settled.n[1:50]; held.n[1:4]; held.n[5]], nothing)
        @test all3.converged
        together = pool_totals(trio, all3.n, nothing)
        @test pools.eggs / pools.heads < together.eggs / together.heads < on_weed.eggs / on_weed.heads
        @test sum(all3.n[1:50]) < 0.7 * sum(settled.n[1:50])                  # apparent competition
        # but in the year a weed stand comes within reach, the same weevils spread their eggs over more heads
        first_year = pool_totals(trio, [settled.n[1:50]; held.n[1:4]; settled.n[51]], nothing)
        @test first_year.eggs / first_year.heads < 0.2 * pools.eggs / pools.heads

        # can the native thistle return when rare? Yes if eggs are spread unevenly, as observed; no if evenly
        @test invasion_growth_rate(native, biocontrol, held.n, nothing) > 1.2
        @test invasion_growth_rate(platte_thistle(cells = 50, aggregated = false), biocontrol, held.n, nothing) < 0.7
        @test invasion_growth_rate(weed, pair, settled.n, nothing) > 1.5

        # whole plants and weevils
        n0 = [round.(alone.n); round.(stand.n); 20.0]
        random = project(trio, n0, fill(nothing, 10); rng = MersenneTwister(2))
        @test all(isinteger, random.u) && all(>=(0), random.u)
        @test abundance(random, :weevil)[end] > 0
    end
end
