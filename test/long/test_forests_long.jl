# Forest trees sharing basal area on the 700-cell mesh of matreex (vignette 15): where each species settles, the long
# runs against matreex 0.4.0 with its default settings (test/reference/R/run_forest_reference.R), and invasion.
using LinearAlgebra

@testset "forests on 700 cells, long runs (vignette 15)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    kinds = (:spruce, :fir, :beech)
    number(x) = x == "NA" ? NaN : parse(Float64, x)
    growth_rate(A) = maximum(abs, eigvals(A))
    tree_cells(sp, kind; cells = 700) = findall(in(stage_names(forest_stage(kind; cells))), stage_names(sp, kind))
    held(sol, kind) = [pool_totals(subcommunity(sol.community, kind), row, nothing).basal_area for row in eachrow(sol[kind])]
    # a run of matreex: basal areas of spruce, fir and beech every five years, and where it settled
    function matreex_run(run; file = "forest_runs_700.csv")
        rows = [r for r in long_reference(file) if r[1] == run]
        recorded = [r for r in rows if r[2] != "-1"]
        return (years = [parse(Int, r[2]) for r in recorded], basal_area = [number(r[j]) for r in recorded, j in 3:5],
                settled = [number(r[j]) for r in rows if r[2] == "-1" for j in 3:5])
    end
    settles(sp, kind; cells = 700) =
        equilibrium(sp, :basal_area => (30.0, 45.0), forest_stand(kind; cells, basal_area = 30.0), nothing)
    trees = (spruce = forest_tree(:spruce), fir = forest_tree(:fir), beech = forest_tree(:beech))

    # where each species settles alone, without iterating the dynamics
    settled = map(k -> settles(trees[k], k), kinds)
    stands = NamedTuple{kinds}(Tuple(s.n for s in settled))
    basal(kind) = pool_totals(trees[kind], stands[kind], nothing).basal_area

    @testset "where each species settles alone" begin
        # vignette 15 prints 36.670, 42.604 and 42.546 m²/ha (matreex 36.669, 42.603, 42.546); 251, 205 and 306 trees
        # and 40, 60 and 137 seedlings and saplings per ha; a year moves the stand by 3e-14, 6e-15 and 6e-15.
        # Achieved (2026-11-18): relative difference from matreex 2.3e-5, 2.7e-5 and 3.7e-6, which matreex's
        # interpolation and quadrature set; the bound is the vignette's own claim, "three parts in 100,000".
        printed = (spruce = (36.670, 251, 40), fir = (42.604, 205, 60), beech = (42.546, 306, 137))
        for (s, kind) in zip(settled, kinds)
            @test s.converged
            n, sp = stands[kind], trees[kind]
            standing = sum(n[tree_cells(sp, kind)])
            theirs = only(filter(!isnan, matreex_run(String(kind)).settled))
            moved = maximum(abs, advance(sp, n, nothing) .- n) / maximum(n)
            println("forest, $kind settles at ", basal(kind), " m²/ha (matreex ", theirs, "); ", standing, " trees and ",
                    sum(n) - standing, " seedlings and saplings per ha; a year moves the stand by ", moved)
            @test basal(kind) ≈ printed[kind][1] atol = 5e-4
            @test abs(basal(kind) / theirs - 1) < 3e-5
            @test standing ≈ printed[kind][2] atol = 0.5
            @test sum(n) - standing ≈ printed[kind][3] atol = 0.5
            @test moved < 1e-12                       # achieved 3e-14, 6e-15 and 6e-15 (2026-11-18)
        end
    end

    @testset "spruce alone over 500 years against matreex" begin
        # vignette 15: largest relative difference 4.8e-4 (year 50), a peak of 51.3 m²/ha in year 174; the main suite
        # has the first 100 years
        alone = project(trees.spruce, forest_stand(:spruce), fill(nothing, 500))
        ref = matreex_run("spruce")
        apart = abs.(alone.pools[ref.years .+ 1, 1] ./ ref.basal_area[:, 1] .- 1)
        println("forest, spruce alone, 500 years: largest relative difference from matreex = ", maximum(apart),
                " in year ", ref.years[argmax(apart)])
        @test maximum(apart) < 7e-4                   # achieved 4.8e-4 (2026-11-18), matreex's coarse rule and interpolation
        @test ref.years[argmax(apart)] == 50
        @test round(maximum(alone.pools[1:500, 1]); digits = 1) == 51.3
        @test argmax(alone.pools[1:500, 1]) - 1 == 174
    end

    @testset "three species over 400 years against matreex" begin
        forest = community(trees...)
        young = reduce(vcat, [forest_stand(k) for k in kinds])          # 1 m²/ha of each
        together = project(forest, young, fill(nothing, 400))
        ref = matreex_run("trio")
        shown = ref.years .<= 400
        # vignette 15: 2.2e-4, 7.3e-4 and 1.5e-4
        for (j, (kind, tol)) in enumerate(zip(kinds, (3.5e-4, 1.1e-3, 2.5e-4)))
            mine = held(together, kind)[ref.years[shown] .+ 1]
            worst = maximum(abs.(mine ./ ref.basal_area[shown, j] .- 1))
            println("forest, three species, 400 years, $kind: largest relative difference from matreex = ", worst)
            @test worst < tol                         # achieved 2.2e-4, 7.3e-4 and 1.5e-4 (2026-11-18)
        end
        # each species frozen at matreex's equilibrium of the three just replaces itself: vignette 15 prints
        # 1.0000001, 1.0000002 and 1.0000000
        there = [number(r[4]) for r in long_reference("forest_states_700.csv") if r[1] == "trio"]
        A = community_matrix(forest, there, nothing)
        ends = cumsum([0; [length(sp) for sp in trees]])
        for (j, printed) in enumerate((1.0000001, 1.0000002, 1.0000000))
            block = (ends[j] + 1):ends[j + 1]
            @test growth_rate(A[block, block]) ≈ printed atol = 5e-8
        end
    end

    @testset "who can invade whom" begin
        # growth rate of an invader holding 0.1 m²/ha in the settled stand of each other species, as vignette 15
        # prints it (rows invader, columns resident)
        printed = Dict((:spruce, :fir) => 1.0219, (:spruce, :beech) => 1.0220, (:fir, :spruce) => 1.0501,
                       (:fir, :beech) => 1.0442, (:beech, :spruce) => 1.0343, (:beech, :fir) => 1.0305)
        for ((invader, resident), rate) in printed
            @test invasion_growth_rate(forest_tree(invader; conspecific = 0.1), trees[resident], stands[resident], nothing) ≈
                  rate atol = 5e-5
        end
        for (b, rate) in ((0.01, 1.0975), (1.0, 1.0227), (10.0, 1.0080))
            @test invasion_growth_rate(forest_tree(:fir; conspecific = b), trees.spruce, stands.spruce, nothing) ≈ rate atol = 5e-5
        end
        # fir added to the settled spruce stand as 0.1 m²/ha of small trees: 0.61, 2.57 and 12.75 m²/ha after 50, 100
        # and 200 years, 1.0291 a year between years 50 and 100
        arrival = project(community(trees.spruce, trees.fir), [stands.spruce; forest_stand(:fir; basal_area = 0.1)],
                          fill(nothing, 200))
        fir = held(arrival, :fir)
        @test round.(fir[[51, 101, 201]]; digits = 2) == [0.61, 2.57, 12.75]
        @test (fir[101] / fir[51])^(1 / 50) ≈ 1.0291 atol = 5e-5
    end

    @testset "the mesh" begin
        # where each species settles on coarser meshes, as vignette 15 prints it (the main suite has spruce on 100,
        # 200 and 400 cells)
        for (cells, printed) in ((100, (37.51, 43.96, 40.51)), (200, (36.90, 43.36, 44.32)), (350, (36.69, 42.76, 43.53)))
            found = [settles(forest_tree(k; cells), k; cells).total for k in kinds]
            println("forest, settled basal areas on $cells cells: ", found)
            @test all(abs.(found .- collect(printed)) .<= 5e-3)
        end
        @test all(abs.([basal(k) for k in kinds] .- [36.67, 42.60, 42.55]) .<= 5e-3)
    end
    BLAS.set_num_threads(threads)
end
