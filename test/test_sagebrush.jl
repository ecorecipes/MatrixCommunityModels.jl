# Four perennial plants of the sagebrush steppe, the multispecies IPM of Adler, Ellner & Levine (2010), compared
# with the R package IPMdoit run on the same parameters (test/reference/R/run_sagebrush_reference.R) and with the
# recruitment coefficients published by Chu & Adler (2015, Appendix H).
using LinearAlgebra

sagebrush_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]

@testset "sagebrush steppe (Adler et al. 2010; IPMdoit)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    kinds = (:artr, :heco, :pose, :pssp)
    codes = Dict("ARTR" => :artr, "HECO" => :heco, "POSE" => :pose, "PSSP" => :pssp)
    number(x) = parse(Float64, x)
    relative(a, b) = abs(a - b) / max(abs(b), 1e-300)
    # a state saved by the R script: genets per m² in each cell, for the species of `run`
    states = sagebrush_reference("sagebrush_state.csv")
    saved(run, k) = [number(r[4]) for r in states if r[1] == run && codes[r[2]] == k]
    covers(run) = [(year = parse(Int, r[2]), cover = number.(r[3:6])) for r in sagebrush_reference("sagebrush_runs.csv") if r[1] == run]
    community_state(c, run) = state_vector(c, (k => saved(run, k) for k in species_names(c))...)
    cover_vector(c, n, env = nothing) = [getproperty(pool_totals(c, n, env), Symbol(k, :_cover)) for k in kinds]

    @testset "parameters" begin
        @test collect(keys(SAGEBRUSH_PARAMETERS)) == collect(kinds)
        @test [SAGEBRUSH_PARAMETERS[k].cells for k in kinds] == [75, 75, 50, 50]
        @test all(length(SAGEBRUSH_PARAMETERS[k].recruitment.intercept_year) == 22 for k in kinds)
        @test all(length(SAGEBRUSH_PARAMETERS[k].growth.slope_year) == 22 for k in kinds)
        # Chu & Adler (2015), Appendix H, Idaho, recruitment: effect of the column species on the row species,
        # printed to four decimals or four figures
        published = [-0.4084 -0.5915 -0.4426 -0.2884; -0.03108 -1.7234 -0.2683 -0.3503;
                     -0.05736 -0.1888 -2.2718 -0.1254; -0.01905 -0.3080 -0.2088 -1.9332]
        for (i, k) in enumerate(kinds)
            @test all(abs.(SAGEBRUSH_PARAMETERS[k].recruitment.crowding .- published[i, :]) .<= 5e-5)
        end
        stage = sagebrush_stage(:artr)
        @test length(stage) == 75
        @test meshpoints(stage)[1] ≈ log(0.2) + step_size(stage) / 2
        @test sagebrush_state(:heco; cover = 0.02) isa Vector{Float64}
        @test pool_totals(sagebrush_plant(:heco), sagebrush_state(:heco; cover = 0.02), nothing).heco_cover ≈ 0.02
        @test_throws ArgumentError sagebrush_plant(:sage)
        @test_throws ArgumentError sagebrush_plant(:pose; niche_differences = (:seeds,))
        # listing every vital rate keeps all niche differences; listing none removes them all
        n = sagebrush_state(:pose; cover = 0.01)
        frozen(niche) = community_matrix(sagebrush_plant(:pose; niche_differences = niche), n, nothing)
        @test frozen((:growth, :survival, :recruitment)) == frozen(true)
        @test frozen(()) == frozen(false)
        @test_throws ArgumentError advance(sagebrush_plant(:pose), sagebrush_state(:pose), 23)
    end

    com = sagebrush_community()
    n30 = community_state(com, "community_30")

    @testset "crowding and kernels against IPMdoit" begin
        pools = pool_totals(com, n30, nothing)
        worst = 0.0
        rows = sagebrush_reference("sagebrush_crowding.csv")
        for k in kinds, (j, other) in enumerate(kinds)
            W = sagebrush_crowding(k, saved("community_30", k), pools)
            mine = [r for r in rows if codes[r[1]] == k && codes[r[2]] == other]
            for (column, M) in ((4, W.growth), (5, W.survival))
                # differences relative to the largest crowding of this pair of species, over the sizes of genet
                scale = maximum(number(r[column]) for r in mine)
                worst = max(worst, maximum(abs(M[parse(Int, r[3]), j] - number(r[column])) for r in mine) / scale)
            end
        end
        # IPMdoit integrates the crowding by conspecifics with R's integrate at its default tolerance (1.2e-4)
        println("sagebrush crowding: largest difference from IPMdoit, relative to the largest crowding = ", worst)
        @test worst < 3e-5
        # the largest of those differences, crowding of Pseudoroegneria in cell 33 by its own species (growth), is
        # R's: integrate with rel.tol = 1e-13 gives 0.273076048328859, where its default gives 0.273079825617989
        W = sagebrush_crowding(:pssp, saved("community_30", :pssp), pools)
        @test W.growth[33, 4] ≈ 0.273076048328859 rtol = 1e-12

        rows = sagebrush_reference("sagebrush_kernels.csv")
        entries, rates = 0.0, 0.0
        for niche in (true, false), year in (0, 3)
            c = niche ? com : sagebrush_community(; niche_differences = false)
            A = community_matrix(c, n30, year == 0 ? nothing : year)
            ends = cumsum([0; [length(stage_names(c, k)) for k in kinds]])
            for (i, k) in enumerate(kinds)
                K = A[(ends[i] + 1):ends[i + 1], (ends[i] + 1):ends[i + 1]]
                mine = [r for r in rows if codes[r[1]] == k && r[2] == (niche ? "niche" : "none") && parse(Int, r[3]) == year]
                entries = max(entries, maximum(abs(K[parse(Int, r[4]), parse(Int, r[5])] - number(r[6])) for r in mine if r[4] != "0"))
                rates = max(rates, relative(maximum(abs, eigvals(K)), number(only(r[6] for r in mine if r[4] == "0"))))
            end
        end
        println("sagebrush kernels: largest difference from IPMdoit = ", entries, "; growth rates ", rates)
        @test entries < 3e-8
        @test rates < 3e-7
    end

    @testset "runs against IPMdoit" begin
        # the community from IPMdoit's starting state, 30 mean years
        start = state_vector(com, (k => sagebrush_state(k) for k in kinds)...)
        sol = project(com, start, fill(nothing, 30))
        ref = covers("community_30")
        worst = maximum(maximum(relative.(cover_vector(com, sol.u[r.year + 1, :]), r.cover)) for r in ref)
        println("sagebrush, 30 mean years: largest relative difference in cover = ", worst)
        @test worst < 1e-6
        # 20 years drawn from the 22 years of the data, in the order the R script drew them
        years = [parse(Int, r[2]) for r in sagebrush_reference("sagebrush_years.csv")]
        sol = project(com, start, years[1:20])
        ref = [r for r in covers("community_years") if r.year <= 20]
        worst = maximum(maximum(relative.(cover_vector(com, sol.u[r.year + 1, :]), r.cover)) for r in ref)
        println("sagebrush, 20 fluctuating years: largest relative difference in cover = ", worst)
        @test worst < 1e-6
    end

    @testset "equilibrium and invasion against IPMdoit" begin
        # IPMdoit's state after 1500 mean years is an equilibrium here too, to within the error of its integration
        eq = community_state(com, "community")
        @test maximum(abs, advance(com, eq, nothing) .- eq) / maximum(eq) < 1e-7
        @test cover_vector(com, eq) ≈ last(covers("community")).cover rtol = 1e-12
        # each species invades the other three at their equilibrium; without niche differences none does
        for r in sagebrush_reference("sagebrush_invasion.csv")
            niche, k = r[1] == "niche", codes[r[2]]
            others = Tuple(x for x in kinds if x != k)
            resident = sagebrush_community(; niche_differences = niche, kinds = others)
            n = community_state(resident, "residents_" * r[1] * "_without_" * r[2])
            rate = invasion_growth_rate(sagebrush_plant(k; niche_differences = niche), resident, n, nothing)
            @test rate ≈ number(r[3]) rtol = 1e-6
            @test niche ? rate > 1 : rate < 1
        end
    end

    @testset "exclusion without niche differences, mean and drawn years" begin
        # From IPMdoit's equilibrium of the full model, the first year each grass is below a hundredth of its
        # starting cover. IPMdoit, run from the same state with the niche differences removed as in
        # run_sagebrush_reference.R, gives 80, 16 and 78 years in the mean year and 43, 15 and 36 along the R
        # script's drawn years (sagebrush_years.csv, set.seed(2010)). Each crossing clears the threshold by at
        # least 1%.
        flat = sagebrush_community(; niche_differences = false)
        n0 = community_state(flat, "community")
        c0 = cover_vector(flat, n0)
        function excluded(envs)
            sol = project(flat, n0, envs)
            C = reduce(hcat, cover_vector(flat, row) for row in eachrow(sol.u))
            return [findfirst(<(0.01 * c0[j]), C[j, :]) - 1 for j in 2:4]
        end
        years = [parse(Int, r[2]) for r in sagebrush_reference("sagebrush_years.csv")]
        mean_years, drawn = excluded(fill(nothing, 82)), excluded(years[1:45])
        @test mean_years == [80, 16, 78]
        @test drawn == [43, 15, 36]
        # drawing the years brings the exclusion forward, here by about half for the two slower grasses
        @test all(drawn .<= mean_years)
        @test all(drawn[[1, 3]] .< 0.6 .* mean_years[[1, 3]])
    end
    BLAS.set_num_threads(threads)
end
