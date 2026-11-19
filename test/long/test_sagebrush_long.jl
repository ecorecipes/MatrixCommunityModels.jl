# Long runs of the sagebrush steppe IPM of Adler, Ellner & Levine (2010) against the R package IPMdoit
# (test/reference/R/run_sagebrush_reference.R), and the exclusion of the grasses in drawn years (vignette 17).
using LinearAlgebra, Random

@testset "sagebrush steppe, long runs (vignette 17)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    kinds = (:artr, :heco, :pose, :pssp)
    codes = Dict("ARTR" => :artr, "HECO" => :heco, "POSE" => :pose, "PSSP" => :pssp)
    number(x) = parse(Float64, x)
    states = long_reference("sagebrush_state.csv")
    saved(run, k) = [number(r[4]) for r in states if r[1] == run && codes[r[2]] == k]
    covers(run) = [(year = parse(Int, r[2]), cover = number.(r[3:6])) for r in long_reference("sagebrush_runs.csv") if r[1] == run]
    community_state(c, run) = state_vector(c, (k => saved(run, k) for k in species_names(c))...)
    cover(c, n, env = nothing) = [getproperty(pool_totals(c, n, env), Symbol(k, :_cover)) for k in kinds]
    years = [parse(Int, r[2]) for r in long_reference("sagebrush_years.csv")]
    four = sagebrush_community()
    start = state_vector(four, (k => sagebrush_state(k) for k in kinds)...)

    @testset "1000 mean years and 300 drawn years against IPMdoit" begin
        # vignette 17: below 2e-6 (1.19e-6 and 1.40e-6); the differences are IPMdoit's integration of the crowding
        sol = project(four, start, fill(nothing, 1000))
        worst = maximum(maximum(abs.(cover(four, sol.u[r.year + 1, :]) ./ r.cover .- 1)) for r in covers("community") if r.year <= 1000)
        println("sagebrush, 1000 mean years: largest relative difference from IPMdoit = ", worst)
        @test worst < 2e-6                  # achieved 1.19e-6 (2026-11-18)
        wild = project(four, start, years[1:300])
        C = reduce(vcat, permutedims(cover(four, row)) for row in eachrow(wild.u))
        worst = maximum(maximum(abs.(C[r.year + 1, :] ./ r.cover .- 1)) for r in covers("community_years"))
        println("sagebrush, 300 drawn years: largest relative difference from IPMdoit = ", worst)
        @test worst < 2e-6                  # achieved 1.40e-6 (2026-11-18)
        # mean cover (%) after year 100, as vignette 17 prints it
        @test all(abs.(100 .* vec(sum(C[102:end, :]; dims = 1)) ./ 200 .- [6.23, 1.53, 0.80, 1.43]) .<= 5e-3)
    end

    @testset "each species alone over 1500 years against IPMdoit" begin
        # the whole of IPMdoit's runs (every fifth year); the main suite has only where they end
        worst = 0.0
        for (j, k) in enumerate(kinds)
            ref = covers("alone_" * uppercase(String(k)))
            sol = project(sagebrush_plant(k), sagebrush_state(k), fill(nothing, 1500))
            mine = [getproperty(pool_totals(sagebrush_plant(k), row, nothing), Symbol(k, :_cover)) for row in eachrow(sol.u)]
            worst = max(worst, maximum(abs(mine[r.year + 1] / r.cover[j] - 1) for r in ref))
        end
        println("sagebrush, each species alone, 1500 mean years: largest relative difference from IPMdoit = ", worst)
        @test worst < 3e-6                  # achieved 1.5e-6 (2026-11-18), IPMdoit's integration (its end: 1.5e-6 too)
    end

    @testset "stochastic invasion rates against IPMdoit" begin
        # the long-run growth rate of each rare invader along the R script's 600 drawn years, the first 100 left out
        # (the main suite has the shrub). Vignette 17 prints 1.092135, 1.309995, 1.521909 and 1.245614, equal to
        # IPMdoit's to the six decimals printed.
        worst = 0.0
        for r in long_reference("sagebrush_stochastic_invasion.csv")
            k = codes[r[1]]
            others = Tuple(x for x in kinds if x != k)
            resident = sagebrush_community(; kinds = others)
            rate = stochastic_invasion_rate(sagebrush_plant(k), resident, community_state(resident, "residents_niche_without_" * r[1]),
                                            years; burnin = parse(Int, r[2]))
            worst = max(worst, abs(rate / number(r[4]) - 1))
            @test rate ≈ number(r[4]) atol = 5e-7
        end
        println("sagebrush, stochastic invasion rates: largest relative difference from IPMdoit = ", worst)
        @test worst < 1e-7                  # achieved 3.3e-8 (2026-11-18), IPMdoit's integration
    end

    @testset "exclusion without niche differences in drawn years" begin
        # From IPMdoit's equilibrium of the full model, with the niche differences removed, the first year each grass
        # is below a hundredth (and a ten-thousandth) of its starting cover, as in vignette 17.
        flat = sagebrush_community(; niche_differences = false)
        together = community_state(flat, "community")
        c0 = cover(flat, together)
        function exclusion_years(envs)
            n = copy(together)
            when = zeros(Int, 3, 2)
            for (t, env) in enumerate(envs)
                n = advance(flat, n, env)
                c = cover(flat, n)
                for j in 1:3, (s, share) in enumerate((1e-2, 1e-4))
                    when[j, s] == 0 && c[j + 1] < share * c0[j + 1] && (when[j, s] = t)
                end
                all(>(0), when) && break
            end
            return when
        end
        in_mean_years = exclusion_years(fill(nothing, 600))
        @test in_mean_years == [80 109; 16 46; 78 115]
        drawn = [exclusion_years(rand(Xoshiro(seed), 1:22, 1500)) for seed in 1:200]
        @test all(d -> all(>(0), d), drawn)
        # Vignette 17, seeds 1 to 40: medians (ranges) 35.0 (13 to 119), 11.0 (5 to 25), 36.5 (17 to 81) years to a
        # hundredth; 50.0 (21 to 137), 25.0 (10 to 52), 56.0 (25 to 124) to a ten-thousandth; in 8 of the 40 some grass
        # lasts longer than in mean years.
        first40 = drawn[1:40]
        times(ds, j, s) = [d[j, s] for d in ds]
        @test [long_median(times(first40, j, 1)) for j in 1:3] == [35.0, 11.0, 36.5]
        @test [extrema(times(first40, j, 1)) for j in 1:3] == [(13, 119), (5, 25), (17, 81)]
        @test [long_median(times(first40, j, 2)) for j in 1:3] == [50.0, 25.0, 56.0]
        @test [extrema(times(first40, j, 2)) for j in 1:3] == [(21, 137), (10, 52), (25, 124)]
        @test count(d -> any(d[:, 1] .> in_mean_years[:, 1]), first40) == 8
        # Seeds 1 to 200 (vignette 17, and IPMdoit on the same sequences, which gives the same year in every case):
        # medians 39, 11 and 41 years to a hundredth, against 80, 16 and 78 in mean years; one sequence keeps
        # H. comata for 905 years, and in the others none lasts beyond 270.
        medians = [long_median(times(drawn, j, 1)) for j in 1:3]
        println("sagebrush, exclusion over 200 drawn sequences: medians ", medians, " (a ten-thousandth: ",
                [long_median(times(drawn, j, 2)) for j in 1:3], "); longest ", sort(vec(reduce(hcat, times(drawn, j, 1) for j in 1:3)))[(end - 2):end])
        @test medians == [39.0, 11.0, 41.0]
        @test [long_median(times(drawn, j, 2)) for j in 1:3] == [52.0, 26.0, 59.0]     # to a ten-thousandth (mean years 109, 46, 115)
        # the 95% ranges of the lab notes (the 5th and 195th of the 200, sorted): 16 to 117, 5 to 21 and 21 to 94
        @test [sort(times(drawn, j, 1))[[5, 195]] for j in 1:3] == [[16, 117], [5, 21], [21, 94]]
        @test all(medians .< in_mean_years[:, 1])
        longest = sort(vec(reduce(hcat, times(drawn, j, 1) for j in 1:3)); rev = true)
        @test longest[1:2] == [905, 270]
        @test drawn[187][1, 1] == 905
    end
    BLAS.set_num_threads(threads)
end
