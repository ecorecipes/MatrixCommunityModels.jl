# Long runs of the grassland IPMs of Chu & Adler (2015) as run by Tredennick et al. (2017), vignette 20: the five
# sites over 2500 drawn years, and the Kansas runs of 500 to 2500 years that the main suite leaves out.
using LinearAlgebra, Random

@testset "grasslands, long runs (vignette 20)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    kinds = (:bocu, :bohi, :scsc)
    number(x) = parse(Float64, x)
    cover(c, n, ks = kinds) = [getproperty(pool_totals(c, n, nothing), Symbol(k, :_cover)) for k in ks]
    covers(c, sol, ks = kinds) = reduce(vcat, permutedims(cover(c, row, ks)) for row in eachrow(sol.u))
    growth_rates(X) = log.(X[2:end, :] ./ X[1:(end - 1), :])
    start(c, site = :kansas, ks = kinds) = state_vector(c, (k => grassland_state(site, k) for k in ks)...)
    years = [parse(Int, r[2]) for r in long_reference("grassland_kansas_years.csv")]

    @testset "five sites, synchrony against Tredennick et al. (2017, Table S1)" begin
        # 2500 years drawn with Xoshiro(2017), without interspecific interactions, the first 500 left out, as in
        # vignette 20. Achieved (2026-11-18): 0.84156, 0.81952, 0.62292, 0.47840, 0.41058, equal to the three decimals
        # the vignette prints. Table S1 comes from the authors' own sequences of years: the largest gap is 0.022
        # (Montana), and the vignette claims only Kansas (0.62) and the order of the sites.
        published = (newmexico = 0.86, arizona = 0.81, kansas = 0.62, montana = 0.50, idaho = 0.42)
        printed = (newmexico = 0.842, arizona = 0.820, kansas = 0.623, montana = 0.478, idaho = 0.411)
        here = Dict{Symbol, Float64}()
        for site in keys(published)
            ks = Tuple(keys(GRASSLAND_PARAMETERS[site]))
            T = length(first(GRASSLAND_PARAMETERS[site]).recruitment.intercept_year)
            c = grassland_community(site; interactions = false)
            sol = project(c, start(c, site, ks), rand(Xoshiro(2017), 1:T, 2499))
            here[site] = synchrony(growth_rates(covers(c, sol, ks)[501:end, :]))
            println("grasslands, $site: synchrony of growth rates over 2000 drawn years = ", here[site],
                    " (Table S1: ", published[site], ")")
            @test here[site] ≈ printed[site] atol = 5e-4
            @test here[site] ≈ published[site] atol = 0.03
        end
        @test here[:kansas] ≈ published.kansas atol = 5e-3
        # the order of the sites, from the deserts to the sagebrush steppe
        @test sortperm([here[s] for s in keys(published)]; rev = true) == sortperm(collect(published); rev = true) == 1:5
    end

    @testset "Kansas, 500 mean years against the authors' code" begin
        # vignette 20: 1.5e-8 with interactions and 9.9e-8 without, which the R script's integration sets
        runs = long_reference("grassland_kansas_runs.csv")
        for (label, inter, tol, printed) in (("mean_500", true, 3e-8, [1.7180, 0.4682, 5.3648]),
                                             ("mean_500_nointer", false, 2e-7, nothing))
            c = grassland_community(:kansas; interactions = inter)
            sol = project(c, start(c), fill(nothing, 500))
            rows = [r for r in runs if r[1] == label]
            worst = maximum(maximum(abs.(cover(c, sol.u[parse(Int, r[2]) + 1, :]) ./ number.(r[3:5]) .- 1)) for r in rows)
            println("grasslands, Kansas, $label: largest relative difference in cover = ", worst)
            @test worst < tol               # achieved 1.5e-8 and 9.9e-8 (2026-11-18)
            # cover (%) at year 500, as vignette 20 prints it
            printed === nothing || @test all(abs.(100 .* cover(c, sol.u[end, :]) .- printed) .<= 5e-5)
        end
    end

    @testset "Kansas, the authors' 2500 drawn years" begin
        saved = long_reference("grassland_kansas_saved.csv")
        # synchrony over 100 windows of 50 years drawn with Xoshiro(1), as in the authors' script for the figures
        function window_synchrony(C; rng = Xoshiro(1))
            S = C[501:end, :]
            return sum(begin
                b = rand(rng, 1:(size(S, 1) - 50))
                synchrony(growth_rates(S[b:(b + 50), :]))
            end for _ in 1:100) / 100
        end
        # (run, interactions, largest difference from the authors' saved covers, vignette 20's synchrony over 2000
        # years and over windows, and its mean covers in %). The main suite has the run without interactions.
        for (run, inter, tol, sync, windows, mean_cover) in
            (("ENVINTER", true, 7e-7, 0.611, 0.614, [2.183, 0.034, 8.924]),
             ("ENVNOINTER", false, 6e-7, 0.618, 0.619, [1.557, 2.477, 10.010]))
            c = grassland_community(:kansas; interactions = inter)
            C = covers(c, project(c, start(c), years[2:2500]))
            theirs = [(parse(Int, r[2]), number.(r[3:5])) for r in saved if r[1] == run]
            worst = maximum(maximum(abs.(C[year, :] ./ x .- 1)) for (year, x) in theirs)
            println("grasslands, Kansas, $run: largest relative difference from the authors' saved covers = ", worst)
            @test worst < tol                # achieved 3.5e-7 and 2.9e-7 (2026-11-18), the R script's integration
            @test synchrony(growth_rates(C[501:end, :])) ≈ sync atol = 5e-4
            @test window_synchrony(C) ≈ windows atol = 5e-4
            @test all(abs.(100 .* vec(sum(C[501:end, :]; dims = 1)) ./ 2000 .- mean_cover) .<= 5e-4)
        end
    end

    @testset "Kansas, growth rate when rare in the authors' drawn years" begin
        # residents at their mean-year equilibrium, then the authors' years 2 to 1001, the first 100 left out;
        # vignette 20 prints 1.5974, 0.9851 and 1.6291 (mean year 1.6060, 1.0923, 1.4386)
        printed = [1.5974, 0.9851, 1.6291]
        for (j, k) in enumerate(kinds)
            others = Tuple(x for x in kinds if x != k)
            resident = grassland_community(:kansas; kinds = others)
            n0 = state_vector(resident, (x => grassland_state(:kansas, x) for x in others)...)
            eq = equilibrium(resident, project(resident, n0, fill(nothing, 500)).u[end, :], nothing; tol = 1e-9, maxiter = 5000)
            @test eq.converged
            rate = stochastic_invasion_rate(grassland_plant(:kansas, k), resident, eq.n, years[2:1001]; burnin = 100)
            println("grasslands, Kansas, $k: growth rate when rare in drawn years = ", rate)
            @test rate ≈ printed[j] atol = 5e-5
        end
    end
    BLAS.set_num_threads(threads)
end
