# The annual plant model of Godoy, Kraft & Levine (2014), checked against the paper's printed values and against
# the paper's one-stage equations run in R on the same data (test/reference/R/run_annualplants_reference.R).

using LinearAlgebra

@testset "Annual plants with a seed bank (Godoy et al. 2014)" begin
    p = ANNUAL_PLANTS_GODOY_2014
    ref = joinpath(@__DIR__, "reference")
    csv(file) = [split(l, ",") for l in readlines(joinpath(ref, file))[2:end]]
    num(x) = x == "NA" ? NaN : parse(Float64, x)
    reldiff(x, y) = abs(x - y) / max(abs(y), 1.0)

    @testset "data" begin
        @test length(p.codes) == 18 && p.codes[1] == :aghe && p.codes[end] == :saca
        @test count(isnan, p.alpha) == 66 && !any(isnan, diag(p.alpha))
        # Appendix S2, Table 3: g to three decimals and s to two agree with the figshare data
        g3 = [0.082, 0.082, 0.047, 0.002, 0.036, 0.004, 0.003, 0.011, 0.030, 0.046, 0.015, 0.004, 0.025, 0.085,
              0.028, 0.005, 0.115, 0.077]
        s3 = [0.56, 0.56, 0.71, 0.15, 0.29, 0.24, 0.40, 1.00, 0.19, 0.13, 1.00, 0.26, 0.39, 0.26, 0.55, 0.11, 0.62, 0.34]
        @test round.(p.g; digits = 3) == g3
        @test round.(p.s; digits = 2) == s3
        # lambda does not: Table 3 prints 784.8 for AGHE and 14606.58 for CEME, the data 765.92 and 10962.99
        λ3 = [784.8, 1946.32, 289.91, 1146.23, 14606.58, 835.17, 295.73, 341.23, 1092.95, 60.92, 991.28, 3689.6,
              481.63, 1484.84, 287.72, 2814.96, 1393.34, 3234.13]
        @test round.(p.lambda; digits = 2) != λ3
        @test round(maximum(abs.(p.lambda .- λ3) ./ λ3); digits = 3) == 0.249      # CEME
        # a missing coefficient is r_i e_j (Appendix S4): ANAR on AGHE, alpha[1, 4], is 0.368 * 0.057
        A = annual_plant_alphas()
        @test isnan(p.alpha[1, 4]) && A[1, 4] == 0.368 * 0.057
        @test A[.!isnan.(p.alpha)] == p.alpha[.!isnan.(p.alpha)]
        @test isnan(annual_plant_alphas(; fill = false)[1, 4])
    end

    @testset "the two-stage net" begin
        plant = annual_plant(:a; g = 0.1, s = 0.5, lambda = 100, alpha = [:a => 0.2])
        @test stage_names(plant, :a) == [:seeds, :plants]
        @test pool_names(plant) == [:a_crowding]
        # with N = 50 seeds before germination, 5 germinate: N' = 0.9 * 0.5 * 50 + 0.1 * 100 * 50 / (1 + 0.2 * 5)
        n = advance(plant, [45.0, 5.0], nothing)
        @test sum(n) ≈ 0.9 * 0.5 * 50 + 0.1 * 100 * 50 / 2 rtol = 1e-14
        @test n ≈ [0.9, 0.1] .* sum(n) rtol = 1e-14
        @test annual_plant_state(:aghe, 100) ≈ [100 * (1 - p.g[1]), 100 * p.g[1]] rtol = 1e-15
        c = annual_plant_community(; codes = (:ceme, :heco))
        @test species_names(c) == [:ceme, :heco] && Set(pool_names(c)) == Set([:ceme_crowding, :heco_crowding])
        @test_throws ArgumentError annual_plant_community(; codes = (:aghe, :anar), fill = false)
        @test_throws ArgumentError annual_plant_state(:xxxx, 1)
    end

    @testset "runs against the paper's equations in R" begin
        rows = csv("annualplants_runs.csv")
        # largest relative difference (floor 1) of the seeds N = seeds + plants from R's one-stage recursion, 100
        # years from 100 seeds of each species; tolerance about 100 times the agreement of 2026-11-19
        for (label, tol) in (("ceme_heco", 1e-13), ("aghe_agre", 1e-13), ("aghe_erci", 1e-13), ("all", 1e-13))
            # achieved: 1.1e-15 for all 18 species, the pairs no more
            mine = [r for r in rows if r[1] == label]
            codes = unique(Symbol.(lowercase.(getindex.(mine, 2))))
            c = annual_plant_community(; codes)
            sol = project(c, reduce(vcat, [annual_plant_state(k, 100) for k in codes]), fill(nothing, 100))
            worst = 0.0
            for r in mine
                N = sum(sol[Symbol(lowercase(r[2]))][parse(Int, r[3]) + 1, :])
                worst = max(worst, reldiff(N, num(r[4])))
            end
            @test worst < tol
        end
    end

    @testset "niche and fitness differences of the 153 pairs" begin
        pairs = annual_plant_pairs()
        rows = csv("annualplants_pairs.csv")
        @test length(pairs) == length(rows) == 153
        @test [(x.i, x.j) for x in pairs] == [(Symbol(lowercase(r[1])), Symbol(lowercase(r[2]))) for r in rows]
        worst = maximum(maximum(reldiff(getproperty(x, f), num(r[k]))
                                for (f, k) in ((:rho, 4), (:fitness_ratio, 5), (:invasion_i, 6), (:invasion_j, 7)))
                        for (x, r) in zip(pairs, rows))
        @test worst < 1e-13                    # achieved 0.0
        @test [String(x.outcome) for x in pairs] == [r[8] for r in rows]
        @test [x.winner === nothing ? "" : uppercase(String(x.winner)) for x in pairs] == [r[9] for r in rows]
        @test [x.filled for x in pairs] == [parse(Int, r[3]) for r in rows]
        # the paper: 12 of 153 pairs predicted to coexist; 18 under the relaxed criterion 0.9 rho < kappa_i / kappa_j
        # (j the fitter); interspecific above intraspecific effects (rho > 1) in 15% of pairs. These data give 10,
        # 11 and 37% (vignette 29).
        outcome(o) = count(x -> x.outcome === o, pairs)
        @test outcome(:coexistence) == 10 && outcome(:priority_effect) == 10 && outcome(:exclusion) == 133
        @test count(x -> 0.9x.rho < 1 / max(x.fitness_ratio, 1 / x.fitness_ratio), pairs) == 11
        @test count(x -> x.rho > 1, pairs) == 56
        # the condition of eqn. 5 is the mutual invasion of the model
        @test all((x.outcome === :coexistence) == (x.rho < x.fitness_ratio < 1 / x.rho) for x in pairs)
        # every pair has both coefficients from the pairwise fits only without filling; 97 pairs have all four
        @test count(x -> x.filled == 0, pairs) == 97
    end

    @testset "invasion and equilibria from the composed nets" begin
        A, η = annual_plant_alphas(), annual_plant_eta()
        alone(k) = (i = findfirst(==(k), p.codes); N = (η[i] - 1) / A[i, i] / p.g[i]; annual_plant_state(k, N))
        pairs = annual_plant_pairs()
        rows = csv("annualplants_pairs.csv")
        worst = 0.0
        for x in pairs[1:17:end]      # nine pairs, from every first species
            c = annual_plant_community(; codes = (x.i, x.j))
            ci, cj = subcommunity(c, x.i), subcommunity(c, x.j)
            # each species alone at its equilibrium, (eta - 1) / alpha_ii germinants
            @test advance(cj, alone(x.j), nothing) ≈ alone(x.j) rtol = 1e-13
            worst = max(worst, reldiff(invasion_growth_rate(ci, cj, alone(x.j), nothing), x.invasion_i),
                        reldiff(invasion_growth_rate(cj, ci, alone(x.i), nothing), x.invasion_j))
        end
        @test worst < 1e-13                    # achieved 2.7e-16
        # the coexistence equilibria: the germinants solve 1 + alpha P = eta for both species; stable
        for (x, r) in zip(pairs, rows)
            x.outcome === :coexistence || continue
            c = annual_plant_community(; codes = (x.i, x.j))
            Pi, Pj = num(r[10]), num(r[11])
            start = vcat(annual_plant_state(x.i, Pi / p.g[findfirst(==(x.i), p.codes)]),
                         annual_plant_state(x.j, Pj / p.g[findfirst(==(x.j), p.codes)]))
            e = periodic_point(c, start .* 1.01, nothing)
            @test e.converged && maximum(abs.(e.multipliers)) < 1
            @test maximum(reldiff.(e.n, start)) < 6e-12        # achieved 6.0e-14 (Newton's tolerance)
        end
    end

    @testset "simulated outcomes" begin
        steps(c, n, k) = (for _ in 1:k; n = advance(c, n, nothing); end; n)
        start(codes, N) = reduce(vcat, [annual_plant_state(k, n) for (k, n) in zip(codes, N)])
        total(n) = [n[1] + n[2], n[3] + n[4]]
        # CEME and HECO coexist from either species rare
        c = annual_plant_community(; codes = (:ceme, :heco))
        @test all(minimum(total(steps(c, start((:ceme, :heco), N), 3000))) > 100 for N in ((1.0, 1000.0), (1000.0, 1.0)))
        # AGRE excludes AGHE
        c = annual_plant_community(; codes = (:aghe, :agre))
        @test total(steps(c, start((:aghe, :agre), (1000.0, 1.0)), 500))[1] < 1e-6
        # AGHE and ERCI: whichever starts commoner wins (a priority effect)
        c = annual_plant_community(; codes = (:aghe, :erci))
        n1 = total(steps(c, start((:aghe, :erci), (1000.0, 10.0)), 3000))
        n2 = total(steps(c, start((:aghe, :erci), (10.0, 1000.0)), 3000))
        @test n1[2] < 1e-3 * n1[1] && n2[1] < 1e-3 * n2[2]
    end

    @testset "competitive ranks (Fig. 2)" begin
        # the R reference ranks species by the paper's sequential simulations (20000 years a round); in its first
        # round HECO and PLER survive (the paper: PLER alone at rank 1, HECO tied third). Checked here without the
        # 20000 years: the HECO-PLER equilibrium resists invasion by each of the other 16 species.
        ranks = Dict(Symbol(lowercase(r[1])) => parse(Int, r[2]) for r in csv("annualplants_ranks.csv"))
        @test sort([k for (k, v) in ranks if v == 1]) == [:heco, :pler]
        @test ranks[:agre] == 2 && sort([k for (k, v) in ranks if v == 3]) == [:ceme, :saca]
        c = annual_plant_community()
        r = only(r for r in csv("annualplants_pairs.csv") if r[1] == "HECO" && r[2] == "PLER")
        g(k) = p.g[findfirst(==(k), p.codes)]
        resident = subcommunity(c, :heco, :pler)
        n = vcat(annual_plant_state(:heco, num(r[10]) / g(:heco)), annual_plant_state(:pler, num(r[11]) / g(:pler)))
        @test advance(resident, n, nothing) ≈ n rtol = 1e-12
        @test all(invasion_growth_rate(subcommunity(c, k), resident, n, nothing) < 1 for k in p.codes if !(k in (:heco, :pler)))
    end
end
