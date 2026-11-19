# The Bd IPM of Wilber et al. (2022) over the whole 529 weeks of the authors' runs at Pennsylvania, Louisiana and
# Tennessee (the main suite runs these three for 120 weeks and Vermont for 529), against the authors' Python code
# (test/reference/python/run_bd_reference.py), and the seeded stochastic runs of vignette 24.
using Random

@testset "chytrid fungus in leopard frogs, long runs (vignette 24)" begin
    number(x) = parse(Float64, x)
    y = meshpoints(chytrid_load_stage())
    runs, states = long_reference("chytrid_runs.csv"), long_reference("chytrid_states.csv")
    start(c) = state_vector(c, :frog => vcat([0.0, 70.0], zeros(300)), :zoospores => [10_000.0])
    deterministic = Dict{Symbol, Any}()

    @testset "529 weeks against the authors' code" begin
        # Achieved (2026-11-18), largest over the 530 saved steps (totals, mean log load) or the seven saved states:
        # - totals (tadpoles, susceptible, infected, zoospores), relative: PA 2.2e-12, LA 7.9e-12, TN 6.6e-12
        # - whole state relative to its largest entry: PA 1.0e-13, LA 2.0e-13, TN 3.6e-13 (at 120 weeks LA 3.5e-14, TN 1.0e-14)
        # - mean log load of the infected frogs, relative: PA 2.6e-13, LA 3.7e-14, TN 4.8e-13
        # - water temperature, absolute: 3.6e-15
        # Tolerances are about 100 times these, the totals' kept at the main suite's 1e-10.
        for (site, state_tol, load_tol) in ((:PA, 1e-11, 3e-11), (:LA, 2e-11, 4e-12), (:TN, 4e-11, 5e-11))
            c = chytrid_community(site; area = 500.0)
            env = chytrid_environment(site, 529)
            sol = project(c, start(c), env)
            deterministic[site] = (community = c, sol = sol)
            frog, z = sol[:frog], sol[:zoospores]
            worst, load, temperature = 0.0, 0.0, 0.0
            for r in runs
                r[1] == String(site) || continue
                k = parse(Int, r[2]) + 1
                mine = (frog[k, 1], frog[k, 2], sum(frog[k, 3:end]), z[k, 1])
                worst = max(worst, maximum(abs(m - t) / max(abs(t), 1e-300) for (m, t) in zip(mine, number.(r[4:7]))))
                infected = frog[k, 3:end]
                load = max(load, sum(infected) > 0 ? abs(sum(infected .* y) / sum(infected) / number(r[8]) - 1) :
                                 r[8] == "nan" ? 0.0 : Inf)
                k <= 529 && (temperature = max(temperature, abs(env[k].temperature - number(r[3]))))
            end
            whole = 0.0
            for r in states
                r[1] == String(site) || continue
                k = parse(Int, r[2]) + 1
                v = number.(r[3:end])
                whole = max(whole, maximum(abs.(vcat(frog[k, :], z[k, 1]) .- v)) / maximum(abs, v))
            end
            println("chytrid, $site, 529 weeks: totals ", worst, ", whole state ", whole, ", mean load ", load,
                    ", temperature ", temperature)
            @test count(r -> r[1] == String(site), runs) == 530
            @test count(r -> r[1] == String(site), states) == 7
            @test worst < 1e-10
            @test whole < state_tol
            @test load < load_tol
            @test temperature < 4e-13
        end
    end

    @testset "Bd kept for ten years in stochastic runs" begin
        # Vignette 24: 30 runs at each location from the end of the deterministic run, whole frogs and zoospores,
        # 522 weeks, drawn with Xoshiro(2022) in the order LA, TN, PA. It prints 77%, 27% and 0% (the authors: 76%,
        # 22% and 0% of 500 runs). Seeded, so the counts are exact on any platform that draws the same numbers.
        function persists(site, rng)
            c = deterministic[site].community
            settled = deterministic[site].sol.u[end, :]
            frog = settled[c.layout[:frog]]
            infected = frog[3:end]
            total = round(Int, sum(infected))
            drawn = zeros(length(infected))
            for _ in 1:total                   # a multinomial draw of the infected frogs over the loads
                u, j = rand(rng) * sum(infected), 1
                while u > infected[j] && j < length(infected)
                    u -= infected[j]; j += 1
                end
                drawn[j] += 1
            end
            n0 = state_vector(c, :frog => vcat(ceil.(frog[1:2]), drawn), :zoospores => [ceil(settled[c.layout[:zoospores]][1])])
            sol = project(c, n0, chytrid_environment(site, 522; start = 7 * 529); rng)
            return sum(sol[:frog][end, 3:end]) > 0 || sol[:zoospores][end, 1] > 0
        end
        rng = Xoshiro(2022)
        kept = [count(persists(site, rng) for _ in 1:30) for site in (:LA, :TN, :PA)]
        println("chytrid, Bd kept for ten years in ", kept, " of 30 runs (LA, TN, PA)")
        @test kept == [23, 8, 0]
        # within about two standard errors of the authors' shares (0.76, 0.224, 0)
        @test all(abs.(kept ./ 30 .- [0.76, 0.224, 0.0]) .<= 0.15)
    end
end
