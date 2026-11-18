# The full Bd IPM of Wilber et al. (2022) against the authors' Python code (test/reference/python/run_bd_reference.py).
using Test, MatrixCommunityModels

@testset "chytrid fungus in leopard frogs" begin
    ref = joinpath(@__DIR__, "reference")
    rows(f) = [split(l, ",") for l in readlines(joinpath(ref, f))[2:end]]
    num(x) = parse(Float64, x)

    @testset "parameters as in the authors' runs" begin
        keys_ = (:s0, :sI, :a, :a_temp, :b, :sigma, :mu_l, :sigma_l, :beta, :beta_temp, :shedding, :fecundity, :K,
                 :larval_period)
        for r in rows("chytrid_parameters.csv")
            s = CHYTRID_SITES[Symbol(r[1])]
            @test Symbol(r[3]) == s.model
            @test [getfield(s, k) for k in keys_] == num.(r[5:18])
            @test Float64(s.aquatic) == num(r[19])
            @test Float64.(s.breeding) == (num(r[20]), num(r[21]))
            @test Float64(s.hibernation) == num(r[22]) && Float64(s.hibernation_aquatic) == num(r[23])
            @test Float64.(s.hibernation_period) == (num(r[24]), num(r[25]))
        end
        z = only(rows("chytrid_zoospore_survival.csv"))
        @test (CHYTRID_ZOOSPORE_SURVIVAL.ld50, CHYTRID_ZOOSPORE_SURVIVAL.s, CHYTRID_ZOOSPORE_SURVIVAL.maxsurv) == Tuple(num.(z))
        # the curve is a least-squares fit to Woodhams et al.'s weekly survival: no nearby curve fits better
        T = [4.0, 14.5, 23.0, 28.0]
        obs = exp.(-[0.0027, 0.0035, 0.0160, 0.0410] .* 24 .* 7)
        ss(l, s) = sum((obs .- maximum(obs) ./ (1 .+ exp.((T .- l) ./ s))) .^ 2)
        best = ss(CHYTRID_ZOOSPORE_SURVIVAL.ld50, CHYTRID_ZOOSPORE_SURVIVAL.s)
        @test all(ss(CHYTRID_ZOOSPORE_SURVIVAL.ld50 + d1, CHYTRID_ZOOSPORE_SURVIVAL.s + d2) >= best - 1e-12
                  for d1 in (-0.01, 0.0, 0.01), d2 in (-0.01, 0.0, 0.01))
    end

    @testset "seasons" begin
        @test chytrid_temperature(:PA, 0) == 4.0
        @test chytrid_temperature(:LA, 182) ≈ 30.0
        envs = chytrid_environment(:TN, 3; start = 7)
        @test [e.day for e in envs] == [7, 14, 21]
        # the breeding pulse falls in one step a year: the middle of the breeding season
        pulse(site) = [e.day for e in chytrid_environment(site, 52) if MatrixCommunityModels._chytrid_breeding_now(CHYTRID_SITES[site], e.day)]
        @test pulse(:PA) == [119] && pulse(:LA) == [77] && pulse(:TN) == [91] && pulse(:VT) == [133]
        water(site, day) = MatrixCommunityModels._chytrid_in_water(CHYTRID_SITES[site], day)
        @test water(:PA, 100) && !water(:PA, 200) && water(:PA, 330) && water(:PA, 10)   # breeding, summer, hibernation
        @test water(:LA, 100) && !water(:LA, 200) && !water(:LA, 330)                    # LA frogs do not hibernate in water
    end

    @testset "the net" begin
        c = chytrid_community(:PA; cells = 300)
        @test species_names(c) == [:frog, :zoospores]
        @test Set(pool_names(c)) == Set([:frogs, :shedding, :zoospores])
        frog = c.system[:frog]
        @test stage_names(frog)[1:3] == [:tadpole, :susceptible, :load1] && length(stage_names(frog)) == 302
        @test Set(transition_names(frog)) ==
              Set([:developing, :metamorphosing, :escaping, :infecting, :growing, :clearing, :breeding])
        # two hosts with their own crowding pools share the zoospores
        two = community(leopard_frog(:PA; name = :pipiens, crowding = :pipiens_frogs, cells = 30),
                        leopard_frog(:TN; name = :sphenocephala, crowding = :sphenocephala_frogs, cells = 30),
                        chytrid_zoospores())
        @test Set(pool_names(two)) == Set([:pipiens_frogs, :sphenocephala_frogs, :shedding, :zoospores])
        @test_throws ArgumentError leopard_frog(:XX)
        @test_throws ArgumentError leopard_frog(:PA; area = 0)
    end

    @testset "runs against the authors' code" begin
        runs, states = rows("chytrid_runs.csv"), rows("chytrid_states.csv")
        y = meshpoints(chytrid_load_stage())
        # Three sites for 120 weekly steps (saved steps 0, 1, 26, 52, 104), and Vermont for the whole 529 (all
        # seven saved steps), where vignette 24 finds the largest difference; 529 steps at every site would add
        # about 10 s. Tolerances, measured 2026-11-17:
        # - totals (tadpoles, susceptible, infected, zoospores), relative to each value: achieved PA 2.2e-12,
        #   LA 7.9e-12, TN 6.2e-12, VT 1.7e-12. 100 times that would loosen the earlier 1e-10, which is kept.
        # - whole state at the saved steps, relative to its largest entry: achieved PA 7.1e-14, LA 3.5e-14,
        #   TN 1.0e-14, VT 1.6e-12 (at step 529; 2.5e-13 up to step 260).
        # - mean log load of the infected frogs, relative: achieved PA 2.6e-13, LA 3.7e-14, TN 4.8e-13, VT 3.1e-13.
        # - water temperature of the environment for the next step, absolute: achieved 3.6e-15 (4 to 30 °C).
        for (site, steps, state_tol) in ((:PA, 120, 7e-12), (:LA, 120, 4e-12), (:TN, 120, 1e-12), (:VT, 529, 2e-10))
            c = chytrid_community(site; area = 500.0)
            n0 = state_vector(c, :frog => vcat([0.0, 70.0], zeros(300)), :zoospores => [10000.0])
            env = chytrid_environment(site, steps)
            sol = project(c, n0, env)
            frog, z = sol[:frog], sol[:zoospores]
            worst, load, temperature = 0.0, 0.0, 0.0
            for r in runs
                (r[1] == String(site) && parse(Int, r[2]) <= steps) || continue
                k = parse(Int, r[2]) + 1
                mine = (frog[k, 1], frog[k, 2], sum(frog[k, 3:end]), z[k, 1])
                worst = max(worst, maximum(abs(m - t) / max(abs(t), 1e-300) for (m, t) in zip(mine, num.(r[4:7]))))
                infected = frog[k, 3:end]
                # Python writes nan for the mean load when no frog is infected (step 0)
                load = max(load, sum(infected) > 0 ? abs(sum(infected .* y) / sum(infected) / num(r[8]) - 1) :
                                 r[8] == "nan" ? 0.0 : Inf)
                k <= steps && (temperature = max(temperature, abs(env[k].temperature - num(r[3]))))
            end
            @test all(env[k + 1].day == 7k for k in 0:(steps - 1))
            @test worst < 1e-10
            @test load < 5e-11
            @test temperature < 4e-13
            for r in states
                (r[1] == String(site) && parse(Int, r[2]) <= steps) || continue
                k = parse(Int, r[2]) + 1
                v = num.(r[3:end])
                @test maximum(abs.(vcat(frog[k, :], z[k, 1]) .- v)) <= state_tol * maximum(abs, v)
            end
        end
    end
end
