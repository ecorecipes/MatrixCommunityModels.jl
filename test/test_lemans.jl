# The length-based multispecies fish model LeMans (Hall et al. 2006) for 21 species of the North Sea, compared with
# the R package LeMaRns 0.1.2 run on its own data (test/reference/R/run_lemans_reference.R).
using LinearAlgebra

lemans_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]

@testset "LeMans, North Sea (LeMaRns)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    number(x) = parse(Float64, x)
    relative(a, b) = b == 0 ? abs(a) : abs(a - b) / abs(b)
    p = lemans_parameters()
    named = Dict(n => s for (n, s) in zip(p.names, p.species))
    classes = lemans_reference("lemans_classes.csv")

    @testset "data and derived quantities" begin
        @test length(p.species) == 21 && size(LEMANS_NORTH_SEA.diet) == (21, 21)
        @test p.species[[1, 4, 20]] == [:sprat, :poor_cod, :cod]
        @test [p.fish[s].classes for s in p.species] == [4, 6, 5, 5, 6, 7, 8, 6, 8, 10, 9, 12, 10, 10, 14, 14, 15, 20, 23, 32, 26]
        @test length(p.predators) == 250 && p.gears == [:industrial, :otter, :beam, :pelagic]
        start = Dict(lemans_state(p))
        worst = Dict{Symbol, Float64}()
        for r in classes
            f, j = p.fish[named[r[1]]], parse(Int, r[2])
            for (name, mine, column) in ((:lower, p.lower[j], 3), (:mid, p.mid[j], 5), (:weight, f.weight[j], 6),
                                         (:leaving, f.leaving[j], 7), (:ration, f.ration[j], 8), (:mature, f.mature[j], 9),
                                         (:natural, f.natural[j], 10), (:catchability, f.catchability[j], 11),
                                         (:start, start[named[r[1]]][j], 12))
                worst[name] = max(get(worst, name, 0.0), relative(mine, number(r[column])))
            end
        end
        @test length(classes) == 250
        @test all(<(5e-13), values(worst))      # achieved 4.8e-15 (rations; 2026-11-17)
        # suitabilities of 1680 predator-prey pairs (whiting, cod, saithe as predators)
        far = maximum(lemans_reference("lemans_suitability.csv")) do r
            row = findfirst(==((named[r[1]], parse(Int, r[2]))), p.predators)
            relative(p.fish[named[r[3]]].suitability[row, parse(Int, r[4])], number(r[5]))
        end
        @test far < 1e-13                       # achieved 1.0e-15 (2026-11-17)
        # each predator's suitabilities sum to one over its prey, except in the smallest class, which eats nothing
        total = reduce(+, [vec(sum(p.fish[s].suitability; dims = 2)) for s in p.species])
        @test all(x -> isapprox(x, 1; atol = 5e-14) || x == 0, total)    # achieved 4.4e-16
        @test all(total[r] == 0 for (r, (sp, j)) in enumerate(p.predators) if j == 1)
        @test_throws ArgumentError lemans_fish(p, :tuna)
    end

    c = lemans_community(p)
    n0 = state_vector(c, lemans_state(p)...)

    @testset "pools and predation" begin
        @test pool_names(c) == [:ration, :prey]
        pools = pool_totals(c, n0, (recruiting = false, effort = (;)))
        @test pools.ration isa Vector{Float64} && length(pools.prey) == 250
        M2 = lemans_predation(p, c, n0, (recruiting = false, effort = (;)))
        far = maximum(relative(M2[named[r[1]]][parse(Int, r[2])], number(r[13])) for r in classes)
        @test far < 1e-13                       # achieved 8.2e-16 (2026-11-17)
        # the recruits of the first step of a year are prey in that step
        with = pool_totals(c, n0, (recruiting = true, effort = (;)))
        @test all(with.prey .>= pools.prey) && any(with.prey .> pools.prey)
        @test with.ration == pools.ration
    end

    @testset "calendar" begin
        steps = findall(e -> e.recruiting, lemans_environment(100))
        @test steps == [parse(Int, r[1]) for r in lemans_reference("lemans_recruitment_steps.csv")]
        envs = lemans_environment(2; effort = [(otter = 1.0,), (otter = 2.0,)])
        @test length(envs) == 20 && envs[10].effort.otter == 1.0 && envs[11].effort.otter == 2.0
        @test envs[11].year == 2 && envs[11].recruiting && !envs[12].recruiting
    end

    runs = lemans_reference("lemans_runs.csv")
    order = species_names(c)
    # the largest relative difference from LeMaRns in biomass and in spawning stock biomass, over the saved steps up to `last`
    function apart(run, sol; last = typemax(Int))
        B, S = lemans_biomass(p, sol), lemans_biomass(p, sol; spawning = true)
        worst = 0.0
        for r in runs
            r[1] == run || continue
            step = parse(Int, r[2])
            step <= min(last, size(B, 1) - 1) || continue
            i = findfirst(==(named[r[3]]), order)
            worst = max(worst, relative(B[step + 1, i], number(r[4])), relative(S[step + 1, i], number(r[5])))
        end
        return worst
    end
    states = lemans_reference("lemans_state.csv")
    saved(run) = state_vector(c, (sp => [number(r[4]) for r in states if r[1] == run && named[r[2]] == sp] for sp in order)...)

    # The runs of vignette 18, 100 unfished years and 50 fished ones (a step costs about a millisecond), with tolerances
    # of about a hundred times the differences found on 2026-11-17.
    @testset "runs against LeMaRns" begin
        century = lemans_environment(100)
        unfished = project(c, n0, century)
        @test apart("unfished", unfished) < 3e-13                 # achieved 2.7e-15
        @test size(unfished.pool_vectors[:prey]) == (1001, 250) && all(isnan, unfished.pools)
        @test all(isnan, unfished.pool_vectors[:ration][end, :]) && !any(isnan, unfished.pool_vectors[:ration][1:1000, :])
        # numbers in every class at the end, and the predation mortality of the last step, which LeMaRns computes from
        # the state before it
        state_apart(run, n) = maximum(relative(a, b) for (a, b) in zip(n, saved(run)))
        M2_apart(run, n, env) = (M2 = lemans_predation(p, c, n, env);
                                 maximum(relative(M2[named[r[2]]][parse(Int, r[3])], number(r[5])) for r in states if r[1] == run))
        @test state_apart("unfished", unfished.u[end, :]) < 3e-13                                          # achieved 2.9e-15
        @test M2_apart("unfished", unfished.u[end - 1, :], century[end]) < 1.5e-13     # achieved 1.3e-15
        # by year 100 the community repeats itself from year to year
        S = lemans_biomass(p, unfished; spawning = true)
        @test maximum(abs.(S[end, :] ./ S[end - 10, :] .- 1)) < 1e-13                                     # achieved 2.2e-16
        start = saved("unfished")
        for (run, effort) in (("fished", (industrial = 1.0, otter = 1.0, beam = 1.0, pelagic = 1.0)), ("otter", (otter = 1.5,)))
            envs = lemans_environment(50; effort)
            sol = project(c, start, envs)
            @test apart(run, sol) < 1.5e-12                       # achieved 1.3e-14 (fished), 7.4e-15 (otter)
            @test state_apart(run, sol.u[end, :]) < 1.5e-12       # achieved 1.3e-14, 7.8e-15
            @test M2_apart(run, sol.u[end - 1, :], envs[end]) < 1.5e-13   # achieved 1.4e-15, 1.3e-15
        end
        # selectivity at the midpoints of the classes, rather than LeMaRns's upper bounds, changes the catchabilities
        mid = lemans_parameters(; selectivity_at = :mid)
        @test mid.fish.cod.catchability != p.fish.cod.catchability && mid.fish.cod.weight == p.fish.cod.weight
        @test_throws ArgumentError lemans_parameters(; selectivity_at = :lower)
        # the rare-invader analysis reads vector pools from the trajectory
        rest = lemans_community(p; species = [s for s in p.species if s != :sprat])
        rate = stochastic_invasion_rate(lemans_fish(p, :sprat), rest,
            state_vector(rest, lemans_state(p; species = species_names(rest))...), lemans_environment(2))
        @test isfinite(rate) && rate > 0
        # both solvers give the same steps
        short = lemans_environment(1)
        @test project(c, n0, short; solver = :ipm).u ≈ project(c, n0, short; solver = :mpm).u rtol = 1e-12
    end
    BLAS.set_num_threads(threads)
end
