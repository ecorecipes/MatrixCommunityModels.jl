# The designed flow regimes of Tonkin et al. (2021): three invertebrates, riparian plants and fish under four flow
# regimes, checked against the authors' scripts (test/reference/R/run_flow_tradeoffs_reference.R).
# Tolerances are about 100 times the agreement measured on 2026-11-17 (macOS, Julia 1.12), never below 1e-13
# relative for runs of 200 years or 30 years of days.

using LinearAlgebra, Random

@testset "Designed flow regimes (Tonkin et al. 2021)" begin
    average(x) = sum(x) / length(x)
    csv(f) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", f))[2:end]]
    code = Dict(:plants => "RIP_OPT", :fish => "FSH_OPT", :invertebrates => "INV_OPT", :natural => "NAT")
    taxon = Dict("Mayfly" => :mayfly, "Dragonfly" => :dragonfly, "Ostracod" => :ostracod)
    inverts = csv("flow_inverts.csv")
    riparian = csv("flow_riparian.csv")
    fish_initial = csv("flow_fish_initial.csv")
    fish_mean = csv("flow_fish_mean.csv")
    targets = Dict(r[1] => parse.(Float64, r[2:4]) for r in csv("flow_targets.csv"))
    rip = riparian_model()
    bugs = flow_invertebrates()
    fish = verde_fish_model()
    @test isempty(pool_names(bugs))
    nonseed(sol, g) = [100 * 350 * sum(sol[g][k, s] / RIPARIAN_GUILDS[g].den[s] for s in 2:6) / rip.K for k in 1:200]
    # largest absolute difference of the invertebrates from R, which integrates numerically (the error is R's):
    # achieved 3.9e-12, 1.1e-12, 1.8e-13, 6.0e-13 (2026-11-17)
    invert_tol = Dict(:plants => 4e-10, :fish => 1e-10, :invertebrates => 2e-11, :natural => 6e-11)
    T = Dict{Symbol, NamedTuple}()
    for s in TONKIN_SCENARIOS
        sc = tonkin_scenario(s)
        @test length(sc.riparian) == 200 && length(sc.fish) == 200 && length(sc.invertebrates) == 30 * 365
        # invertebrates, daily, against R's closed form with numerical integration
        sol = project(bugs, sc.invertebrate_n0, sc.invertebrates)
        for (name, x) in taxon
            rows = [r for r in inverts if r[1] == name && r[2] == code[s]]
            days = [parse(Int, r[3]) for r in rows]
            @test maximum(abs.(sol[x][days .+ 1, 1] .- [parse(Float64, r[5]) for r in rows])) < invert_tol[s]
            @test [sc.invertebrates[d].Kpre[x] for d in days] == [parse(Float64, r[4]) for r in rows]
        end
        invertebrate_target = average((sol[:mayfly][3652:end, 1] .+ sol[:dragonfly][3652:end, 1]) ./ 2)
        @test invertebrate_target ≈ targets[code[s]][3] rtol = 1e-13    # achieved <= 4.4e-16 (2026-11-17)
        # riparian plants, yearly
        r = project(rip.community, sc.riparian_n0, sc.riparian[1:199])
        rows = [x for x in riparian if x[1] == code[s]]
        for (j, g) in enumerate((:C, :T, :W, :S, :M))
            R = [parse(Float64, x[2 + j]) for x in rows]
            # achieved <= 5.5e-16 over the five guilds and four regimes (2026-11-17)
            @test maximum(abs.(nonseed(r, g) .- R) ./ max.(abs.(R), 1e-3)) < 1e-13
        end
        plant_target = average(nonseed(r, :C)[11:200])
        @test plant_target ≈ targets[code[s]][1] rtol = 1e-13           # achieved <= 2.2e-16 (2026-11-17)
        # fish, replaying the authors' 100 random starts of this regime
        acc = zeros(200, 7)
        for r in fish_initial
            r[1] == code[s] || continue
            f = project(fish.community, [Symbol(sp) => [0.0, 0.0, parse(Float64, r[2 + j])] for (j, sp) in enumerate(VERDE_SPECIES)],
                        sc.fish[1:199])
            for (j, sp) in enumerate(VERDE_SPECIES)
                acc[:, j] .+= (f[Symbol(sp)][:, 2] .+ f[Symbol(sp)][:, 3]) ./ 100
            end
        end
        M = permutedims(reduce(hcat, [parse.(Float64, r[3:9]) for r in fish_mean if r[1] == code[s]]))
        # achieved 3.0e-15, 1.4e-15, 1.5e-15, 1.3e-15 (2026-11-17; the vignette, dividing by 100 once, gets 0.0)
        @test maximum(abs.(acc .- M) ./ max.(M, 1e-6)) < 3e-13
        fish_target = average(100 .* (acc[11:200, 1] .+ acc[11:200, 2] .+ acc[11:200, 4]) ./ 47660)
        @test fish_target ≈ targets[code[s]][2] rtol = 1e-13            # achieved <= 4.4e-16 (2026-11-17)
        T[s] = (plants = plant_target, fish = fish_target, invertebrates = invertebrate_target)
    end
    # the targets as vignette 22 prints them (equal to R to two decimals)
    @test [round.(values(T[s]); digits = 2) for s in (:plants, :fish, :invertebrates, :natural)] ==
          [(49.33, 2.65, 40.87), (2.34, 70.94, 46.58), (3.19, 3.51, 72.87), (32.76, 52.77, 52.21)]
    # the ratios against the paper, in per cent: natural / designed for each group, then the cost of each
    # prescription to the others. All within a point of the paper but fish (74.4%, published 77%: the
    # authors' runs drew other random starting biomasses)
    nat = T[:natural]
    ratios = 100 .* [nat.plants / T[:plants].plants, nat.fish / T[:fish].fish, nat.invertebrates / T[:invertebrates].invertebrates,
                     T[:fish].plants / nat.plants, T[:plants].fish / nat.fish, T[:invertebrates].plants / nat.plants,
                     T[:invertebrates].fish / nat.fish]
    @test round.(ratios; digits = 1) == [66.4, 74.4, 71.6, 7.2, 5.0, 9.7, 6.6]
    @test all(abs.(ratios .- [66, 77, 72, 7, 5, 10, 6])[[1, 3, 4, 5, 6, 7]] .<= 1)
    # the frequency search of the authors' code (dod-3taxa-inverts-findopt.R): ten years of floods every 91, 121,
    # 182, 365 or 730 days, mean and geometric mean after the first year. Achieved: means 3.2e-14, geometric
    # means 1.1e-13 relative (both largest for dragonflies flooded every 91 days, at 1.7e-7 and 4.7e-28; 2026-11-17)
    gm(x) = exp(sum(log.(x[x .> 0])) / length(x))
    found = csv("flow_findopt.csv")
    for (name, x, Q) in (("Dragonfly", :dragonfly, 1000), ("Mayfly", :mayfly, 1000), ("Ostracod", :ostracod, 30))
        for (k, every) in enumerate((365 / 4, 365 / 3, 365 / 2, 365, 730))
            n = project(flood_invertebrate(x), [50.0], invertebrate_flows(flood_pulses(every, Q, 3650)...))[x][367:end, 1]
            r = only(r for r in found if r[1] == name && r[2] == "freq$k")
            @test average(n) ≈ parse(Float64, r[3]) rtol = 3e-12
            @test gm(n) ≈ parse(Float64, r[4]) rtol = 1e-11
        end
    end
    # the random starts of the authors' code can be drawn here too: mean CAIN adults 34331 over 2000 draws
    # (standard error 660, so 5% is 2.6 standard errors)
    starts = [tonkin_fish_start(Random.Xoshiro(k)) for k in 1:2000]
    @test average(s[4].second[3] for s in starts) ≈ TONKIN_FISH_INITIAL.CAIN[2] rtol = 0.05
    # the corrected carrying capacity relaxes towards Kb between floods
    c = invertebrate_flows(flood_pulses(91.25, 50, 3650)...; corrected = true)
    @test c[end].Kpre.dragonfly < invertebrate_flows(flood_pulses(91.25, 50, 3650)...)[end].Kpre.dragonfly
    # and lowers the invertebrate target of the invertebrate prescription from 72.9 to 52.8 (as coded / corrected:
    # plants 40.9 / 40.9, fish 46.6 / 45.6, natural 52.2 / 47.3), as vignette 22 prints them
    corrected_target(s) = (sc = tonkin_scenario(s; corrected = true); sol = project(bugs, sc.invertebrate_n0, sc.invertebrates);
                           average((sol[:mayfly][3652:end, 1] .+ sol[:dragonfly][3652:end, 1]) ./ 2))
    @test [round(corrected_target(s); digits = 1) for s in (:plants, :fish, :invertebrates, :natural)] == [40.9, 45.6, 52.8, 47.3]
end
