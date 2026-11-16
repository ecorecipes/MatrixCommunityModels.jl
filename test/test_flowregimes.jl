# The designed flow regimes of Tonkin et al. (2021): three invertebrates, riparian plants and fish under four flow
# regimes, checked against the authors' scripts (test/reference/R/run_flow_tradeoffs_reference.R).

using LinearAlgebra, Random

@testset "Designed flow regimes (Tonkin et al. 2021)" begin
    average(x) = sum(x) / length(x)
    csv(f) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", f))[2:end]]
    code = Dict(:plants => "RIP_OPT", :fish => "FSH_OPT", :invertebrates => "INV_OPT", :natural => "NAT")
    taxon = Dict("Mayfly" => :mayfly, "Dragonfly" => :dragonfly, "Ostracod" => :ostracod)
    inverts = csv("flow_inverts.csv")
    riparian = csv("flow_riparian.csv")
    targets = Dict(r[1] => parse.(Float64, r[2:4]) for r in csv("flow_targets.csv"))
    rip = riparian_model()
    bugs = flow_invertebrates()
    @test isempty(pool_names(bugs))
    nonseed(sol, g) = [100 * 350 * sum(sol[g][k, s] / RIPARIAN_GUILDS[g].den[s] for s in 2:6) / rip.K for k in 1:200]
    for s in TONKIN_SCENARIOS
        sc = tonkin_scenario(s)
        @test length(sc.riparian) == 200 && length(sc.fish) == 200 && length(sc.invertebrates) == 30 * 365
        # invertebrates, daily, against R's closed form with numerical integration
        sol = project(bugs, sc.invertebrate_n0, sc.invertebrates)
        for (name, x) in taxon
            rows = [r for r in inverts if r[1] == name && r[2] == code[s]]
            days = [parse(Int, r[3]) for r in rows]
            @test maximum(abs.(sol[x][days .+ 1, 1] .- [parse(Float64, r[5]) for r in rows])) < 1e-9
            @test [sc.invertebrates[d].Kpre[x] for d in days] == [parse(Float64, r[4]) for r in rows]
        end
        target = average((sol[:mayfly][3652:end, 1] .+ sol[:dragonfly][3652:end, 1]) ./ 2)
        @test target ≈ targets[code[s]][3] rtol = 1e-10
        # riparian plants, yearly
        r = project(rip.community, sc.riparian_n0, sc.riparian[1:199])
        rows = [x for x in riparian if x[1] == code[s]]
        for (j, g) in enumerate((:C, :T, :W, :S, :M))
            R = [parse(Float64, x[2 + j]) for x in rows]
            @test maximum(abs.(nonseed(r, g) .- R) ./ max.(abs.(R), 1e-3)) < 1e-12
        end
        @test average(nonseed(r, :C)[11:200]) ≈ targets[code[s]][1] rtol = 1e-12
    end
    # fish under the natural regime, replaying the authors' 100 random starts
    fish = verde_fish_model()
    sc = tonkin_scenario(:natural)
    acc = zeros(200, 7)
    for r in csv("flow_fish_initial.csv")
        r[1] == "NAT" || continue
        f = project(fish.community, [Symbol(sp) => [0.0, 0.0, parse(Float64, r[2 + j])] for (j, sp) in enumerate(VERDE_SPECIES)],
                    sc.fish[1:199])
        for (j, sp) in enumerate(VERDE_SPECIES)
            acc[:, j] .+= (f[Symbol(sp)][:, 2] .+ f[Symbol(sp)][:, 3]) ./ 100
        end
    end
    M = permutedims(reduce(hcat, [parse.(Float64, r[3:9]) for r in csv("flow_fish_mean.csv") if r[1] == "NAT"]))
    @test maximum(abs.(acc .- M) ./ max.(M, 1e-6)) < 1e-12
    @test average(100 .* (acc[11:200, 1] .+ acc[11:200, 2] .+ acc[11:200, 4]) ./ 47660) ≈ targets["NAT"][2] rtol = 1e-12
    # the random starts of the authors' code can be drawn here too
    starts = [tonkin_fish_start(Random.Xoshiro(k)) for k in 1:2000]
    @test average(s[4].second[3] for s in starts) ≈ TONKIN_FISH_INITIAL.CAIN[2] rtol = 0.05
    # the corrected carrying capacity relaxes towards Kb between floods
    c = invertebrate_flows(pulses(91.25, 50, 3650)...; corrected = true)
    @test c[end].Kpre.dragonfly < invertebrate_flows(pulses(91.25, 50, 3650)...)[end].Kpre.dragonfly
end
