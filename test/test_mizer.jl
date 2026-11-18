# The size-spectrum model of mizer 3.4.0 (North Sea, NS_params) against mizer run in R
# (test/reference/R/run_mizer_reference.R).
using Test, MatrixCommunityModels, LinearAlgebra

@testset "mizer size spectrum" begin
    ref = joinpath(@__DIR__, "reference")
    table(f) = (ls = [split(l, ",") for l in readlines(joinpath(ref, f))]; (String.(ls[1]), ls[2:end]))
    num(x) = parse(Float64, x)
    rel(a, b) = b == 0 ? abs(a) : abs(a / b - 1)
    p = mizer_parameters()
    named = Dict(n => s for (n, s) in zip(p.names, p.species))

    @testset "data and derived quantities" begin
        @test length(p.species) == 12 && p.species[3] == :n_pout
        @test length(p.w) == 100 && length(p.w_full) == 226 && p.offset == 126
        @test MIZER_NORTH_SEA.effort == (industrial = 0.0, pelagic = 1.0, beam = 0.5, otter = 0.5)
        # grid, derived arrays and kernel: achieved <= 4.4e-16 relative, kernel 8.4e-15 absolute (2026-11-17);
        # 100 times that would loosen the tolerances below, which are kept
        _, grid = table("mizer_grid.csv")
        @test maximum(rel(p.w_full[i], num(r[2])) for (i, r) in enumerate(grid)) < 1e-14
        @test maximum(rel(p.dw_full[i], num(r[3])) for (i, r) in enumerate(grid)) < 1e-14
        @test maximum(rel(p.resource.rate[i], num(r[4])) for (i, r) in enumerate(grid)) < 1e-14
        @test maximum(rel(p.resource.capacity[i], num(r[5])) for (i, r) in enumerate(grid)) < 1e-14
        h, arrays = table("mizer_arrays.csv")
        for name in ("search_vol", "intake_max", "metab", "psi", "mu_b")
            j = findfirst(==(name), h)
            @test maximum(rel(getfield(p.fish[named[r[1]]], Symbol(name))[parse(Int, r[2])], num(r[j])) for r in arrays) < 1e-14
        end
        j = findfirst(==("f_mort"), h)
        @test all(get(MIZER_NORTH_SEA.effort, p.fish[named[r[1]]].gear, 0.0) * p.fish[named[r[1]]].catchability *
                  p.fish[named[r[1]]].selectivity[parse(Int, r[2])] == num(r[j]) for r in arrays)
        _, kernel = table("mizer_kernel.csv")
        @test maximum(abs(p.kernel[findfirst(==(named[r[1]]), p.species), parse(Int, r[2]) + 1] - num(r[3])) for r in kernel) < 1e-13
    end

    @testset "the step is mizer's implicit upwind step" begin
        # T solves the bidiagonal system of mizer in numbers: B T = I with B[j, j] = b_j, B[j, j-1] = -g dt / dw
        g = rand(100) .* p.w .* 3
        mort = rand(100)
        T = MatrixCommunityModels._mizer_step_matrix(p, g, mort)
        B = zeros(100, 100)
        for j in 1:100
            B[j, j] = 1 + g[j] * p.dt / p.dw[j] + mort[j] * p.dt
            j > 1 && (B[j, j - 1] = -g[j - 1] * p.dt / p.dw[j - 1])
        end
        @test B * T ≈ I atol = 1e-13    # norm of B T - I achieved 6e-16 for one draw (2026-11-17)
        @test all(T .>= 0) && istril(T)
    end

    c = mizer_community(p)
    n0 = state_vector(c, mizer_state(p)...)
    env = (effort = MIZER_NORTH_SEA.effort,)

    @testset "the community" begin
        @test species_names(c) == vcat(p.species, :resource)
        @test Set(pool_names(c)) == Set([:prey, :consumers])
        pools = pool_totals(c, n0, env)
        @test length(pools.prey) == 1200 && length(pools.consumers) == 1200
        @test Set(transition_names(c.system[:cod])) == Set([:staying, :growing, :spawning])
        @test_throws ArgumentError mizer_fish(p, :tuna)
    end

    @testset "rates at the stored state" begin
        R = mizer_rates(p, c, n0, env)
        # achieved (2026-11-17): encounter 2.9e-14, feeding level 1.3e-14, e_repro 3.5e-14, e_growth 3.8e-14, mort
        # 2.3e-14, R_di 1.6e-14, R_dd 9.0e-15 relative; 100 times would loosen 1e-12, which is kept. Predation and
        # resource mortality, relative to their largest values: 5.3e-15 and 4.6e-15
        h, rates = table("mizer_rates.csv")
        for name in ("encounter", "feeding_level", "e_repro", "e_growth", "mort")
            j = findfirst(==(name), h)
            @test maximum(rel(getfield(R[named[r[1]]], Symbol(name))[parse(Int, r[2])], num(r[j])) for r in rates) < 1e-12
        end
        # predation mortality is zero or rounding noise of mizer's Fourier transforms in the largest classes
        j = findfirst(==("pred_mort"), h)
        largest = maximum(num(r[j]) for r in rates)
        @test maximum(abs(R[named[r[1]]].pred_mort[parse(Int, r[2])] - num(r[j])) for r in rates) < 5e-13 * largest
        _, repro = table("mizer_reproduction.csv")
        @test maximum(rel(R[named[r[1]]].rdi, num(r[2])) for r in repro) < 1e-12
        @test maximum(rel(R[named[r[1]]].rdd, num(r[3])) for r in repro) < 1e-12
        _, res = table("mizer_resource_rates.csv")
        largest = maximum(num(r[2]) for r in res)
        @test maximum(abs(R[:resource][i] - num(r[2])) for (i, r) in enumerate(res)) < 5e-13 * largest
    end

    @testset "runs" begin
        # Tolerances are about 100 times the agreement measured on 2026-11-17 (macOS, Julia 1.12), or the earlier,
        # tighter tolerance where 100 times would loosen it.
        _, runs = table("mizer_runs.csv")
        _, states = table("mizer_states.csv")
        _, steady = table("mizer_steady.csv")
        settled = Dict(r[1] => round(Int, 10 * num(r[2])) for r in steady)    # steps findSteadyState took
        @test settled == Dict("stored_effort" => 15, "unfished" => 555) && all(r[3] == "TRUE" for r in steady)
        # yearly biomass of every species up to `years`, relative
        biomass_gap(name, sol, years) = (B = mizer_biomass(p, sol);
            maximum(rel(B[10 * parse(Int, r[2]) + 1, findfirst(==(named[r[3]]), p.species)], num(r[4]))
                    for r in runs if r[1] == name && parse(Int, r[2]) <= years))
        # the whole state, relative, in the classes holding more than 1e-6 of the species' largest
        function state_gap(run, n)
            worst = 0.0
            for r in states
                r[1] == run || continue
                sp = r[2] == "Resource" ? :resource : named[r[2]]
                j = parse(Int, r[3])
                theirs = num(r[4]) * (sp === :resource ? p.dw_full[j] : p.dw[j])
                theirs > 1e-6 * maximum(n[c.layout[sp]]) || continue
                worst = max(worst, rel(n[c.layout[sp]][j], theirs))
            end
            return worst
        end
        # stored effort, 20 years: its first step, its 15th (where findSteadyState stopped) and its end
        sol = project(c, n0, mizer_environment(20.0))
        _, step = table("mizer_step.csv")
        # one step: achieved 1.8e-14 (100 times would loosen 1e-12)
        @test maximum(rel(sol[named[r[1]]][2, parse(Int, r[2])] / p.dw[parse(Int, r[2])], num(r[3])) for r in step) < 1e-12
        _, res = table("mizer_resource_rates.csv")
        @test maximum(rel(sol[:resource][2, i] / p.dw_full[i], num(r[4])) for (i, r) in enumerate(res)) < 2.5e-13  # achieved 2.4e-15
        @test biomass_gap("stored_effort", sol, 1) < 2e-12                        # achieved 1.7e-14 to year 1
        @test biomass_gap("stored_effort", sol, 20) < 6e-12                       # achieved 6.0e-14
        @test state_gap("stored_effort", sol.u[end, :]) < 1e-11                   # achieved 1.0e-13
        @test state_gap("steady_stored_effort", sol.u[settled["stored_effort"] + 1, :]) < 7e-12   # achieved 7.2e-14
        # otter trawling at three times the stored effort, 30 years
        otter = project(c, n0, mizer_environment(30.0; effort = merge(MIZER_NORTH_SEA.effort, (otter = 1.5,))))
        @test biomass_gap("otter_1.5", otter, 30) < 2.5e-11                       # achieved 2.5e-13
        @test state_gap("otter_1.5", otter.u[end, :]) < 4e-11                     # achieved 3.8e-13
        # unfished: one run of 555 steps is both mizer's 30-year run (its first 300 steps) and findSteadyState's.
        # Zero effort for every gear gives the same states as the vignette's empty effort `(;)` (checked
        # 2026-11-17) and saves compiling a second environment type, about 4 s.
        unfished = project(c, n0, mizer_environment(settled["unfished"] / 10; effort = map(_ -> 0.0, MIZER_NORTH_SEA.effort)))
        @test size(unfished.u, 1) == 556
        @test biomass_gap("unfished", unfished, 30) < 3e-11                       # achieved 2.7e-13
        @test state_gap("unfished", unfished.u[301, :]) < 2e-10                   # achieved 1.9e-12
        @test state_gap("steady_unfished", unfished.u[end, :]) < 1.5e-10          # achieved 1.4e-12
        # biomass relative to unfished, as vignette 25 prints it: steady state at stored effort, and otter x3 after
        # 30 years
        final(B, k) = mizer_biomass(p, B)[k, :]
        @test round.(final(sol, 16) ./ final(unfished, 556); digits = 3) ==
              [1.290, 1.003, 0.961, 0.776, 0.622, 0.904, 0.811, 0.604, 0.831, 1.016, 0.497, 0.816]
        @test round.(final(otter, 301) ./ final(unfished, 301); digits = 3) ==
              [0.584, 0.767, 0.937, 0.751, 0.436, 0.786, 0.636, 0.068, 1.515, 0.985, 0.225, 0.627]
    end
end
