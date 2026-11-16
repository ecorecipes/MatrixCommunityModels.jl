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
        @test B * T ≈ I atol = 1e-12
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
        h, rates = table("mizer_rates.csv")
        for name in ("encounter", "feeding_level", "e_repro", "e_growth", "mort")
            j = findfirst(==(name), h)
            @test maximum(rel(getfield(R[named[r[1]]], Symbol(name))[parse(Int, r[2])], num(r[j])) for r in rates) < 1e-12
        end
        # predation mortality is zero or rounding noise of mizer's Fourier transforms in the largest classes
        j = findfirst(==("pred_mort"), h)
        largest = maximum(num(r[j]) for r in rates)
        @test maximum(abs(R[named[r[1]]].pred_mort[parse(Int, r[2])] - num(r[j])) for r in rates) < 1e-12 * largest
        _, repro = table("mizer_reproduction.csv")
        @test maximum(rel(R[named[r[1]]].rdi, num(r[2])) for r in repro) < 1e-12
        @test maximum(rel(R[named[r[1]]].rdd, num(r[3])) for r in repro) < 1e-12
        _, res = table("mizer_resource_rates.csv")
        largest = maximum(num(r[2]) for r in res)
        @test maximum(abs(R[:resource][i] - num(r[2])) for (i, r) in enumerate(res)) < 1e-12 * largest
    end

    @testset "runs" begin
        sol = project(c, n0, mizer_environment(1.0))
        _, step = table("mizer_step.csv")
        @test maximum(rel(sol[named[r[1]]][2, parse(Int, r[2])] / p.dw[parse(Int, r[2])], num(r[3])) for r in step) < 1e-12
        _, res = table("mizer_resource_rates.csv")
        @test maximum(rel(sol[:resource][2, i] / p.dw_full[i], num(r[4])) for (i, r) in enumerate(res)) < 1e-12
        B = mizer_biomass(p, sol)
        _, runs = table("mizer_runs.csv")
        year1 = [r for r in runs if r[1] == "stored_effort" && r[2] == "1"]
        @test maximum(rel(B[11, findfirst(==(named[r[3]]), p.species)], num(r[4])) for r in year1) < 1e-11
    end
end
