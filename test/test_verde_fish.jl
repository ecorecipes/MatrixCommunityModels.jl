# The Verde River fish MCM of Rogosch et al. (2019), rebuilt by composition and compared
# with the output of the original R script (test/reference/R/run_fish_reference.R).

@testset "Verde River fish (Rogosch et al. 2019)" begin
    m = verde_fish_model()
    fish = m.community
    @test species_names(fish) == Symbol.(VERDE_SPECIES)
    @test pool_names(fish) == [:biomass]
    @test length(m.envs) == 54 && m.years[1] == 1964

    sol = project(fish, m.n0, m.envs)

    rows = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", "fish_det_biomass.csv"))[2:end]]
    ref = permutedims(reduce(hcat, [parse.(Float64, r[2:end]) for r in rows]))
    @test size(ref) == (54, 21)
    got = sol.u[1:54, :]
    err = maximum(abs.(got .- ref) ./ max.(abs.(ref), 1e-300))
    println("fish: max relative difference from R = ", err, "; identical entries = ", count(got .== ref), " of ", length(ref))
    @test err < 1e-13

    # state after the 54th update, printed by the R driver
    @test sol[:CACL][end, :] ≈ [11.886757092113552, 62.746436211955050, 153.259534732547678] rtol = 1e-13
    @test sol[:CAIN][end, :] ≈ [563.86632120713205, 3838.72013864657811, 7816.54696811462236] rtol = 1e-13

    # the pool is total biomass, and the frozen community matrix reproduces one step
    @test sol.pools[1, 1] == sum(sol.u[1, :])
    A = community_matrix(fish, sol.u[1, :], m.envs[1])
    @test A * sol.u[1, :] ≈ sol.u[2, :] rtol = 1e-14
end
