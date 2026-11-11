# The two-species Cypripedium community of adapt3's comproj vignette, rebuilt by
# composition and compared with project3() (test/reference/R/run_adapt3_reference.R).

function adapt3_reference(tag)
    rows = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", "adapt3_$(tag)_N.csv"))[2:end]]
    return permutedims(reduce(hcat, [parse.(Float64, r[2:4]) for r in rows]))
end

@testset "adapt3 project3 (Cypripedium)" begin
    for (tag, equivalence, integeronly) in (("continuous", true, false), ("noequiv", false, false),
                                            ("vignette", true, true))
        m = cypripedium_model(; equivalence, integeronly)
        @test length.(stage_names.(Ref(m.community), species_names(m.community))) == [9, 9]
        sol = project(m.community, m.n0, m.envs)
        ref = adapt3_reference(tag)
        got = hcat(sum(sol[:candidum]; dims = 2), sum(sol[:parviflorum]; dims = 2))
        err = maximum(abs.(got .- ref[:, 1:2]) ./ max.(abs.(ref[:, 1:2]), 1.0))
        println("adapt3 $tag: max scaled difference from project3() = ", err)
        @test err < 1e-10
        # aggregate density used for each step
        @test sol.pools[1:100, 1] ≈ ref[1:100, 3] rtol = 1e-10
    end
end
