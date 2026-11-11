# The riparian guild MCM of Lytle et al. (2017), rebuilt by composition and compared with
# the output of the original R script (test/reference/R/run_riparian_reference.R).

function riparian_reference(file)
    rows = [parse.(Float64, split(l, ",")) for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]
    M = permutedims(reduce(hcat, rows))
    return (y = Int.(M[:, 2]), N = M[:, 5:34])
end

@testset "Riparian guilds (Lytle et al. 2017)" begin
    m = riparian_model()
    @test length(m.envs) == 83
    @test sum(e.flood for e in m.envs) == 47 && sum(e.drought for e in m.envs) == 17
    @test m.K == 41_770_400
    @test pool_names(m.community) == [:space, :willow_adults]

    for file in ("rip_seq_traj.csv", "rip_seed_traj.csv")
        ref = riparian_reference(file)
        sol = project(m.community, m.n0, m.envs[ref.y[1:249]])
        err = maximum(abs.(sol.u .- ref.N) ./ max.(abs.(ref.N), 1.0))
        println("riparian $file: max scaled difference from R = ", err,
                "; identical entries = ", count(sol.u .== ref.N), " of ", length(ref.N))
        @test err < 1e-12
    end

    @testset "corrected variant" begin
        c = riparian_model(; corrected = true)
        @test pool_names(c.community) == [:space]          # no hidden willow-meadow wire
        ref = riparian_reference("rip_seq_traj.csv")
        envs = c.envs[ref.y[1:249]]
        orig = project(m.community, m.n0, envs)
        corr = project(c.community, c.n0, envs)
        # cottonwood and tamarisk start identically; the variants then diverge
        @test corr.u[2, 1:12] == orig.u[2, 1:12]
        @test corr[:M] != orig[:M]
        @test all(isfinite, corr.u) && all(>=(0), corr.u)
        # stage-5 sagebrush that does not mature is retained only in the corrected model
        env = c.envs[1]
        ctx = Context(ones(6), (space = 0.0,), env)
        Sc = only(s for s in MatrixCommunityModels._species(c.community) if s.name === :S)
        @test projection_matrix(Sc, ctx)[5, 5] > 0
    end
end
