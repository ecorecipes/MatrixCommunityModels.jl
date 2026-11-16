# The predator-prey matrix model of Barraquand & Gimenez (2019), checked against the authors' difference
# equations (test/reference/R/run_predator_prey_reference.R).

using LinearAlgebra, Random

@testset "Predator and prey (Barraquand & Gimenez 2019)" begin
    rows = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", "predprey_trajectories.csv"))[2:end]]
    A, B = BARRAQUAND_GIMENEZ_2019, BARRAQUAND_GIMENEZ_2019_SET2
    with7(p, a7) = merge(p, (alpha = Base.setindex(p.alpha, a7, 7),))
    runs = ["set1_table5" => A, "set1_appendix" => with7(A, 1.75), "set2" => B, "set1_high" => with7(A, 5.0),
            "set2_high" => with7(B, 4.7), "no_intra" => merge(A, (alpha = (0.5, 0.0, 0.5, -0.025, 0.0, 0.004, 1.75, 0.0),))]
    for (label, p) in runs
        R = permutedims(reduce(hcat, [parse.(Float64, r[3:6]) for r in rows if r[1] == label]))
        m = predator_prey(p)
        sol = project(m.community, m.n0, fill(nothing, size(R, 1) - 1))
        @test maximum(abs.(sol.u .- R) ./ max.(abs.(R), 1.0)) < 1e-11
    end
    m = predator_prey()
    @test pool_names(m.community) == [:juvenile_prey, :adult_predators]
    # set 1 settles to a stable fixed point
    e = equilibrium(m.community, state_vector(m.community, m.n0...), nothing)
    @test e.converged && maximum(abs.(eigvals(jacobian(m.community, e.n, nothing)))) < 1
    # survival probabilities and marked births allow demographic stochasticity
    sol = project(m.community, m.n0, fill(nothing, 10); rng = Random.Xoshiro(1))
    @test all(isinteger, sol.u) && all(>=(0), sol.u)
end
