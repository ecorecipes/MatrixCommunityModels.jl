# The predator-prey matrix model of Barraquand & Gimenez (2019), checked against the authors' difference
# equations (test/reference/R/run_predator_prey_reference.R).

using LinearAlgebra, Random

@testset "Predator and prey (Barraquand & Gimenez 2019)" begin
    rows = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", "predprey_trajectories.csv"))[2:end]]
    A, B = BARRAQUAND_GIMENEZ_2019, BARRAQUAND_GIMENEZ_2019_SET2
    with7(p, a7) = merge(p, (alpha = Base.setindex(p.alpha, a7, 7),))
    # largest relative difference (floor 1) from the reference, with its tolerance: about 100 times the agreement
    # measured on 2026-11-17, not below 1e-13
    runs = ["set1_table5" => (A, 1e-13),                    # 100 years, achieved 5.0e-16
            "set1_appendix" => (with7(A, 1.75), 1e-13),     # 100 years, achieved 2.0e-16
            "set2" => (B, 1e-13),                           # 100 years, achieved 0.0
            "set1_high" => (with7(A, 5.0), 8e-12),          # 500 years of cycles, achieved 7.7e-14
            "set2_high" => (with7(B, 4.7), 1e-11),          # 500 years of cycles, achieved 2.4e-13
            "no_intra" => (merge(A, (alpha = (0.5, 0.0, 0.5, -0.025, 0.0, 0.004, 1.75, 0.0),)), 1e-13)]  # achieved 4.0e-16
    for (label, (p, tol)) in runs
        R = permutedims(reduce(hcat, [parse.(Float64, r[3:6]) for r in rows if r[1] == label]))
        m = predator_prey(p)
        sol = project(m.community, m.n0, fill(nothing, size(R, 1) - 1))
        @test maximum(abs.(sol.u .- R) ./ max.(abs.(R), 1.0)) < tol
    end
    m = predator_prey()
    @test pool_names(m.community) == [:juvenile_prey, :adult_predators]
    # both sets settle to a stable fixed point, as vignette 23 prints them: set 1 prey (113.7, 170.5), predators
    # (13.9, 32.4), largest |eigenvalue| 0.883; set 2 prey (5.0, 7.5), predators (22.4, 52.2), 0.962
    for (p, n, λ) in ((A, [113.7, 170.5, 13.9, 32.4], 0.883), (B, [5.0, 7.5, 22.4, 52.2], 0.962))
        mm = predator_prey(p)
        e = equilibrium(mm.community, state_vector(mm.community, mm.n0...), nothing)
        @test e.converged && round.(e.n; digits = 1) == n
        @test round(maximum(abs.(eigvals(jacobian(mm.community, e.n, nothing)))); digits = 3) == λ
    end
    # survival probabilities and marked births allow demographic stochasticity
    sol = project(m.community, m.n0, fill(nothing, 10); rng = Random.Xoshiro(1))
    @test all(isinteger, sol.u) && all(>=(0), sol.u)
    # with demographic stochasticity (vignette: 500 runs of 100 years from Xoshiro(2019)) set 1 never loses its
    # prey; set 2 (about 13 prey at equilibrium) loses them in 0.0, 0.1 and 0.672 of runs by years 10, 30 and 100;
    # neither loses its predators. Set 1 with 200 runs here; set 2 with the vignette's 500, standard errors 0.013
    # and 0.021 at years 30 and 100, so the tolerances are about 3.5 standard errors.
    lost(p, replicates) = (mm = predator_prey(p);
                           ens = ensemble(mm.community, state_vector(mm.community, mm.n0...), fill(nothing, 100);
                                          replicates, rng = Xoshiro(2019));
                           (prey = [count(abundance(ens, :prey)[:, t + 1] .== 0) / replicates for t in (10, 30, 100)],
                            predators = count(abundance(ens, :predator)[:, end] .== 0) / replicates))
    one, two = lost(A, 200), lost(B, 500)
    @test one.prey == [0, 0, 0] && one.predators == 0
    @test two.prey[1] == 0 && two.predators == 0
    @test two.prey[2] ≈ 0.1 atol = 0.05
    @test two.prey[3] ≈ 0.672 atol = 0.07
end
