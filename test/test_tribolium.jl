# The competition LPA model of two flour beetles (Edmunds et al. 2003; Desharnais et al. 2005), checked against
# the values printed in the papers.

using LinearAlgebra, Random

@testset "Flour beetles (competition LPA)" begin
    @testset "Desharnais et al. (2005)" begin
        com = flour_beetle_competition(TRIBOLIUM_DESHARNAIS_2005)
        @test species_names(com) == [:castaneum, :confusum]
        @test Set(pool_names(com)) == Set([:castaneum_eggs, :castaneum_pupae, :confusum_eggs, :confusum_pupae])
        # the vertical-axis equilibrium (0, 0, 0, 114.04, 91.23, 148.14)
        v = periodic_point(com, [0, 0, 0, 100.0, 100, 100], nothing)
        @test v.converged
        @test round.(v.n[4:6]; digits = 2) == [114.04, 91.23, 148.14]
        # the unstable horizontal-axis equilibrium and its eigenvalues
        h = periodic_point(com, [98.0, 78, 40, 0, 0, 0], nothing)
        @test round.(h.n[1:3]; digits = 2) == [98.36, 78.69, 40.93]
        @test real(sort(h.multipliers; by = real)[1]) ≈ -1.1220 atol = 1e-4
        # the paper's eigenvalues are those of the Jacobian at the equilibrium rounded as printed
        λ = sort(eigvals(jacobian(com, [98.36, 78.69, 40.93, 0, 0, 0], nothing)); by = real)
        @test round(real(λ[1]); digits = 4) == -1.1219
        @test round(real(λ[2]); digits = 4) == -0.1586
        @test round.(real.(λ[3:4]); digits = 4) == [-0.0026, -0.0026]
        @test round.(abs.(imag.(λ[3:4])); digits = 5) == [0.06754, 0.06754]
        @test round.(real.(λ[5:6]); digits = 5) == [0.59352, 0.86528]
        # the two-cycle on the horizontal axis
        c = periodic_point(com, [9.0, 247, 45, 0, 0, 0], nothing; period = 2)
        @test round.(c.n[1:3]; digits = 2) == [9.36, 247.55, 45.48]
        @test round.(advance(com, c.n, nothing)[1:3]; digits = 2) == [309.43, 7.49, 52.59]
    end

    @testset "Edmunds et al. (2003)" begin
        P = TRIBOLIUM_EDMUNDS_2003
        beetles(c) = flour_beetle_competition(merge(P, (c_el = [0.01 c; c 0.01],)))
        com = beetles(0.018)
        # each species alone settles at about 160 animals (Fig. 2)
        a1 = periodic_point(com, [30.0, 30, 100, 0, 0, 0], nothing)
        a2 = periodic_point(com, [0, 0, 0, 30.0, 30, 100], nothing)
        @test 155 < sum(a1.n) < 170 && 155 < sum(a2.n) < 170
        transverse(c, a) = maximum(abs.(periodic_point(beetles(c), a.n, nothing).multipliers))
        # region A: neither axis is stable; B: castaneum's is; D: both are
        @test transverse(0.018, a1) > 1 && transverse(0.018, a2) > 1
        @test transverse(0.022, a1) > 1 && transverse(0.022, a2) < 1
        @test transverse(0.028, a1) < 1 && transverse(0.028, a2) < 1
        # region A: a stable coexistence equilibrium
        e = periodic_point(com, [8.0, 6, 34, 21, 17, 69], nothing)
        @test e.converged && all(>(0), e.n) && maximum(abs.(e.multipliers)) < 1
        # region E: a stable coexistence two-cycle, larvae of both species high and low together
        cyc = periodic_point(beetles(0.032), [0.62, 16.68, 40.14, 3.05, 32.19, 67.78], nothing; period = 2)
        @test cyc.converged && all(>(0), cyc.n) && maximum(abs.(cyc.multipliers)) < 1
        next = advance(beetles(0.032), cyc.n, nothing)
        @test next[1] > 10cyc.n[1] && next[4] > 10cyc.n[4]
        # the adult equation adds surviving adults: with no pupae, adults fall by mu_a
        @test advance(com, [0, 0, 50.0, 0, 0, 0], nothing)[3] ≈ 50 * (1 - 0.02)
    end
end
