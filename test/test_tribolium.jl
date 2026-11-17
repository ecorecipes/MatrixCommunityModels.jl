# The competition LPA model of two flour beetles (Edmunds et al. 2003; Desharnais et al. 2005), checked against
# the values printed in the papers and the further values vignette 21 reports.

using LinearAlgebra, Random

@testset "Flour beetles (competition LPA)" begin
    steps(com, n, k) = (for _ in 1:k; n = advance(com, n, nothing); end; n)
    stable(x) = maximum(abs.(x.multipliers)) < 1

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
        # at the exact equilibrium, as vignette 21 prints them (-1.1220203 here, 2026-11-17)
        μ = sort(h.multipliers; by = real)
        @test round.(real.(μ[[1, 2, 3, 5, 6]]); digits = 5) == [-1.12202, -0.15859, -0.00264, 0.59356, 0.86528]
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
        # the invasion criterion of T. castaneum at the T. confusum equilibrium, printed as 0.29 (0.28861 here)
        p = TRIBOLIUM_DESHARNAIS_2005
        L2, A2 = v.n[4], v.n[6]
        @test round(p.b[1] * (1 - p.mu_l[1]) / p.mu_a[1] *
                    exp(-p.c_el[1, 2] * L2 - (p.c_ea[1, 2] + p.c_pa[1, 2]) * A2); digits = 2) == 0.29
        # varying c_el[2, 1] from four adults of each (paper: T. confusum alone below 0.061, invariant
        # double loop 0.061 to 0.092, two-cycle 0.092 to 0.17)
        with21(x) = flour_beetle_competition(merge(p, (c_el = [0.012 0.01; x 0.003],)))
        outcome = map((0.05, 0.07, 0.10, 0.15)) do x
            com21 = with21(x)
            n = steps(com21, [0, 0, 4.0, 0, 0, 4], 5000)
            e = periodic_point(com21, n, nothing; period = 2)
            sum(n[1:3]) < 1e-6 ? :castaneum_lost :
            e.converged && maximum(abs.(e.n .- n)) < 1e-6 && stable(e) ? :two_cycle : :not_on_a_cycle
        end
        @test outcome == (:castaneum_lost, :not_on_a_cycle, :two_cycle, :two_cycle)
        # demographic stochasticity on the double loop, 150 steps: paper 30% both persist, 15% T. castaneum
        # lost, 55% T. confusum lost. 500 runs here (vignette: 2000), so a standard error of about 0.02 a
        # share; this seed gives 0.278, 0.146, 0.576 (2026-11-17). atol 0.07 is about 3.5 standard errors.
        ens = ensemble(with21(0.062), [0, 0, 4.0, 0, 0, 4], fill(nothing, 150); replicates = 500, rng = Xoshiro(2005))
        lost1 = abundance(ens, :castaneum)[:, end] .== 0
        lost2 = abundance(ens, :confusum)[:, end] .== 0
        share(b) = count(b) / length(b)
        @test share(.!lost1 .& .!lost2) ≈ 0.30 atol = 0.07
        @test share(lost1) ≈ 0.15 atol = 0.07
        @test share(lost2) ≈ 0.55 atol = 0.07
        @test !any(lost1 .& lost2)
    end

    @testset "Edmunds et al. (2003)" begin
        P = TRIBOLIUM_EDMUNDS_2003
        beetles(c) = flour_beetle_competition(merge(P, (c_el = [0.01 c; c 0.01],)))
        # region A test point: c = 0.016, the vignette's first point, well below the period-doubling bubble
        # at 0.0184 to 0.019 (the leading multiplier is 0.994 here, 0.997 at 0.018 and -0.998 at 0.0182)
        com = beetles(0.016)
        # each species alone settles at 163.7 (T. confusum) and 158.3 (T. castaneum) animals, whatever c
        # (about 160 in Fig. 2)
        a1 = periodic_point(com, [30.0, 30, 100, 0, 0, 0], nothing)
        a2 = periodic_point(com, [0, 0, 0, 30.0, 30, 100], nothing)
        @test round(sum(a1.n); digits = 1) == 163.7 && round(sum(a2.n); digits = 1) == 158.3
        transverse(c, a) = maximum(abs.(periodic_point(beetles(c), a.n, nothing).multipliers))
        # region A: neither axis is stable; B: castaneum's is; C: castaneum's only; D: both are
        @test transverse(0.016, a1) > 1 && transverse(0.016, a2) > 1
        @test transverse(0.022, a1) > 1 && transverse(0.022, a2) < 1
        @test transverse(0.0248, a1) > 1 && transverse(0.0248, a2) < 1
        @test transverse(0.028, a1) < 1 && transverse(0.028, a2) < 1
        # region A: a stable coexistence equilibrium with 61.5 T. confusum and 97.3 T. castaneum (about 58
        # and 95 read from Fig. 2)
        e = periodic_point(com, [8.0, 6, 34, 21, 17, 69], nothing)
        @test e.converged && all(>(0), e.n) && stable(e)
        @test round(sum(e.n[1:3]); digits = 1) == 61.5 && round(sum(e.n[4:6]); digits = 1) == 97.3
        # the boundaries, as vignette 21 prints them (0.0201, 0.0243, 0.0252, 0.0289; Fig. 2: 0.0200, about
        # 0.0243, 0.0252, 0.0288): each lies inside the interval that rounds to the printed value. Exact
        # values by bisection (2026-11-17): 0.0200514, 0.0242723, 0.0252493, 0.0289294.
        resists(c, alone) = stable(periodic_point(beetles(c), alone.n, nothing))
        @test !resists(0.02005, a2) && resists(0.02015, a2)          # A-B: T. castaneum resists invasion
        @test !resists(0.02515, a1) && resists(0.02525, a1)          # C-D: T. confusum resists invasion
        function coexistence_equilibrium(c)      # reached from T. confusum alone with a few T. castaneum adults
            com = beetles(c)
            x = steps(com, a1.n .+ [0, 0, 0, 0, 0, 0.5], 20000)
            min(sum(x[1:3]), sum(x[4:6])) < 1e-3 && return false
            q = periodic_point(com, x, nothing)
            q.converged && all(>(0), q.n) && stable(q)
        end
        function coexistence_cycle(c)            # reached from 16 adults of each species
            com = beetles(c)
            x = steps(com, [0, 0, 16.0, 0, 0, 16], 20000)
            min(sum(x[1:3]), sum(x[4:6])) < 1e-3 && return false
            q = periodic_point(com, x, nothing; period = 2)
            q.converged && all(>(0), q.n) && stable(q) && maximum(abs.(advance(com, q.n, nothing) .- q.n)) > 1e-6
        end
        @test !coexistence_equilibrium(0.02425) && coexistence_equilibrium(0.02435)   # B-C
        @test coexistence_equilibrium(0.0248)    # region C: beside the stable T. castaneum state (above)
        @test !coexistence_cycle(0.02885) && coexistence_cycle(0.02895)               # D-E
        # the period-doubling bubble inside region A, from 0.0184 to 0.019 (exact 0.0183656 and 0.0190118):
        # the leading multiplier is -1.00156 at 0.0186
        lead(c) = (q = periodic_point(beetles(c), [8.0, 6, 34, 21, 17, 69], nothing); q.multipliers[argmax(abs.(q.multipliers))])
        @test abs(lead(0.01835)) < 1 && abs(lead(0.01845)) > 1 && abs(lead(0.01895)) > 1 && abs(lead(0.01905)) < 1
        @test round(real(lead(0.0186)); digits = 5) == -1.00156
        # and the small two-cycle there: adults (31.06, 71.81) and (30.95, 71.6)
        b = beetles(0.0186)
        q = periodic_point(b, steps(b, periodic_point(b, [8.0, 6, 34, 21, 17, 69], nothing).n .+ 0.01, 10000), nothing; period = 2)
        @test q.converged && stable(q)
        @test round.(q.n[[3, 6]]; digits = 2) == [31.06, 71.81]
        @test round.(advance(b, q.n, nothing)[[3, 6]]; digits = 2) == [30.95, 71.6]
        # region E: a stable coexistence two-cycle, larvae of both species high and low together:
        # larvae (0.62, 3.05) then (20.85, 40.23), adults (40.1, 67.8) then (40.9, 69.0)
        cyc = periodic_point(beetles(0.032), [0.62, 16.68, 40.14, 3.05, 32.19, 67.78], nothing; period = 2)
        @test cyc.converged && all(>(0), cyc.n) && stable(cyc)
        next = advance(beetles(0.032), cyc.n, nothing)
        @test next[1] > 10cyc.n[1] && next[4] > 10cyc.n[4]
        @test round.(cyc.n[[1, 4]]; digits = 2) == [0.62, 3.05] && round.(next[[1, 4]]; digits = 2) == [20.85, 40.23]
        @test round.(cyc.n[[3, 6]]; digits = 1) == [40.1, 67.8] && round.(next[[3, 6]]; digits = 1) == [40.9, 69.0]
        # basins from adults only, 5000 steps, along 16 T. castaneum adults and 2 to 128 T. confusum
        # (1: T. confusum wins, 2: T. castaneum wins, 3: both persist): at 0.028 T. castaneum wins unless T.
        # confusum starts with several times as many; at 0.032 half to four times as many leads to coexistence
        outcomes(c) = [(x = steps(beetles(c), [0, 0, a, 0, 0, 16.0], 5000);
                        sum(x[4:6]) < 1e-3 ? 1 : sum(x[1:3]) < 1e-3 ? 2 : 3) for a in 2.0 .^ (1:7)]
        @test outcomes(0.028) == [2, 2, 2, 2, 2, 2, 1]
        @test outcomes(0.032) == [2, 2, 3, 3, 3, 3, 1]
        # the coexistence culture's start (32 T. confusum, 4 T. castaneum adults) loses T. castaneum
        @test all(sum(steps(beetles(c), [0, 0, 32.0, 0, 0, 4], 20000)[4:6]) < 1e-3 for c in (0.028, 0.030, 0.032))
        # the adult equation adds surviving adults: with no pupae, adults fall by mu_a
        @test advance(com, [0, 0, 50.0, 0, 0, 0], nothing)[3] ≈ 50 * (1 - 0.02)
    end
end
