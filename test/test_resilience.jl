# Resilience metrics through ResilienceCore.jl: linearise, periodic_operator and the community interface (MCMState)
using LinearAlgebra
import ResilienceCore

@testset "ResilienceCore interface" begin
    RC = ResilienceCore
    bykey(z) = (real(z), imag(z))

    @testset "one namespace with ResilienceCore" begin
        # names the two packages share are the same functions, so loading both is unambiguous
        common = intersect(names(MatrixCommunityModels), names(RC))
        @test :jacobian in common && :linearise in common
        @test all(name -> getfield(MatrixCommunityModels, name) === getfield(RC, name), common)
        together = Module()
        Core.eval(together, :(using MatrixCommunityModels, ResilienceCore))
        @test Core.eval(together, :jacobian) === RC.jacobian
        @test Core.eval(together, :linearise) === RC.linearise
    end

    @testset "linearise at fixed points and cycles" begin
        com = flour_beetle_competition(TRIBOLIUM_DESHARNAIS_2005)
        # the unstable horizontal-axis equilibrium: the eigenvalues are periodic_point's multipliers
        h = periodic_point(com, [98.0, 78, 40, 0, 0, 0], nothing)
        J = linearise(com, h.n, nothing)
        @test J isa RC.DiscreteJacobian
        @test RC.operator(J) == jacobian(com, h.n, nothing)
        @test sort(eigvals(RC.operator(J)); by = bykey) ≈ sort(h.multipliers; by = bykey)
        @test !RC.isstable(J)
        # the two-cycle on the same axis: the Jacobian of the two-step map
        c = periodic_point(com, [9.0, 247, 45, 0, 0, 0], nothing; period = 2)
        J2 = linearise(com, c.n, nothing; period = 2)
        @test sort(eigvals(RC.operator(J2)); by = bykey) ≈ sort(c.multipliers; by = bykey)
        @test_throws ArgumentError linearise(com, c.n, nothing)                 # not a fixed point of one step
        @test_throws ArgumentError linearise(com, c.n, nothing; period = 0)
        @test linearise(com, c.n, nothing; check = false) isa RC.DiscreteJacobian
        # a stable coexistence equilibrium (Edmunds et al. 2003, region A)
        beetles = flour_beetle_competition(merge(TRIBOLIUM_EDMUNDS_2003, (c_el = [0.01 0.018; 0.018 0.01],)))
        e = periodic_point(beetles, [8.0, 6, 34, 21, 17, 69], nothing)
        Je = linearise(beetles, e.n, nothing)
        @test RC.isstable(Je) && RC.resilience(Je) > 0
        @test Je.equilibrium == e.n
    end

    @testset "periodic_operator" begin
        # an annual plant that germinates in season 1 and sets seed in season 2 (as in test_analysis.jl)
        life = ValuedProjectionNet([:seed, :plant],
            :dormancy => [(:seed => :seed) => 0.9], :germination => [(:seed => :plant) => 0.5],
            :growth => [(:plant => :plant) => 0.9], :seed_set => [(:plant => :seed) => 10.0])
        on(season) = ctx -> ctx.env == season ? 1.0 : 0.0
        sp = species(:x, respond(life, :germination => on(1), :seed_set => on(2),
                                 :germination => ctx -> exp(-ctx.pool.space), :growth => ctx -> ctx.env == 2 ? 0.0 : 1.0);
                     occupies = [:space => [:plant => 1.0]])
        empty = (space = 0.0,)
        P = periodic_operator(sp, empty, [1, 2])
        @test P isa RC.PeriodicOperator && RC.period(P) == 2
        @test RC.operator(RC.monodromy(P)) ≈ periodic_matrix(sp, empty, [1, 2])
        @test RC.growth_rate(P) ≈ sqrt(maximum(abs, eigvals(periodic_matrix(sp, empty, [1, 2]))))
        @test_throws ArgumentError periodic_operator(sp, empty, Int[])
        # a continuous stage carries its mesh weights
        plant = ContinuousStage(:plant, 0.0, 10.0, 50)
        normpdf(x, m, s) = exp(-0.5 * ((x - m) / s)^2) / (s * sqrt(2pi))
        ipm = species(:plant, kernel_net(plant, :growing => (y, x) -> 0.8 * normpdf(y, 1.2 + 0.85x, 0.7),
                                         :seeding => (y, x) -> 0.2 * normpdf(y, 1.0, 0.5)))
        Q = periodic_operator(ipm, NamedTuple(), [nothing])
        @test RC.weights(only(Q.parts)) ≈ fill(step_size(plant), 50)
    end

    @testset "MCMState: presses and removals" begin
        beetles = flour_beetle_competition(merge(TRIBOLIUM_EDMUNDS_2003, (c_el = [0.01 0.018; 0.018 0.01],)))
        e = periodic_point(beetles, [8.0, 6, 34, 21, 17, 69], nothing)
        s = MCMState(beetles, e.n, nothing)
        @test s.converged
        @test s.n ≈ e.n rtol = 1e-8
        @test RC.richness(s) == 2 && RC.present_species(s) == [1, 2]
        @test RC.state(s) ≈ [sum(e.n[1:3]), sum(e.n[4:6])]
        @test RC.abundances(s) == RC.state(s)
        # the interface's Jacobian is the community's, and linearise gives a DiscreteJacobian
        @test RC.jacobian(s) ≈ jacobian(beetles, s.n, nothing) rtol = 1e-6
        @test linearise(s) isa RC.DiscreteJacobian && RC.isstable(linearise(s))
        # baseline mortality factors leave the community where it is
        @test RC.state(RC.with_mortality(s, [1.0, 1.0])) ≈ RC.state(s) rtol = 1e-8
        # a small press on one species lowers it and raises its competitor
        p = RC.with_mortality(s, [1.02, 1.0])
        @test p.converged && RC.richness(p) == 2
        @test RC.state(p)[1] < RC.state(s)[1] && RC.state(p)[2] > RC.state(s)[2]
        # a large one excludes it, and the competitor settles where it lives alone
        q = RC.with_mortality(s, [1.5, 1.0])
        @test RC.present_species(q) == [2] && RC.state(q)[1] == 0
        alone = MCMState(beetles, [0, 0, 0, 30.0, 30, 100], nothing)
        @test RC.state(q)[2] ≈ RC.state(alone)[2] rtol = 1e-6
        # removal: the other species is released from competition
        r = RC.without_species(s, 1)
        @test RC.present_species(r) == [2] && RC.state(r)[2] > RC.state(s)[2]
        @test_throws ArgumentError RC.without_species(r, 1)
        # ResilienceCore's experiments run unchanged
        rel, l1, lost = RC.press_response(s, [1.02, 1.0])
        @test lost == 0 && l1 > 0
        # the local tolerance is the press that first excludes the species
        tol = first(RC.local_mortality_tolerances(s))
        @test RC.press_response(s, [1 + tol, 1.0])[3] == 1
        @test RC.press_response(s, [1 + tol - 0.1, 1.0])[3] == 0
        @test [x.survivors for x in RC.removal_responses(s)] == [1, 1]
        @test RC.global_mortality_tolerance(s) > 1
        @test all(>(0), RC.local_mortality_tolerances(s))
        @test_throws ArgumentError RC.with_mortality(s, [1.0])
        # a state that cycles is not accepted as a steady state
        cyc = flour_beetle_competition(merge(TRIBOLIUM_EDMUNDS_2003, (c_el = [0.01 0.032; 0.032 0.01],)))
        @test_throws ArgumentError MCMState(cyc, [0.62, 16.68, 40.14, 3.05, 32.19, 67.78], nothing; newton = false,
                                            maxiter = 2_000)
    end
end
