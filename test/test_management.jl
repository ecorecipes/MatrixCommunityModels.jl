# Managed communities (src/management.jl): parameters, actions and outcomes, against loops written by hand on
# plain matrices. The DecisionAnalysis.jl adapter built on them is tested in test/decision_analysis.

using Random: Xoshiro
const StructuredPopulationCore = MatrixCommunityModels.StructuredPopulationCore

@testset "Managed communities" begin
    # nothing here needs DecisionAnalysis.jl: in the package's test environment its extension is not loaded
    @test Base.get_extension(MatrixCommunityModels, :MatrixCommunityModelsDecisionAnalysisExt) === nothing

    # two species with constant matrices; the fecundity of `a` is multiplied by `env.fecundity`
    U = [0.0 0.0; 0.5 0.8]
    F = [0.0 2.0; 0.0 0.0]
    function pair_of(survival_b)
        a = births(respond(matrix_net([:juvenile, :adult], :survival => U, :breeding => F),
                           :breeding => ctx -> ctx.env.fecundity), :breeding)
        b = matrix_net([:only], :survival => fill(survival_b, 1, 1))
        return community(species(:a, a), species(:b, b))
    end
    c = pair_of(0.9)
    n0 = [:a => [10.0, 20.0], :b => [50.0]]
    envs = fill((fecundity = 1.0,), 6)
    step(n, fecundity, sb) = [fecundity * 2.0 * n[2], 0.5 * n[1] + 0.8 * n[2], sb * n[3]]   # the model, by hand

    @testset "no action is the plain projection" begin
        m = ManagedCommunity(c, n0, envs)
        run = project(m)
        plain = project(c, n0, envs)
        @test run isa ManagedTrajectory
        @test run.trajectory.u == plain.u
        @test run.trajectory.t == plain.t
        @test all(iszero, run.removed) && all(iszero, run.added)
        @test run.problem === m
        @test m.outcomes(run) == (a = sum(plain.u[end, 1:2]), b = plain.u[end, 3])
        @test project(ManagedCommunity(m; action = nothing)).trajectory.u == plain.u
        @test occursin("ManagedCommunity", sprint(show, m)) && occursin("ManagedTrajectory", sprint(show, run))
    end

    @testset "parameters reach the model" begin
        # merged into the environment of every step
        m = ManagedCommunity(c, n0, envs; parameters = (fecundity = 1.5,))
        n = [10.0, 20.0, 50.0]
        for _ in 1:6
            n = step(n, 1.5, 0.9)
        end
        @test project(m).trajectory.u[end, :] ≈ n rtol = 1e-12
        # through a build function, with a parameter that is not in the environment
        build(θ) = (community = pair_of(θ.survival_b), n0 = n0, envs = fill((fecundity = θ.fecundity,), 4), note = "ignored")
        mb = ManagedCommunity(build; parameters = (survival_b = 0.7, fecundity = 1.2))
        n = [10.0, 20.0, 50.0]
        for _ in 1:4
            n = step(n, 1.2, 0.7)
        end
        @test project(mb).trajectory.u[end, :] ≈ n rtol = 1e-12
        # an action replaces a parameter, by value or by a function of the value in force
        for (change, value) in ((0.5, 0.5), (x -> x / 2, 0.35))
            n = [10.0, 20.0, 50.0]
            for _ in 1:4
                n = step(n, 1.2, value)
            end
            acted = ManagedCommunity(mb; action = (; parameters = (; survival_b = change)))
            @test project(acted).trajectory.u[end, :] ≈ n rtol = 1e-12
        end
        @test mb.parameters == (survival_b = 0.7, fecundity = 1.2)              # the problem is not changed
        @test_throws ArgumentError project(ManagedCommunity(mb; action = (; parameters = (; absent = x -> x))))
        @test_throws ArgumentError project(ManagedCommunity(c, n0, fill(nothing, 2); parameters = (fecundity = 1.0,)))
    end

    @testset "actions on the state, against a loop by hand" begin
        function by_hand(; remove = nothing, take = nothing, add = nothing, steps = 1:6, fecundity = t -> 1.0)
            n = [10.0, 20.0, 50.0]
            removed, added = zeros(7, 3), zeros(7, 3)
            for t in 1:6
                n = step(n, fecundity(t), 0.9)
                t in steps || continue
                if remove !== nothing
                    i, f = remove
                    removed[t + 1, i] += f * n[i]
                    n[i] -= f * n[i]
                end
                if take !== nothing               # from both stages of `a`, in proportion
                    share = min(take / (n[1] + n[2]), 1.0)
                    removed[t + 1, 1:2] .+= share .* n[1:2]
                    n[1:2] .-= share .* n[1:2]
                end
                if add !== nothing
                    i, k = add
                    added[t + 1, i] += k
                    n[i] += k
                end
            end
            return n, removed, added
        end
        m = ManagedCommunity(c, n0, envs)
        function check(action; kw...)
            run = project(ManagedCommunity(m; action))
            n, removed, added = by_hand(; kw...)
            @test run.trajectory.u[end, :] ≈ n rtol = 1e-12
            @test run.removed ≈ removed rtol = 1e-12
            @test run.added ≈ added rtol = 1e-12
            return run
        end
        check((; remove = (:a, :adult) => 0.25); remove = (2, 0.25))
        check((; remove = :b => 0.1); remove = (3, 0.1))
        check((; remove = [(:a, :adult) => 0.25], steps = 3:4); remove = (2, 0.25), steps = 3:4)
        check((; take = :a => 7); take = 7)
        check((; take = (:a, [:juvenile, :adult]) => 1e6); take = 1e6)                # more than there are: all of them
        check((; add = (:a, :juvenile) => 5); add = (1, 5))
        check((; add = :b => [3.0], steps = [2, 5]); add = (3, 3.0), steps = [2, 5])
        check((; remove = (:a, :adult) => 0.25, take = :a => 7, add = (:a, :juvenile) => 5);
              remove = (2, 0.25), take = 7, add = (1, 5))                              # in that order
        run = check((; environment = (t, env) -> t <= 2 ? (fecundity = 0.0,) : env); fecundity = t -> t <= 2 ? 0.0 : 1.0)
        @test run.trajectory.u[2, 1] == 0
        # a custom outcome reads the record of the action
        harvest = ManagedCommunity(m; action = (; remove = (:a, :adult) => 0.25),
                                   outcomes = r -> (; harvested = sum(r.removed), adults = r.trajectory[:a][end, 2]))
        out = harvest.outcomes(project(harvest))
        n, removed, _ = by_hand(; remove = (2, 0.25))
        @test out.harvested ≈ sum(removed) && out.adults ≈ n[2]
    end

    @testset "invalid actions" begin
        m = ManagedCommunity(c, n0, envs)
        @test_throws ArgumentError ManagedCommunity(m; action = (; cull = :a => 0.1))
        @test_throws ArgumentError ManagedCommunity(m; action = :none)
        @test_throws ArgumentError ManagedCommunity(m; action = (; parameters = 3))
        @test_throws ArgumentError ManagedCommunity(c, n0, envs; parameters = [1.0])
        for action in ((; remove = :zebra => 0.1), (; remove = (:a, :pupa) => 0.1), (; remove = :a => 1.5),
                       (; remove = "a" => 0.1), (; take = :a => -1), (; add = :a => 5), (; add = :a => [1.0]),
                       (; add = :b => -2))
            @test_throws ArgumentError project(ManagedCommunity(m; action))
        end
        @test_throws ArgumentError project(ManagedCommunity(θ -> (community = c, n0 = n0)))
        @test_throws ArgumentError project(ManagedCommunity(θ -> (community = 1, n0 = n0, envs = envs)))
    end

    @testset "demographic stochasticity and the two streams" begin
        whole = [:a => [10.0, 20.0], :b => [50.0]]
        m = ManagedCommunity(c, whole, envs; demographic = true)
        @test_throws ArgumentError project(m)
        @test_throws ArgumentError project(ManagedCommunity(m; action = (; take = :a => 2.5)); rng = Xoshiro(1))
        @test_throws ArgumentError project(ManagedCommunity(m; action = (; add = :b => 0.5)); rng = Xoshiro(1))
        # by hand: the first draw of the generator seeds the stream of the action, and one more draw seeds each step
        function by_hand(seed, fraction, taken)
            rng = Xoshiro(seed)
            action_rng = Xoshiro(rand(rng, UInt64))
            n = [10.0, 20.0, 50.0]
            removed = 0.0
            for t in 1:6
                n = advance(c, n, envs[t]; rng = Xoshiro(rand(rng, UInt64)))
                x = StructuredPopulationCore.rand_binomial(action_rng, Int(n[2]), fraction)
                n[2] -= x
                y = min(taken, n[3])
                n[3] -= y
                removed += x + y
            end
            return n, removed
        end
        for seed in 1:5
            run = project(ManagedCommunity(m; action = (; remove = (:a, :adult) => 0.3, take = :b => 4)); rng = Xoshiro(seed))
            n, removed = by_hand(seed, 0.3, 4)
            @test run.trajectory.u[end, :] == n
            @test sum(run.removed) == removed
            @test all(isinteger, run.trajectory.u)
        end
        # an action that starts at step 4 leaves the first three steps as they were (common random numbers)
        late = project(ManagedCommunity(m; action = (; remove = :a => 0.5, steps = 4:6)); rng = Xoshiro(7))
        none = project(m; rng = Xoshiro(7))
        @test late.trajectory.u[1:4, :] == none.trajectory.u[1:4, :]
        @test late.trajectory.u[5, 1:2] != none.trajectory.u[5, 1:2]
        # an action on b leaves a, which is stepped before b, individual for individual as it was
        cull = project(ManagedCommunity(m; action = (; remove = :b => 0.5)); rng = Xoshiro(7))
        @test cull.trajectory.u[:, 1:2] == none.trajectory.u[:, 1:2] && cull.trajectory.u[end, 3] != none.trajectory.u[end, 3]
        # `take` from several states draws whole individuals without replacement and never more than there are
        for seed in 1:20
            run = project(ManagedCommunity(m; action = (; take = :a => 6)); rng = Xoshiro(seed))
            @test all(isinteger, run.removed) && all(>=(0), run.trajectory.u)
            @test all(sum(run.removed[t, 1:2]) == min(6, sum(run.trajectory.u[t, 1:2]) + sum(run.removed[t, 1:2])) for t in 2:7)
        end
        # a number taken from a single state leaves whole numbers, whatever the numbers (no k / total * n rounding)
        for total in 41:99, k in (7, 40)
            n, removed = [Float64(total)], zeros(1)
            MatrixCommunityModels._intervene!(n, removed, zeros(1), [], [([1], k)], [], Xoshiro(1))
            @test n == [total - k] && removed == [k]
        end
        quota = project(ManagedCommunity(m; action = (; take = (:a, :adult) => 3)); rng = Xoshiro(3))
        @test all(isinteger, quota.trajectory.u) && all(x -> x in (0, 3), quota.removed[2:end, 2])
        # the share taken from each state is in proportion to its numbers, on average (hypergeometric: mean 300 and
        # 100, standard deviation of the mean over 200 draws 0.6)
        taken = zeros(2)
        for seed in 1:200
            n, removed = [3000.0, 1000.0], zeros(2)
            MatrixCommunityModels._intervene!(n, removed, zeros(2), [], [([1, 2], 400)], [], Xoshiro(seed))
            @test sum(removed) == 400 && n .+ removed == [3000.0, 1000.0]
            taken .+= removed ./ 200
        end
        @test abs(taken[1] - 300) < 3
    end
end
