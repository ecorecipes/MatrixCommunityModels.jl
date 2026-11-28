# The DecisionAnalysis.jl adapter (ext/MatrixCommunityModelsDecisionAnalysisExt.jl), in the environment of
# this directory. The oracles are loops written by hand on plain matrices; they share no code with
# src/management.jl or the extension.

using Test
using MatrixCommunityModels
using DecisionAnalysis
using Random: Xoshiro
using Statistics: mean, median

const SPC = MatrixCommunityModels.StructuredPopulationCore

@testset "DecisionAnalysis adapter" begin
    @test Base.get_extension(MatrixCommunityModels,
                             :MatrixCommunityModelsDecisionAnalysisExt) !== nothing

    # A harvested species `a` (juveniles and adults) beside a second species `b`; the uncertain parameters are
    # the fecundity of `a` and the survival of `b`.
    function build(θ)
        a = births(matrix_net([:juvenile, :adult], :survival => [0.0 0.0; 0.5 0.8],
                              :breeding => [0.0 θ.fecundity; 0.0 0.0]),
                   :breeding)
        b = matrix_net([:only], :survival => fill(θ.survival_b, 1, 1))
        return (community=community(species(:a, a), species(:b, b)),
                n0=[:a => [10.0, 20.0], :b => [50.0]],
                envs=fill(nothing, 8))
    end
    outs(run) = (; harvested=sum(run.removed), adults=run.trajectory[:a][end, 2],
                 b=run.trajectory[:b][end, 1],
                 adults_at_3=run.trajectory[:a][4, 2])
    m = ManagedCommunity(build; parameters=(fecundity=1.0, survival_b=0.9), outcomes=outs)
    actions = [:none => (;), :adults_20 => (; remove=(:a, :adult) => 0.2),
               :adults_40 => (; remove=(:a, :adult) => 0.4),
               :late_40 => (; remove=(:a, :adult) => 0.4, steps=4:8),
               :cull_b => (; take=:b => 5)]
    draws = (; fecundity=[0.2, 0.6, 1.4, 1.8, 0.2, 1.0],
             survival_b=[0.5, 0.9, 0.7, 0.8, 0.95, 0.6])

    # the model and the actions, by hand
    function by_hand(fecundity, survival_b, label)
        n = [10.0, 20.0, 50.0]
        harvested, adults_at_3 = 0.0, 0.0
        for t in 1:8
            n = [fecundity * n[2], 0.5 * n[1] + 0.8 * n[2], survival_b * n[3]]
            fraction = label === :adults_20 ? 0.2 :
                       label === :adults_40 || (label === :late_40 && t >= 4) ? 0.4 : 0.0
            harvested += fraction * n[2]
            n[2] -= fraction * n[2]
            if label === :cull_b
                x = min(5.0, n[3])
                harvested += x
                n[3] -= x
            end
            t == 3 && (adults_at_3 = n[2])
        end
        return (; harvested, adults=n[2], b=n[3], adults_at_3)
    end

    @testset "the four verbs" begin
        θ = (fecundity=1.4, survival_b=0.7)
        p = with_parameters(m, θ)
        @test p isa ManagedCommunity && p.parameters == θ
        @test m.parameters == (fecundity=1.0, survival_b=0.9)                    # not mutated
        @test with_parameters(m, (; fecundity=2.0)).parameters ==
              (fecundity=2.0, survival_b=0.9)
        q = apply_action(p, (; remove=(:a, :adult) => 0.2))
        @test q.action == (; remove=(:a, :adult) => 0.2) && isempty(p.action)
        @test isempty(apply_action(q, nothing).action) &&
              isempty(apply_action(q, (;)).action)
        @test_throws ArgumentError apply_action(p, (; harvest=0.2))
        run = simulate(q; rng=Xoshiro(1))
        @test run isa ManagedTrajectory
        o = outcome(run)
        expected = by_hand(1.4, 0.7, :adults_20)
        @test keys(o) == keys(expected)
        @test all(isapprox(o[k], expected[k]; rtol=1e-12) for k in keys(o))
        # the default outcome: the final abundance of each species
        d = outcome(simulate(ManagedCommunity(build; parameters=θ); rng=Xoshiro(1)))
        none = by_hand(1.4, 0.7, :none)
        @test keys(d) == (:a, :b) && d.b ≈ none.b
    end

    table = outcomes(m, actions, draws; replicates=3, seed=11)

    @testset "the driver against the loop by hand" begin
        @test action_labels(table) == first.(actions)
        @test outcome_names(table) == [:harvested, :adults, :b, :adults_at_3]
        @test run_counts(table) == (runs=5 * 6 * 3, ok=5 * 6 * 3, failed=0, absent=0)
        @test parameter_draws(table) == draws
        for (a, (label, _)) in enumerate(actions), d in 1:6
            expected = by_hand(draws.fecundity[d], draws.survival_b[d], label)
            for (k, name) in enumerate(outcome_names(table)), r in 1:3
                @test table.values[a, d, r, k] ≈ expected[name] rtol = 1e-12
            end
        end
        # a SampleBelief with unequal weights is carried through
        belief = SampleBelief(draws; weights=[1, 2, 3, 1, 2, 1])
        @test draw_weights(outcomes(m, actions, belief; seed=1)) == belief.weights
    end

    @testset "a deterministic model has no Monte Carlo error" begin
        payoff = payoff_table(table, IndividualsUtility(:harvested); baseline=:none)
        # zero up to the rounding of a mean of equal numbers
        @test payoff.standard_errors !== nothing && all(<(1e-12), payoff.standard_errors)
        @test all(==(3), payoff.replicates) && !single_replicate(payoff)
        @test all(table.values[:, :, 1, :] .== table.values[:, :, 2, :] .==
                  table.values[:, :, 3, :])
        @test outcomes(m, actions, draws; replicates=3, seed=99).values == table.values     # the seed is not used
        @test single_replicate(payoff_table(outcomes(m, actions, draws; seed=1),
                                            IndividualsUtility(:harvested)))
    end

    @testset "the full pipeline" begin
        utility = IndividualsUtility(:harvested)
        payoff = payoff_table(table, utility)
        Q = [by_hand(draws.fecundity[d], draws.survival_b[d], label).harvested
             for (label, _) in actions, d in 1:6]
        @test Float64.(payoff.Q) ≈ Q rtol = 1e-12
        problem = DecisionProblem(table, utility; baseline=:none)
        means = vec(sum(Q; dims=2)) ./ 6
        result = decide(problem, MaxExpectedUtility())
        @test first_choice(result) == first.(actions)[argmax(means)]
        @test first_choice(decide(problem, Maximin())) ==
              first.(actions)[argmax(vec(minimum(Q; dims=2)))]
        # EVPI by hand: the mean over draws of the best payoff, less the best mean payoff
        by_hand_evpi = sum(maximum(Q; dims=1)) / 6 - maximum(means)
        @test Float64(evpi(problem)) ≈ by_hand_evpi rtol = 1e-10
        @test by_hand_evpi > 1                                    # the best action depends on the draw
        gain = Float64.(gain_over_baseline(problem).Q)
        @test gain ≈ Q .- Q[1:1, :] rtol = 1e-10
        # EVPPI by a partition of the draws, by hand: low and high fecundity
        cells = [[1, 2, 5], [3, 4, 6]]
        partition = [d in cells[1] ? :low : :high for d in 1:6]
        by_hand_evppi = sum(maximum(sum(Q[:, cell]; dims=2)) for cell in cells) / 6 -
                        maximum(means)
        @test Float64(evppi(problem, partition)) ≈ by_hand_evppi rtol = 1e-10
        # and by the sample-based method, from the parameter draws that the table carries
        grouped = evppi(problem, :fecundity, Grouping(; bins=2))
        @test grouped isa VOIEstimate
        @test grouped.value ≈ by_hand_evppi rtol = 1e-10
        @test evppi(problem, :survival_b, Grouping(; bins=2)).value <=
              Float64(evpi(problem)) + 1e-12
    end

    @testset "common random numbers" begin
        whole = ManagedCommunity(build; parameters=(fecundity=1.0, survival_b=0.9),
                                 outcomes=outs, demographic=true)
        stochastic = outcomes(whole, actions, draws; replicates=4, seed=11)
        k3, kb = 4, 3
        none, late, cull = 1, 4, 5
        # before step 4 the late harvest is no harvest: the same seed gives the same draws to both actions
        @test stochastic.values[late, :, :, k3] == stochastic.values[none, :, :, k3]
        @test stochastic.values[late, :, :, 2] != stochastic.values[none, :, :, 2]
        # culling b does not touch a, which is stepped before it: a's run is the same individual for individual
        @test stochastic.values[cull, :, :, 2] == stochastic.values[none, :, :, 2]
        @test stochastic.values[cull, :, :, kb] != stochastic.values[none, :, :, kb]
        # replicates differ from one another, and another seed gives other runs
        @test stochastic.values[none, :, 1, 2] != stochastic.values[none, :, 2, 2]
        @test outcomes(whole, actions, draws; replicates=4, seed=12).values !=
              stochastic.values
        @test outcomes(whole, actions, draws; replicates=4, seed=11).values ==
              stochastic.values
        @test all(isinteger, stochastic.values)
        # one cell by hand: the verbs with a known generator (the driver's seed for a cell is not public)
        d = 3
        θ = (fecundity=draws.fecundity[d], survival_b=draws.survival_b[d])
        c = build(θ).community
        function by_hand_stochastic(rng)
            action_rng = Xoshiro(rand(rng, UInt64))
            n = [10.0, 20.0, 50.0]
            harvested = 0.0
            for t in 1:8
                n = advance(c, n, nothing; rng=Xoshiro(rand(rng, UInt64)))
                x = SPC.rand_binomial(action_rng, Int(n[2]), 0.4)
                n[2] -= x
                harvested += x
            end
            return harvested, n[2]
        end
        direct(rng) = outcome(simulate(apply_action(with_parameters(whole, θ),
                                                    (; remove=(:a, :adult) => 0.4)); rng))
        for seed in 1:10
            o = direct(Xoshiro(seed))
            harvested, adults = by_hand_stochastic(Xoshiro(seed))
            @test o.harvested == harvested && o.adults == adults
        end
        # the paired differences are less noisy than independent runs would be
        utility = IndividualsUtility(:adults)
        @test crn_correlation(stochastic, utility) > 0.2
        payoff = payoff_table(stochastic, utility)
        @test payoff.standard_errors !== nothing && any(>(0), payoff.standard_errors)
    end

    @testset "runs that fail are kept and counted" begin
        bad = [:none => (;), :zebra => (; remove=:zebra => 0.1)]
        t = outcomes(m, bad, draws; seed=1)
        @test run_counts(t).failed == 6
        @test all(f -> f.action === :zebra && f.status == "error:ArgumentError",
                  failed_runs(t))
    end

    @testset "Verde River fish, Tonkin flow regimes, ResilienceCore & Track 3/4" begin
        # Live Verde River 7-species fish community under Tonkin et al. (2021) flow regimes
        # with uncertain carrying capacity K and drought frequency, plus ResilienceCore metrics.
        const_nat = tonkin_scenario(:natural; years=30)
        const_fsh = tonkin_scenario(:fish; years=30)
        const_plt = tonkin_scenario(:plants; years=30)
        const_inv = tonkin_scenario(:invertebrates; years=30)

        function build_verde(θ)
            vm = verde_fish_model(; K=47660.0 * θ.K_scale)
            envs = map(const_nat.fish) do e
                if e.nonevent == 1 && θ.extra_drought > 0.5
                    (SP_highflood=0, SU_highflood=0, medflood=0, drought=1, nonevent=0)
                else
                    e
                end
            end
            return (community=vm.community, n0=vm.n0, envs=envs)
        end

        function verde_outcomes(run)
            tr = run.trajectory
            native = mean(sum(tr[sp][11:end, 2] .+ tr[sp][11:end, 3]
                              for sp in (:CACL, :GIRO, :CAIN)))
            nonnative = mean(sum(tr[sp][11:end, 2] .+ tr[sp][11:end, 3]
                                 for sp in (:LECY, :MIDO, :CYLU, :AMNA)))
            share = native / max(native + nonnative, 1e-12)
            # ResilienceCore local linearisation of the final community state under a non-event year
            J = linearise(tr.community, tr.u[end, :],
                          (SP_highflood=0, SU_highflood=0, medflood=0, drought=0,
                           nonevent=1);
                          check=false)
            react = MatrixCommunityModels.ResilienceCore.reactivity(J)
            return (; native, nonnative, native_share=share, reactivity=react)
        end

        mv = ManagedCommunity(build_verde; parameters=(K_scale=1.0, extra_drought=0.0),
                              outcomes=verde_outcomes)
        regime_actions = [:natural => (;),
                          :fish_flow => (; environment=(t, _) -> const_fsh.fish[t]),
                          :plant_flow => (; environment=(t, _) -> const_plt.fish[t]),
                          :invert_flow => (; environment=(t, _) -> const_inv.fish[t]),
                          :natural_cull => (; remove=[:MIDO => 0.3, :LECY => 0.3])]
        rng = Xoshiro(310)
        n_draws = 16
        vdraws = (; K_scale=0.7 .+ 0.6 .* rand(rng, n_draws),
                  extra_drought=rand(rng, n_draws))

        vtable = outcomes(mv, regime_actions, vdraws; replicates=1, seed=310)
        @test run_counts(vtable) == (runs=5 * n_draws, ok=5 * n_draws, failed=0, absent=0)
        @test outcome_names(vtable) == [:native, :nonnative, :native_share, :reactivity]

        # Direct oracle check on one cell against plain project(verde_fish_model(...))
        θ1 = (K_scale=vdraws.K_scale[1], extra_drought=vdraws.extra_drought[1])
        vm1 = verde_fish_model(; K=47660.0 * θ1.K_scale)
        sol_fsh = project(vm1.community, vm1.n0, const_fsh.fish)
        expected_native = mean(sum(sol_fsh[sp][11:end, 2] .+ sol_fsh[sp][11:end, 3]
                                   for sp in (:CACL, :GIRO, :CAIN)))
        @test vtable.values[2, 1, 1, 1] ≈ expected_native rtol = 1e-12
        @test all(isfinite, vtable.values[:, :, 1, 4])

        # Multi-objective Pareto front & Track 3 rules / CE-VOI / PRIM
        p_native = DecisionProblem(vtable, IndividualsUtility(:native); baseline=:natural)
        p_share = DecisionProblem(vtable, IndividualsUtility(:native_share);
                                  baseline=:natural)
        @test first_choice(decide(p_native, MaxExpectedUtility())) === :fish_flow

        pf = decide((; native=p_native.payoff, native_share=p_share.payoff),
                    p_native.belief, ParetoFront())
        @test :fish_flow in pf.chosen

        # Stochastic dominance (second-order) & certainty-equivalent EVPI
        sd = decide(p_native, StochasticDominance(; order=2))
        @test :fish_flow in sd.chosen
        ce_val = ce_evpi(p_native, ExponentialUtility(:native; risk_aversion=1e-4))
        @test isfinite(Float64(ce_val.value)) && Float64(ce_val.value) >= -1e-9

        # PRIM scenario discovery on native fish shortfall under natural flows
        nat_idx = findfirst(==(:natural), action_labels(vtable))
        nat_vals = vec(vtable.values[nat_idx, :, 1, 1])
        prim = prim_boxes(p_native.belief, nat_vals .< median(nat_vals); peel_alpha=0.25,
                          min_mass=0.2)
        @test !isempty(prim.frontier) && maximum(b -> b.density, prim.frontier) >= 0.5
    end
end
