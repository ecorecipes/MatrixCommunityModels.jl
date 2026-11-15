# Sea fans and aspergillosis: an integral projection model stratified by infection status (Bruno et al. 2011)
using LinearAlgebra
using Random

@testset "sea fans and aspergillosis" begin
    P, phases = SEA_FAN_PARAMETERS, SEA_FAN_PHASES
    cells = 50
    fan = sea_fan_stage(; cells)
    z, h = meshpoints(fan), step_size(fan)
    # the positions of a stratum's cells, found by name
    where(sp, stratum) = Int.(indexin(stratum_stage.(stage_names(fan), stratum), stage_names(sp, :sea_fan)))
    # growth rate, stable state and reproductive values of a species held in one phase, and the size of the second
    # eigenvalue relative to the first
    function steady(sp, phase)
        A = community_matrix(sp, zeros(length(sp)), phase)
        values, vectors = eigen(A)
        order = sortperm(abs.(values); rev = true)
        w = abs.(real.(vectors[:, order[1]]))
        left = eigen(Matrix(A'))
        v = abs.(real.(left.vectors[:, argmax(abs.(left.values))]))
        rate = abs(values[order[1]])
        net = frozen_net(sp, :sea_fan, w, phase)
        # share of the elasticity of the growth rate that belongs to each named transition
        shares = Dict(t => 100 * dot(v, transition_matrix(net, t) * w) / (rate * dot(v, w)) for t in transition_names(net))
        return (rate = rate, state = w ./ sum(w), second = abs(values[order[2]]) / rate, shares = shares)
    end
    # the paper's fitted functions, written out (p. 130, Eqs. 5 to 7 and B.2)
    logistic(u) = 1 / (1 + exp(-u))
    normal(y, m, s) = exp(-0.5 * ((y - m) / s)^2) / (s * sqrt(2pi))
    height(x) = 1.3 * x^1.5
    storms(x) = 0.093 * logistic(2.3 * (x - 8.7))
    survival(x) = 1 - logistic(-0.55 - 0.62x) - storms(x)
    infected_survival(x) = 1 - logistic(1.3 - 0.67x) - storms(x)
    risk(x) = x <= 5.6 ? 0.011 : 0.053
    growth(y, x) = normal(y, 0.682 + 0.945x, 0.53)
    change(shift) = (y, x) -> 0.95 * normal(y, x + shift, 0.57) + 0.05 * normal(y, x - 3.2, 0.52)
    recruits(x) = height(x) >= 20 ? P.fecundity * x^3 : 0.0
    recruit_size(y) = normal(y, P.recruit_size[1], P.recruit_size[2])
    two_by_two = Dict(:escaping => [1 0; 0 0], :encountering => [0 0; 1 0], :recovering => [0 1; 0 0],
                      :lingering => [0 0; 0 1], :recruiting => [1 0; 0 0])
    # the model built by hand from those functions; `pressure` multiplies the risk of infection
    function by_hand(pressure)
        life = kernel_net(fan, :escaping => growth, :encountering => growth, :recovering => change(0.27),
                          :lingering => change(-0.33), :recruiting => (y, x) -> recruits(x) * recruit_size(y))
        fates = respond(life,
            :escaping => ctx -> by_trait(fan, x -> survival(x) * (1 - pressure(ctx) * risk(x))),
            :encountering => ctx -> by_trait(fan, x -> survival(x) * pressure(ctx) * risk(x)),
            :recovering => ctx -> by_trait(fan, x -> infected_survival(x) * 0.38),
            :lingering => ctx -> by_trait(fan, x -> infected_survival(x) * 0.62))
        return stratify(births(fates, :recruiting), [:healthy, :infected]; coupling = two_by_two)
    end

    @testset "strata" begin
        reef = sea_fan(; cells)
        net = system(reef)[:sea_fan]
        @test stage_names(reef, :sea_fan) == vcat(stratum_stage.(stage_names(fan), :healthy), stratum_stage.(stage_names(fan), :infected))
        @test Set(transition_names(net)) == Set([:escaping_healthy, :encountering_healthy_to_infected,
                                                 :recovering_infected_to_healthy, :lingering_infected, :recruiting_healthy])
        @test net.fecundity == [:recruiting_healthy]
        @test pool_names(reef) == [:tissue, :infectious]
        healthy, infected = where(reef, :healthy), where(reef, :infected)
        @test healthy == 1:cells && infected == (cells + 1):(2cells)

        # the four blocks of Eq. 1 of the paper, written out
        block(f) = [h * f(y, x) for y in z, x in z]
        expected = zeros(2cells, 2cells)
        expected[healthy, healthy] = block((y, x) -> survival(x) * (1 - risk(x)) * growth(y, x) + recruits(x) * recruit_size(y))
        expected[infected, healthy] = block((y, x) -> survival(x) * risk(x) * growth(y, x))
        expected[healthy, infected] = block((y, x) -> infected_survival(x) * 0.38 * change(0.27)(y, x))
        expected[infected, infected] = block((y, x) -> infected_survival(x) * 0.62 * change(-0.33)(y, x))
        A = community_matrix(reef, zeros(2cells), phases.endemic)
        @test maximum(abs, A - expected) < 1e-14
        @test community_matrix(reef, zeros(2cells), nothing) == A              # without a phase, the endemic one
        # the same matrix from the net built by hand, and from two continuous stages with no strata
        hand = by_hand(ctx -> 1.0)
        @test maximum(abs, A - to_matrix(hand, context(hand, zeros(2cells), NamedTuple(), nothing))) < 1e-14
        well, sick = ContinuousStage(:well, 0.0, 20.0, cells), ContinuousStage(:sick, 0.0, 20.0, cells)
        blocks = kernel_net([well, sick],
            :escaping => [(:well => :well) => (y, x) -> survival(x) * (1 - risk(x)) * growth(y, x)],
            :encountering => [(:well => :sick) => (y, x) -> survival(x) * risk(x) * growth(y, x)],
            :recovering => [(:sick => :well) => (y, x) -> infected_survival(x) * 0.38 * change(0.27)(y, x)],
            :lingering => [(:sick => :sick) => (y, x) -> infected_survival(x) * 0.62 * change(-0.33)(y, x)],
            :recruiting => [(:well => :well) => (y, x) -> recruits(x) * recruit_size(y)])
        @test maximum(abs, A - to_matrix(blocks)) < 1e-14

        # without infection the healthy stratum is the life cycle of a fan that is never infected, and nothing
        # enters the infected stratum
        never = kernel_net(fan, :growing => (y, x) -> survival(x) * growth(y, x), :recruiting => (y, x) -> recruits(x) * recruit_size(y))
        clean = community_matrix(reef, zeros(2cells), merge(phases.endemic, (small_fan_risk = 0.0, large_fan_risk = 0.0)))
        @test maximum(abs, clean[healthy, healthy] - to_matrix(never)) < 1e-14
        @test all(iszero, clean[infected, healthy])
        @test clean[:, infected] == A[:, infected]
        # resistant fans alone are that life cycle too
        alone = sea_fan(; cells, resistant = :only)
        @test stage_names(alone, :sea_fan) == stratum_stage.(stage_names(fan), :resistant)
        @test maximum(abs, community_matrix(alone, zeros(cells), phases.endemic) - to_matrix(never)) < 1e-14
        @test pool_names(alone) == [:tissue]
        # and stay so if infection is passed among fans: they have no infected tissue to pass
        @test community_matrix(sea_fan(; cells, resistant = :only, transmission = 5.0), ones(cells), phases.endemic) ==
              community_matrix(alone, ones(cells), phases.endemic)

        # the state read by name: a stratum, a stage in a stratum, the cells of the stage whatever the stratum,
        # and its cells in one stratum
        n = collect(1.0:(2cells))
        ctx = context(net, n, (tissue = 0.0, infectious = 0.0), phases.endemic)
        @test ctx.n[:healthy] == sum(1:cells) && ctx.n[:infected] == sum((cells + 1):(2cells))
        @test ctx.n[:fan_infected] == ctx.n[:infected] && ctx.n[:fan] == sum(n)
        @test ctx.n[fan] == n[healthy] .+ n[infected]
        @test ctx.n[fan, :healthy] == n[healthy] && ctx.n[fan, :infected] == n[infected]
        @test ctx.n[:fan7_infected] == n[cells + 7] && ctx.n[:fan7] == n[7] + n[cells + 7]
        pools = pool_totals(reef, n, phases.endemic)
        @test pools.tissue ≈ sum(z .^ 3 .* ctx.n[fan]) && pools.infectious ≈ sum(z .^ 3 .* ctx.n[fan, :infected])
        # survival and growth are probabilities, and only healthy fans recruit
        U, F = survival_fecundity_matrices(net, context(net, n, (tissue = 0.0, infectious = 0.0), phases.peak))
        @test maximum(sum(U; dims = 1)) < 1 && all(>=(0), U)
        @test all(iszero, F[infected, :]) && all(iszero, F[:, infected]) && sum(F) > 0
        @test_throws ArgumentError sea_fan(; cells, resistant = :some)
        @test_throws ArgumentError sea_fan(; cells, transmission = -1.0)
    end

    @testset "values in Bruno et al. (2011)" begin
        reef = sea_fan(; cells)
        infected = where(reef, :infected)
        # 22% of recruits die in their first year (p. 128): this sets the mean size of a recruit
        @test sum(h * recruit_size(x) * logistic(-0.55 - 0.62x) for x in z) ≈ 0.22 atol = 2e-3
        # a fan of size 5.6 is about 17.5 cm tall; sizes 7 and 10 are heights of about 25 and 40 cm (Appendix B)
        @test height(5.6) ≈ 17.5 atol = 0.5
        @test height(7.0) ≈ 25 atol = 1.5
        @test height(10.0) ≈ 40 atol = 1.5

        pre, peak, endemic = steady(reef, phases.pre_epizootic), steady(reef, phases.peak), steady(reef, phases.endemic)
        # growth rates read from Figs. E2 and E3: 0.95 before the epizootic, 0.851 at its peak, 1.009 afterwards, and
        # 0.935 at the peak if recruitment were that of 2003
        @test pre.rate ≈ 0.95 atol = 4e-3
        @test peak.rate ≈ 0.851 atol = 4e-3
        @test endemic.rate ≈ 1.009 atol = 5e-3
        @test steady(reef, merge(phases.peak, (recruitment = 1.0,))).rate ≈ 0.935 atol = 4e-3
        # without infection 1.01, and the second eigenvalue is 0.86 of the first (p. 136); 0.84 in the endemic phase (p. 132)
        none = steady(reef, merge(phases.pre_epizootic, (recruitment = 1.0,)))
        @test none.rate ≈ 1.01 atol = 6e-3
        @test none.second ≈ 0.86 atol = 0.015
        @test endemic.second ≈ 0.84 atol = 0.01
        # prevalence in the endemic phase: 7.05% (Fig. E3D)
        @test 100 * sum(endemic.state[infected]) ≈ 7.05 atol = 0.6
        @test sum(pre.state[infected]) < 1e-12

        # Table E1: the elasticity of the growth rate, by kind of transition and by infection status
        @test sum(values(endemic.shares)) ≈ 100
        @test endemic.shares[:recruiting_healthy] ≈ 3.5 atol = 0.3
        @test endemic.shares[:escaping_healthy] + endemic.shares[:recruiting_healthy] ≈ 91.1 atol = 0.6
        @test endemic.shares[:encountering_healthy_to_infected] ≈ 2.8 atol = 0.2
        @test endemic.shares[:recovering_infected_to_healthy] ≈ endemic.shares[:encountering_healthy_to_infected] rtol = 1e-6
        @test endemic.shares[:lingering_infected] ≈ 3.3 atol = 0.2
        @test pre.shares[:recruiting_healthy] ≈ 2.1 atol = 0.25
        @test pre.shares[:escaping_healthy] + pre.shares[:recruiting_healthy] ≈ 100
        @test peak.shares[:recruiting_healthy] ≈ 3.3 atol = 0.2
        @test peak.shares[:escaping_healthy] + peak.shares[:recruiting_healthy] ≈ 98.6 atol = 0.2
        # the table prints 0.05, 0.05 and 0.03 for the other three, which leaves the column summing to 98.7
        @test peak.shares[:encountering_healthy_to_infected] ≈ 0.5 atol = 0.1
        @test peak.shares[:lingering_infected] ≈ 0.3 atol = 0.1

        # prevalence by size at the peak (Fig. 5D): 0.36, 0.55 and 0.62 at heights of 50, 80 and 100 cm
        function infected_by_size(steady_state, phase)
            n = context(system(reef)[:sea_fan], steady_state.state, (tissue = 0.0, infectious = 0.0), phase).n
            return n[fan, :infected] ./ n[fan]
        end
        at(curve, height_cm) = curve[argmin(abs.(height.(z) .- height_cm))]
        at_peak = infected_by_size(peak, phases.peak)
        @test at(at_peak, 50) ≈ 0.36 atol = 0.05
        @test at(at_peak, 80) ≈ 0.55 atol = 0.05
        @test at(at_peak, 100) ≈ 0.62 atol = 0.06
        @test issorted(at_peak[20 .< height.(z) .< 100])                # prevalence rises with size among large fans
        # and under the risks of the endemic phase (Fig. 5B, the asymptotic curve): 0.14, 0.105, 0.095 and 0.09 at 22, 40, 60
        # and 75 cm. The curve is matched with the recruitment of the epizootic years, which the paper does not state
        few_recruits = merge(phases.endemic, (recruitment = 1 / 11,))
        at_few = infected_by_size(steady(reef, few_recruits), few_recruits)
        @test at(at_few, 22) ≈ 0.14 atol = 0.015
        @test at(at_few, 40) ≈ 0.105 atol = 0.01
        @test at(at_few, 60) ≈ 0.095 atol = 0.01
        @test at(at_few, 75) ≈ 0.09 atol = 0.01
        @test at(infected_by_size(endemic, phases.endemic), 22) < 0.115  # with the recruits of 2003 the hump is lower

        # the recruitment constant as printed, 2 x 10^-4: growth rates lower by up to 0.01
        printed = sea_fan(; cells, parameters = merge(P, (fecundity = 2e-4,)))
        @test 0.004 < endemic.rate - steady(printed, phases.endemic).rate < 0.008
        # the mesh does not matter
        fine = steady(sea_fan(), phases.endemic)
        @test fine.rate ≈ endemic.rate rtol = 1e-3
        @test sum(fine.state[101:200]) ≈ sum(endemic.state[infected]) rtol = 5e-3
    end

    @testset "resistant fans" begin
        reef, three = sea_fan(; cells), sea_fan(; cells, resistant = true)
        net = system(three)[:sea_fan]
        @test stage_names(three, :sea_fan)[[1, cells + 1, 2cells + 1]] == [:fan1_healthy, :fan1_infected, :fan1_resistant]
        @test Set(transition_names(net)) == Set([:escaping_healthy, :encountering_healthy_to_infected, :recovering_infected_to_healthy,
            :lingering_infected, :recruiting_healthy, :escaping_resistant, :encountering_resistant, :recruiting_resistant])
        @test Set(net.fecundity) == Set([:recruiting_healthy, :recruiting_resistant])
        # a resistant fan that meets the fungus stays as it is: the stratum is closed, and is a fan that is never infected
        A = community_matrix(three, zeros(3cells), phases.peak)
        healthy, infected, resistant = where(three, :healthy), where(three, :infected), where(three, :resistant)
        @test all(iszero, A[vcat(healthy, infected), resistant]) && all(iszero, A[resistant, vcat(healthy, infected)])
        @test A[vcat(healthy, infected), vcat(healthy, infected)] == community_matrix(reef, zeros(2cells), phases.peak)
        @test A[resistant, resistant] ≈ community_matrix(sea_fan(; cells, resistant = :only), zeros(cells), phases.peak)

        # Fig. 8A: 15% of fans resistant in 1994; the peak phase from then on, with the risk of the others raised by
        # 15% and recruitment rising 11-fold between 1995 and 2003 (Appendix C)
        before = steady(reef, phases.pre_epizootic).state[where(reef, :healthy)]
        n0 = [0.85 .* before; zeros(cells); 0.15 .* before]
        raised = merge(phases.peak, (small_fan_risk = 1.15 * phases.peak.small_fan_risk,
                                     large_fan_risk = 1.15 * phases.peak.large_fan_risk, risk_per_cm = 1.15 * phases.peak.risk_per_cm))
        envs = [merge(raised, (recruitment = min(1.0, exp(-0.3 * (2003 - year))),)) for year in 1994:2007]
        sol = project(three, n0, envs)
        fans = abundance(sol, :sea_fan)
        infected_share = 100 .* stage_total(sol, :sea_fan, :infected) ./ fans
        resistant_share = 100 .* stage_total(sol, :sea_fan, :resistant) ./ fans
        published = [18, 20, 19, 17, 14.5, 12.5, 11, 9.5, 8.5, 7.5, 6.5, 6, 5.5, 5]      # % infected, 1995 to 2008, read from the figure
        @test maximum(abs, infected_share[2:end] .- published) < 1.5
        @test argmax(infected_share) == 3                                              # the peak is in 1996
        # in the first year nothing tells the resistant fans from the others: those newly infected are still alive
        @test resistant_share[1] ≈ 15 && resistant_share[2] ≈ 15
        @test 50 < resistant_share[end] < 60                                           # about 53% in 2008
        @test issorted(resistant_share[2:end])

        # resistant fans as a species of their own give the same numbers
        both = community(reef, sea_fan(; cells, resistant = :only, name = :resistant_fans))
        @test species_names(both) == [:sea_fan, :resistant_fans] && pool_names(both) == [:tissue, :infectious]
        @test project(both, n0, envs).u ≈ sol.u rtol = 1e-12
    end

    @testset "infection passed among fans" begin
        reef = sea_fan(; cells)
        healthy, infected = where(reef, :healthy), where(reef, :infected)
        endemic = steady(reef, phases.endemic)
        pools = pool_totals(reef, endemic.state, phases.endemic)
        share = pools.infectious / pools.tissue                       # of tissue that is infected, in the published model
        @test 0.08 < share < 0.10
        # the same dependence in a response that reads the reef's own state, written before the strata exist: the tissue
        # in each cell whatever its stratum (`ctx.n[fan]`), and that of its infected fans (`ctx.n[fan, :infected]`)
        area = z .^ 3
        own(q) = ctx -> 1 - q + q / share * sum(area .* ctx.n[fan, :infected]) / sum(area .* ctx.n[fan])
        for q in (0.5, 1.0)
            passed = sea_fan(; cells, outside = 1 - q, transmission = q / share)
            @test pool_names(passed) == [:tissue, :infectious]
            # at the published stable state the risk is the published one, so that state persists
            @test advance(passed, endemic.state, phases.endemic) ≈ advance(reef, endemic.state, phases.endemic)
            hand = species(:sea_fan, by_hand(own(q)))
            start = endemic.state .* (1 .+ 0.5 .* sin.(1:(2cells)))   # away from the stable state
            @test project(hand, start, fill(nothing, 30)).u ≈ project(passed, start, fill(phases.endemic, 30)).u rtol = 1e-10
        end
        # with all infection passed among fans there is a threshold: at nine tenths of the published risk the infection
        # dies out, where infection from outside only falls in proportion
        scaled(m) = merge(phases.endemic, (small_fan_risk = m * phases.endemic.small_fan_risk, large_fan_risk = m * phases.endemic.large_fan_risk))
        prevalence(sp, m; years = 400) = (n = project(sp, endemic.state, fill(scaled(m), years)).u[end, :]; sum(n[infected]) / sum(n))
        among = sea_fan(; cells, outside = 0.0, transmission = 1 / share)
        @test prevalence(among, 1.0) ≈ sum(endemic.state[infected]) rtol = 1e-6
        @test prevalence(among, 0.9) < 1e-3
        @test prevalence(among, 1.1) > 1.8 * sum(endemic.state[infected])
        @test 0.85 < prevalence(reef, 0.9) / sum(endemic.state[infected]) < 0.95
        # resistant fans, a species of their own, dilute infection passed among fans and leave infection from outside alone
        resistant = sea_fan(; cells, resistant = :only, name = :resistant_fans)
        n0 = [endemic.state; (endemic.state[healthy] .+ endemic.state[infected]) ./ 3]      # a quarter of all fans
        later(sp) = (sol = project(community(sp, resistant), n0, fill(phases.endemic, 100));
                     stage_total(sol, :sea_fan, :infected)[end] / abundance(sol, :sea_fan)[end])
        @test later(reef) ≈ sum(endemic.state[infected]) rtol = 1e-6
        @test later(among) < 1e-3
        # the same as a third stratum
        three = sea_fan(; cells, resistant = true, outside = 0.0, transmission = 1 / share)
        @test project(three, n0, fill(phases.endemic, 20)).u ≈ project(community(among, resistant), n0, fill(phases.endemic, 20)).u rtol = 1e-10
    end

    @testset "whole fans" begin
        reef = sea_fan(; cells)
        start = round.(2000 .* steady(reef, phases.endemic).state)
        random = project(reef, start, fill(phases.endemic, 10); rng = MersenneTwister(16))
        @test all(isinteger, random.u) && all(>=(0), random.u)
        @test stage_total(random, :sea_fan, :infected)[end] > 0
    end
end
