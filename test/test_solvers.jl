# A community with a continuous stage is solved by IntegralProjectionModels.jl, any other by
# MatrixProjectionModels.jl; both solvers give the same projection, deterministic or with whole individuals.
using Random

@testset "solvers" begin
    IPM = MatrixCommunityModels.IntegralProjectionModels
    CPD = MatrixCommunityModels.CategoricalPopulationDynamics
    meshes = CPD.meshes
    filled(c, value) = reduce(vcat, [fill(value, length(stage_names(c, sp))) for sp in species_names(c)])
    function same(c, n0, envs; rng_seed = nothing)
        a = project(c, n0, envs; solver = :ipm)
        b = project(c, n0, envs; solver = :mpm)
        @test a.u == b.u
        @test isequal(a.pools, b.pools)
        @test project(c, n0, envs).u == a.u                                     # :auto
        if rng_seed !== nothing
            whole = round.(n0)
            x = project(c, whole, envs; rng = Xoshiro(rng_seed), solver = :ipm)
            y = project(c, whole, envs; rng = Xoshiro(rng_seed), solver = :mpm)
            @test x.u == y.u && all(isinteger, x.u)
            @test isequal(x.pools, y.pools)
            @test project(c, whole, envs; rng = Xoshiro(rng_seed)).u == x.u
        end
        return a
    end

    @testset "routing" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 20)
        life = kernel_net(plant, :growing => (y, x) -> 0.7 * exp(-(y - x - 0.5)^2 / 2) / sqrt(2pi),
                          :seeding => (y, x) -> 0.05x * exp(-(y - 1)^2 / 0.5) / sqrt(0.5pi))
        @test length(meshes(life)) == 1 && only(meshes(life)).cells == stage_names(plant)
        crowded = species(:plant, births(respond(life, :seeding => ctx -> max(0.0, 1 - ctx.pool.space / 500)), :seeding);
                          occupies = [:space => trait_weights(plant, z -> 0.1z^2)])
        @test !isempty(meshes(system(crowded)[:plant]))
        @test MatrixCommunityModels._continuous(crowded)
        @test_throws ArgumentError project(crowded, fill(1.0, 20), fill(nothing, 2); solver = :fast)
        # a matrix species with it is a component of discrete classes
        shrub_life = ValuedProjectionNet([:seedling, :adult], :maturing => [(:seedling => :adult) => 0.3],
                                         :persisting => [(:adult => :adult) => 0.9], :seeding => [(:adult => :seedling) => 2.0])
        shrub = species(:shrub, births(respond(shrub_life, :seeding => ctx -> max(0.0, 1 - ctx.pool.space / 500)), :seeding);
                        occupies = [:space => [0.1, 1.0]])
        com = community(crowded, shrub)
        layout = CPD.state_layout(system(com)[:shrub])
        @test keys(layout) == (:seedling, :adult)
        @test keys(CPD.state_layout(system(com)[:plant])) == (:plant,)
        same(com, [fill(3.0, 20); 4.0; 2.0], fill(nothing, 25); rng_seed = 3)
        # without a continuous stage, :auto is MatrixProjectionModels.jl, and :ipm gives the same
        @test !MatrixCommunityModels._continuous(shrub)
        same(shrub, [4.0, 2.0], fill(nothing, 25); rng_seed = 4)
        # a species frozen as an IPMProblem has a continuous and discrete states
        prob = ipm_problem(com, :plant, [fill(3.0, 20); 4.0; 2.0], nothing)
        @test prob.kernel isa IPM.DiscretizedKernel && prob.domain === prob.kernel.states
        @test keys(prob.domain) == (:plant,) && prob.domain.plant == plant.domain
        @test IPM.materialize(prob.kernel) == community_matrix(com, [fill(3.0, 20); 4.0; 2.0], nothing)[1:20, 1:20]
    end

    @testset "the published IPMs" begin
        # vignette 14's thistles and weevil
        trio = community(platte_thistle(), musk_thistle(), flowerhead_weevil())
        same(trio, filled(trio, 2.0), fill(nothing, 15); rng_seed = 14)
        # the sea fan, stratified by infection, and with a third stratum
        reef = sea_fan(; cells = 40)
        endemic = SEA_FAN_PHASES.endemic
        start = [fill(20.0, 40); fill(2.0, 40)]
        same(reef, start, fill(endemic, 15); rng_seed = 16)
        three = sea_fan(; cells = 40, resistant = true)
        @test keys(CPD.state_layout(system(three)[:sea_fan])) ==
              (:fan_healthy, :fan_infected, :fan_resistant)
        same(three, [start; fill(5.0, 40)], fill(endemic, 10); rng_seed = 17)
        # two trees on a tabulated kernel
        cells = 40
        pair = community(forest_tree(:spruce; cells), forest_tree(:fir; cells))
        stand = [forest_stand(:spruce; cells, basal_area = 15.0); forest_stand(:fir; cells, basal_area = 15.0)]
        same(pair, stand, fill(nothing, 12); rng_seed = 15)
        # the sagebrush steppe in fluctuating years
        steppe = sagebrush_community()
        state = state_vector(steppe, (k => sagebrush_state(k) for k in species_names(steppe))...)
        same(steppe, state, [nothing, 1, 2, nothing, 5, 7])
    end

    @testset "kernels of IntegralProjectionModels.jl" begin
        domain = ContinuousDomain(0.0, 10.0, 30)
        plant = ContinuousStage(:plant, domain)
        z = meshpoints(plant)
        P = IPM.PKernel(IPM.LinearSurvival(-0.5, 0.4), IPM.NormalGrowth(1.2, 0.85, 0.7), domain)
        F = IPM.FKernel(IPM.FecundityRate(-1.0, 0.5, 1.0, 0.5, 0.2), domain)
        life = kernel_net(plant, :growing => P, :seeding => F)
        @test transition_matrix(life.base, :growing) == IPM.materialize(P)
        @test to_matrix(life) == IPM.materialize(P) .+ IPM.materialize(F)
        # with an eviction correction
        corrected = IPM.PKernel(IPM.LinearSurvival(-0.5, 0.4), IPM.NormalGrowth(1.2, 0.85, 0.7), domain;
                                eviction = IPM.DiscreteExtrema)
        @test to_matrix(kernel_net(plant, :growing => corrected)) == IPM.materialize(corrected)
        # between the continuous stages of a larger net, and returned by a responsive kernel
        net = kernel_net([:seed, plant], :growing => [(:plant => :plant) => P],
                         :sprouting => [(:seed => :plant) => y -> 0.1 * exp(-(y - 1)^2)],
                         :shading => [(:plant => :plant) => responsive(ctx -> P)])
        ctx = context(net, ones(31), NamedTuple(), nothing)
        @test survival_fecundity_matrices(net, ctx)[1][2:end, 2:end] ≈ 2 .* IPM.materialize(P)
        # the domain must be the mesh, and the rule the midpoint rule
        @test_throws ArgumentError kernel_net(ContinuousStage(:plant, 0.0, 10.0, 20), :growing => P)
        @test_throws ArgumentError kernel_net(plant, :growing => P; rule = :simpson)
        @test_throws ArgumentError kernel_net([:seed, plant], :growing => [(:seed => :plant) => P])
        # a community of such a species is solved by IntegralProjectionModels.jl, like any IPM
        sp = species(:plant, births(respond(life, :seeding => ctx -> max(0.0, 1 - ctx.pool.space / 300)), :seeding);
                     occupies = [:space => trait_weights(plant, x -> 0.02x^2)])
        same(sp, fill(2.0, 30), fill(nothing, 20); rng_seed = 6)
    end
end
