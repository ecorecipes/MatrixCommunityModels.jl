# Integral projection models as nets: continuous stages, kernels, responses and pool weights by trait
using LinearAlgebra
using Random

@testset "kernel nets" begin
    normpdf(x, m, s) = exp(-0.5 * ((x - m) / s)^2) / (s * sqrt(2pi))
    survival(z) = 1 / (1 + exp(0.5 - 0.4z))
    P(y, x) = survival(x) * normpdf(y, 1.2 + 0.85x, 0.7)
    F(y, x) = exp(-1.0 + 0.5x) * 0.2 * normpdf(y, 1.0, 0.5)
    dominant(A) = maximum(abs, eigvals(A))

    @testset "continuous stages" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 50)
        @test length(plant) == 50
        @test stage_names(plant)[[1, 50]] == [:plant1, :plant50]
        @test step_size(plant) ≈ 0.2
        @test meshpoints(plant)[[1, 50]] ≈ [0.1, 9.9]
        @test meshpoints(ContinuousStage(:plant, ContinuousDomain(0.0, 10.0, 50))) == meshpoints(plant)
        # Names shared with the rest of the stack are the same functions there, so they stay usable when the
        # other packages are loaded alongside
        upstream = MatrixCommunityModels.CategoricalPopulationDynamics
        @test meshpoints === upstream.meshpoints && step_size === upstream.step_size
        for other in (upstream, MatrixCommunityModels.StructuredPopulationCore, MatrixCommunityModels.MatrixProjectionModels,
                      MatrixCommunityModels.IntegralProjectionModels)
            common = intersect(names(MatrixCommunityModels), names(other))
            @test all(name -> getfield(MatrixCommunityModels, name) === getfield(other, name), common)
        end
        together = Module()
        Core.eval(together, :(using MatrixCommunityModels, MatrixCommunityModels.CategoricalPopulationDynamics,
                              MatrixCommunityModels.MatrixProjectionModels, MatrixCommunityModels.IntegralProjectionModels))
        @test Core.eval(together, :meshpoints) === meshpoints && Core.eval(together, :step_size) === step_size
        @test Core.eval(together, :stratify) === stratify && Core.eval(together, :quasi_extinction) === quasi_extinction
    end

    @testset "one continuous stage" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 50)
        z, h = meshpoints(plant), step_size(plant)
        life = kernel_net(plant, :growing => P, :seeding => F)
        @test life isa ResponsiveProjectionNet
        @test stage_names(life) == stage_names(plant)
        @test transition_names(life) == [:growing, :seeding]
        A = to_matrix(life)
        @test A ≈ [h * (P(y, x) + F(y, x)) for y in z, x in z]
        # the name of the stage is the total of its cells
        ctx = context(life, collect(1.0:50.0), NamedTuple(), nothing)
        @test ctx.n[:plant] == sum(1:50)
        # the same net from ready-made matrices
        frozen = matrix_net(stage_names(plant), :growing => [h * P(y, x) for y in z, x in z],
                            :seeding => [h * F(y, x) for y in z, x in z])
        @test to_matrix(frozen) ≈ A
        # the growth rate does not depend on the mesh
        rate(cells) = dominant(to_matrix(kernel_net(ContinuousStage(:plant, 0.0, 10.0, cells), :growing => P, :seeding => F)))
        @test rate(50) ≈ rate(200) rtol = 1e-4
        # other quadrature rules agree
        @test dominant(to_matrix(kernel_net(plant, :growing => P, :seeding => F; rule = :simpson))) ≈ rate(200) rtol = 1e-3
    end

    @testset "ordinary and continuous stages together" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 40)
        z, h = meshpoints(plant), step_size(plant)
        net = kernel_net([:seed, plant],
            :dormancy    => [(:seed => :seed) => 0.5],
            :germinating => [(:seed => :plant) => y -> 0.1 * normpdf(y, 1.0, 0.5)],
            :growing     => [(:plant => :plant) => P],
            :seeding     => [(:plant => :seed) => x -> exp(-1.0 + 0.5x)])
        @test stage_names(net) == vcat(:seed, stage_names(plant))
        A = to_matrix(net)
        @test A[1, 1] == 0.5
        @test A[2:end, 1] ≈ [h * 0.1 * normpdf(y, 1.0, 0.5) for y in z]
        @test A[1, 2:end] ≈ [exp(-1.0 + 0.5x) for x in z]
        @test A[2:end, 2:end] ≈ [h * P(y, x) for y in z, x in z]
        ctx = context(net, vcat(7.0, ones(40)), NamedTuple(), nothing)
        @test ctx.n[:seed] == 7 && ctx.n[:plant] == 40
        # the numbers in the cells, found by name wherever the stage sits in the net
        counted = context(net, vcat(7.0, collect(1.0:40.0)), NamedTuple(), nothing)
        @test counted.n[plant] == collect(1.0:40.0) && counted.n[:plant] == sum(1:40)

        # two traits on different meshes: the state is individuals per cell, so a kernel is weighted by
        # the width of the destination cell
        small, large = ContinuousStage(:small, 0.0, 2.0, 10), ContinuousStage(:large, 0.0, 10.0, 25)
        two = kernel_net([small, large], :maturing => [(:small => :large) => (y, x) -> 0.3 * normpdf(y, 3.0 + x, 0.8)])
        B = to_matrix(two)
        @test B[11:35, 1:10] ≈ [step_size(large) * 0.3 * normpdf(y, 3.0 + x, 0.8) for y in meshpoints(large), x in meshpoints(small)]
        @test all(isapprox.(sum(B[11:35, 1:10]; dims = 1), 0.3; rtol = 1e-3))     # 30% of every small cell matures
    end

    @testset "responsive kernels" begin
        K = 1000.0
        slowed(pool) = (y, x) -> survival(x) * normpdf(y, 1.2 * (1 - 0.5 * min(1.0, pool / K)) + 0.85x, 0.7)
        function perennial(cells)
            plant = ContinuousStage(:plant, 0.0, 10.0, cells)
            net = kernel_net(plant, :growing => responsive(ctx -> slowed(ctx.pool.space)), :seeding => F)
            net = births(respond(net, :seeding => ctx -> max(0.0, 1 - ctx.pool.space / K)), :seeding)
            return plant, species(:perennial, net; occupies = [:space => trait_weights(plant, z -> 0.02z^2)])
        end
        plant, sp = perennial(50)
        z, h = meshpoints(plant), step_size(plant)
        net = system(sp)[:perennial]
        empty = to_matrix(net, context(net, zeros(50), (space = 0.0,), nothing))
        @test empty ≈ [h * (P(y, x) + F(y, x)) for y in z, x in z]
        half = context(net, zeros(50), (space = 500.0,), nothing)
        U, B = survival_fecundity_matrices(net, half)
        @test U ≈ [h * slowed(500.0)(y, x) for y in z, x in z]          # the growth kernel has changed shape
        @test B ≈ 0.5 .* [h * F(y, x) for y in z, x in z]               # the seeding kernel has only been scaled
        # the net can still be refined: replicated in two strata, each stratum has the kernel of the moment
        two = stratify(net, [:north, :south])
        @test stage_names(two)[[1, 51]] == [:plant1_north, :plant1_south]
        U2, B2 = survival_fecundity_matrices(two, context(two, zeros(100), (space = 500.0,), nothing))
        @test U2 == kron(Matrix(1.0I, 2, 2), U) && B2 == kron(Matrix(1.0I, 2, 2), B)
        # equilibrium and invasion analysis do not depend on the mesh
        function at_equilibrium(cells)
            plant, sp = perennial(cells)
            n0 = zeros(cells)
            n0[meshpoints(plant) .< 2] .= 20 / count(meshpoints(plant) .< 2)
            e = equilibrium(sp, n0, nothing)
            return e, sp
        end
        (e50, sp50), (e100, sp100) = at_equilibrium(50), at_equilibrium(100)
        @test e50.converged && e100.converged
        @test sum(e50.n) ≈ sum(e100.n) rtol = 2e-3
        @test pool_totals(sp50, e50.n, nothing).space ≈ pool_totals(sp100, e100.n, nothing).space rtol = 2e-3
        # the same equilibrium without iterating: the pool value at which the growth rate is one
        direct = equilibrium(sp50, :space => (1.0, K), e50.n ./ 7, nothing)
        @test direct.converged
        @test direct.total ≈ pool_totals(sp50, e50.n, nothing).space rtol = 1e-5
        @test pool_totals(sp50, direct.n, nothing).space ≈ direct.total
        @test direct.n ≈ e50.n rtol = 1e-4
        @test advance(sp50, direct.n, nothing) ≈ direct.n rtol = 1e-8
        @test_throws ArgumentError equilibrium(sp50, :space => (1.0, 2.0), e50.n, nothing)       # not bracketed
        @test_throws ArgumentError equilibrium(sp50, :space => (1.0, K), zeros(50), nothing)

        # composed with a matrix model through the pool
        shrub_life = ValuedProjectionNet([:seedling, :juvenile, :adult],
            :establishing => [(:seedling => :juvenile) => 0.30], :maturing => [(:juvenile => :adult) => 0.50],
            :persisting => [(:adult => :adult) => 0.95], :seeding => [(:adult => :seedling) => 3.0])
        shrub = species(:shrub, births(respond(shrub_life, :seeding => ctx -> max(0.0, 1 - ctx.pool.space / K)), :seeding);
                        occupies = [:space => [0.1, 1.0, 4.0]])
        @test invasion_growth_rate(shrub, sp50, e50.n, nothing) ≈ invasion_growth_rate(shrub, sp100, e100.n, nothing) rtol = 1e-3
        @test invasion_growth_rate(shrub, sp50, e50.n, nothing) > 1
        com = community(sp50, shrub)
        sol = project(com, [e50.n; 5.0; 0.0; 1.0], fill(nothing, 50))
        @test size(sol[:perennial]) == (51, 50) && size(sol[:shrub]) == (51, 3)
        @test stage_total(sol, :perennial, :plant) ≈ abundance(sol, :perennial)
        @test abundance(sol, :shrub)[end] > abundance(sol, :shrub)[1]
        # whole individuals: demographic stochasticity over the mesh cells
        random = project(com, [round.(e50.n); 5.0; 0.0; 1.0], fill(nothing, 20); rng = MersenneTwister(1))
        @test all(isinteger, random.u) && all(>=(0), random.u)
    end

    @testset "responses and pool weights by trait" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 20)
        z, h = meshpoints(plant), step_size(plant)
        net = kernel_net([:seed, plant], :sprouting => [(:seed => :plant) => y -> normpdf(y, 1.0, 0.5)],
                         :seeding => [(:plant => :seed) => x -> 2.0])
        damaged = respond(net, :seeding => ctx -> by_trait(plant, x -> exp(-0.1 * ctx.pool.eggs * x)),
                          :sprouting => ctx -> by_trait(plant, x -> 0.0))       # leaves entries from :seed alone
        A = to_matrix(damaged, context(damaged, zeros(21), (eggs = 3.0,), nothing))
        @test A[1, 2:end] ≈ [2.0 * exp(-0.3x) for x in z]
        @test A[2:end, 1] ≈ [h * normpdf(y, 1.0, 0.5) for y in z]
        # the factors can be given as a vector, one for each cell
        halved = respond(net, :seeding => ctx -> by_trait(plant, fill(0.5, 20)))
        @test to_matrix(halved, context(halved, zeros(21), NamedTuple(), nothing))[1, 2:end] == fill(1.0, 20)
        @test_throws DimensionMismatch by_trait(plant, fill(0.5, 19))
        weights = trait_weights(plant, x -> x^2)
        @test first.(weights) == stage_names(plant) && last.(weights) ≈ z .^ 2
        sp = species(:thistle, damaged; occupies = [:heads => weights], reads = [:eggs])
        n = vcat(4.0, collect(1.0:20.0))
        @test pool_totals(system(sp), Dict(:thistle => n), nothing).heads ≈ sum(z .^ 2 .* (1:20))
        # both keep working when the net is refined: the cells of every stratum take the factor, and the
        # weight, of the cell they descend from
        strata = stratify(damaged, [:north, :south])
        @test stage_names(strata)[[1, 2, 22, 23]] == [:seed_north, :plant1_north, :seed_south, :plant1_south]
        @test to_matrix(strata, context(strata, zeros(42), (eggs = 3.0,), nothing)) == kron(Matrix(1.0I, 2, 2), A)
        both = species(:thistle, strata; occupies = [:heads => weights], reads = [:eggs])
        n2 = vcat(4.0, collect(1.0:20.0), 5.0, 2 .* collect(1.0:20.0))
        @test pool_totals(system(both), Dict(:thistle => n2), nothing).heads ≈ 3 * sum(z .^ 2 .* (1:20))
        @test context(strata, n2, (eggs = 3.0,), nothing).n[plant] == 3 .* collect(1.0:20.0)   # each cell, over both strata
        # the cells of one stratum, and pool weights for one stratum
        @test context(strata, n2, (eggs = 3.0,), nothing).n[plant, :north] == collect(1.0:20.0)
        @test context(strata, n2, (eggs = 3.0,), nothing).n[plant, :south] == 2 .* collect(1.0:20.0)
        @test_throws KeyError context(strata, n2, (eggs = 3.0,), nothing).n[plant, :west]
        southern = trait_weights(plant, x -> x^2; stratum = :south)
        @test first.(southern) == stratum_stage.(stage_names(plant), :south) && last.(southern) == last.(weights)
        south = species(:thistle, strata; occupies = [:heads => southern], reads = [:eggs])
        @test pool_totals(system(south), Dict(:thistle => n2), nothing).heads ≈ 2 * sum(z .^ 2 .* (1:20))
        # a response that reads the stage's own numbers by name sees those totals
        crowded = respond(net, :seeding => ctx -> by_trait(plant, 1 ./ (1 .+ ctx.n[plant] ./ 10)))
        A1 = to_matrix(crowded, context(crowded, vcat(4.0, 3 .* collect(1.0:20.0)), NamedTuple(), nothing))
        @test A1[1, 2:end] ≈ [2.0 / (1 + 3i / 10) for i in 1:20]
        split = stratify(crowded, [:north, :south])
        @test to_matrix(split, context(split, n2, NamedTuple(), nothing)) == kron(Matrix(1.0I, 2, 2), A1)
    end

    @testset "responses by trait in one stratum, attached after stratify" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 6)
        z, h = meshpoints(plant), step_size(plant)
        life = kernel_net([:seed, plant], :sprouting => [(:seed => :plant) => y -> normpdf(y, 1.0, 0.5)],
                          :growing => [(:plant => :plant) => P], :seeding => [(:plant => :seed) => x -> 2.0])
        # growing plants stay healthy or are infected; coupling matrices are [to, from]
        two = stratify(life, [:healthy, :infected]; coupling = Dict(:growing => [1 0; 1 1]))
        @test Set(transition_names(two)) == Set([:sprouting_healthy, :sprouting_infected, :growing_healthy,
            :growing_healthy_to_infected, :growing_infected, :seeding_healthy, :seeding_infected])
        risk(x) = 0.01x
        fewer(x) = exp(-0.2x)
        sick = respond(two,
            :growing_healthy => ctx -> by_trait(plant, x -> 1 - risk(x); stratum = :healthy),
            :growing_healthy_to_infected => ctx -> by_trait(plant, risk; stratum = :healthy),
            :seeding_infected => ctx -> by_trait(plant, fewer; stratum = :infected),
            :sprouting_infected => ctx -> by_trait(plant, x -> 0.0; stratum = :infected))   # leaves entries from :seed alone
        A = to_matrix(sick, context(sick, zeros(14), NamedTuple(), nothing))
        at = Dict(s => i for (i, s) in enumerate(stage_names(sick)))
        cell(i, stratum) = at[stratum_stage(stage_names(plant)[i], stratum)]
        seed(stratum) = at[stratum_stage(:seed, stratum)]
        expected = zeros(14, 14)
        for stratum in (:healthy, :infected), i in 1:6
            expected[cell(i, stratum), seed(stratum)] = h * normpdf(z[i], 1.0, 0.5)
            expected[seed(stratum), cell(i, stratum)] = stratum === :healthy ? 2.0 : 2.0 * fewer(z[i])
        end
        for i in 1:6, j in 1:6
            expected[cell(i, :healthy), cell(j, :healthy)] = h * P(z[i], z[j]) * (1 - risk(z[j]))
            expected[cell(i, :infected), cell(j, :healthy)] = h * P(z[i], z[j]) * risk(z[j])
            expected[cell(i, :infected), cell(j, :infected)] = h * P(z[i], z[j])
        end
        @test all(isapprox(A[i, j], expected[i, j]; rtol = 1e-14) for i in 1:14, j in 1:14)
        @test iszero.(A) == iszero.(expected)
        # called directly: entries leaving the cells of the stratum, and the rest
        g = by_trait(plant, fewer; stratum = :infected)
        @test g(:plant3_infected, :seed_infected) == fewer(z[3])
        @test g(:seed_infected, :plant3_infected) == 1.0
        @test by_trait(plant, fill(0.5, 6); stratum = :infected)(:plant6_infected, :seed_infected) == 0.5
        @test_throws DimensionMismatch by_trait(plant, fill(0.5, 5); stratum = :infected)
        # the response stays with its cells when the net is refined again
        patches = stratify(sick, [:north, :south])
        @test to_matrix(patches, context(patches, zeros(28), NamedTuple(), nothing)) == kron(Matrix(1.0I, 2, 2), A)

        # without a stratum nothing changes: before stratify a response applies in every stratum, and after it a
        # response written for the unstratified cells matches none of them
        base = to_matrix(two)
        before = stratify(respond(life, :seeding => ctx -> by_trait(plant, fewer)), [:healthy, :infected];
                          coupling = Dict(:growing => [1 0; 1 1]))
        B = to_matrix(before, context(before, zeros(14), NamedTuple(), nothing))
        @test B[seed(:healthy), cell(4, :healthy)] == B[seed(:infected), cell(4, :infected)] == 2.0 * fewer(z[4])
        after = respond(two, :seeding_infected => ctx -> by_trait(plant, fewer))
        @test to_matrix(after, context(after, zeros(14), NamedTuple(), nothing)) == base

        # an unknown stratum, the wrong one, or one on a net without strata is an error, raised when the
        # response is evaluated, rather than a factor that silently does nothing
        evaluated(net, k) = to_matrix(net, context(net, zeros(k), NamedTuple(), nothing))
        @test_throws ArgumentError evaluated(respond(two, :seeding_infected => ctx -> by_trait(plant, fewer; stratum = :sick)), 14)
        @test_throws ArgumentError evaluated(respond(two, :seeding_infected => ctx -> by_trait(plant, fewer; stratum = :healthy)), 14)
        @test_throws ArgumentError evaluated(respond(life, :seeding => ctx -> by_trait(plant, fewer; stratum = :healthy)), 7)
        # stage names that only look alike are other stages, and left alone
        @test by_trait(plant, fewer; stratum = :sick)(:plantation_sick, :seed_sick) == 1.0
        @test by_trait(plant, fewer; stratum = :sick)(:seed_infected, :seed_sick) == 1.0
    end

    @testset "errors" begin
        plant = ContinuousStage(:plant, 0.0, 10.0, 10)
        @test_throws ArgumentError kernel_net([:seed, plant], :growing => P)                          # which stages?
        @test_throws ArgumentError kernel_net(plant, :growing => x -> x)                              # needs (z′, z)
        @test_throws ArgumentError kernel_net([:seed, plant], :sprouting => [(:seed => :plant) => P])  # needs z′ only
        @test_throws ArgumentError kernel_net([:seed, plant], :dormancy => [(:seed => :seed) => x -> x])
        @test_throws ArgumentError kernel_net([:seed, plant], :sprouting => [(:seed => :tree) => 0.5])
        @test_throws ArgumentError kernel_net([:plant1, plant], :growing => [(:plant => :plant) => P]) # name clash
        @test_throws ArgumentError kernel_net([plant, plant], :growing => [(:plant => :plant) => P])
    end
end
