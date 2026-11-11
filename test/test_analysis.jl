# Composition laws and the analysis layer, on a small two-species community with a
# known equilibrium.

using LinearAlgebra
using Catlab: @relation
import CategoricalPopulationDynamics

const TOY_K = 1000.0
toy_crowding(ctx) = max(0.0, 1 - ctx.pool.space / TOY_K)

toy_grass(; f = 6.0) = Species(:grass, [:seedling, :adult],
                               [survival(:establish, :seedling => :adult, 0.30),
                                survival(:persist, :adult => :adult, 0.50),
                                fecundity(:seed, :adult => :seedling, ctx -> f * toy_crowding(ctx))];
                               occupancy = [:space => [0.1, 1.0]])
toy_shrub(; f = 3.0) = Species(:shrub, [:seedling, :juvenile, :adult],
                               [survival(:establish, :seedling => :juvenile, 0.30),
                                survival(:mature, :juvenile => :adult, 0.50),
                                survival(:persist, :adult => :adult, 0.95),
                                fecundity(:seed, :adult => :seedling, ctx -> f * toy_crowding(ctx))];
                               occupancy = [:space => [0.1, 1.0, 4.0]])

@testset "composition" begin
    grass, shrub = toy_grass(), toy_shrub()
    com = community(grass, shrub)
    n0 = [50.0, 10.0, 5.0, 0.0, 1.0]
    envs = fill(nothing, 60)
    sol = project(com, n0, envs)

    # an explicit wiring diagram gives the same model
    d = @relation (space,) begin
        grass(space)
        shrub(space)
    end
    @test project(compose(d, [grass, shrub]), n0, envs).u == sol.u

    # composites nest, and grouping does not matter
    herb = Species(:herb, [:plant], [survival(:persist, :plant => :plant, 0.2),
                                     fecundity(:seed, :plant => :plant, ctx -> 2.0 * toy_crowding(ctx))];
                   occupancy = [:space => [0.5]])
    outer = @relation (space,) begin
        a(space)
        b(space)
    end
    n3 = [n0; 20.0]
    flat = project(community(grass, shrub, herb), n3, envs).u
    @test project(compose(outer, [com, herb]), n3, envs).u ≈ flat rtol = 1e-13
    regrouped = project(compose(outer, [grass, community(shrub, herb)]), n3, envs).u
    @test regrouped ≈ flat rtol = 1e-13

    # a pool that is not exposed is still shared inside
    closed = @relation () begin
        grass(space)
        shrub(space)
    end
    cc = compose(closed, [grass, shrub])
    @test isempty(pool_names(cc))
    @test project(cc, n0, envs).u == sol.u

    # errors
    @test_throws DimensionMismatch compose(d, [grass])
    @test_throws ArgumentError community(grass, grass)
    @test_throws ArgumentError Species(:x, [:a], [survival(:s, :a => :b, 1.0)])
    @test_throws DimensionMismatch project(com, [1.0], envs)
end

@testset "species representations" begin
    grass = toy_grass()
    ctx = Context([0.0, 0.0], (space = 400.0,), nothing)
    U, F = vital_matrices(grass, ctx)
    @test U == [0 0; 0.3 0.5] && F ≈ [0 3.6; 0 0]
    mpm = lower(grass, ctx)
    @test mpm.U == U && mpm.F == F && mpm.A == U + F
    @test CategoricalPopulationDynamics.to_matrix(valued_net(grass, ctx)) ≈ U + F

    # the same species given as whole matrices behaves identically
    mgrass = MatrixSpecies(:grass, [:seedling, :adult],
                           c -> ([0 0; 0.3 0.5], [0 6.0 * toy_crowding(c); 0 0]); occupancy = [:space => [0.1, 1.0]])
    shrub = toy_shrub()
    n0 = [50.0, 10.0, 5.0, 0.0, 1.0]
    @test project(community(mgrass, shrub), n0, fill(nothing, 40)).u ≈
          project(community(grass, shrub), n0, fill(nothing, 40)).u rtol = 1e-14
    @test species(community(mgrass, shrub), :grass) === mgrass
end

@testset "equilibrium, stability, invasion" begin
    grass, shrub = toy_grass(), toy_shrub()
    com = community(grass, shrub)
    eq = equilibrium(com, [:grass => [50.0, 10.0], :shrub => [5.0, 0.0, 1.0]], nothing)
    @test eq.converged
    # the shrub replaces itself when crowding = 1/9, the grass needs 1/3.6: shrub wins
    @test pool_totals(com, eq.n, nothing).space ≈ TOY_K * (1 - 1 / 9) rtol = 1e-6
    @test all(eq.n[1:2] .< 1e-6)
    @test advance(com, eq.n, nothing) ≈ eq.n rtol = 1e-8

    J = jacobian(com, eq.n, nothing)
    @test maximum(abs, eigvals(J)) < 1
    # the frozen community matrix has the shrub at exactly replacement; the Jacobian adds
    # the feedback through the pool
    A = community_matrix(com, eq.n, nothing)
    @test maximum(abs, eigvals(A[3:5, 3:5])) ≈ 1 rtol = 1e-6
    @test J[1:2, 3:5] ≈ zeros(2, 3) atol = 1e-6      # grass is absent, so nothing acts on it
    @test J[3, 5] < A[3, 5]                           # more adults, fewer seedlings per adult

    # mutual invasibility: shrub invades the grass equilibrium, not the reverse
    geq = equilibrium(open_community(grass), [50.0, 10.0], nothing)
    @test geq.converged
    @test pool_totals(open_community(grass), geq.n, nothing).space ≈ TOY_K * (1 - 0.5 / 1.8) rtol = 1e-6
    @test invasion_growth_rate(shrub, open_community(grass), geq.n, nothing) > 1
    @test invasion_growth_rate(grass, com, eq.n, nothing) < 1
    @test invasion_growth_rate(grass, (space = 0.0,), nothing) ≈ (0.5 + sqrt(0.25 + 4 * 0.3 * 6)) / 2
    @test_throws ArgumentError invasion_growth_rate(grass, (biomass = 0.0,), nothing)
end

@testset "sensitivity" begin
    # at equilibrium total space is K(1 - 1/R0) with R0 = f * 0.3 * 0.5 / 0.05 = 3f for the shrub
    build(f) = (community = community(toy_shrub(; f)), n0 = [5.0, 0.0, 1.0])
    values = 2.9:0.02:3.1
    s = abundance_sensitivity(build, collect(values), fill(nothing, 3000); burnin = 2500)
    @test size(s.means) == (length(values), 3)
    @test all(s.slope .> 0) && all(s.r2 .> 0.99)
    w = [0.1, 1.0, 4.0]
    @test dot(w, s.slope) ≈ TOY_K / (3 * 3.0^2) rtol = 2e-2     # d/df of K(1 - 1/(3f))
    @test mean_abundance(build(3.0).community, [5.0, 0.0, 1.0], fill(nothing, 3000); burnin = 2500) ≈
          s.means[6, :] rtol = 1e-12
    @test_throws ArgumentError abundance_sensitivity(build, [3.0], fill(nothing, 10))

    g = with_rate(RIPARIAN_GUILDS, :C, :S, 6, 0.05)
    @test g.C.S[6] == 0.05 && RIPARIAN_GUILDS.C.S[6] == 0.02 && g.T === RIPARIAN_GUILDS.T
    @test riparian_model(; guilds = g).community isa OpenCommunity
end
