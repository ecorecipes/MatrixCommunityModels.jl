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

@testset "equilibrium sensitivity" begin
    # shrub alone: space at equilibrium is K(1 - 1/(3f)), so d(space)/df = K/(3f^2)
    build(f) = (community = community(toy_shrub(; f)), n0 = [5.0, 0.0, 1.0])
    s = equilibrium_sensitivity(build, 3.0, nothing)
    @test s.converged
    @test dot([0.1, 1.0, 4.0], s.dn) ≈ TOY_K / (3 * 3.0^2) rtol = 1e-5
    # agrees with the regression estimate
    r = abundance_sensitivity(build, collect(2.98:0.01:3.02), fill(nothing, 4000); burnin = 3500)
    @test s.dn ≈ r.slope rtol = 1e-3
end

@testset "stochastic invasion rate" begin
    grass, shrub = toy_grass(), toy_shrub()
    alone = open_community(grass)
    geq = equilibrium(alone, [50.0, 10.0], nothing)
    # constant environment: equals the deterministic invasion growth rate
    @test stochastic_invasion_rate(shrub, alone, geq.n, fill(nothing, 400); burnin = 200) ≈
          invasion_growth_rate(shrub, alone, geq.n, nothing) rtol = 1e-8
    # a resident cannot invade itself at its own equilibrium
    @test stochastic_invasion_rate(grass, alone, geq.n, fill(nothing, 400); burnin = 200) ≈ 1 rtol = 1e-6
    @test_throws ArgumentError stochastic_invasion_rate(shrub, alone, geq.n, fill(nothing, 5); burnin = 5)
end

@testset "demographic stochasticity" begin
    using Random
    grass, shrub = toy_grass(), toy_shrub()
    com = community(grass, shrub)
    n0 = [50.0, 10.0, 5.0, 0.0, 1.0]
    envs = fill(nothing, 30)
    a = project(com, n0, envs; rng = Xoshiro(1))
    @test all(isinteger, a.u) && all(>=(0), a.u)
    @test a.u == project(com, n0, envs; rng = Xoshiro(1)).u          # reproducible
    @test a.u != project(com, n0, envs; rng = Xoshiro(2)).u
    @test project(com, n0, envs).u == project(com, n0, envs; rng = nothing).u
    # the mean of many one-step draws approaches the deterministic step
    big = [5000.0, 1000.0, 500.0, 100.0, 100.0]
    rng = Xoshiro(3)
    draws = reduce(hcat, [advance(com, big, nothing; rng) for _ in 1:4000])
    @test vec(sum(draws; dims = 2)) ./ 4000 ≈ advance(com, big, nothing) rtol = 0.02
    # requirements
    @test_throws ArgumentError project(com, [50.5, 10.0, 5.0, 0.0, 1.0], envs; rng = Xoshiro(1))
    fish = verde_fish_model()
    n0f = round.(reduce(vcat, last.(fish.n0)))
    @test_throws ArgumentError project(fish.community, n0f, fish.envs[1:2]; rng = Xoshiro(1))
end

@testset "interaction network" begin
    effect = [-1.0 0.8 0.0; 0.6 -2.0 0.05; 0.0 0.3 -1.0]
    net = InteractionNetwork(effect, ["a", "b", "c"]; threshold = 0.1)
    @test count(net.adjacency) == 3
    @test connectance(net) == 3 / 6 && link_density(net) == 1.0
    @test reciprocity(net) == 2 / 3
    @test out_degree(net) == [1, 2, 0] && in_degree(net) == [1, 1, 1]
    l = links(net)
    @test l[1] == (from = "b", to = "a", effect = 0.8) && length(l) == 3
    @test_throws DimensionMismatch InteractionNetwork(effect, ["a", "b"]; threshold = 0.1)
end

@testset "periodic matrix" begin
    # a species whose rates differ between two seasons
    sp = Species(:x, [:seed, :plant],
                 [survival(:germinate, :seed => :plant, ctx -> ctx.env == 1 ? 0.5 * toy_crowding(ctx) : 0.0),
                  survival(:dormant, :seed => :seed, ctx -> ctx.env == 1 ? 0.4 : 0.9),
                  fecundity(:set_seed, :plant => :seed, ctx -> ctx.env == 2 ? 10.0 : 0.0)];
                 occupancy = [:space => [0.0, 1.0]])
    A(env, space) = projection_matrix(sp, Context(zeros(2), (space = space,), env))
    @test periodic_matrix(sp, (space = 0.0,), [1, 2]) == A(2, 0.0) * A(1, 0.0)
    @test periodic_matrix(sp, (space = 300.0,), [1]) == A(1, 300.0)
    # projecting through the seasons with nobody else present multiplies by it
    n = [10.0, 0.0]
    empty = (space = 0.0,)
    @test periodic_matrix(sp, empty, [1, 2, 1, 2]) ≈ periodic_matrix(sp, empty, [1, 2])^2
end

@testset "ensembles" begin
    using Random
    grass, shrub = toy_grass(), toy_shrub()
    com = community(grass, shrub)
    n0 = [50.0, 10.0, 5.0, 0.0, 1.0]
    envs = fill(nothing, 40)
    ens = ensemble(com, n0, envs; replicates = 25, rng = Xoshiro(4))
    @test size(ens.u) == (25, 41, 5) && ens.t == collect(0:40)
    @test all(ens.u[:, 1, :] .== n0')
    N = abundance(ens, :shrub)
    @test size(N) == (25, 41) && N[3, :] == vec(sum(ens.u[3, :, 3:5]; dims = 2))
    sol = project(com, n0, envs)
    @test abundance(sol, :grass) == vec(sum(sol[:grass]; dims = 2))
    # without demographic noise and with fixed environments, replicates are identical
    det = ensemble(com, n0, envs; replicates = 3, demographic = false)
    @test det.u[1, :, :] == det.u[3, :, :] == sol.u
    # per-replicate environments
    varied = ensemble(com, n0, r -> fill(nothing, 40); replicates = 2, demographic = false)
    @test varied.u[2, :, :] == sol.u
    # extinction summaries on a doomed population
    doomed = community(Species(:d, [:a], [survival(:s, :a => :a, 0.5)]; occupancy = [:space => [1.0]]))
    e = ensemble(doomed, [20.0], fill(nothing, 30); replicates = 50, rng = Xoshiro(5))
    times = extinction_times(e, :d)
    @test all(!isnothing, times) && all(t -> 1 <= t <= 30, times)
    q = quasi_extinction(e, :d)
    @test q[1] == 0 && q[end] == 1 && issorted(q)
    @test quasi_extinction(e, :d; threshold = 21)[1] == 1
    @test extinction_times(det, :shrub) == fill(nothing, 3)
    @test_throws ArgumentError ensemble(com, n0, envs; replicates = 2)
    @test_throws ArgumentError ensemble(com, n0, envs; replicates = 0, rng = Xoshiro(1))
end

@testset "Graphs extension" begin
    import Graphs
    effect = [-1.0 0.8 0.0; 0.6 -2.0 0.05; 0.0 0.3 -1.0]
    net = InteractionNetwork(effect, ["a", "b", "c"]; threshold = 0.1)
    g = Graphs.SimpleDiGraph(net)
    @test Graphs.nv(g) == 3 && Graphs.ne(g) == 3
    @test Graphs.has_edge(g, 2, 1) && Graphs.has_edge(g, 1, 2) && Graphs.has_edge(g, 2, 3)
    @test Graphs.outdegree(g) == out_degree(net) && Graphs.indegree(g) == in_degree(net)
end
