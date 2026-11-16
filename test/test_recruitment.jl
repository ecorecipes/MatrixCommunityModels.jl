# The lottery for space, against the hand-written versions it replaces, and the
# sensitivity of invasion growth rates, against finite differences.

import StructuredPopulationCore

# ---------------------------------------------------------------- the lottery

# the annual plant of the seasons vignette, with its lottery written by hand ...
function annual_by_hand(name; germinate, seed, K = 1000.0, flood = 1, g = 0.5)
    life = ValuedProjectionNet([:seed, :plant],
        :dormancy => [(:seed => :seed) => 0.9], :germination => [(:seed => :plant) => g],
        :growth => [(:plant => :plant) => 0.9], :seed_set => [(:plant => :seed) => 20.0])
    during(season) = ctx -> ctx.env == season ? 1.0 : 0.0
    establishment(ctx) = min(1.0, max(K - ctx.pool.space, 0.0) / max(ctx.pool.demand, 1e-12))
    plant = respond(life, :dormancy => ctx -> ctx.env == germinate ? 1 - g : 1.0,
        :germination => during(germinate), :germination => establishment, :seed_set => during(seed),
        :growth => ctx -> ctx.env == seed ? 0.0 : (ctx.env == flood ? 0.1 / 0.9 : 1.0))
    species(name, births(plant, :seed_set);
            occupies = [:space => [:plant => 1.0], :demand => (n, env) -> env == germinate ? g * n[:seed] : 0.0])
end

# ... and with `lottery`
function annual_by_lottery(name; germinate, seed, K = 1000.0, flood = 1, g = 0.5)
    life = ValuedProjectionNet([:seed, :plant],
        :dormancy => [(:seed => :seed) => 0.9], :germination => [(:seed => :plant) => g],
        :growth => [(:plant => :plant) => 0.9], :seed_set => [(:plant => :seed) => 20.0])
    during(season) = ctx -> ctx.env == season ? 1.0 : 0.0
    plant = respond(life, :dormancy => ctx -> ctx.env == germinate ? 1 - g : 1.0,
        :germination => during(germinate), :seed_set => during(seed),
        :growth => ctx -> ctx.env == seed ? 0.0 : (ctx.env == flood ? 0.1 / 0.9 : 1.0))
    plant = lottery(plant, :germination; capacity = K)
    species(name, births(plant, :seed_set);
            occupies = [:space => [:plant => 1.0], :demand => lottery_demand(plant, :germination)])
end

# musk thistle with its open-ground response from `lottery` (the rest as in src/thistles.jl)
function musk_thistle_by_lottery(; p = MUSK_THISTLE_PARAMETERS, weevil_survival = WEEVIL_PARAMETERS.survival)
    rosettes = [:small, :medium, :large]
    bolters = [p.survival[j] * p.bolting[j] for j in 1:3]
    staying = [p.survival[j] * (1 - p.bolting[j]) for j in 1:3]
    to_medium, to_large = p.growth_small
    sized(value) = [(rosettes[k], value * p.seedling_size[k]) for k in 1:3 if p.seedling_size[k] > 0]
    life = ValuedProjectionNet([:seed_bank, :small, :medium, :large],
        :dormancy => [(:seed_bank => :seed_bank) => p.seed_survival],
        :emerging => [(:seed_bank => to) => v for (to, v) in sized(p.emergence)],
        :growing => [(:small => :small) => staying[1] * (1 - to_medium - to_large),
                     (:small => :medium) => staying[1] * to_medium, (:small => :large) => staying[1] * to_large,
                     (:medium => :medium) => staying[2] * (1 - p.growth_medium),
                     (:medium => :large) => staying[2] * p.growth_medium,
                     (:large => :large) => staying[3]],
        :seed_rain => [(rosettes[j] => :seed_bank) => bolters[j] * p.seeds[j] * p.to_seed_bank for j in 1:3],
        :establishing => [(rosettes[j] => to) => v for j in 1:3 for (to, v) in sized(bolters[j] * p.seeds[j] * p.establishment)])
    function escaping(ctx)
        (ctx.pool.heads > 0 && ctx.pool.rearing > 0) || return 1.0
        eggs_per_head = p.attraction * ctx.pool.eggs / ctx.pool.heads
        larvae = weevil_survival * eggs_per_head / (1 + weevil_survival * ctx.pool.eggs / ctx.pool.rearing)
        return 1 - p.exposed * min(1.0, larvae * p.seeds_per_larva / p.seeds_per_head)
    end
    net = respond(life, :seed_rain => escaping, :establishing => escaping)
    net = lottery(net, :emerging, :establishing; capacity = p.sites, space = rosettes, demand = :own, rule = :saturating)
    net = births(net, :seed_rain, :establishing)
    early = [bolters[j] * p.exposed * p.seeds[j] / p.seeds_per_head for j in 1:3]
    occupies = Pair{Symbol, Any}[:heads => [rosettes[j] => p.attraction * early[j] for j in 1:3],
                                 :rearing => [rosettes[j] => p.weevils_per_head * early[j] for j in 1:3]]
    return species(:musk_thistle, net; occupies, reads = [:eggs]), net
end

@testset "lottery" begin
    # the seasons vignette: identical trajectories and pools for its three pairs of annuals
    envs = repeat([1, 2, 3, 4], 100)
    for (ga, sa, gb, sb) in [(2, 3, 4, 1), (2, 4, 3, 1), (2, 3, 2, 3)]
        hand = community(annual_by_hand(:a; germinate = ga, seed = sa), annual_by_hand(:b; germinate = gb, seed = sb))
        drawn = community(annual_by_lottery(:a; germinate = ga, seed = sa), annual_by_lottery(:b; germinate = gb, seed = sb))
        h, l = project(hand, [100.0, 0, 100, 0], envs), project(drawn, [100.0, 0, 100, 0], envs)
        @test h.u == l.u
        @test isequal(h.pools, l.pools)
        for k in (1, 2, 3, 50, 101), env in 1:4
            @test community_matrix(hand, h.u[k, :], env) == community_matrix(drawn, h.u[k, :], env)
        end
    end
    # a pulse never overshoots: plants never exceed the space, however many seeds germinate
    alone = annual_by_lottery(:a; germinate = 2, seed = 3)
    pulse = project(alone, [1e7, 0.0], repeat([1, 2, 3, 4], 5))
    @test maximum(pulse[:a][:, 2]) <= 1000.0 * (1 + 1e-12)
    @test pulse[:a][3, 2] ≈ 1000.0                                  # the pulse fills the space exactly

    # musk thistle: the same model within rounding, alone and with the weevil
    weed, net = musk_thistle_by_lottery()
    @test transition_names(net) == transition_names(system(musk_thistle())[:musk_thistle])
    n0 = [1000.0, 10.0, 5.0, 1.0]
    years = fill(nothing, 60)
    @test project(weed, n0, years).u ≈ project(musk_thistle(), n0, years).u rtol = 1e-12
    weevil = flowerhead_weevil()
    with_lottery, by_hand = community(weed, weevil), community(musk_thistle(), weevil)
    sol = project(by_hand, [n0; 5.0], years)
    @test project(with_lottery, [n0; 5.0], years).u ≈ sol.u rtol = 1e-12
    for k in (1, 10, 61)
        @test community_matrix(with_lottery, sol.u[k, :], nothing) ≈ community_matrix(by_hand, sol.u[k, :], nothing) rtol = 1e-13
    end

    # the rules, read from the species' own stages
    life = ValuedProjectionNet([:seed, :adult], :recruiting => [(:seed => :adult) => 0.5],
                               :staying => [(:adult => :adult) => 0.8])
    for (rule, expected) in [(:lottery, min(1, 40 / 50)), (:saturating, 40 / 90), ((f, d) -> f / (f + 2d), 40 / 140)]
        c = species(:p, lottery(life, :recruiting; capacity = 100.0, space = [:adult => 2.0], demand = :own, rule))
        @test community_matrix(c, [100.0, 30.0], nothing)[2, 1] ≈ 0.5 * expected     # demand 50, free 100 - 60
    end
    # no demand: a recruit establishes if there is any free space
    own = species(:p, lottery(life, :recruiting; capacity = 100.0, space = [:adult], demand = :own))
    @test community_matrix(own, [0.0, 30.0], nothing)[2, 1] == 0.5
    @test community_matrix(own, [0.0, 100.0], nothing)[2, 1] == 0.0
    @test community_matrix(own, [10.0, 120.0], nothing)[2, 1] == 0.0          # more held than the capacity

    # misuse
    @test_throws ArgumentError lottery(life, :flowering; capacity = 1.0)
    @test_throws ArgumentError lottery(life; capacity = 1.0)
    @test_throws ArgumentError lottery(life, :recruiting; capacity = 1.0, rule = :fair)
    @test_throws ArgumentError lottery(life, :recruiting; capacity = 1.0, space = 3)
    @test_throws ArgumentError lottery(life, :recruiting; capacity = 1.0, demand = 3)
    @test_throws ArgumentError lottery_demand(net, :establishing)([1.0, 1.0, 1.0, 1.0], nothing)   # reads a pool
    @test lottery_demand(net, :emerging)([10.0, 0.0, 0.0, 0.0], nothing) ≈ 10 * MUSK_THISTLE_PARAMETERS.emergence
end

# ---------------------------------------------------------------- sensitivity of invasion

@testset "invasion sensitivity" begin
    grass = toy_grass()
    pool = (space = 900.0,)
    s = invasion_sensitivity(grass, pool, nothing)
    A = MatrixCommunityModels._rare_matrix(grass, pool, nothing)
    λ(B) = StructuredPopulationCore.lambda(B)
    @test s.lambda ≈ invasion_growth_rate(grass, pool, nothing)
    @test s.lambda < 1                                          # excluded at this much space held
    @test sum(s.elasticity) ≈ 1
    @test sum(x.elasticity for x in s.transitions) ≈ 1

    # every entry, against central differences
    δ = 1e-6
    for i in 1:2, j in 1:2
        up, dn = copy(A), copy(A)
        up[i, j] += δ
        dn[i, j] -= δ
        @test s.sensitivity[i, j] ≈ (λ(up) - λ(dn)) / 2δ rtol = 1e-6
    end
    # every named transition: scale all of its entries by 1 ± δ
    for row in s.transitions
        At = zeros(2, 2)
        stages = Dict(:seedling => 1, :adult => 2)
        for x in s.entries
            x.transition === row.transition && (At[stages[x.to], stages[x.from]] += x.value)
        end
        @test row.sensitivity ≈ (λ(A + δ * At) - λ(A - δ * At)) / 2δ rtol = 1e-6
        @test row.elasticity ≈ row.sensitivity / s.lambda
    end
    seeding = only(x for x in s.entries if x.transition === :seeding)
    @test (seeding.from, seeding.to) == (:adult, :seedling)
    @test seeding.value ≈ 6.0 * (1 - 900 / 1000)
    @test seeding.sensitivity == s.sensitivity[1, 2]

    # at a resident's state
    com = community(toy_grass(), toy_shrub())
    eq = equilibrium(com, [50.0, 10.0, 5.0, 0.0, 1.0], nothing)
    @test invasion_sensitivity(grass, com, eq.n, nothing).lambda ≈ invasion_growth_rate(grass, com, eq.n, nothing)

    # a parameter through a builder: seed production f, analytic λ = (0.5 + √(0.25 + 1.2 f c)) / 2
    c = 1 - 900 / 1000
    exact(f) = 0.3 * c / sqrt(0.25 + 1.2 * f * c)
    by_pool = invasion_sensitivity(f -> (invader = toy_grass(; f), pool = pool, env = nothing), 6.0)
    @test by_pool.lambda ≈ s.lambda
    @test by_pool.derivative ≈ exact(6.0) rtol = 1e-6
    @test by_pool.derivative ≈ only(x for x in s.transitions if x.transition === :seeding).sensitivity / 6 rtol = 1e-6
    @test by_pool.elasticity ≈ 6 / s.lambda * by_pool.derivative
    by_resident = invasion_sensitivity(f -> (invader = toy_grass(; f), resident = com, n = eq.n, env = nothing), 6.0)
    @test by_resident.lambda ≈ invasion_growth_rate(grass, com, eq.n, nothing)
end
