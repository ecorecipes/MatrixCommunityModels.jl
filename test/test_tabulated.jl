# Kernels tabulated at values of a pool and interpolated between them, against kernels rebuilt every step
using LinearAlgebra

@testset "tabulated kernels" begin
    normpdf(x, m, s) = exp(-0.5 * ((x - m) / s)^2) / (s * sqrt(2pi))
    survival(z) = 1 / (1 + exp(0.5 - 0.4z))
    F(y, x) = exp(-1.0 + 0.5x) * 0.2 * normpdf(y, 1.0, 0.5)
    K = 1000.0
    slowed(pool) = (y, x) -> survival(x) * normpdf(y, 1.2 * (1 - 0.5 * min(1.0, pool / K)) + 0.85x, 0.7)
    plant = ContinuousStage(:plant, 0.0, 10.0, 50)
    z, h = meshpoints(plant), step_size(plant)
    built = Ref(0)                                    # kernels built for the tables
    counted(pool) = (built[] += 1; slowed(pool))
    grid = 0:100:1000
    make(growth) = kernel_net(plant, :growing => growth, :seeding => F)
    exact = make(responsive(ctx -> slowed(ctx.pool.space)))
    table = make(tabulated(counted; pool = :space, grid))
    at(net, space) = to_matrix(net, context(net, zeros(50), (space = space,), nothing))

    @testset "at and between the grid values" begin
        @test built[] == 0                             # nothing is built before it is needed
        for space in (0.0, 300.0, 1000.0)
            @test at(table, space) ≈ at(exact, space) rtol = 1e-14
        end
        @test built[] == 3
        at(table, 300.0)
        @test built[] == 3                             # kept, not built again
        # between grid values: the linear interpolation of the two neighbouring kernels
        mixed = at(table, 330.0)
        @test built[] == 4
        @test mixed ≈ 0.7 .* at(exact, 300.0) .+ 0.3 .* at(exact, 400.0) rtol = 1e-13
        # and close to the kernel rebuilt there: the interpolation error, relative to the largest entry
        worst = maximum(maximum(abs, at(table, s) .- at(exact, s)) / maximum(at(exact, s)) for s in 5.0:10.0:995.0)
        println("tabulated kernel, grid every 100 of 1000: largest error between grid values, relative to the ",
                "largest entry = ", worst)
        @test worst < 1e-4
        finer = make(tabulated(slowed; pool = :space, grid = 0:10:1000))
        fine = maximum(maximum(abs, at(finer, s) .- at(exact, s)) / maximum(at(exact, s)) for s in 5.0:10.0:995.0)
        @test fine < worst / 50                       # second order in the spacing
    end

    @testset "outside the grid" begin
        @test_throws DomainError at(table, 1200.0)
        @test_throws DomainError at(table, -1.0)
        @test_throws DomainError at(make(tabulated(slowed; pool = :space, grid, outside = :exact)), NaN)
        clamped = make(tabulated(slowed; pool = :space, grid, outside = :clamp))
        @test at(clamped, 1200.0) ≈ at(exact, 1000.0) rtol = 1e-14
        @test at(clamped, -5.0) ≈ at(exact, 0.0) rtol = 1e-14
        rebuilt = make(tabulated(slowed; pool = :space, grid, outside = :exact))
        @test at(rebuilt, 1200.0) ≈ at(exact, 1200.0) rtol = 1e-14
        @test at(rebuilt, 250.0) ≈ (at(exact, 200.0) .+ at(exact, 300.0)) ./ 2 rtol = 1e-13
        net = make(tabulated(slowed; pool = :crowding, grid))
        @test_throws ArgumentError at(net, 10.0)        # the pool is not in the context
        @test_throws ArgumentError tabulated(slowed; pool = :space, grid = [0.0])
        @test_throws ArgumentError tabulated(slowed; pool = :space, grid = [0.0, 2.0, 1.0])
        @test_throws ArgumentError tabulated(slowed; pool = :space, grid = [0.0, Inf])
        @test_throws ArgumentError tabulated(slowed; pool = :space, grid, outside = :nearest)
    end

    @testset "several stages and kernels" begin
        # a seed bank: a tabulated kernel from an ordinary stage, a tabulated and a responsive kernel in one
        # transition, and a tabulated kernel next to fixed ones
        recruits(pool) = y -> 0.1 * (1 - pool / 2K) * normpdf(y, 1.0, 0.5)
        tabbed = kernel_net([:seed, plant],
            :dormancy => [(:seed => :seed) => 0.5],
            :germinating => [(:seed => :plant) => tabulated(recruits; pool = :space, grid)],
            :growing => [(:plant => :plant) => tabulated(slowed; pool = :space, grid),
                         (:plant => :seed) => responsive(ctx -> (x -> 2.0 * survival(x) * (1 - ctx.pool.space / 2K)))])
        rebuilt = kernel_net([:seed, plant],
            :dormancy => [(:seed => :seed) => 0.5],
            :germinating => [(:seed => :plant) => responsive(ctx -> recruits(ctx.pool.space))],
            :growing => [(:plant => :plant) => responsive(ctx -> slowed(ctx.pool.space)),
                         (:plant => :seed) => responsive(ctx -> (x -> 2.0 * survival(x) * (1 - ctx.pool.space / 2K)))])
        seedbank(net, space) = to_matrix(net, context(net, zeros(51), (space = space,), nothing))
        @test seedbank(tabbed, 700.0) ≈ seedbank(rebuilt, 700.0) rtol = 1e-14
        A, B, C = seedbank(tabbed, 750.0), seedbank(rebuilt, 700.0), seedbank(rebuilt, 800.0)
        @test A[2:end, 1:end] ≈ (B[2:end, 1:end] .+ C[2:end, 1:end]) ./ 2 rtol = 1e-13
        @test A[1, :] ≈ seedbank(rebuilt, 750.0)[1, :] rtol = 1e-14      # seeding is rebuilt, not tabulated
        # the net can be refined, and the tables follow
        two = stratify(table, [:north, :south])
        M2 = to_matrix(two, context(two, zeros(100), (space = 330.0,), nothing))
        @test M2 == kron(Matrix(1.0I, 2, 2), at(table, 330.0))
        # in a species, the pool is read under the species' own name for it
        sp = species(:perennial, table; occupies = [:space => trait_weights(plant, z -> 0.02z^2)])
        n = fill(10.0, 50)
        total = pool_totals(sp, n, nothing).space
        @test community_matrix(sp, n, nothing) ≈ at(exact, total) rtol = 2e-3
    end

    @testset "banded storage" begin
        B = [abs(i - j - 3) <= 2 ? Float64(i + j) : 0.0 for i in 1:30, j in 1:20]
        B[:, 7] .= 0.0
        e = MatrixCommunityModels._Envelope(B)
        @test [e[i, j] for i in 1:30, j in 1:20] == B
        @test length(e.values) == count(!iszero, B)
    end

    @testset "forest trees" begin
        threads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        frozen(sp, stand) = periodic_matrix(sp, (basal_area = stand,), [nothing])
        cells = 100
        spruce = forest_tree(:spruce; cells)
        tabbed = forest_tree(:spruce; cells, tabulate = 0:1:80)
        @test frozen(tabbed, 30.0) ≈ frozen(spruce, 30.0) rtol = 1e-14
        @test frozen(tabbed, 30.5) ≈ (frozen(spruce, 30.0) .+ frozen(spruce, 31.0)) ./ 2 rtol = 1e-12
        @test frozen(tabbed, 95.0) ≈ frozen(spruce, 95.0) rtol = 1e-14      # outside the grid: built for it
        # runs: alone on 100 and 700 cells, and three species on 100 cells
        trees(kinds, cells, tabulate) = [forest_tree(k; cells, tabulate) for k in kinds]
        stand(kinds, cells) = reduce(vcat, [forest_stand(k; cells) for k in kinds])
        for (kinds, cells, years) in (([:spruce], 100, 400), ([:spruce], 700, 100), ([:spruce, :fir, :beech], 100, 300))
            a = project(community(trees(kinds, cells, nothing)...), stand(kinds, cells), fill(nothing, years))
            b = project(community(trees(kinds, cells, 0:1:80)...), stand(kinds, cells), fill(nothing, years))
            difference = maximum(abs.(b.pools[1:years, 1] ./ a.pools[1:years, 1] .- 1))
            println("forest, ", join(kinds, ", "), ", $cells cells, $years years: largest relative difference in basal ",
                    "area between kernels tabulated at whole basal areas and rebuilt each year = ", difference)
            @test difference < 2e-4
        end
        BLAS.set_num_threads(threads)
    end
end
