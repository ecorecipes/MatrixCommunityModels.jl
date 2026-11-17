# Trees as integral projection models sharing the basal area of the stand, compared with the R package
# matreex 0.4.0 (test/reference/R/run_forest_reference.R). The reference on 100 mesh cells was made with
# finer numerical settings than matreex's defaults, so that the two programs can be compared closely; the
# reference on 700 cells, the mesh of matreex, was made with its defaults.
using LinearAlgebra
using Random

forest_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]

@testset "forest trees sharing basal area (matreex)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)        # the eigenvalues of matrices of a few hundred rows are found fastest on one thread
    kinds = [:spruce, :fir, :beech]
    number(x) = parse(Float64, x)
    growth_rate(A) = maximum(abs, eigvals(A))
    # the matrix of a species with the pool held at `stand`; a species without trees has no recruits
    frozen(sp, stand) = periodic_matrix(sp, (basal_area = stand,), [nothing])
    # the rows and columns of the tree cells in the matrix of a species
    trees(sp, kind, cells) = findall(in(stage_names(forest_stage(kind; cells))), stage_names(sp, kind))
    # the basal area a species holds through a projection
    held(sol, kind) = [pool_totals(subcommunity(sol.community, kind), row, nothing).basal_area for row in eachrow(sol[kind])]
    # a reference run: the basal area of `kind` in the years recorded, and where it settled
    function reference_run(file, run, kind)
        rows = [r for r in forest_reference(file) if r[1] == run]
        column = 2 + findfirst(==(kind), kinds)
        settled = [number(r[column]) for r in rows if r[2] == "-1"]
        return (years = [parse(Int, r[2]) for r in rows if r[2] != "-1"],
                basal_area = [number(r[column]) for r in rows if r[2] != "-1"],
                settled = isempty(settled) ? NaN : only(settled))
    end
    reference_state(file, run, kind) = [number(r[4]) for r in forest_reference(file) if r[1] == run && Symbol(r[2]) == kind]
    # the basal area at which a species alone replaces itself: its growth rate there is one (secant method)
    function settles(sp, kind, cells; a = 30.0, b = 45.0)
        excess(x) = growth_rate(community_matrix(sp, forest_stand(kind; cells, basal_area = x), nothing)) - 1
        fa, fb = excess(a), excess(b)
        for _ in 1:30
            c = b - fb * (b - a) / (fb - fa)
            a, fa, b, fb = b, fb, c, excess(c)
            abs(b - a) < 1e-10 && break
        end
        return b
    end

    @testset "fitted vital rates" begin
        @test collect(keys(FOREST_PARAMETERS)) == kinds
        @test [FOREST_PARAMETERS[k].delay for k in kinds] == [14, 28, 48]
        @test [FOREST_PARAMETERS[k].largest for k in kinds] == [1095.0, 1420.0, 1649.0]
        @test FOREST_CLIMATE.sgdd ≈ 1444.67 atol = 0.01
        # every value against the functions that matreex builds from the same coefficients
        worst = 0.0
        for r in forest_reference("forest_vital_rates.csv")
            rates = forest_vital_rates(Symbol(r[1]))
            a, b, value = number(r[3]), number(r[4]), number(r[5])
            got = r[2] == "log_increment" ? rates.log_increment(a, b) : r[2] == "death" ? rates.death(a, b) : rates.recruits(a, b)
            worst = max(worst, abs(got - value) / abs(value))
        end
        println("forest vital rates: largest relative difference from matreex = ", worst)
        @test worst < 1e-11
        # trees grow more slowly in a crowded stand, and a species has fewer recruits among others
        for kind in kinds
            rates = forest_vital_rates(kind)
            @test rates.log_increment(300.0, 40.0) < rates.log_increment(300.0, 0.0)
            @test rates.recruits(10.0, 30.0) < rates.recruits(10.0, 0.0)
            @test 0 < rates.death(300.0, 30.0) < 0.02
        end
        # the climate matters: a warmer, drier site
        @test forest_vital_rates(:spruce; climate = (sgdd = 2185.0, wai = 0.03)).log_increment(300.0, 20.0) !=
              forest_vital_rates(:spruce).log_increment(300.0, 20.0)
        @test_throws ArgumentError forest_vital_rates(:oak)
        odd = (spruce = merge(FOREST_PARAMETERS.spruce, (growth = ["size:BATOTcomp" => 1.0],)),)
        @test_throws ArgumentError forest_vital_rates(:spruce; parameters = odd)
        # the normal distribution function, against R's pnorm
        for (x, p) in ((-8.5, 9.4795348222033192e-18), (-6.0, 9.8658764503769809e-10), (-3.3, 0.00048342414238377744),
                       (-1.0, 0.15865525393145705), (-0.25, 0.4012936743170763), (0.0, 0.5), (0.6, 0.72574688224992645),
                       (2.0, 0.97724986805182079), (4.2, 0.9999866542509841), (7.5, 0.99999999999996814))
            @test abs(MatrixCommunityModels._normal_cdf(x) - p) < 2e-16
        end
    end

    @testset "a species as a net" begin
        cells = 50
        spruce = forest_tree(:spruce; cells)
        names = stage_names(spruce, :spruce)
        @test length(names) == 14 + cells == length(spruce)
        @test names[1:3] == [:seedling, :sapling1, :sapling2] && names[14:15] == [:sapling13, :tree1]
        @test names[15:end] == stage_names(forest_stage(:spruce; cells))
        @test pool_names(spruce) == [:basal_area]
        net = system(spruce)[:spruce]
        @test Set(transition_names(net)) == Set([:growing, :recruiting, :establishing, :ingrowth, :sapling_ageing])
        # the mesh of matreex: 90 mm to 1.1 times the largest tree in the data
        fir = forest_stage(:fir)
        @test length(fir) == 700 && step_size(fir) ≈ (1.1 * 1420 - 90) / 700
        @test meshpoints(fir)[1] ≈ 90 + step_size(fir) / 2

        # a young stand, and its first year
        stand = forest_stand(:spruce; cells, basal_area = 12.5)
        @test length(stand) == length(spruce) && all(==(0), stand[1:14])
        @test pool_totals(spruce, stand, nothing).basal_area ≈ 12.5
        @test all(d -> 90 <= d < 150, meshpoints(forest_stage(:spruce; cells))[stand[15:end] .> 0])
        next = advance(spruce, stand, nothing)
        recruits = forest_vital_rates(:spruce).recruits(12.5, 0.0)
        @test next[1] ≈ recruits / 2 && next[2] ≈ recruits / 2          # half of them a year ahead
        @test all(==(0), next[3:14])
        @test 0.98 * sum(stand) < sum(next[15:end]) < sum(stand)        # some trees die or are felled, none appear
        @test all(==(0), advance(spruce, zeros(length(spruce)), nothing))
        # recruits are not lost on the way: after 13 years all of them so far are still seedlings or saplings,
        # and a year later the first have become trees
        sol = project(spruce, stand, fill(nothing, 15))
        waiting = stage_total(sol, :spruce, [:seedling, :sapling])
        recruited = cumsum([forest_vital_rates(:spruce).recruits(sol.pools[year, 1], 0.0) for year in 1:14])
        @test waiting[14] ≈ recruited[13]
        @test waiting[15] ≈ recruited[14] - recruited[1] / 2
        # trees do not shrink, and survival is a probability
        K = frozen(spruce, 20.0)
        idx = trees(spruce, :spruce, cells)
        @test istril(K[idx, idx]) && all(>=(0), K)
        U, F = survival_fecundity_matrices(net, context(net, stand, (basal_area = 12.5,), nothing))
        @test maximum(sum(U; dims = 1)) <= 1
        @test sum(F) ≈ recruits / 12.5 * sum(pi .* (meshpoints(forest_stage(:spruce; cells)) ./ 2000) .^ 2)
        # growth is integrated over the cell: the mean growth of a tree is that of the log-normal distribution
        rates = forest_vital_rates(:spruce)
        fine = forest_tree(:spruce; cells = 400, harvest = 0.0)
        z = meshpoints(forest_stage(:spruce; cells = 400))
        column = frozen(fine, 20.0)[trees(fine, :spruce, 400), trees(fine, :spruce, 400)][:, 100]
        @test sum(column) ≈ 1 - rates.death(z[100], 20.0) rtol = 2e-8         # achieved 8.4e-9 (2026-11-17), growth cut off at 6 sd
        @test sum((z .- z[100]) .* column) / sum(column) ≈ exp(rates.log_increment(z[100], 20.0) + rates.sigma^2 / 2) rtol = 2.5e-4  # achieved 1.0e-4

        # options
        short = forest_tree(:spruce; cells, delay = 2)
        @test stage_names(short, :spruce)[1:3] == [:seedling, :sapling, :tree1]
        @test length(forest_stand(:spruce; cells, delay = 2)) == length(short)
        @test forest_tree(:beech; cells, name = :fagus) |> species_names == [:fagus]
        unfelled = forest_tree(:spruce; cells, harvest = 0.0)
        @test sum(advance(unfelled, stand, nothing)[15:end]) ≈ sum(next[15:end]) / (1 - 0.006)
        @test_throws ArgumentError forest_tree(:oak)
        @test_throws ArgumentError forest_tree(:spruce; cells, delay = 1)
        @test_throws ArgumentError forest_tree(:spruce; cells, harvest = 1.5)
        @test_throws ArgumentError forest_tree(:spruce; cells, conspecific = 0.0)
        @test_throws ArgumentError forest_stand(:spruce; cells, diameters = (10, 20))
    end

    @testset "kernels against matreex" begin
        # (cells, species): the entries matreex integrates by Gauss-Legendre quadrature, and the share of the trees
        # of a cell that stay on the mesh alive. With its default settings (700 cells) matreex uses a coarse rule for
        # growth above 50 mm a year, which is common only in an empty stand: the shares differ there. Achieved
        # (2026-11-17): 100 cells, entries 3.3e-10 and shares 2.9e-8; 700 cells, entries 2.0e-14 (spruce), 1.5e-11 (fir)
        # and 8.6e-10 (beech), shares 2.2e-4 in an empty stand and 7.8e-6 at 30 m²/ha. The differences in the entries are
        # matreex's quadrature, not rounding.
        entry700 = Dict(:spruce => 1e-12, :fir => 5e-11, :beech => 2e-9)
        for (cells, checked, entry, share) in ((100, kinds, k -> 7e-10, ba -> 6e-8), (700, kinds, k -> entry700[k], ba -> ba == 0 ? 3e-4 : 1.2e-5))
            ref = forest_reference("forest_kernel_$cells.csv")
            for kind in checked
                sp = forest_tree(kind; cells, harvest = 0.0)
                idx = trees(sp, kind, cells)
                for ba in (0.0, 30.0)
                    K = frozen(sp, ba)[idx, idx]
                    rows = [r for r in ref if Symbol(r[1]) == kind && number(r[2]) == ba]
                    entries = maximum(abs(K[parse(Int, r[3]), parse(Int, r[4])] - number(r[5])) for r in rows if r[3] != "0")
                    shares = maximum(abs(sum(K[:, parse(Int, r[4])]) - number(r[5])) for r in rows if r[3] == "0")
                    println("forest kernel, $kind, $cells cells, basal area $ba: largest difference from matreex in ",
                            count(r -> r[3] != "0", rows), " entries = ", entries, "; in the share kept = ", shares)
                    @test entries < entry(kind)
                    @test shares < share(ba)
                end
            end
        end
    end

    @testset "one species (100 cells)" begin
        cells = 100
        runs, states = "forest_runs_100.csv", "forest_states_100.csv"
        spruce = forest_tree(:spruce; cells)
        ref = reference_run(runs, "spruce", :spruce)
        sol = project(spruce, forest_stand(:spruce; cells), fill(nothing, maximum(ref.years) + 1))
        mine = sol.pools[ref.years .+ 1, 1]
        println("forest, spruce alone, 100 cells: largest relative difference from matreex over ", maximum(ref.years),
                " years = ", maximum(abs.(mine ./ ref.basal_area .- 1)))
        @test maximum(abs.(mine ./ ref.basal_area .- 1)) < 1.2e-5       # achieved 6.0e-6 (2026-11-17), matreex's interpolation
        @test held(sol, :spruce)[1:(end - 1)] ≈ sol.pools[1:(end - 1), 1]
        @test maximum(mine) > 1.2 * ref.settled                     # the stand overshoots, then cycles slowly
        # What remains is matreex's interpolation: it tabulates its kernels (here every 0.25 m² per hectare) and
        # uses a mixture of the two nearest tables. With the kernels of this model used in the same way, and the
        # recruits added to the first two classes as matreex does, the run is that of matreex.
        tables = Dict{Float64, Matrix{Float64}}()
        tabulated(b) = get!(() -> frozen(spruce, b), tables, b)      # growth, ageing and felling in a stand of b
        recruits = forest_vital_rates(:spruce).recruits
        n = forest_stand(:spruce; cells)
        interpolated = Float64[]
        for year in 0:maximum(ref.years)
            b = pool_totals(spruce, n, nothing).basal_area
            push!(interpolated, b)
            lower = 0.25 * floor(b / 0.25)
            share = (b - lower) / 0.25
            n = (1 - share) .* (tabulated(lower) * n) .+ share .* (tabulated(lower + 0.25) * n)
            n[1:2] .+= recruits(b, 0.0) / 2
        end
        println("forest, spruce alone, 100 cells, kernels interpolated as in matreex: largest relative difference = ",
                maximum(abs.(interpolated[ref.years .+ 1] ./ ref.basal_area .- 1)))
        @test maximum(abs.(interpolated[ref.years .+ 1] ./ ref.basal_area .- 1)) < 5e-8   # achieved 2.5e-8 (2026-11-17)
        # the equilibrium matreex reaches is one of this model, and the nearest
        there = reference_state(states, "spruce", :spruce)
        @test pool_totals(spruce, there, nothing).basal_area ≈ ref.settled rtol = 1e-10            # achieved 1.4e-12 (2026-11-17)
        @test maximum(abs, advance(spruce, there, nothing) .- there) < 6e-8 * maximum(there)    # achieved 2.7e-8
        eq = equilibrium(spruce, there, nothing)
        @test eq.converged
        @test pool_totals(spruce, eq.n, nothing).basal_area ≈ ref.settled rtol = 2.5e-7           # achieved 1.0e-7
        @test maximum(abs, eq.n .- there) < 2e-7 * maximum(there)   # the same size distribution; achieved 9.0e-8
        @test growth_rate(community_matrix(spruce, eq.n, nothing)) ≈ 1 atol = 3e-9                # achieved 2.6e-11
        # the settled basal area without iterating: where the growth rate of the frozen species is one
        for kind in kinds
            sp = kind == :spruce ? spruce : forest_tree(kind; cells)
            @test settles(sp, kind, cells) ≈ reference_run(runs, String(kind), kind).settled rtol = 5e-7  # achieved 2.4e-7 (2026-11-17)
        end
        # without felling the stand holds more, and matreex's run is followed for 3000 years of cycles
        free = reference_run(runs, "unharvested", :spruce)
        unfelled = forest_tree(:spruce; cells, harvest = 0.0)
        long = project(unfelled, forest_stand(:spruce; cells), fill(nothing, maximum(free.years) + 1))
        @test maximum(abs.(long.pools[free.years .+ 1, 1] ./ free.basal_area .- 1)) < 5e-5   # achieved 2.6e-5 (2026-11-17)
        @test minimum(long.pools[2001:2991, 1]) < 40 < 50 < maximum(long.pools[2001:2991, 1])
    end

    @testset "two and three species (100 cells)" begin
        cells = 100
        runs, states = "forest_runs_100.csv", "forest_states_100.csv"
        spruce, fir, beech = (forest_tree(k; cells) for k in kinds)
        for (run, members) in (("pair", [:spruce, :fir]), ("trio", kinds))
            forest = community((forest_tree(k; cells) for k in members)...)
            @test species_names(forest) == members && pool_names(forest) == [:basal_area]
            years = reference_run(runs, run, :spruce).years
            sol = project(forest, reduce(vcat, [forest_stand(k; cells) for k in members]), fill(nothing, maximum(years) + 1))
            for kind in members
                ref = reference_run(runs, run, kind)
                mine = held(sol, kind)[years .+ 1]
                println("forest, $run, 100 cells, $kind: largest relative difference from matreex over ", maximum(years),
                        " years = ", maximum(abs.(mine ./ ref.basal_area .- 1)))
                @test maximum(abs.(mine ./ ref.basal_area .- 1)) < 1.4e-5   # achieved 6.7e-6 (2026-11-17)
            end
            @test sol.pools[1:(end - 1), 1] ≈ sum(held(sol, k) for k in members)[1:(end - 1)]
        end
        # the equilibrium of the three that matreex reaches is one of this model
        trio = community(spruce, fir, beech)
        there = reduce(vcat, [reference_state(states, "trio", k) for k in kinds])
        @test length(there) == length(trio)
        @test maximum(abs, advance(trio, there, nothing) .- there) < 5e-7 * maximum(there)    # achieved 2.1e-7 (2026-11-17)
        totals = [reference_run(runs, "trio", k).settled for k in kinds]
        @test pool_totals(trio, there, nothing).basal_area ≈ sum(totals) rtol = 3e-12          # achieved 2.1e-14
        @test totals[2] > totals[1] > totals[3] > 5                 # fir, then spruce, then beech: all three persist
        # each species frozen there just replaces itself
        A = community_matrix(trio, there, nothing)
        edges = cumsum([0; [length(sp) for sp in (spruce, fir, beech)]])
        for j in 1:3
            block = (edges[j] + 1):edges[j + 1]
            @test growth_rate(A[block, block]) ≈ 1 atol = 3e-8                                # achieved 1.3e-8 (2026-11-17)
        end
    end

    @testset "rare invaders (100 cells)" begin
        cells = 100
        # growth rate and net reproductive rate of a species holding 0.1 m² per hectare of a stand, against the same
        # quantities from matreex's kernels and recruitment functions
        for r in forest_reference("forest_invasion_100.csv")
            kind, stand = Symbol(r[1]), number(r[2])
            invader = forest_tree(kind; cells, conspecific = number(r[3]))
            @test invasion_growth_rate(invader, (basal_area = stand,), nothing) ≈ number(r[4]) rtol = 1e-9   # achieved 1.1e-10 (2026-11-17)
            net = system(invader)[kind]
            U, F = survival_fecundity_matrices(net, context(net, zeros(length(invader)), (basal_area = stand,), nothing))
            @test growth_rate(F / (I - U)) ≈ number(r[5]) rtol = 1e-7                                       # achieved 4.3e-8
        end
        spruce = forest_tree(:spruce; cells)
        there = reference_state("forest_states_100.csv", "spruce", :spruce)
        stand = pool_totals(spruce, there, nothing).basal_area
        # a species without trees has no recruits, so the rarity of an invader has to be stated
        @test invasion_growth_rate(forest_tree(:fir; cells), spruce, there, nothing) < 1
        rates = [invasion_growth_rate(forest_tree(:fir; cells, conspecific = c), spruce, there, nothing) for c in (0.01, 0.1, 1.0)]
        @test rates[1] > rates[2] > rates[3] > 1                     # the rarer, the faster; fir invades spruce
        @test invasion_growth_rate(forest_tree(:beech; cells, conspecific = 0.1), spruce, there, nothing) > 1
        # a species that holds the whole stand at its equilibrium just replaces itself
        @test invasion_growth_rate(forest_tree(:spruce; cells, conspecific = stand), (basal_area = stand,), nothing) ≈ 1 atol = 1e-8  # achieved 2.7e-9
        # an invader does invade: fir added to the spruce stand increases
        pair = community(spruce, forest_tree(:fir; cells))
        sol = project(pair, [there; forest_stand(:fir; cells, basal_area = 0.1)], fill(nothing, 150))
        @test held(sol, :fir)[end] > 10 * held(sol, :fir)[1]
    end

    @testset "the mesh" begin
        # the settled basal area of spruce alone on coarser meshes, against matreex on 700 cells
        fine = reference_run("forest_runs_700.csv", "spruce", :spruce).settled
        found = [settles(forest_tree(:spruce; cells), :spruce, cells) for cells in (100, 200, 400)]
        println("forest, spruce alone settles at ", found, " m²/ha on 100, 200 and 400 cells; matreex on 700 cells: ", fine)
        @test abs(found[1] / fine - 1) < 0.025      # achieved 0.0229 (2026-11-17)
        @test abs(found[2] / fine - 1) < 0.007      # achieved 0.0063
        @test abs(found[3] / fine - 1) < 5e-4       # achieved 4.0e-4
        @test found[1:2] ≈ [37.51, 36.90] atol = 0.005                     # as vignette 15 prints them
        @test abs(found[3] - fine) < abs(found[2] - fine) < abs(found[1] - fine)
    end

    @testset "the mesh of matreex (700 cells)" begin
        cells = 700
        spruce = forest_tree(:spruce; cells)
        @test length(spruce) == 714
        ref = reference_run("forest_runs_700.csv", "spruce", :spruce)
        # the first century of matreex's run, with its default settings: it interpolates its kernels between whole
        # basal areas, which shows most while the young stand grows fast
        years = ref.years[ref.years .<= 100]
        sol = project(spruce, forest_stand(:spruce; cells), fill(nothing, 101))
        println("forest, spruce alone, 700 cells: largest relative difference from matreex over 100 years = ",
                maximum(abs.(sol.pools[years .+ 1, 1] ./ ref.basal_area[1:length(years)] .- 1)))
        @test maximum(abs.(sol.pools[years .+ 1, 1] ./ ref.basal_area[1:length(years)] .- 1)) < 7e-4   # achieved 4.8e-4 (2026-11-17)
        # matreex's equilibrium is all but an equilibrium here: its kernels are interpolated between whole basal areas
        there = reference_state("forest_states_700.csv", "spruce", :spruce)
        @test pool_totals(spruce, there, nothing).basal_area ≈ ref.settled rtol = 1e-10            # achieved 1.2e-12
        @test maximum(abs, advance(spruce, there, nothing) .- there) < 7e-5 * maximum(there)    # achieved 3.3e-5
        # and the three species together
        trio = community(spruce, forest_tree(:fir; cells), forest_tree(:beech; cells))
        @test length(trio) == 714 + 728 + 748
        all3 = reduce(vcat, [reference_state("forest_states_700.csv", "trio", k) for k in kinds])
        @test maximum(abs, advance(trio, all3, nothing) .- all3) < 1.5e-5 * maximum(all3)     # achieved 6.3e-6 (2026-11-17)
    end

    @testset "whole trees" begin
        cells = 50
        pair = community(forest_tree(:spruce; cells), forest_tree(:fir; cells))
        start = round.([forest_stand(:spruce; cells, basal_area = 15.0); forest_stand(:fir; cells, basal_area = 15.0)])
        random = project(pair, start, fill(nothing, 30); rng = MersenneTwister(15))
        @test all(isinteger, random.u) && all(>=(0), random.u)
        @test stage_total(random, :spruce, :seedling)[end] + stage_total(random, :spruce, :sapling)[end] > 0
        @test 20 < random.pools[30, 1] < 60
    end
    BLAS.set_num_threads(threads)
end
