# Grasslands of five sites, the multispecies IPMs of Chu & Adler (2015) as run by Tredennick et al. (2017),
# compared with the authors' simulation code run in R (test/reference/R/run_grassland_reference.R) and with the
# coefficients printed in Tredennick et al. (2017, Tables S11-S13) for Kansas.
using LinearAlgebra

grassland_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]

@testset "grasslands (Chu & Adler 2015; Tredennick et al. 2017)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    kinds = (:bocu, :bohi, :scsc)
    number(x) = parse(Float64, x)
    cover(c, n) = [getproperty(pool_totals(c, n, nothing), Symbol(k, :_cover)) for k in kinds]

    @testset "parameters" begin
        @test collect(keys(GRASSLAND_PARAMETERS)) == [:arizona, :idaho, :kansas, :montana, :newmexico]
        @test collect(keys(GRASSLAND_PARAMETERS.kansas)) == collect(kinds)
        @test [length(GRASSLAND_PARAMETERS[s][1].recruitment.intercept_year) for s in keys(GRASSLAND_PARAMETERS)] ==
              [18, 22, 34, 13, 32]
        @test [GRASSLAND_PARAMETERS.kansas[k].cells for k in kinds] == [75, 50, 75]
        k = GRASSLAND_PARAMETERS.kansas
        # Tredennick et al. (2017), Tables S11 (growth) and S12 (survival): row species i, column neighbour j;
        # Table S13 (recruitment) is printed the other way round, row = the neighbour
        growth = [-0.0017 0.0011 0.0010; 0.0002 -0.0046 -0.0078; -0.0003 0.0007 -0.0124]
        survival = [-0.1084 0.0006 0.0010; 0.0021 -0.0403 -0.0153; -0.0064 -0.0066 -0.0271]
        recruitment = [-1.1544 -0.0312 -0.0643; 0.3410 -0.8681 0.0041; 0.0619 -0.0336 -0.9611]
        for (i, s) in enumerate(kinds)
            @test all(abs.(k[s].growth.crowding .- growth[i, :]) .<= 5e-5)
            @test all(abs.(k[s].survival.crowding .- survival[i, :]) .<= 5e-5)
            @test all(abs.(k[s].recruitment.crowding .- recruitment[:, i]) .<= 5e-5)
        end
        @test pool_totals(grassland_plant(:kansas, :scsc), grassland_state(:kansas, :scsc; cover = 0.03), nothing).scsc_cover ≈ 0.03
        @test_throws ArgumentError grassland_plant(:utah, :scsc)
        @test_throws ArgumentError grassland_plant(:kansas, :artr)
        @test_throws ArgumentError advance(grassland_plant(:kansas, :bocu), grassland_state(:kansas, :bocu), 35)
        # without interactions only the diagonal coefficients remain
        n = grassland_state(:kansas, :bohi; cover = 0.01)
        alone(inter) = community_matrix(grassland_plant(:kansas, :bohi; interactions = inter), n, nothing)
        @test alone(true) == alone(false)       # a species alone has no neighbours to remove
        @test synchrony([1.0 2.0; 2.0 4.0; 3.0 6.0]) ≈ 1
        @test synchrony([1.0 -1.0; -1.0 1.0; 1.0 -1.0]) ≈ 0 atol = 1e-12
    end

    @testset "runs against the authors' code" begin
        years = [parse(Int, r[2]) for r in grassland_reference("grassland_kansas_years.csv")]
        runs = grassland_reference("grassland_kansas_runs.csv")
        for (label, inter, drawn) in (("mean_30", true, false), ("mean_30_nointer", false, false),
                                      ("years_30", true, true), ("years_30_nointer", false, true))
            c = grassland_community(:kansas; interactions = inter)
            n0 = state_vector(c, (k => grassland_state(:kansas, k) for k in kinds)...)
            rows = [r for r in runs if r[1] == label]
            T = length(rows) - 1
            sol = project(c, n0, drawn ? years[2:(T + 1)] : fill(nothing, T))
            worst = maximum(maximum(abs.(cover(c, sol.u[parse(Int, r[2]) + 1, :]) ./ number.(r[3:5]) .- 1)) for r in rows)
            println("grasslands, Kansas, $label: largest relative difference in cover = ", worst)
            # the authors integrate the crowding by conspecifics with R's integrate at its default tolerance
            @test worst < 1e-7
        end
        # their state after 500 mean years is an equilibrium here too
        states = grassland_reference("grassland_kansas_state.csv")
        for (label, inter) in (("mean_500", true), ("mean_500_nointer", false))
            c = grassland_community(:kansas; interactions = inter)
            n = state_vector(c, (k => [number(r[4]) for r in states if r[1] == label && r[2] == uppercase(String(k))] for k in kinds)...)
            @test maximum(abs, advance(c, n, nothing) .- n) / maximum(n) < 1e-5
            last_row = only(r for r in runs if r[1] == label && r[2] == "500")
            @test cover(c, n) ≈ number.(last_row[3:5]) rtol = 1e-12
        end
    end
    BLAS.set_num_threads(threads)
end
