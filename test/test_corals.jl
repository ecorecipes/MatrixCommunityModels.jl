# Three coral genera of Mo'orea, the open-population IPM of Kayal et al. (2018), compared with the authors' kernel
# functions run in R (test/reference/R/run_coral_reference.R) and with the published simulations (Fig. 3, Fig. S2
# and Fig. S4 of the paper, read from the figures).
using LinearAlgebra

coral_reference(file) = [split(l, ",") for l in readlines(joinpath(@__DIR__, "reference", file))[2:end]]

@testset "corals of Mo'orea (Kayal et al. 2018)" begin
    threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    taxa = (:acropora, :pocillopora, :porites)
    codes = Dict("Acr" => :acropora, "Poc" => :pocillopora, "Por" => :porites)
    reefs = Dict("1" => :lter1, "2" => :lter2, "4" => :lter4, "5" => :lter5)
    number(x) = parse(Float64, x)

    @testset "parameters" begin
        @test collect(keys(CORAL_PARAMETERS)) == [:lter1, :lter2, :lter4, :lter5]
        # Table S1, a few entries
        @test CORAL_PARAMETERS.lter1.pocillopora.recruitment.intercept == 2.688
        @test CORAL_PARAMETERS.lter4.porites.survival.intercept == -1.417
        @test CORAL_PARAMETERS.lter5.acropora.growth.sd == 0.199
        @test all(CORAL_PARAMETERS[r].acropora.recruitment.slope == 0.040 for r in keys(CORAL_PARAMETERS))
        @test CORAL_PARAMETERS.lter2.pocillopora.recruitment.crowding === :assemblage_log_size
        stage = coral_stage()
        @test length(stage) == 100 && step_size(stage) ≈ 0.07 && meshpoints(stage)[1] ≈ 0.035
        @test coral_recruits(:pocillopora, 0.0) ≈ 0.4 * 10^2.688
        @test coral_recruits(:pocillopora, 0.0; reading = :as_written) ≈ exp(2.688)
        @test coral_recruits(:acropora, 1000.0) == CORAL_RECRUITMENT_CAP
        @test_throws ArgumentError coral(:montipora)
        @test_throws ArgumentError coral(:acropora; reef = :lter3)
        @test_throws ArgumentError coral(:acropora; reading = :guessed)
    end

    @testset "kernels and runs against the R reference" begin
        # survival-growth kernels from the authors' P_z1z functions
        worst = 0.0
        kernels = Dict{Tuple{String, String}, Matrix{Float64}}()
        for r in coral_reference("coral_kernels.csv")
            # the kernel does not depend on the state: freeze it at the empty reef, once per species and reef
            K = get!(kernels, (r[1], r[2])) do
                c = coral(codes[r[2]]; reef = reefs[r[1]])
                community_matrix(c, zeros(length(c)), nothing)
            end
            worst = max(worst, abs(K[parse(Int, r[3]), parse(Int, r[4])] - number(r[5])))
        end
        @test worst < 1e-12
        rows = coral_reference("coral_runs.csv")
        for reading in (:calibrated, :as_written), (code, reef) in reefs
            c = coral_community(reef; reading)
            sol = project(c, zeros(length(c)), fill(nothing, 40))
            mine = [r for r in rows if r[1] == String(reading) && r[2] == code]
            worst = 0.0
            for r in mine
                p = pool_totals(c, sol.u[parse(Int, r[3]) + 1, :], nothing)
                k = codes[r[4]]
                theirs = number(r[5])
                theirs > 0 && (worst = max(worst, abs(getproperty(p, Symbol(k, :_colonies)) / theirs - 1),
                                           abs(getproperty(p, Symbol(k, :_log_size)) / number(r[6]) - 1)))
            end
            @test worst < 1e-10
        end
        # the size distributions after four steps (2013)
        states = coral_reference("coral_state.csv")
        c = coral_community(:lter1)
        sol = project(c, zeros(length(c)), fill(nothing, 4))
        for k in taxa
            theirs = [number(r[6]) for r in states if r[1] == "calibrated" && r[2] == "1" && r[3] == "4" && codes[r[4]] == k]
            @test sol[k][end, :] ≈ theirs rtol = 1e-10
        end
    end

    @testset "against the published simulations" begin
        # Composition of the assemblage by summed log size, in %, printed on Fig. 3a-d (2013, step 4) and Fig. S4a-d
        # (2015, step 6), and colonies per 10 m² of Pocillopora read from Fig. S2a (2013)
        fig3 = Dict(:lter1 => [3, 96, 2], :lter2 => [2, 97, 1], :lter4 => [6, 93, 1], :lter5 => [2, 95, 3])
        figS4 = Dict(:lter1 => [4, 93, 3], :lter2 => [3, 96, 1], :lter4 => [6, 93, 0], :lter5 => [2, 94, 4])
        figS2 = Dict(:lter1 => 612, :lter2 => 428, :lter4 => 82, :lter5 => 128)
        shares(c, n) = (p = pool_totals(c, n, nothing); x = [getproperty(p, Symbol(k, :_log_size)) for k in taxa]; 100 .* x ./ sum(x))
        for reef in keys(fig3)
            c = coral_community(reef)
            sol = project(c, zeros(length(c)), fill(nothing, 6))
            @test maximum(abs.(shares(c, sol.u[5, :]) .- fig3[reef])) < 1.5
            @test maximum(abs.(shares(c, sol.u[7, :]) .- figS4[reef])) < 3.5
            @test abs(pool_totals(c, sol.u[5, :], nothing).pocillopora_colonies / figS2[reef] - 1) < 0.3
            # the cap of 348 recruits is not reached in these runs (Appendix S3)
            @test project(coral_community(reef; cap = Inf), zeros(length(c)), fill(nothing, 6)).u ≈ sol.u
            # read as written (natural logarithm, per 10 m²), Acropora and Porites hold far more of the assemblage
            w = coral_community(reef; reading = :as_written)
            @test shares(w, project(w, zeros(length(w)), fill(nothing, 4)).u[end, :])[2] < fig3[reef][2] - 15
        end
    end
    BLAS.set_num_threads(threads)
end
