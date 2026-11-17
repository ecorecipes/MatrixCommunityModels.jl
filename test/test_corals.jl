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
        @test coral_recruits(:pocillopora, 0.0) ≈ exp(2.688)
        @test coral_recruits(:pocillopora, 0.0; reading = :calibrated) ≈ 0.4 * 10^2.688
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
        @test worst < 5e-14                 # achieved 3.6e-16 (2026-11-17)
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
            @test worst < 1e-12             # achieved 7.5e-15 (2026-11-17), 40 years
            # in these runs from empty reef, Acropora recruitment does not reach the cap of 348 within 40 years
            if reading === :calibrated
                most = maximum(coral_recruits(:acropora, pool_totals(c, n, nothing).acropora_colonies; reef, reading) for n in eachrow(sol.u))
                @test most < CORAL_RECRUITMENT_CAP
            end
        end
        # the size distributions after four steps (2013)
        states = coral_reference("coral_state.csv")
        c = coral_community(:lter1; reading = :calibrated)
        sol = project(c, zeros(length(c)), fill(nothing, 4))
        for k in taxa
            theirs = [number(r[6]) for r in states if r[1] == "calibrated" && r[2] == "1" && r[3] == "4" && codes[r[4]] == k]
            @test sol[k][end, :] ≈ theirs rtol = 1e-13   # achieved 7.1e-16 (2026-11-17)
        end
    end

    @testset "against the published simulations" begin
        # Composition of the assemblage by summed log size, in %, printed on Fig. 3a-d (2013, step 4) and Fig. S4a-d
        # (2015, step 6), and colonies per 10 m² of each genus read from Fig. S2a (2013). Vignette 19: every share of
        # 2013 within a point of the published one (allowing for the rounding of both), every share of 2015 within
        # three, Pocillopora colonies within 22% and those of Acropora and Porites within 10 colonies. The bounds are
        # a little above the largest differences found on 2026-11-17, in the comments.
        fig3 = Dict(:lter1 => [3, 96, 2], :lter2 => [2, 97, 1], :lter4 => [6, 93, 1], :lter5 => [2, 95, 3])
        figS4 = Dict(:lter1 => [4, 93, 3], :lter2 => [3, 96, 1], :lter4 => [6, 93, 0], :lter5 => [2, 94, 4])
        figS2 = Dict(:lter1 => [24, 612, 14], :lter2 => [16, 428, 4], :lter4 => [8, 82, 4], :lter5 => [4, 128, 8])
        shares(c, n) = (p = pool_totals(c, n, nothing); x = [getproperty(p, Symbol(k, :_log_size)) for k in taxa]; 100 .* x ./ sum(x))
        colonies(c, n) = (p = pool_totals(c, n, nothing); [getproperty(p, Symbol(k, :_colonies)) for k in taxa])
        for reef in keys(fig3)
            c = coral_community(reef; reading = :calibrated)
            sol = project(c, zeros(length(c)), fill(nothing, 6))
            @test maximum(abs.(shares(c, sol.u[5, :]) .- fig3[reef])) < 1.2          # 1.13
            @test maximum(abs.(shares(c, sol.u[7, :]) .- figS4[reef])) < 3.0         # 2.83
            here = colonies(c, sol.u[5, :])
            @test abs(here[2] / figS2[reef][2] - 1) < 0.22                          # 0.219
            @test maximum(abs.(here[[1, 3]] .- figS2[reef][[1, 3]])) < 10           # 9.5
            # the cap of 348 recruits is not reached in these runs (Appendix S3)
            @test project(coral_community(reef; reading = :calibrated, cap = Inf), zeros(length(c)), fill(nothing, 6)).u ≈ sol.u
            # by default, read as written (natural logarithm, per 10 m²), Acropora and Porites hold far more of the
            # assemblage: Pocillopora 19 to 24 points less than published (18.9 to 24.0)
            w = coral_community(reef)
            @test 18.5 < fig3[reef][2] - shares(w, project(w, zeros(length(w)), fill(nothing, 4)).u[end, :])[2] < 24.5
        end
    end

    @testset "elasticities (Table S2, LTER 1)" begin
        # The elasticity of the ratio of summed log size in 2013 (step 4) to that in 2010 (step 1) to each parameter,
        # by central differences, as in vignette 19. Table S2 of the paper, LTER 1:
        table_s2 = [0.24 0.31 0.11; 0.02 0.05 0.19; 0.82 0.78 0.82; 0.24 0.79 0.57; 0.00 0.00 0.00;
                    0.44 1.56 0.02; 0.29 0.28 0.01; 0.78 0.99 0.61; 0.02 0.01 0.02]
        items = [(:survival, :intercept), (:survival, :slope), (:growth, :intercept), (:growth, :slope), (:growth, :sd),
                 (:recruitment, :intercept), (:recruitment, :slope), (:recruit_size, :mean), (:recruit_size, :sd)]
        P = CORAL_PARAMETERS
        function surface_ratio(taxon, field, which, factor)
            sizes = NamedTuple{taxa}(Tuple(t === taxon && field === :recruit_size ?
                (which === :mean ? (0.6 * factor, 0.2) : (0.6, 0.2 * factor)) : (0.6, 0.2) for t in taxa))
            pars = P
            if field !== :recruit_size
                p = P.lter1[taxon]
                rate = merge(getproperty(p, field), NamedTuple{(which,)}((getproperty(getproperty(p, field), which) * factor,)))
                pars = merge(P, (lter1 = merge(P.lter1, NamedTuple{(taxon,)}((merge(p, NamedTuple{(field,)}((rate,))),))),))
            end
            c = coral_community(:lter1; parameters = pars, recruit_size = sizes, reading = :calibrated)
            u = project(c, zeros(length(c)), fill(nothing, 4)).u
            s(n) = getproperty(pool_totals(c, n, nothing), Symbol(taxon, :_log_size))
            return s(u[end, :]) / s(u[2, :])
        end
        elasticity(t, f, w) = abs((surface_ratio(t, f, w, 1 + 1e-4) - surface_ratio(t, f, w, 1 - 1e-4)) / (2e-4 * surface_ratio(t, f, w, 1.0)))
        E = [elasticity(t, f, w) for (f, w) in items, t in taxa]
        # as vignette 19 prints them, to two decimals
        printed = [0.25 0.28 0.11; 0.02 0.04 0.19; 0.87 0.76 0.82; 0.26 0.68 0.57; 0.00 0.00 0.00;
                   0.09 0.48 0.01; 0.04 0.08 0.00; 0.83 0.78 0.61; 0.01 0.01 0.02]
        @test maximum(abs.(E .- printed)) <= 0.005 + 1e-9
        # survival, growth and the size of recruits within 0.05 of Table S2 (largest difference 0.049), except for
        # Pocillopora's growth slope and recruit mean size (0.68 against 0.79, 0.78 against 0.99)
        close = trues(9, 3)
        close[[6, 7], :] .= false
        close[4, 2] = close[8, 2] = false
        @test maximum(abs.(E .- table_s2)[close]) < 0.05
        # the elasticities to the intercept and slope of recruitment are three to seven times smaller than published
        # (3.3 to 7.2)
        @test all(3 .< table_s2[6:7, :] ./ E[6:7, :] .< 7.5)
    end
    BLAS.set_num_threads(threads)
end
