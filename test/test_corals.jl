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

    @testset "from the observed assemblages (Fig. 4, Tables S3 and S4)" begin
        # The colonies mapped by the MCR LTER on the fore reefs in 2013 and 2015 (data/corals/mcr_lter_observed_colonies.csv,
        # from EDI knb-lter-mcr.4009.5), as colonies per 10 m²: Acropora mapped over 20 m², the other two over 4 m².
        counts = Dict((:lter1, 2013) => [41, 289, 132], (:lter2, 2013) => [29, 189, 42], (:lter4, 2013) => [19, 87, 15],
                      (:lter5, 2013) => [14, 73, 22], (:lter1, 2015) => [50, 215, 103], (:lter2, 2015) => [32, 128, 38],
                      (:lter4, 2015) => [41, 70, 7], (:lter5, 2015) => [48, 87, 23])
        for ((reef, year), n) in counts
            s = coral_observed_state(reef, year)
            @test first.(s) == collect(taxa) && all(length(last(p)) == 100 && all(>=(0), last(p)) for p in s)
            @test sum.(last.(s)) ≈ n .* [0.5, 2.5, 2.5] rtol = 1e-14
        end
        s = coral_observed_state(:lter1, 2013)
        last(s[1])[1] = 1e6
        @test last(coral_observed_state(:lter1, 2013)[1])[1] == 0           # a copy
        @test_throws ArgumentError coral_observed_state(:lter1, 2014)
        @test_throws ArgumentError coral_observed_state(:lter3, 2013)

        # Shares of the summed log size printed on Fig. 3e-h (observed 2013) and Fig. 4e-h (observed 2015), Acr Poc Por %
        fig3obs = Dict(:lter1 => [2, 67, 30], :lter2 => [3, 78, 19], :lter4 => [5, 78, 17], :lter5 => [4, 75, 21])
        fig4obs = Dict(:lter1 => [3, 68, 29], :lter2 => [4, 75, 20], :lter4 => [11, 83, 6], :lter5 => [9, 70, 21])
        fig4pred = Dict(:lter1 => [3, 67, 30], :lter2 => [4, 82, 14], :lter4 => [5, 87, 7], :lter5 => [3, 77, 20])
        fig4proj = Dict(:lter1 => [27, 50, 24], :lter2 => [11, 74, 15], :lter4 => [14, 84, 2], :lter5 => [8, 63, 29])
        shares(c, n) = (p = pool_totals(c, n, nothing); x = [getproperty(p, Symbol(k, :_log_size)) for k in taxa]; 100 .* x ./ sum(x))
        c0 = coral_community(:lter1)
        observed_share(reef, year) = shares(c0, state_vector(c0, coral_observed_state(reef, year)))
        # 2013: every share within 1.6 points (LTER 4 Pocillopora 76.4 against 78); 2015: within 1 point (0.93) except
        # LTER 4, where the data give Pocillopora 4.2 points less and Porites 3.6 more than printed
        @test maximum(maximum(abs.(observed_share(r, 2013) .- fig3obs[r])) for r in keys(fig3obs)) < 1.7      # 1.59
        @test maximum(maximum(abs.(observed_share(r, 2015) .- fig4obs[r])) for r in (:lter1, :lter2, :lter5)) < 1.0   # 0.93
        @test -4.3 < observed_share(:lter4, 2015)[2] - fig4obs[:lter4][2] < -4.1                             # -4.22

        # The updated runs: Acropora recruits of mean 1.2 and sd 0.4, from 2013 two steps to 2015, from 2015 five to 2020
        size4 = (acropora = (1.2, 0.4), pocillopora = (0.6, 0.2), porites = (0.6, 0.2))
        rows = coral_reference("coral_observed_runs.csv")
        for reading in (:as_written, :calibrated)
            worst, d15, d20, most = 0.0, 0.0, Dict{Symbol, Vector{Float64}}(), 0.0
            for (code, reef) in reefs
                c = coral_community(reef; reading, recruit_size = size4)
                for (start, steps) in ((2013, 2), (2015, 5))
                    u = project(c, coral_observed_state(reef, start), fill(nothing, steps)).u
                    for r in rows
                        (r[1] == String(reading) && r[2] == code && parse(Int, r[3]) == start) || continue
                        p = pool_totals(c, u[parse(Int, r[4]) + 1, :], nothing)
                        k = codes[r[5]]
                        worst = max(worst, abs(getproperty(p, Symbol(k, :_colonies)) / number(r[6]) - 1),
                                    abs(getproperty(p, Symbol(k, :_log_size)) / number(r[7]) - 1))
                    end
                    if start == 2013
                        d15 = max(d15, maximum(abs.(shares(c, u[end, :]) .- fig4pred[reef])))
                    else
                        d20[reef] = shares(c, u[end, :]) .- fig4proj[reef]
                        most = max(most, maximum(coral_recruits(:acropora, pool_totals(c, n, nothing).acropora_colonies;
                                                                reef, reading, cap = Inf) for n in eachrow(u)))
                    end
                end
            end
            @test worst < 2e-13             # against the R reference; achieved 1.6e-15 and 1.8e-15 (2026-11-18)
            # Acropora recruitment stays below the cap of 348 everywhere until 2020 (the paper: reached on the north shore)
            @test most < CORAL_RECRUITMENT_CAP
            if reading === :as_written
                @test d15 < 1.9                                         # Fig. 4a-d within 1.83 points
                @test maximum(abs.(d20[:lter2])) < 3.4 && maximum(abs.(d20[:lter4])) < 1.8   # Fig. 4j, 4k: 3.36, 1.78
                @test -19.9 < d20[:lter1][1] < -19.7                    # Fig. 4i: Acropora 7.2% against 27% (-19.79)
                @test 150 < most < 160                                  # largest Acropora recruitment, LTER 1: 158.3
            else
                @test d15 < 3.0                                         # Fig. 4a-d within 2.97 points
                @test -24.1 < d20[:lter1][1] < -23.9                    # Fig. 4i: Acropora 3.0% against 27% (-23.99)
                @test most < 6                                          # 5.5
            end
        end

        # Elasticities of the summed log size of 2015 (Table S3) and of 2020 (Table S4) to each parameter, LTER 1, as written
        items = [(:survival, :intercept), (:survival, :slope), (:growth, :intercept), (:growth, :slope), (:growth, :sd),
                 (:recruitment, :intercept), (:recruitment, :slope), (:recruit_size, :mean), (:recruit_size, :sd)]
        P = CORAL_PARAMETERS
        function final_size(taxon, field, which, factor, start, steps)
            sizes = size4
            if field === :recruit_size
                m, sd = size4[taxon]
                sizes = merge(size4, NamedTuple{(taxon,)}(((which === :mean ? (m * factor, sd) : (m, sd * factor)),)))
            end
            pars = P
            if field !== :recruit_size
                p = P.lter1[taxon]
                rate = merge(getproperty(p, field), NamedTuple{(which,)}((getproperty(getproperty(p, field), which) * factor,)))
                pars = merge(P, (lter1 = merge(P.lter1, NamedTuple{(taxon,)}((merge(p, NamedTuple{(field,)}((rate,))),))),))
            end
            c = coral_community(:lter1; parameters = pars, recruit_size = sizes)
            u = project(c, coral_observed_state(:lter1, start), fill(nothing, steps)).u
            return getproperty(pool_totals(c, u[end, :], nothing), Symbol(taxon, :_log_size))
        end
        elasticity(t, f, w, start, steps) = abs((final_size(t, f, w, 1 + 1e-4, start, steps) - final_size(t, f, w, 1 - 1e-4, start, steps)) /
                                                (2e-4 * final_size(t, f, w, 1.0, start, steps)))
        table_s3 = [0.17 0.22 0.03; 0.03 0.06 0.10; 0.68 0.67 0.74; 0.29 0.98 0.83; 0.00 0.00 0.00;
                    0.82 0.05 0.05; 0.53 0.03 0.03; 0.28 0.00 0.02; 0.00 0.00 0.00]
        table_s4 = [0.50 0.41 0.04; 0.08 0.15 0.14; 0.49 0.91 0.97; 0.17 2.08 1.45; 0.00 0.00 0.00;
                    3.70 0.13 0.16; 4.35 0.10 0.08; 0.58 0.00 0.03; 0.03 0.00 0.00]
        E3 = [elasticity(t, f, w, 2013, 2) for (f, w) in items, t in taxa]
        E4 = [elasticity(t, f, w, 2015, 5) for (f, w) in items, t in taxa]
        # as vignette 19 prints them, to two decimals
        printed3 = [0.18 0.23 0.03; 0.03 0.06 0.11; 0.76 0.67 0.75; 0.34 0.98 0.83; 0.00 0.00 0.00;
                    0.33 0.00 0.02; 0.31 0.00 0.02; 0.18 0.00 0.01; 0.00 0.00 0.00]
        printed4 = [0.70 0.42 0.03; 0.12 0.16 0.13; 0.76 0.93 0.98; 0.33 2.13 1.48; 0.00 0.00 0.00;
                    1.96 0.00 0.06; 3.19 0.00 0.05; 0.34 0.00 0.02; 0.02 0.00 0.00]
        @test maximum(abs.(E3 .- printed3)) <= 0.005 + 1e-9 && maximum(abs.(E4 .- printed4)) <= 0.005 + 1e-9
        rest = [1, 2, 3, 4, 5, 8, 9]                                    # all but the intercept and slope of recruitment
        # Table S3: survival, growth and recruit size within 0.11 (0.10, Acropora recruit mean size 0.18 against 0.28)
        @test maximum(abs.(E3 .- table_s3)[rest, :]) < 0.11
        # Table S4: within 0.06 for Pocillopora and Porites (0.052), within 0.28 for Acropora (0.27, growth intercept)
        @test maximum(abs.(E4 .- table_s4)[rest, 2:3]) < 0.06 && maximum(abs.(E4 .- table_s4)[rest, 1]) < 0.28
        # recruitment: Acropora's elasticities are 1.4 to 2.5 times smaller than published (0.33, 0.31 against 0.82,
        # 0.53 in 2013-2015; 1.96, 3.19 against 3.70, 4.35 in 2015-2020), but its slope is, as published, the largest
        # elasticity of all in 2015-2020; Pocillopora's are near zero (0.001 to 0.005 against 0.03 to 0.13)
        @test all(1.3 .< table_s3[6:7, 1] ./ E3[6:7, 1] .< 2.6) && all(1.3 .< table_s4[6:7, 1] ./ E4[6:7, 1] .< 2.6)
        @test argmax(E4) == argmax(table_s4) == CartesianIndex(7, 1)
        @test maximum(E3[6:7, 2]) < 0.01 && maximum(E4[6:7, 2]) < 0.01
    end
    BLAS.set_num_threads(threads)
end
