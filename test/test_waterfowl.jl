# Canvasbacks and redheads (Gibson et al. 2025), checked against the authors' process equations and LTRE script
# (test/reference/R/run_waterfowl_reference.R), against their stored posterior draws, and against Appendix S1.

@testset "Canvasbacks and redheads (Gibson et al. 2025)" begin
    ref = joinpath(@__DIR__, "reference")
    rows = [split(l, ",") for l in readlines(joinpath(ref, "waterfowl_projection.csv"))[2:end]]
    # canvasback years, then redhead years; columns sy_f, asy_f, sy_m, asy_m
    table(label) = permutedims(reduce(hcat, [parse.(Float64, r[4:7]) for r in rows if r[1] == label]))
    both(sol) = vcat(sol[:canvasback], sol[:redhead])
    relerr(a, b) = maximum(abs.(a ./ b .- 1))
    p = WATERFOWL_GIBSON_2025
    env = waterfowl_environment()
    @test length(env) == 59 && env[1].year == 1961 && env[end].canvasback.ss_f_next == p.rates.canvasback.ss_f[end]

    # largest relative difference, with its tolerance: about 100 times the agreement measured on 2026-11-19, not
    # below 1e-13

    # (a) the published model at the posterior means, against the authors' equations evaluated in R
    m = waterfowl_community()
    sol = project(m.community, m.n0, env)
    @test relerr(both(sol), table("posterior_mean")) < 4e-13          # achieved 4.0e-15
    @test isempty(pool_names(m.community))

    # three single draws: rates and the two constants of each draw, against the authors' equations and against
    # the abundances the authors stored for that draw
    dh = Symbol.(split(readline(joinpath(ref, "waterfowl_draws.csv")), ","))
    drows = [split(l, ",") for l in readlines(joinpath(ref, "waterfowl_draws.csv"))[2:end]]
    col(name) = findfirst(==(name), dh)
    for j in ("1", "4000", "8000")
        dj = [r for r in drows if r[1] == j]
        constant(sp, name) = parse(Float64, first(r for r in dj if r[2] == sp)[col(name)])
        q = merge(p, (rates = MatrixCommunityModels._waterfowl_rates(dh, dj),
                      omega = (canvasback = constant("canvasback", :omega), redhead = constant("redhead", :omega)),
                      indirect_harvest = (canvasback = constant("canvasback", :indirect_harvest),
                                          redhead = constant("redhead", :indirect_harvest))))
        stored = table("stored_$j")
        q = merge(q, (n0 = (canvasback = stored[1, :], redhead = stored[61, :]),))
        mj = waterfowl_community(q)
        s = both(project(mj.community, mj.n0, waterfowl_environment(q)))
        @test relerr(s, table("draw_$j")) < 2e-13                    # achieved 1.8e-15 (draw 4000)
        @test relerr(s, stored) < 4e-12                              # achieved 3.4e-14 (draw 8000)
    end

    # (b) the feedback extension: fecundity responds to the pools of breeding females
    mf = waterfowl_community(; feedback = true)
    @test pool_names(mf.community) == [:canvasback_females, :redhead_females]
    solf = project(mf.community, mf.n0, env)
    @test relerr(both(solf), table("feedback")) < 2e-11                # achieved 2.1e-13
    # with both slopes zero the feedback leaves the published model
    flat = merge(p, (beta_density = (canvasback = 0.0, redhead = 0.0), beta_ratio = (canvasback = 0.0, redhead = 0.0)))
    m0 = waterfowl_community(flat; feedback = true)
    @test both(project(m0.community, m0.n0, env)) == both(sol)

    # the authors' transient LTRE (LTRE_Code.R) at the posterior means, against the script evaluated in R
    lrows = [split(l, ",") for l in readlines(joinpath(ref, "waterfowl_ltre.csv"))[2:end]]
    lt = waterfowl_ltre()
    for sp in (:canvasback, :redhead)
        r = getproperty(lt, sp)
        val(q, name) = parse(Float64, first(x for x in lrows if x[1] == q && x[2] == name && x[3] == String(sp))[4])
        @test maximum(abs.(r.sensitivity .- [val("sensitivity", String(n)) for n in r.names])) < 1e-13  # achieved 6.7e-16
        @test maximum(abs.(r.contribution .- [val("contribution", String(n)) for n in r.names])) / r.total < 1e-13  # achieved 9.8e-16
        @test abs(r.total / val("total", "all") - 1) < 2e-13           # achieved 2.1e-15
        @test abs(r.lambda / val("lambda", "total") - 1) < 1e-13       # achieved 0.0
    end

    # printed values: Appendix S1 Table S3 (omega 0.81, 0.80) and Table S7, 60-year means of the female rates
    # (natural winter survival of hatch-years 0.82 and 0.67, of adults 0.88 and 0.80; harvest survival of
    # hatch-years 0.81 and 0.86, of adults 0.92 and 0.96; fecundity of older females 1.24 and 1.35), and
    # lambda 0.998 and 1.005 (Table S7) (geometric mean growth of the posterior-mean breeding numbers, 1961 to 2021)
    mean60(x) = sum(x) / length(x)
    @test round.((p.omega.canvasback, p.omega.redhead); digits = 2) == (0.81, 0.80)
    for (sp, sw_hy, sw_ahy, sh_hy, sh_ahy, f) in ((:canvasback, 0.82, 0.88, 0.81, 0.92, 1.24),
                                                  (:redhead, 0.67, 0.80, 0.86, 0.96, 1.35))
        r = getproperty(p.rates, sp)
        @test round(mean60(r.sw_hy_f); digits = 2) == sw_hy
        @test round(mean60(r.sw_ahy_f); digits = 2) == sw_ahy
        @test round(mean60(1 .- r.kappa_hy_f); digits = 2) == sh_hy
        @test round(mean60(1 .- r.kappa_ahy_f); digits = 2) == sh_ahy
        @test round(mean60(r.fecundity); digits = 2) == f
    end
    @test round(lt.canvasback.lambda; digits = 3) == 0.998
    @test round(lt.redhead.lambda; digits = 3) == 1.005
end
