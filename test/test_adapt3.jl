# The two-species Cypripedium community of adapt3's comproj vignette, rebuilt by
# composition and compared with project3() (test/reference/R/run_adapt3_reference.R).

const A3REF = joinpath(@__DIR__, "reference")

a3csv(path) = [split(replace(l, "\"" => ""), ",") for l in readlines(path)[2:end]]

function adapt3_matrices()
    mats = Dict{Int,Dict{Int,Matrix{Float64}}}()
    for row in a3csv(joinpath(A3REF, "adapt3_matrices.csv"))
        sp, k, i, j = parse(Int, row[1]), parse(Int, row[2]), parse(Int, row[4]), parse(Int, row[5])
        M = get!(get!(mats, sp, Dict{Int,Matrix{Float64}}()), k, zeros(9, 9))
        M[i, j] = parse(Float64, row[6])
    end
    return [[mats[sp][k] for k in 1:length(mats[sp])] for sp in 1:2]
end

# One species: every nonzero entry of the annual matrices is a named transition whose rate
# looks up the matrix of the current year; Ricker density dependence on reproduction into
# SD (alpha 0.5, beta 1) and P1 (alpha 1.2, beta 2).
function cypripedium(name, mats, stages, weights; integeronly = false)
    reproductive = (:XSm, :Sm, :Md, :Lg, :XLg)
    year(ctx) = mod(ctx.env, length(mats)) + 1
    transitions = Transition[]
    for (j, from) in enumerate(stages), (i, to) in enumerate(stages)
        any(M[i, j] != 0 for M in mats) || continue
        if to === :SD && from in reproductive
            rate = ctx -> mats[year(ctx)][i, j] * 0.5 * exp(-1.0 * ctx.pool.density)
            push!(transitions, fecundity(Symbol(from, :_to_, to), from => to, rate))
        elseif to === :P1 && from in reproductive
            rate = ctx -> mats[year(ctx)][i, j] * 1.2 * exp(-2.0 * ctx.pool.density)
            push!(transitions, survival(Symbol(from, :_to_, to), from => to, rate))
        else
            push!(transitions, survival(Symbol(from, :_to_, to), from => to, ctx -> mats[year(ctx)][i, j]))
        end
    end
    return Species(name, stages, transitions; occupancy = [:density => weights],
                   threshold = integeronly ? (n -> floor.(n)) : nothing)
end

function adapt3_reference(tag)
    rows = a3csv(joinpath(A3REF, "adapt3_$(tag)_N.csv"))
    return permutedims(reduce(hcat, [parse.(Float64, r[2:4]) for r in rows]))
end

@testset "adapt3 project3 (Cypripedium)" begin
    stages = Symbol.(readlines(joinpath(A3REF, "adapt3_stages.txt")))
    mats = adapt3_matrices()
    @test length.(mats) == [5, 17]
    n0 = [:candidum => [500.0, 100, 200, 0, 0, 0, 0, 0, 0], :parviflorum => [5000.0, 1000, 2000, 0, 0, 0, 0, 0, 0]]
    nonseed = [0.0; ones(8)]

    for (tag, weights, integeronly) in (("continuous", nonseed, false), ("noequiv", ones(9), false),
                                        ("vignette", nonseed, true))
        com = community(cypripedium(:candidum, mats[1], stages, weights; integeronly),
                        cypripedium(:parviflorum, mats[2], stages, weights; integeronly))
        sol = project(com, n0, 0:99)
        ref = adapt3_reference(tag)
        got = hcat(sum(sol[:candidum]; dims = 2), sum(sol[:parviflorum]; dims = 2))
        err = maximum(abs.(got .- ref[:, 1:2]) ./ max.(abs.(ref[:, 1:2]), 1.0))
        println("adapt3 $tag: max scaled difference from project3() = ", err)
        @test err < 1e-10
        # aggregate density used for each step
        @test sol.pools[1:100, 1] ≈ ref[1:100, 3] rtol = 1e-10
    end
end
