# Analysis of composed communities: equilibria, local stability, invasion, and the
# cross-species sensitivity analysis of Lytle & Tonkin (2023).

"""
    equilibrium(c, n0, env; maxiter = 100_000, tol = 1e-10)

Iterate a closed community in a constant environment `env` from `n0` until the largest
change in any state, relative to `max(|n|, 1)`, falls below `tol`. Returns
`(n, converged, iterations)`. A stable equilibrium is found only if the iteration
converges to one; cycles and chaos return `converged = false`.
"""
function equilibrium(c::OpenCommunity, n0::AbstractVector{<:Real}, env; maxiter::Integer = 100_000, tol::Real = 1e-10)
    n = collect(Float64, n0)
    for k in 1:maxiter
        n1 = advance(c, n, env)
        if maximum(abs.(n1 .- n) ./ max.(abs.(n), 1.0)) < tol
            return (n = n1, converged = true, iterations = k)
        end
        n = n1
    end
    return (n = n, converged = false, iterations = maxiter)
end

function equilibrium(c::OpenCommunity, n0::AbstractVector{<:Pair}, env; kwargs...)
    d = Dict(n0)
    return equilibrium(c, reduce(vcat, [collect(Float64, d[s]) for s in species_names(c)]), env; kwargs...)
end

"""
    jacobian(c, n, env; h = 1e-6)

Jacobian of one time step of a closed community at state `n`, by central differences with
step `h * max(|n_j|, 1)`. Unlike [`community_matrix`](@ref), which freezes the pools, this
includes the effect of each state on every species through the pools, so its off-diagonal
blocks are the direct between-species effects. At an equilibrium, the equilibrium is
locally stable if every eigenvalue has modulus below one.
"""
function jacobian(c::OpenCommunity, n::AbstractVector{<:Real}, env; h::Real = 1e-6)
    k = length(c)
    J = zeros(Float64, k, k)
    for j in 1:k
        d = h * max(abs(n[j]), 1.0)
        up = collect(Float64, n)
        dn = collect(Float64, n)
        up[j] += d
        dn[j] -= d
        J[:, j] = (advance(c, up, env) .- advance(c, dn, env)) ./ (2d)
    end
    return J
end

"""
    invasion_growth_rate(invader, pool, env)
    invasion_growth_rate(invader, resident, n, env)

Asymptotic growth rate of a rare `invader` (a species) when the pools are held at the
totals `pool` (a `NamedTuple`), or at the totals of the `resident` community in state `n`.
A value above one means the invader increases when rare. The invader's own contribution
to the pools is neglected, as is any dependence of its rates on its own abundance.
"""
function invasion_growth_rate(sp::AbstractSpecies, pool::NamedTuple, env)
    names = Tuple(p.first for p in sp.occupancy)
    all(q -> haskey(pool, q), names) ||
        throw(ArgumentError("invader :$(sp.name) needs pools $(names), got $(keys(pool))"))
    ctx = Context(zeros(length(sp.stages)), NamedTuple{names}(Tuple(pool[q] for q in names)), env)
    return StructuredPopulationCore.lambda(projection_matrix(sp, ctx))
end

function invasion_growth_rate(sp::AbstractSpecies, resident::OpenCommunity, n::AbstractVector{<:Real}, env)
    return invasion_growth_rate(sp, pool_totals(resident, n, env), env)
end

"""
    mean_abundance(c, n0, envs; burnin = 0)

Mean of each state of a closed community over the time steps after `burnin`, projecting
from `n0` through `envs`.
"""
function mean_abundance(c::OpenCommunity, n0, envs; burnin::Integer = 0)
    sol = project(c, n0, envs)
    rows = (burnin + 2):size(sol.u, 1)
    return vec(sum(sol.u[rows, :]; dims = 1)) ./ length(rows)
end

"""
    abundance_sensitivity(build, values, envs; burnin = 0, replicates = 1)

The simulation-based sensitivity of Lytle & Tonkin (2023). `build(v)` returns a model (a
`NamedTuple` with `community` and `n0`) in which one vital rate has the value `v`. For
each `v` in `values` the model is projected and the mean abundance of every state after
`burnin` is recorded, averaged over `replicates` runs. `envs` is a sequence of
environments, or a function `r -> sequence` giving the sequence for replicate `r` (use the
same sequences for every `v` to remove sampling noise from the comparison).

Each state's mean abundance is then regressed on `values`. Returns `(slope, r2, means)`:
the slope and R² for each state (in layout order) and the matrix of mean abundances
(one row per value). The slope estimates the derivative of mean abundance with respect to
the vital rate; R² says how consistent the response is across the range.
"""
function abundance_sensitivity(build, values::AbstractVector{<:Real}, envs; burnin::Integer = 0, replicates::Integer = 1)
    length(values) >= 2 || throw(ArgumentError("need at least two values to estimate a slope"))
    means = nothing
    for (i, v) in enumerate(values)
        m = build(v)
        acc = zeros(Float64, length(m.community))
        for r in 1:replicates
            seq = envs isa Function ? envs(r) : envs
            acc .+= mean_abundance(m.community, m.n0, seq; burnin)
        end
        means === nothing && (means = zeros(Float64, length(values), length(acc)))
        means[i, :] = acc ./ replicates
    end
    x = collect(Float64, values)
    xc = x .- sum(x) / length(x)
    sxx = sum(abs2, xc)
    slope = zeros(Float64, size(means, 2))
    r2 = zeros(Float64, size(means, 2))
    for j in axes(means, 2)
        y = means[:, j]
        yc = y .- sum(y) / length(y)
        syy = sum(abs2, yc)
        sxy = sum(xc .* yc)
        slope[j] = sxy / sxx
        r2[j] = syy == 0 ? 0.0 : sxy^2 / (sxx * syy)
    end
    return (slope = slope, r2 = r2, means = means)
end
