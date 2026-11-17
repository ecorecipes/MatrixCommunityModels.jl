# Analysis of composed communities: equilibria, local stability, invasion, and the
# cross-species sensitivity analysis of Lytle & Tonkin (2023).

"""
    equilibrium(c, n0, env; maxiter = 100_000, tol = 1e-10)

Iterate a closed community in a constant environment `env` from `n0` until the largest
change in any state, relative to `max(|n|, 1)`, falls below `tol`. Returns
`(n, converged, iterations)`. A stable equilibrium is found only if the iteration
converges to one; cycles and chaos return `converged = false`.
"""
function equilibrium(c::OpenCommunity, n0::AbstractVector, env; maxiter::Integer = 100_000, tol::Real = 1e-10)
    n = state_vector(c, n0)
    for k in 1:maxiter
        n1 = advance(c, n, env)
        if maximum(abs.(n1 .- n) ./ max.(abs.(n), 1.0)) < tol
            return (n = n1, converged = true, iterations = k)
        end
        n = n1
    end
    return (n = n, converged = false, iterations = maxiter)
end

"""
    equilibrium(species, pool => (low, high), n0, env; tol = 1e-10, maxiter = 200)

The equilibrium of a single species that regulates itself through one pool, found without
iterating the dynamics: the value of `pool` between `low` and `high` at which the species,
with its matrix evaluated there, has a growth rate of one. The state returned is the
stable distribution of that matrix, scaled so that the species' contribution to the pool
equals the value found. Returns `(n, total, converged, iterations)`, where `total` is the
pool value.

Use it when iteration is slow, for example a long-lived species on a large mesh, whose
approach to equilibrium takes thousands of steps. It assumes that the matrix depends on
the state only through the pool (and through totals that scale with it), that the
species' contribution to the pool is proportional to its abundance, and that the growth
rate falls below one somewhere between `low` and `high` and exceeds it elsewhere in that
range. `n0` is any state of the species with a positive contribution to the pool; only its
shape is used, to evaluate the matrix. The equilibrium found need not be stable: check
with [`jacobian`](@ref), or by projecting from it.
"""
function equilibrium(c::OpenCommunity, bracket::Pair{Symbol, <:Tuple{Real, Real}}, n0::AbstractVector, env;
        tol::Real = 1e-10, maxiter::Integer = 200)
    length(c.names) == 1 || throw(ArgumentError("equilibrium by pool value needs a single species"))
    pool = bracket.first
    held(n) = getproperty(pool_totals(c, n, env), pool)
    start = state_vector(c, n0)
    held(start) > 0 || throw(ArgumentError("the starting state contributes nothing to the pool :$pool"))
    # the state with the stable shape of the matrix at pool value `x`, and its growth rate
    function at(x)
        n = start .* (x / held(start))
        rate = NaN
        for _ in 1:50
            F = eigen(community_matrix(c, n, env))
            lead = argmax(abs.(F.values))
            shape = abs.(F.vectors[:, lead])
            stable = shape .* (x / held(shape))
            done = abs(abs(F.values[lead]) - rate) <= 1e-13
            n, rate = stable, abs(F.values[lead])
            done && break
        end
        return n, rate
    end
    a, b = float(bracket.second[1]), float(bracket.second[2])
    (na, ra), (nb, rb) = at(a), at(b)
    (ra - 1) * (rb - 1) <= 0 || throw(ArgumentError(
        "the growth rate is $ra at :$pool = $a and $rb at $b: the equilibrium is not between them"))
    fa, fb = ra - 1, rb - 1
    n, x, iterations, converged = nb, b, 0, false
    for k in 1:maxiter                       # regula falsi with the Illinois correction
        iterations = k
        x = fb == fa ? (a + b) / 2 : b - fb * (b - a) / (fb - fa)
        n, r = at(x)
        f = r - 1
        if abs(f) <= tol || abs(b - a) <= tol * max(abs(x), 1.0)
            converged = true
            break
        end
        if f * fb < 0
            a, fa = b, fb
        else
            fa /= 2
        end
        b, fb = x, f
    end
    return (n = n, total = x, converged = converged, iterations = iterations)
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

Asymptotic growth rate of a rare `invader` (a single species, from [`species`](@ref)) when
the pools are held at the totals `pool` (a `NamedTuple`), or at the totals of the
`resident` community in state `n`. A value above one means the invader increases when
rare. The invader's own contribution to the pools is neglected, as is any dependence of
its rates on its own abundance.
"""
function invasion_growth_rate(invader::OpenCommunity, pool::NamedTuple, env)
    return StructuredPopulationCore.lambda(_rare_matrix(invader, pool, env))
end

function invasion_growth_rate(invader::OpenCommunity, resident::OpenCommunity, n::AbstractVector, env)
    return invasion_growth_rate(invader, pool_totals(resident, n, env), env)
end

# the net and response context of a single rare species when the pools are held at `pool`
function _rare_context(invader::OpenCommunity, pool::NamedTuple, env)
    length(invader.names) == 1 || throw(ArgumentError("the invader must be a single species"))
    all(p -> haskey(pool, p), invader.ports) ||
        throw(ArgumentError("invader :$(only(invader.names)) needs pools $(Tuple(invader.ports)), got $(keys(pool))"))
    net = invader.system[only(invader.names)]
    return net, context(net, zeros(length(stage_names(net))), pool, env)
end

# the matrix of a single rare species when the pools are held at `pool`
function _rare_matrix(invader::OpenCommunity, pool::NamedTuple, env)
    net, ctx = _rare_context(invader, pool, env)
    return to_matrix(net, ctx)
end

"""
    invasion_sensitivity(invader, pool, env)
    invasion_sensitivity(invader, resident, n, env)

Sensitivity and elasticity of the [`invasion_growth_rate`](@ref) `λ` of a rare `invader` to
its own vital rates, with the pools held at the totals `pool` (a `NamedTuple`) or at the
totals of the `resident` community in state `n`. It answers the question: which change to
an excluded species (`λ < 1`) would let it persist?

They follow from the invader's matrix `A` at those pools by the eigenvector formula
`∂λ/∂a_ij = v_i w_j / ⟨v, w⟩`, with `w` and `v` the right and left eigenvectors of `λ`.
The pools are held fixed: a small change to a rare invader does not change the
residents. Returns a `NamedTuple` with

- `lambda`: the invasion growth rate;
- `sensitivity` and `elasticity`: matrices over the invader's stages (rows are
  destinations, columns sources), `∂λ/∂a_ij` and `(a_ij/λ) ∂λ/∂a_ij`; the elasticities sum
  to one;
- `transitions`: one row per named transition, `(transition, sensitivity, elasticity)`,
  where `sensitivity` is `∂λ/∂m` for a factor `m` multiplying every value of the
  transition (at `m = 1`), that is the sum of `a_ij ∂λ/∂a_ij` over its entries, and
  `elasticity` is `sensitivity / λ`, the proportional change in `λ` per proportional change
  in the transition;
- `entries`: one row per entry of each transition, `(transition, from, to, value,
  sensitivity, elasticity)`, with `value` the entry evaluated at the pools,
  `sensitivity = ∂λ/∂a_ij` and `elasticity = (value/λ) ∂λ/∂a_ij`.

Two transitions that share an entry each get its contribution in full. For the
sensitivity to a parameter that enters several rates, or the residents too, use the
method below.

    invasion_sensitivity(build, θ; h = 1e-6)

Derivative of the invasion growth rate with respect to a parameter `θ`, by central
differences with step `h * max(|θ|, 1)`. `build(θ)` returns a `NamedTuple`, either
`(invader, pool, env)` or `(invader, resident, n, env)`; whether the residents' state `n`
is recomputed for each `θ` (for example by [`equilibrium`](@ref)) is up to `build`.
Returns `(lambda, derivative, elasticity)` with `elasticity = θ/λ × derivative`.
"""
function invasion_sensitivity(invader::OpenCommunity, pool::NamedTuple, env)
    net, ctx = _rare_context(invader, pool, env)
    A = to_matrix(net, ctx)
    λ = StructuredPopulationCore.lambda(A)
    S = StructuredPopulationCore.sensitivity(A)
    E = (A .* S) ./ λ
    frozen = evaluate(net, ctx)
    idx = Dict(s => i for (i, s) in enumerate(stage_names(net)))
    entries = [(transition = t, from = e.first.first, to = e.first.second, value = Float64(e.second),
                sensitivity = S[idx[e.first.second], idx[e.first.first]],
                elasticity = e.second * S[idx[e.first.second], idx[e.first.first]] / λ)
               for t in transition_names(net) for e in frozen.transition_values[t]]
    transitions = map(transition_names(net)) do t
        s = sum((x.value * x.sensitivity for x in entries if x.transition === t); init = 0.0)
        (transition = t, sensitivity = s, elasticity = s / λ)
    end
    return (lambda = λ, sensitivity = S, elasticity = E, transitions = transitions, entries = entries)
end

function invasion_sensitivity(invader::OpenCommunity, resident::OpenCommunity, n::AbstractVector, env)
    return invasion_sensitivity(invader, pool_totals(resident, n, env), env)
end

function invasion_sensitivity(build, θ::Real; h::Real = 1e-6)
    rate(m) = haskey(m, :pool) ? invasion_growth_rate(m.invader, m.pool, m.env) :
                                 invasion_growth_rate(m.invader, m.resident, m.n, m.env)
    λ = rate(build(θ))
    d = h * max(abs(θ), 1.0)
    derivative = (rate(build(θ + d)) - rate(build(θ - d))) / (2d)
    return (lambda = λ, derivative = derivative, elasticity = θ / λ * derivative)
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

"""
    equilibrium_sensitivity(build, θ, env; h = 1e-6, kwargs...)

Derivative of the equilibrium of a community with respect to a parameter, without
simulation over a range of values. `build(θ)` returns a model (a `NamedTuple` with
`community` and `n0`) for parameter value `θ`. The equilibrium `n*` is found by
[`equilibrium`](@ref) (which receives `kwargs`), and then

    dn*/dθ = (I - J)⁻¹ ∂f/∂θ

where `f` is one time step, `J` its [`jacobian`](@ref) at `n*`, and `∂f/∂θ` is a central
difference with step `h * max(|θ|, 1)` holding the state at `n*`. Returns
`(n, dn, converged)`.

This applies to a constant environment and a stable equilibrium at which the step is
differentiable; thresholds and `max(0, ...)` terms that are active at `n*` break it. For
fluctuating environments use [`abundance_sensitivity`](@ref).
"""
function equilibrium_sensitivity(build, θ::Real, env; h::Real = 1e-6, kwargs...)
    m = build(θ)
    eq = equilibrium(m.community, m.n0, env; kwargs...)
    J = jacobian(m.community, eq.n, env)
    d = h * max(abs(θ), 1.0)
    dfdθ = (advance(build(θ + d).community, eq.n, env) .- advance(build(θ - d).community, eq.n, env)) ./ (2d)
    return (n = eq.n, dn = (I - J) \ dfdθ, converged = eq.converged)
end

"""
    stochastic_invasion_rate(invader, resident, n0, envs; burnin = 0)

Long-run growth rate of a rare `invader` in a fluctuating environment. The `resident`
community is projected from `n0` through `envs`; at each step after `burnin` the invader's
projection matrix is evaluated at the resident's pool totals and that step's environment,
and applied to the invader's stage vector, which is renormalised. Returns the geometric
mean growth factor per step, so a value above one means the invader increases when rare.
In a constant environment at the resident's equilibrium it equals
[`invasion_growth_rate`](@ref).
"""
function stochastic_invasion_rate(invader::OpenCommunity, resident::OpenCommunity, n0, envs; burnin::Integer = 0)
    sol = project(resident, n0, envs)
    names = Tuple(resident.ports)
    # the pools of step t as the solver logged them, numbers or vectors
    logged(t) = NamedTuple{names}(Tuple(haskey(sol.pool_vectors, p) ? sol.pool_vectors[p][t, :] : sol.pools[t, j]
                                        for (j, p) in enumerate(names)))
    k = length(invader)
    v = fill(1.0 / k, k)
    total = 0.0
    steps = 0
    for (t, env) in enumerate(envs)
        v = _rare_matrix(invader, logged(t), env) * v
        g = sum(v)
        g > 0 || return 0.0
        v ./= g
        if t > burnin
            total += log(g)
            steps += 1
        end
    end
    steps > 0 || throw(ArgumentError("burnin leaves no time steps"))
    return exp(total / steps)
end

"""
    InteractionNetwork(effect, labels; threshold)

A directed network read off a matrix of net effects, such as the slopes from
[`abundance_sensitivity`](@ref): `effect[i, j]` is the effect of perturbing node `j` on
node `i`, and there is a link `j → i` when `|effect[i, j]| > threshold`. Self-effects are
not links. Use [`links`](@ref), [`connectance`](@ref), [`link_density`](@ref),
[`reciprocity`](@ref), [`out_degree`](@ref) and [`in_degree`](@ref).
"""
struct InteractionNetwork
    labels::Vector{String}
    effect::Matrix{Float64}
    adjacency::BitMatrix
end

function InteractionNetwork(effect::AbstractMatrix{<:Real}, labels::AbstractVector; threshold::Real)
    k = size(effect, 1)
    (size(effect, 2) == k && length(labels) == k) ||
        throw(DimensionMismatch("effect must be square with one label per node"))
    adjacency = BitMatrix(i != j && abs(effect[i, j]) > threshold for i in 1:k, j in 1:k)
    return InteractionNetwork(String.(string.(labels)), Matrix{Float64}(effect), adjacency)
end

Base.show(io::IO, net::InteractionNetwork) =
    print(io, "InteractionNetwork(", length(net.labels), " nodes, ", count(net.adjacency), " links)")

"""
    links(net)

The links of an [`InteractionNetwork`](@ref) as `(from, to, effect)` tuples, strongest first.
"""
function links(net::InteractionNetwork)
    out = [(from = net.labels[j], to = net.labels[i], effect = net.effect[i, j])
           for j in axes(net.adjacency, 2), i in axes(net.adjacency, 1) if net.adjacency[i, j]]
    return sort(vec(out); by = l -> -abs(l.effect))
end

"""
    connectance(net)

Fraction of the possible directed links (excluding self-links) that are present.
"""
connectance(net::InteractionNetwork) = count(net.adjacency) / (length(net.labels) * (length(net.labels) - 1))

"""
    link_density(net)

Mean number of links per node.
"""
link_density(net::InteractionNetwork) = count(net.adjacency) / length(net.labels)

"""
    reciprocity(net)

Fraction of links whose reverse link is also present.
"""
function reciprocity(net::InteractionNetwork)
    total = count(net.adjacency)
    return total == 0 ? 0.0 : count(net.adjacency .& net.adjacency') / total
end

"""
    out_degree(net)

Number of nodes each node affects, in label order.
"""
out_degree(net::InteractionNetwork) = vec(sum(net.adjacency; dims = 1))

"""
    in_degree(net)

Number of nodes each node is affected by, in label order.
"""
in_degree(net::InteractionNetwork) = vec(sum(net.adjacency; dims = 2))

"""
    periodic_matrix(species, pool, envs)

The projection matrix of `species` over one cycle of environments `envs` (for example the
seasons of a year) with the pools held at `pool` (a `NamedTuple`): the product
`A(envs[end]) ⋯ A(envs[1])`. Its dominant eigenvalue is the growth factor per cycle, so
with `pool` describing an empty community it gives the species' annual growth rate when
rare.
"""
function periodic_matrix(invader::OpenCommunity, pool::NamedTuple, envs)
    k = length(invader)
    P = Matrix{Float64}(I, k, k)
    for env in envs
        P = _rare_matrix(invader, pool, env) * P
    end
    return P
end

"""
    stage_total(sol, species, stages)

Sum through time of some stages of `species`, from a [`CommunityTrajectory`](@ref).
`stages` is a stage or group name (as understood by `StageVector`, so `:pupa` is the
total over its age classes and strata), or a vector of them.
"""
function stage_total(sol::CommunityTrajectory, sp::Symbol, stages)
    groups = sol.community.system[sp].groups
    idx = Int[]
    for s in (stages isa Symbol ? (stages,) : stages)
        haskey(groups, s) || throw(ArgumentError("species :$sp has no stage or group :$s"))
        append!(idx, groups[s])
    end
    return vec(sum(sol[sp][:, unique(idx)]; dims = 2))
end

"""
    abundance(sol, species)
    abundance(ens, species)

Total abundance of `species` (summed over its stages) through time: a vector for a
[`CommunityTrajectory`](@ref), a replicate-by-time matrix for a
[`CommunityEnsemble`](@ref).
"""
abundance(sol::CommunityTrajectory, sp::Symbol) = vec(sum(sol[sp]; dims = 2))

"""
    CommunityEnsemble

Replicate projections of one community, from [`ensemble`](@ref). `u[r, t + 1, :]` is the
state of replicate `r` at time `t`.
"""
struct CommunityEnsemble
    t::Vector{Int}
    u::Array{Float64,3}
    community::OpenCommunity
end

Base.show(io::IO, ens::CommunityEnsemble) =
    print(io, "CommunityEnsemble(", size(ens.u, 1), " replicates, ", length(ens.t) - 1, " steps, ",
          size(ens.u, 3), " states)")

abundance(ens::CommunityEnsemble, sp::Symbol) = dropdims(sum(ens.u[:, :, ens.community.layout[sp]]; dims = 3); dims = 3)

"""
    ensemble(c, n0, envs; replicates, rng, demographic = true)

Project a closed community `replicates` times and collect the results in a
[`CommunityEnsemble`](@ref). `envs` is a sequence of environments used for every
replicate, or a function `r -> sequence` giving replicate `r` its own (environmental
stochasticity). With `demographic = true` each projection has demographic stochasticity
drawn from `rng` (see [`project`](@ref)); with `false`, `rng` is not used and replicates
differ only through `envs`.
"""
function ensemble(c::OpenCommunity, n0, envs; replicates::Integer, rng = nothing, demographic::Bool = true)
    replicates >= 1 || throw(ArgumentError("replicates must be at least 1"))
    demographic && rng === nothing && throw(ArgumentError("demographic = true needs an rng"))
    u = nothing
    t = Int[]
    for r in 1:replicates
        seq = envs isa Function ? envs(r) : envs
        sol = project(c, n0, seq; rng = demographic ? rng : nothing)
        if u === nothing
            u = zeros(Float64, replicates, size(sol.u, 1), size(sol.u, 2))
            t = sol.t
        end
        size(sol.u, 1) == size(u, 2) || throw(DimensionMismatch("every replicate needs the same number of steps"))
        u[r, :, :] = sol.u
    end
    return CommunityEnsemble(t, u, c)
end

"""
    extinction_times(ens, species; threshold = 1)

For each replicate of a [`CommunityEnsemble`](@ref), the first time at which the total
abundance of `species` is below `threshold`, or `nothing` if that never happens.
"""
function extinction_times(ens::CommunityEnsemble, sp::Symbol; threshold::Real = 1)
    N = abundance(ens, sp)
    return [(i = findfirst(<(threshold), view(N, r, :)); i === nothing ? nothing : ens.t[i]) for r in axes(N, 1)]
end

"""
    quasi_extinction(ens, species; threshold = 1)

Cumulative probability, at each time, that the total abundance of `species` has fallen
below `threshold` at or before that time, estimated as the fraction of replicates of a
[`CommunityEnsemble`](@ref).
"""
function quasi_extinction(ens::CommunityEnsemble, sp::Symbol; threshold::Real = 1)
    times = extinction_times(ens, sp; threshold)
    return [count(x -> x !== nothing && x <= t, times) / length(times) for t in ens.t]
end

"""
    periodic_point(c, n0, env; period = 1, tol = 1e-12, maxiter = 100)

A point of a closed community that returns to itself after `period` steps in the constant environment `env`, found
by Newton's method from `n0`: an equilibrium for `period = 1`, a point of a two-cycle for `period = 2`. Unlike
[`equilibrium`](@ref), which iterates the dynamics, this also finds unstable equilibria and cycles. Returns `(n,
converged, iterations, multipliers)`; `multipliers` are the eigenvalues of the Jacobian of the `period`-step map at
`n` (the product of the one-step Jacobians around the cycle), so the point is locally stable if all have modulus
below one. If Newton's method leaves the region where the model is finite, it stops with `converged = false`.
"""
function periodic_point(c::OpenCommunity, n0::AbstractVector, env; period::Integer = 1, tol::Real = 1e-12,
                        maxiter::Integer = 100)
    n = state_vector(c, n0)
    k = length(n)
    function along(x)                       # the period-step map and its Jacobian
        J, y = Matrix{Float64}(I, k, k), x
        for _ in 1:period
            J = jacobian(c, y, env) * J
            y = advance(c, y, env)
        end
        return y, J
    end
    converged, iterations = false, maxiter
    for it in 1:maxiter
        y, J = along(n)
        all(isfinite, J) && all(isfinite, y) || break
        step = (J - I) \ (y .- n)
        all(isfinite, step) || break
        n = n .- step
        if maximum(abs.(step) ./ max.(abs.(n), 1.0)) < tol
            converged, iterations = true, it
            break
        end
    end
    _, J = along(n)
    multipliers = all(isfinite, J) ? eigvals(J) : fill(complex(NaN), k)
    return (n = n, converged = converged, iterations = iterations, multipliers = multipliers)
end
