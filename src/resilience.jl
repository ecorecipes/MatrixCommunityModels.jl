# Resilience of composed communities through ResilienceCore.jl:
#
# - `linearise`: the Jacobian of a closed community at a fixed point (or of the k-step map at a
#   point of a k-cycle) as a discrete-time `EquilibriumJacobian`;
# - `periodic_operator`: the matrices of a rare invader over a cycle of environments as a
#   `PeriodicOperator`;
# - `MCMState`: a community at a fixed point with ResilienceCore's community interface, so that
#   its simulated press, removal and extinction experiments run on matrix community models.
#
# `jacobian` (analysis.jl) and `linearise` are ResilienceCore's functions, extended here.

"""
    linearise(c, n, env; period = 1, check = true, tol = 1e-8, h = 1e-6)

The Jacobian of the closed community `c` at a fixed point `n` of its dynamics in the constant
environment `env`, as a ResilienceCore `DiscreteJacobian` whose equilibrium is `n`. The
metrics of ResilienceCore.jl for discrete time then apply: `isstable`, `resilience`,
`reactivity`, `max_amplification`, the median metrics and the press sensitivities.

With `period = k`, `n` is a point of a `k`-cycle and the result is the Jacobian of the
`k`-step map: the product of the one-step [`jacobian`](@ref)s around the cycle. Its
eigenvalues are the `multipliers` that [`periodic_point`](@ref) returns, and its metrics
describe one pass around the cycle.

`check` verifies that `n` returns to itself after `period` steps, to the relative tolerance
`tol`; find the point first with [`periodic_point`](@ref) (which also finds unstable points)
or [`equilibrium`](@ref). `h` is the relative step of the central differences.
"""
function linearise(c::OpenCommunity, n::AbstractVector, env; period::Integer = 1, check::Bool = true,
                   tol::Real = 1e-8, h::Real = 1e-6)
    period >= 1 || throw(ArgumentError("period must be at least 1, got $period"))
    x = collect(Float64, state_vector(c, n))
    k = length(x)
    J, y = Matrix{Float64}(I, k, k), x
    for _ in 1:period
        J = jacobian(c, y, env; h = h) * J
        y = advance(c, y, env)
    end
    if check
        gap = maximum(abs.(y .- x) ./ max.(abs.(x), 1.0))
        gap <= tol || throw(ArgumentError(
            "the state does not return to itself after $period step(s) (relative change $gap > tol = $tol); " *
            "find a fixed point with `periodic_point` or `equilibrium`, or pass `check = false`"))
    end
    return ResilienceCore.DiscreteJacobian(J, x)
end

# mesh weights of a species' states, in layout order: bin widths of continuous stages, ones for
# discrete ones, so that norms of a discretised kernel approximate integrals
function _state_weights(net)
    layout = CategoricalPopulationDynamics.state_layout(net)
    w = reduce(vcat, (ResilienceCore.mesh_weights(d) for d in values(layout)); init = Float64[])
    return length(w) == length(stage_names(net)) ? w : ones(length(stage_names(net)))
end

"""
    periodic_operator(invader, pool, envs)

The projection matrices of a rare `invader` (a single species, from [`species`](@ref)) over
one cycle of environments `envs`, with the pools held at `pool` (a `NamedTuple`), as a
ResilienceCore `PeriodicOperator` whose parts are applied in the order of `envs`. Its
`monodromy` is [`periodic_matrix`](@ref)`(invader, pool, envs)`, so its `growth_rate` is the
invader's growth rate per step (the per-cycle rate to the power `1/length(envs)`) and its
`resilience` the rate at which its stage structure settles, per step. Continuous stages carry
their mesh weights.
"""
function periodic_operator(invader::OpenCommunity, pool::NamedTuple, envs)
    isempty(envs) && throw(ArgumentError("envs must hold at least one environment"))
    net, _ = _rare_context(invader, pool, first(envs))
    w = _state_weights(net)
    return ResilienceCore.PeriodicOperator([ResilienceCore.ProjectionMatrix(_rare_matrix(invader, pool, env); weights = w)
                                            for env in envs])
end

"""
    MCMState(c, n, env; settle = true, newton = true, extinct = 1e-6, tol = 1e-10, maxiter = 100_000)

The closed community `c` at state `n` in the constant environment `env`, with the community
interface of ResilienceCore.jl. ResilienceCore's simulated experiments then run unchanged:
`press_response`, `removal_responses`, `global_mortality_tolerance`,
`local_mortality_tolerances` and `extinction_tolerance`.

- **Species and abundance.** Species are those of `c`, in the order of
  [`species_names`](@ref). `state(s)` is each species' abundance, the sum over its states
  (as in [`abundance`](@ref)); a species is present while that sum exceeds `extinct`.
- **Mortality.** `with_mortality(s, factors)` multiplies the mortality hazard of every stage
  of species `i` by `factors[i]`: a stage that survives a step with probability `σ` survives
  with probability `σ^factors[i]`. It scales the species' survival transitions (those of
  `survival_fecundity_matrices`) and leaves births alone; `factors = ones(...)` is the model
  itself. Scaling the hazard rather than the probability of dying keeps a small press small
  for stages that almost all die, such as seedlings.
- **Removal.** `without_species(s, i)` sets the states of species `i` to zero and holds them
  there.
- **Steady state.** After a perturbation, and on construction when `settle = true`, the
  community is iterated in `env` to a fixed point, trying Newton's method first when
  `newton = true`; species whose abundance falls below `extinct` are removed. If no fixed
  point is reached within `maxiter` steps (cycles, chaos), the state reached is kept and
  `s.converged` is `false`; on construction this is an error. The keyword `steady` of the
  experiments may be a `NamedTuple` overriding `tol` and `maxiter`.
- **Local stability.** `jacobian(s)` is the Jacobian of one (perturbed) step over the states
  of the species present, and `linearise(s)` the corresponding `DiscreteJacobian`.

The state of every stage is `s.n`, in the layout of `c`.

Loading both this package and EcologicalResilience.jl leaves names that differ between them,
such as `equilibrium` and `subcommunity`: qualify them, as in
`MatrixCommunityModels.equilibrium`.
"""
struct MCMState{C<:OpenCommunity,E}
    community::C
    n::Vector{Float64}
    env::E
    factors::Vector{Float64}
    alive::BitVector
    converged::Bool
    newton::Bool
    extinct::Float64
    tol::Float64
    maxiter::Int
end

function MCMState(c::OpenCommunity, n::AbstractVector, env; settle::Bool = true, newton::Bool = true,
                  extinct::Real = 1e-6, tol::Real = 1e-10, maxiter::Integer = 100_000)
    x = collect(Float64, state_vector(c, n))
    alive = BitVector([sum(x[c.layout[sp]]) > extinct for sp in c.names])
    _zero_absent!(x, c, alive)
    s = MCMState(c, x, env, ones(length(c.names)), alive, true, newton, Float64(extinct), Float64(tol), Int(maxiter))
    settle || return s
    t = _settled(s; newton = newton)
    t.converged || throw(ArgumentError(
        "no fixed point reached from this state in $maxiter steps: the dynamics may cycle or be chaotic. " *
        "Analyse a cycle with `linearise(c, n, env; period = k)`, or pass `settle = false`"))
    return t
end

function _zero_absent!(x, c, alive)
    for (i, sp) in enumerate(c.names)
        alive[i] || (x[c.layout[sp]] .= 0)
    end
    return x
end

_copy(s::MCMState; n = copy(s.n), factors = s.factors, alive = copy(s.alive), converged = s.converged) =
    MCMState(s.community, n, s.env, collect(Float64, factors), alive, converged, s.newton, s.extinct, s.tol, s.maxiter)

# One step of the perturbed community: survival scaled by the mortality factors, absent
# species held at zero.
function _step(s::MCMState, n::AbstractVector)
    c, env = s.community, s.env
    y = advance(c, n, env)
    for (i, sp) in enumerate(c.names)
        r = c.layout[sp]
        if !s.alive[i]
            y[r] .= 0
            continue
        end
        f = s.factors[i]
        f == 1 && continue
        U, _ = survival_fecundity_matrices(c.system[sp], _context(c, sp, n, env))
        # pressed / baseline survival, by stage: the mortality hazard -log(σ) times f, so σ^f / σ
        keep = map(σ -> σ <= 0 ? 0.0 : σ >= 1 ? 1.0 : σ^(f - 1), vec(sum(U; dims = 1)))
        y[r] .= max.(y[r] .- U * ((1 .- keep) .* n[r]), 0.0)
    end
    return y
end

_present_states(s::MCMState) =
    reduce(vcat, (collect(s.community.layout[sp]) for (i, sp) in enumerate(s.community.names) if s.alive[i]);
           init = Int[])

function _step_jacobian(s::MCMState, n::AbstractVector, idx::AbstractVector{<:Integer}; h::Real = 1e-6)
    J = zeros(length(idx), length(idx))
    for (col, j) in enumerate(idx)
        d = h * max(abs(n[j]), 1.0)
        up, dn = copy(n), copy(n)
        up[j] += d
        dn[j] -= d
        J[:, col] = (_step(s, up)[idx] .- _step(s, dn)[idx]) ./ (2d)
    end
    return J
end

# Newton's method for a fixed point of the perturbed step, over the states of the species
# present; accepted only if it converges to a non-negative, locally stable point at which
# every species present stays above the extinction threshold.
function _newton(s::MCMState, x::Vector{Float64}; maxiter::Integer = 30)
    idx = _present_states(s)
    isempty(idx) && return nothing
    n = copy(x)
    for _ in 1:maxiter
        y = _step(s, n)
        J = _step_jacobian(s, n, idx)
        all(isfinite, J) && all(isfinite, y) || return nothing
        δ = try
            (J - I) \ (y[idx] .- n[idx])
        catch
            return nothing
        end
        n[idx] .-= δ
        all(isfinite, n) || return nothing
        if maximum(abs.(δ) ./ max.(abs.(n[idx]), 1.0)) < s.tol
            all(>=(-s.tol), n) || return nothing
            n .= max.(n, 0.0)
            all(sum(n[s.community.layout[sp]]) > s.extinct for (i, sp) in enumerate(s.community.names) if s.alive[i]) ||
                return nothing
            ρ = maximum(abs, eigvals(_step_jacobian(s, n, idx)))
            return ρ < 1 ? n : nothing
        end
    end
    return nothing
end

# Iterate to a fixed point, removing species that fall below the extinction threshold.
function _settled(s::MCMState; newton::Bool = s.newton, tol::Real = s.tol, maxiter::Integer = s.maxiter)
    c = s.community
    alive = copy(s.alive)
    x = _zero_absent!(copy(s.n), c, alive)
    t = _copy(s; n = x, alive = alive)
    if newton
        z = _newton(t, x)
        z === nothing || return _copy(t; n = z, converged = true)
    end
    for _ in 1:maxiter
        y = _step(t, x)
        dead = [alive[i] && sum(y[c.layout[sp]]) <= t.extinct for (i, sp) in enumerate(c.names)]
        if any(dead)
            alive .&= .!dead
            _zero_absent!(y, c, alive)
            t = _copy(t; n = y, alive = alive)
        elseif maximum(abs.(y .- x) ./ max.(abs.(x), 1.0)) < tol
            return _copy(t; n = y, converged = true)
        end
        x = y
    end
    return _copy(t; n = x, converged = false)
end

_steady(s::MCMState, steady) = steady === nothing ? (tol = s.tol, maxiter = s.maxiter) :
                                (tol = get(steady, :tol, s.tol), maxiter = get(steady, :maxiter, s.maxiter))

ResilienceCore.state(s::MCMState) = [sum(s.n[s.community.layout[sp]]) for sp in s.community.names]

ResilienceCore.present_species(s::MCMState) = findall(s.alive)

function ResilienceCore.with_mortality(s::MCMState, factors::AbstractVector; newton::Bool = true, steady = nothing)
    length(factors) == length(s.alive) ||
        throw(ArgumentError("factors has $(length(factors)) entries for $(length(s.alive)) species"))
    all(>=(0), factors) || throw(ArgumentError("mortality factors must be non-negative"))
    return _settled(_copy(s; factors = factors); newton = newton && s.newton, _steady(s, steady)...)
end

function ResilienceCore.without_species(s::MCMState, i::Integer; steady = nothing)
    s.alive[i] || throw(ArgumentError("species $i ($(s.community.names[i])) is not present"))
    alive = copy(s.alive)
    alive[i] = false
    n = _zero_absent!(copy(s.n), s.community, alive)
    return _settled(_copy(s; n = n, alive = alive); _steady(s, steady)...)
end

ResilienceCore.jacobian(s::MCMState; h::Real = 1e-6) = _step_jacobian(s, s.n, _present_states(s); h = h)

linearise(s::MCMState) = ResilienceCore.DiscreteJacobian(jacobian(s), s.n[_present_states(s)])

function Base.show(io::IO, s::MCMState)
    names = s.community.names[s.alive]
    print(io, "MCMState(", length(names), " of ", length(s.alive), " species present: ", join(names, ", "),
          s.converged ? "" : "; not at a fixed point", ")")
end
