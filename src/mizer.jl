# The multispecies size-spectrum model of the R package mizer (Scott, Blanchard & Andersen 2014), North Sea
# version (the MizerParams object NS_params of mizer 3.4.0; Blanchard et al. 2014).
#
# Fish of each species are counted in weight classes of equal width in log weight, shared by all species, from
# the egg weight to the largest asymptotic weight. Below them lie more classes of the same width that only the
# background resource (plankton) occupies. A fish of weight w encounters food at the rate
#
#   E_i(w) = γ_i w^q Σ_classes φ_i(w / w_p) w_p (θ_i,resource N_resource(w_p) + Σ_j θ_ij N_j(w_p)),
#
# with φ_i a lognormal preference for prey a factor β_i smaller (log-sd σ_i) and θ an interaction matrix. It
# eats a share f = E / (E + h w^n) of what it could (the feeding level, a Holling type II response), uses
# α f h w^n - k_s w^p for growth and reproduction, puts a share ψ(w) of that into reproduction, and grows at the
# rate g = (1 - ψ)(that), or not at all if it is negative. Fish die of external mortality, of fishing (knife-edge
# selectivity) and of predation: each fish eats the prey of each class in proportion to φ, so its predators
# remove prey at the rate
#
#   μ_p,j(w_p) = Σ_i θ_ij Σ_classes φ_i(w / w_p) (1 - f_i(w)) γ_i w^q N_i(w).
#
# The eggs produced, 0.5 ε Σ ψ e N / w_egg a year (R_di), are reduced by a Beverton-Holt relation with maximum
# R_max (R_dd) and enter the smallest class. The resource grows back towards its carrying capacity
# (semichemostat) and is eaten as the fish are.
#
# mizer solves the McKendrick-von Foerster equation ∂N/∂t + ∂(gN)/∂w = -μN by its default first-order scheme
# ("euler", upwind flux): every rate is computed from the state at the start of the step, the resource is
# advanced by the exact solution of its equation, and the fish by one implicit upwind step,
#
#   (1 + g_j dt / Δw_j + μ_j dt) n'_j - (g_{j-1} dt / Δw_j) n'_{j-1} = n_j  (+ R_dd dt / Δw_1 in the first class),
#
# a lower-bidiagonal system that mizer solves by forward substitution. Its solution is a matrix applied to the
# state: in numbers N_j = n_j Δw_j, a fish in class k ends the step in class j ≥ k with chance
#
#   T_jk = (1 / b_k) Π_{m=k+1..j} g_{m-1} dt / (Δw_{m-1} b_m),   b_m = 1 + g_m dt / Δw_m + μ_m dt,
#
# so the step is exactly a projection matrix whose entries are recomputed every step from the rates. That is
# how each species is written here: a net over its weight classes whose transitions `staying`, `growing` and
# `spawning` take their values from T. The species share two pools with a value for each predator class of the
# community: the prey biomass each predator class meets (`prey`) and the search volume of each predator class
# (`consumers`), from which every species and the resource compute how fast they are eaten.
#
# mizer computes the two sums over classes by fast Fourier transforms, which are circular convolutions over
# the classes: the sums here are those convolutions written out, with the kernel point-sampled at the grid
# ratios, the own-class term zero, and mizer's rule that values below 1e-18 are zero. Equations, defaults and
# the order of operations were read from the R source of mizer 3.4.0. Weights are in g, time in years, numbers
# are of fish (or resource organisms) in the North Sea. The data are in data/mizer (GPL-3; see
# data/mizer/LICENSE.txt).

const _MIZER_DIR = joinpath(DATA_DIR, "mizer")

_mizer_symbol(name::AbstractString) = Symbol(lowercase(replace(strip(name), "." => "_", " " => "_")))

function _mizer_north_sea(dir)
    table(file) = (lines = [split(l, ",") for l in readlines(joinpath(dir, file)) if !isempty(strip(l))];
        (String.(lines[1]), lines[2:end]))
    num(x) = parse(Float64, x)
    sh, rows = table("north_sea_species.csv")
    col(header, rows, name) = (j = findfirst(==(name), header);
        j === nothing ? throw(ArgumentError("column $name missing")) : [r[j] for r in rows])
    names = String.(col(sh, rows, "species"))
    numbers = Dict(Symbol(c) => num.(col(sh, rows, c)) for c in sh[2:end])
    ih, irows = table("north_sea_interaction.csv")
    ih[2:end] == names && [String(r[1]) for r in irows] == names ||
        throw(ArgumentError("the interaction matrix does not list the species in order"))
    interaction = [num(irows[i][j + 1]) for i in eachindex(names), j in eachindex(names)]
    gh, grows = table("north_sea_gears.csv")
    gear_species = String.(col(gh, grows, "species"))
    gear_species == names || throw(ArgumentError("the gears are not listed one per species, in order"))
    all(==("knife_edge"), col(gh, grows, "sel_func")) || throw(ArgumentError("only knife-edge selectivity is read"))
    gears = _mizer_symbol.(col(gh, grows, "gear"))
    effort = Dict{Symbol, Float64}()
    for (g, e) in zip(gears, num.(col(gh, grows, "initial_effort")))
        effort[g] = e
    end
    order = unique(gears)
    rh, rrows = table("north_sea_resource.csv")
    resource = NamedTuple{Tuple(Symbol.(rh))}(Tuple(num.(only(rrows))))
    _, srows = table("north_sea_initial_state.csv")
    fish = Dict(_mizer_symbol(n) => zeros(Int(resource.no_w)) for n in names)
    plankton = Float64[]
    for r in srows
        if r[1] == "Resource"
            push!(plankton, num(r[3]))
        else
            fish[_mizer_symbol(r[1])][parse(Int, r[2])] = num(r[3])
        end
    end
    return (species = _mizer_symbol.(names), names = names, numbers..., interaction = interaction,
            gear = gears, knife_edge_size = num.(col(gh, grows, "knife_edge_size")),
            catchability = num.(col(gh, grows, "catchability")),
            effort = NamedTuple{Tuple(order)}(Tuple(effort[g] for g in order)), resource = resource,
            initial = (fish = fish, resource = plankton))
end

"""
    MIZER_NORTH_SEA

The North Sea model of the R package mizer 3.4.0 (its `NS_params`, after Blanchard et al.
2014; GPL-3), read from `data/mizer`, for [`mizer_parameters`](@ref). A `NamedTuple`:

- `species` (symbols such as `:cod`, `:n_pout`) and `names` ("Cod", "N.pout"), 12 of them;
- one vector per species parameter of mizer, in mizer's names and units (g, years):
  `w_min` (egg weight), `w_mat`, `w_mat25` (weights at which half and a quarter are
  mature), `w_max`, `w_repro_max`, `beta`, `sigma` (the predation kernel), `alpha`
  (assimilation efficiency), `gamma`, `q` (search volume `γ w^q`), `h`, `n` (maximum intake
  `h w^n`), `ks`, `p`, `k` (metabolism `k_s w^p + k w`), `m` (the exponent of the
  reproductive investment), `z0` (external mortality per year), `erepro` (reproductive
  efficiency), `R_max`, `interaction_resource` and `k_vb`;
- `interaction`: the matrix θ, `[i, j]` the availability of species `j` to predator `i`;
- `gear`, `knife_edge_size`, `catchability`: the fleet that catches each species (one per
  species), the weight (g) from which it is caught, and its catchability;
- `effort`: the fishing effort stored with `NS_params`, by fleet;
- `resource`: the weight grid (`min_w`, `max_w`, `no_w`, `min_w_pp`) and the resource
  parameters `kappa`, `lambda`, `r_pp`, `n` and `w_pp_cutoff`;
- `initial`: the state stored with `NS_params`, numbers per g by class (`fish[:cod]`,
  `resource`).
"""
const MIZER_NORTH_SEA = _mizer_north_sea(_MIZER_DIR)

"""
    mizer_parameters(data = MIZER_NORTH_SEA; dt = 0.1)

The quantities of the mizer model derived from `data` as mizer derives them for
`NS_params`, for a time step of `dt` years (mizer's default). A `NamedTuple` with

- `species`, `names`, `gears` and `dt`;
- `w`, `dw`: the lower edge and width (g) of the fish weight classes, equally spaced in
  log weight; `w_full`, `dw_full`: those of all classes, the resource's below the fish
  classes and then the fish classes; `offset`, the number of resource-only classes;
- `fish`: one `NamedTuple` per species (`p.fish.cod`) with its [`ContinuousStage`](@ref)
  `stage` (cells `weight1`, ..., one per fish class; its mesh points are log10 of the
  geometric middle of each class) and, by class, `search_vol` (`γ w^q`), `intake_max`
  (`h w^n`), `metab`, `psi` (the share of energy put into reproduction: maturity, zero
  below 10⁻⁸, times `(w / w_repro_max)^(m - n)`, at most one, and one above
  `w_repro_max`), `mu_b` (external mortality), `selectivity`, and its `gear`,
  `catchability`, `alpha`, `erepro`, `R_max`, `w_min` and `interaction_resource`;
- `kernel`: one row per species, the predation kernel `φ_i(w_full[1 + m] / w_full[1])` at
  the `m`-th class offset, `exp(-(log(ratio) - log β)² / (2σ²))`, zero at offset zero;
- `interaction`; `mask`, whether each class of `w_full` is below each species' `w_max`;
- `resource`: its [`ContinuousStage`](@ref) `stage` (cells `resource1`, ...), `rate`
  (`r_pp w^(n - 1)`) and `capacity` (`κ w^-λ` below `w_pp_cutoff`, zero above), by class.
"""
function mizer_parameters(data = MIZER_NORTH_SEA; dt::Real = 0.1)
    dt > 0 || throw(ArgumentError("the time step dt must be positive"))
    r = data.resource
    no_w = Int(r.no_w)
    dx = log10(r.max_w / r.min_w) / (no_w - 1)
    x0 = log10(r.min_w)
    w = [10^(x0 + j * dx) for j in 0:(no_w - 1)]
    offset = round(Int, log10(r.min_w / r.min_w_pp) / dx)
    w_full = vcat([10^((x0 - k * dx) - dx) for k in (offset - 1):-1:0], w)
    grow = 10^dx - 1
    dw, dw_full = grow .* w, grow .* w_full
    nsp = length(data.species)
    N = length(w_full)

    function one(i)
        U = log(3) / log(data.w_mat[i] / data.w_mat25[i])
        maturity = [1 / (1 + (x / data.w_mat[i])^(-U)) for x in w]
        maturity[maturity .< 1e-8] .= 0.0
        psi = min.(maturity .* (w ./ data.w_repro_max[i]) .^ (data.m[i] - data.n[i]), 1.0)
        psi[w .> data.w_repro_max[i]] .= 1.0
        stage = ContinuousStage(:weight, x0, x0 + no_w * dx, no_w)
        return (name = data.names[i], stage = stage,
                search_vol = data.gamma[i] .* w .^ data.q[i], intake_max = data.h[i] .* w .^ data.n[i],
                metab = data.ks[i] .* w .^ data.p[i] .+ data.k[i] .* w, psi = psi, mu_b = fill(data.z0[i], no_w),
                selectivity = Float64.(w .>= data.knife_edge_size[i]), gear = data.gear[i],
                catchability = data.catchability[i], alpha = data.alpha[i], erepro = data.erepro[i],
                R_max = data.R_max[i], w_min = data.w_min[i], interaction_resource = data.interaction_resource[i])
    end
    fish = [one(i) for i in 1:nsp]
    all(f -> f.w_min == w[1], fish) || throw(ArgumentError("every species must have its egg weight at the first class"))
    ratio = w_full ./ w_full[1]
    kernel = [m == 1 ? 0.0 : exp(-(log(ratio[m]) - log(data.beta[i]))^2 / (2 * data.sigma[i]^2)) for i in 1:nsp, m in 1:N]
    mask = [w_full[l] < data.w_max[i] for i in 1:nsp, l in 1:N]
    cut = r.w_pp_cutoff
    resource = (stage = ContinuousStage(:resource, x0 - offset * dx, x0 + no_w * dx, N),
                rate = r.r_pp .* w_full .^ (r.n - 1),
                capacity = [x > cut ? 0.0 : r.kappa * x^(-r.lambda) for x in w_full])
    p = (species = data.species, names = data.names, gears = collect(keys(data.effort)), dt = Float64(dt),
         w = w, dw = dw, w_full = w_full, dw_full = dw_full, offset = offset,
         fish = NamedTuple{Tuple(data.species)}(Tuple(fish)), kernel = kernel, interaction = data.interaction,
         mask = mask, resource = resource, cache = Ref{Any}(nothing), lock = ReentrantLock())
    return p
end

_mizer_fish(p, sp::Symbol) = haskey(p.fish, sp) ? p.fish[sp] :
    throw(ArgumentError("unknown species :$sp; the species are $(join(p.species, ", "))"))

# The kernel weight of a predator class `K` (of all classes) on a prey class `l`: mizer's circular convolution.
@inline _mizer_phi(p, i, K, l) = p.kernel[i, mod(K - l, size(p.kernel, 2)) + 1]

# The feeding level of every predator class (rows of the pools: species by species, class by class) and the
# rate at which the predators of each species remove prey of each class of all classes (mizer's pred_rate),
# from the two pools. All species and the resource read the same, so it is kept for the pools it was made from.
function _mizer_predation(p, pool)
    prey, consumers = pool.prey, pool.consumers
    lock(p.lock) do
        c = p.cache[]
        c !== nothing && c[1] === prey && c[2] === consumers && return c[3]
        nsp, no_w, N = length(p.species), length(p.w), length(p.w_full)
        feeding = similar(prey)
        for (i, sp) in enumerate(p.species), k in 1:no_w
            row = (i - 1) * no_w + k
            avail = prey[row] < 1e-18 ? 0.0 : prey[row]
            E = p.fish[sp].search_vol[k] * avail
            feeding[row] = E / (E + p.fish[sp].intake_max[k])
        end
        rate = zeros(nsp, N)
        for i in 1:nsp, l in 1:N
            total = 0.0
            for k in 1:no_w
                row = (i - 1) * no_w + k
                total += _mizer_phi(p, i, p.offset + k, l) * (1 - feeding[row]) * consumers[row]
            end
            rate[i, l] = total < 1e-18 ? 0.0 : total * p.mask[i, l]
        end
        # predation mortality of the fish of each species by class, and of the resource
        fish = transpose(p.interaction) * rate[:, (p.offset + 1):end]
        theta = [p.fish[sp].interaction_resource for sp in p.species]
        resource = vec(transpose(theta) * rate)
        result = (feeding = feeding, rate = rate, fish = fish, resource = resource)
        p.cache[] = (prey, consumers, result)
        return result
    end
end

# The rates of species `sp` (index `i`) in a step, from its numbers `N`, the pools and the fishing effort:
# mizer's encounter, feeding level, e, e_repro, e_growth, mortality, R_di and R_dd.
function _mizer_rates(p, sp::Symbol, N::AbstractVector, pool, effort)
    f = p.fish[sp]
    i = findfirst(==(sp), p.species)
    no_w = length(p.w)
    rows = ((i - 1) * no_w + 1):(i * no_w)
    avail = [x < 1e-18 ? 0.0 : x for x in pool.prey[rows]]
    encounter = f.search_vol .* avail
    feeding = encounter ./ (encounter .+ f.intake_max)
    e = f.alpha .* (1 .- feeding) .* encounter .- f.metab
    e_repro = f.psi .* max.(e, 0.0)
    e_growth = max.(e .- e_repro, 0.0)
    pred_mort = _mizer_predation(p, pool).fish[i, :]
    f_mort = get(effort, f.gear, 0.0) * f.catchability .* f.selectivity
    mort = pred_mort .+ f.mu_b .+ f_mort
    rdi = 0.5 * f.erepro * dot(e_repro, N) / f.w_min
    rdd = rdi / (1 + rdi / f.R_max)
    return (encounter = encounter, feeding_level = feeding, e = e, e_repro = e_repro, e_growth = e_growth,
            pred_mort = pred_mort, f_mort = f_mort, mort = mort, rdi = rdi, rdd = rdd)
end

# The implicit upwind step of mizer as a matrix on numbers: entry [j, k] is the chance that a fish of class k
# ends the step in class j.
function _mizer_step_matrix(p, g, mort)
    dt, dw = p.dt, p.dw
    m = length(dw)
    b = [1 + g[j] * dt / dw[j] + mort[j] * dt for j in 1:m]
    onward = [j == 1 ? 0.0 : g[j - 1] * dt / (dw[j - 1] * b[j]) for j in 1:m]   # passing on into class j
    T = zeros(m, m)
    for k in 1:m
        x = 1 / b[k]
        T[k, k] = x
        for j in (k + 1):m
            x *= onward[j]
            x == 0 && break
            T[j, k] = x
        end
    end
    return T
end

"""
    mizer_fish(p, species)

One species of the mizer model with parameters `p` (see [`mizer_parameters`](@ref)), as an
[`OpenCommunity`](@ref) that meets the others through two pools with a value for each
predator class of the community (the fish classes of every species, species by species):

- `prey`: the prey biomass each predator class meets, `Σ φ(w / w_p) w_p θ N(w_p)` over the
  classes of every species and of the resource (mizer's available energy);
- `consumers`: the search volume of the predators of each class, `γ w^q N(w)`, which
  their feeding level turns into the rate at which they remove prey.

The species adds its fish to both and reads both. Its stages are its weight classes, the
cells of `p.fish[species].stage` (`weight1`, ...), holding numbers of fish, and it has three
transitions, whose values are mizer's implicit upwind step (see `src/mizer.jl`):

- `staying`: the fish of a class that survive the step in their class;
- `growing`: those that survive and end the step in a larger class;
- `spawning`: the eggs, `R_dd dt` a step, which enter the smallest class and die and grow
  in the same step as the fish there do, shared among the spawners in proportion to the
  energy they put into reproduction. Marked as births.

The environment of a step is a `NamedTuple` with `effort`, the fishing effort by fleet
(`(otter = 1.0,)`, missing fleets not fishing); see [`mizer_environment`](@ref).
"""
function mizer_fish(p, sp::Symbol)
    f = _mizer_fish(p, sp)
    weight = f.stage
    cells = stage_names(weight)
    m = length(cells)
    spawners = [k for k in 1:m if f.psi[k] > 0]
    life = ValuedProjectionNet(cells,
        :staying => [(c => c) => 1.0 for c in cells],
        :growing => [(cells[k] => cells[j]) => 1.0 for k in 1:m for j in (k + 1):m],
        :spawning => [(cells[k] => cells[j]) => 1.0 for k in spawners for j in 1:m])
    step(ctx) = (r = _mizer_rates(p, sp, ctx.n[weight], ctx.pool, ctx.env.effort);
                 (rates = r, T = _mizer_step_matrix(p, r.e_growth, r.mort)))
    # eggs from a fish of class k that end the step in class j: R_dd per unit of R_di, times the fish's share of
    # R_di, times dt, times the chance that an egg of the first class ends the step in class j
    function eggs(ctx)
        s = step(ctx)
        share = s.rates.rdi > 0 ? s.rates.rdd / s.rates.rdi : 1.0      # R_dd / R_di tends to one as R_di falls to zero
        per_fish = 0.5 * f.erepro .* s.rates.e_repro ./ f.w_min .* share .* p.dt
        return [s.T[j, 1] * per_fish[k] for j in 1:m, k in 1:m]
    end
    net = births(respond(life,
        :staying => ctx -> step(ctx).T,
        :growing => ctx -> step(ctx).T,
        :spawning => eggs), :spawning)

    # this species' rows among the predator classes, and what it adds to the two pools
    i = findfirst(==(sp), p.species)
    nsp, N = length(p.species), length(p.w_full)
    rows(j) = ((j - 1) * m + 1):(j * m)
    prey = zeros(nsp * m, m)
    for j in 1:nsp, k in 1:m, l in 1:m
        prey[rows(j)[k], l] = p.interaction[j, i] * _mizer_phi(p, j, p.offset + k, p.offset + l) * p.w[l]
    end
    consumers = zeros(nsp * m, m)
    for k in 1:m
        consumers[rows(i)[k], k] = f.search_vol[k]
    end
    return species(sp, net; occupies = [:prey => prey, :consumers => consumers])
end

"""
    mizer_resource(p; name = :resource)

The background resource of the mizer model with parameters `p` (plankton, by weight class
from the smallest class of `p.w_full` to the cutoff), as an [`OpenCommunity`](@ref). Its
stages are the cells of `p.resource.stage` (`resource1`, ...), numbers in each class, and
it follows mizer's semichemostat: in a class with growth rate `r` and capacity `c` (numbers
`c Δw`), eaten at the rate `μ` read from the pools `prey` and `consumers`, a share
`exp(-(r + μ) dt)` of the organisms stays (`renewing`) and `r c Δw (1 - exp(-(r + μ) dt)) /
(r + μ)` are added (an inflow), the exact solution over a step with the rates held fixed.
It adds itself to the pool `prey`.
"""
function mizer_resource(p; name::Symbol = :resource)
    stage = p.resource.stage
    cells = stage_names(stage)
    N = length(cells)
    mort(ctx) = _mizer_predation(p, ctx.pool).resource
    staying(ctx) = exp.(-(p.resource.rate .+ mort(ctx)) .* p.dt)
    function added(ctx, l)
        total = p.resource.rate[l] + mort(ctx)[l]
        steady = p.resource.rate[l] * p.resource.capacity[l] * p.dw_full[l] / total
        return steady * (1 - exp(-total * p.dt))
    end
    net = respond(ValuedProjectionNet(cells, :renewing => [(c => c) => 1.0 for c in cells]),
                  :renewing => ctx -> by_trait(stage, staying(ctx)))
    net = inflow(net, (cells[l] => (ctx -> added(ctx, l)) for l in 1:N)...)
    # the resource on offer to each predator class (rows: species by species, class by class)
    m = length(p.w)
    prey = zeros(length(p.species) * m, N)
    for (j, sp) in enumerate(p.species), k in 1:m, l in 1:N
        prey[(j - 1) * m + k, l] = p.fish[sp].interaction_resource * _mizer_phi(p, j, p.offset + k, l) * p.w_full[l]
    end
    return species(name, net; occupies = [:prey => prey], reads = [:consumers])
end

"""
    mizer_community(p = mizer_parameters(); species = p.species)

The species of the mizer model with parameters `p` ([`mizer_fish`](@ref)) and the
background resource ([`mizer_resource`](@ref)), composed through the pools `prey` and
`consumers`; or the species named in `species` with the resource. A species left out is
neither predator nor prey.
"""
mizer_community(p = mizer_parameters(); species = p.species) =
    community((mizer_fish(p, sp) for sp in species)..., mizer_resource(p))

"""
    mizer_state(p = mizer_parameters(); data = MIZER_NORTH_SEA, species = p.species)

The state stored with mizer's `NS_params` (`data.initial`), as `species => numbers` pairs
for [`project`](@ref) or `state_vector`: the numbers in each weight class (density per g
times the class width), for the fish species named in `species` and the resource.
"""
mizer_state(p = mizer_parameters(); data = MIZER_NORTH_SEA, species = p.species) =
    vcat([sp => data.initial.fish[sp] .* p.dw for sp in species], [:resource => data.initial.resource .* p.dw_full])

"""
    mizer_environment(years; effort = MIZER_NORTH_SEA.effort, dt = 0.1)

The environments of `years` years of the mizer model in steps of `dt` years
(`round(years / dt)` steps): `(effort = ...,)` for each step, with `effort` the fishing
effort by fleet, the same in every step. By default it is the effort stored with
`NS_params` (industrial 0, pelagic 1, beam 0.5, otter 0.5); `(;)` is no fishing.
"""
function mizer_environment(years::Real; effort = MIZER_NORTH_SEA.effort, dt::Real = 0.1)
    steps = round(Int, years / dt)
    return fill((effort = effort,), steps)
end

"""
    mizer_biomass(p, sol)

The biomass in g of each fish species through a run of the mizer model (a
[`CommunityTrajectory`](@ref)): a matrix with one row per time and one column per fish
species of the community, in the order of `species_names(sol.community)` without the
resource.
"""
function mizer_biomass(p, sol::CommunityTrajectory)
    kinds = [sp for sp in species_names(sol.community) if haskey(p.fish, sp)]
    return reduce(hcat, [sol[sp] * p.w for sp in kinds])
end

"""
    mizer_rates(p, c, n, env)

The rates of the fish of the mizer community `c` in state `n` and environment `env`, as
mizer's `getRates` gives them: a `Dict` from species to a `NamedTuple` with `encounter`,
`feeding_level`, `e` (energy for growth and reproduction), `e_repro`, `e_growth`,
`pred_mort`, `f_mort` and `mort` (one value per weight class, per year) and `rdi`, `rdd`
(eggs per year); and, under `:resource`, the mortality of the resource by class.
"""
function mizer_rates(p, c::OpenCommunity, n::AbstractVector, env)
    pool = pool_totals(c, n, env; internal = true)
    out = Dict{Symbol, Any}()
    for sp in species_names(c)
        haskey(p.fish, sp) || continue
        out[sp] = _mizer_rates(p, sp, n[c.layout[sp]], pool, env.effort)
    end
    out[:resource] = copy(_mizer_predation(p, pool).resource)
    return out
end
