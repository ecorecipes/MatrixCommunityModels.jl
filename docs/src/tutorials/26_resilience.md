# Resilience metrics with ResilienceCore.jl: beetles and riparian guilds
Simon Frost

- [Matrix communities in ResilienceCore’s
  terms](#matrix-communities-in-resiliencecores-terms)
- [One beetle: a stable equilibrium that
  amplifies](#one-beetle-a-stable-equilibrium-that-amplifies)
- [Raising adult mortality: cycles and
  chaos](#raising-adult-mortality-cycles-and-chaos)
- [Two beetles: presses and removals](#two-beetles-presses-and-removals)
- [Riparian guilds in a constant flood
  regime](#riparian-guilds-in-a-constant-flood-regime)
- [With EcologicalResilience.jl](#with-ecologicalresiliencejl)
- [References](#references)

``` julia
ENV["GKSwstype"] = "100"
using MatrixCommunityModels
using ResilienceCore
using LinearAlgebra, Printf, Statistics
using Plots
```

## Matrix communities in ResilienceCore’s terms

[ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl)
computes stability and resilience metrics that dispatch on what a matrix
*means*: the Jacobian at a fixed point, a density-independent projection
matrix, a periodic model. This package hands it three kinds of object:

- `linearise(c, n, env)`: the Jacobian of one step of a closed community
  at a fixed point `n`, as a `DiscreteJacobian`. With `period = k` it is
  the Jacobian of `k` steps at a point of a `k`-cycle.
- `periodic_operator(invader, pool, envs)`: the matrices of a rare
  invader over a cycle of environments, as a `PeriodicOperator`.
- `MCMState(c, n, env)`: a community at a fixed point with
  ResilienceCore’s community interface, so that its simulated press,
  removal and extinction experiments run unchanged.

For a discrete-time fixed point, ResilienceCore’s metrics are:

- `isstable`: whether every eigenvalue of `J` has modulus below one;
- `resilience`: `-ln ρ(J)`, the rate at which small perturbations die
  out per step;
- `reactivity`: `ln ‖J‖₂`, positive when some perturbation grows in the
  first step;
- `max_amplification`: the largest growth of any perturbation before it
  decays, and when.

## One beetle: a stable equilibrium that amplifies

The larva-pupa-adult (LPA) beetle of [vignette
21](21_flour_beetles.md), alone, with the rates of
`TRIBOLIUM_EDMUNDS_2003`. The state is `[L, P, A]`, in two-week steps.

``` julia
beetle(; mu_a = 0.02, cpa = 0.012 / 0.4) =
    flour_beetle(:castaneum; b = 10, mu_l = 0.2, mu_a = mu_a,
                 eats = [:castaneum => (larvae_eggs = 0.01 / 0.4, adults_eggs = 0.01 / 0.4, adults_pupae = cpa)])
one = beetle()
e = periodic_point(one, [30.0, 30, 100], nothing)
J = linearise(one, e.n, nothing)
(state = round.(e.n; digits = 1), stable = isstable(J), resilience = resilience(J),
 reactivity = reactivity(J), max_amplification = max_amplification(J))
```

    (state = [36.7, 29.4, 92.2], stable = true, resilience = 0.06950822732957873, reactivity = 0.2977001026387896, max_amplification = (1.3467578381889962, 1))

The equilibrium is stable, but perturbations die out slowly (a
resilience of 0.07 per step, so a perturbation takes about 33 steps to
shrink tenfold), and the equilibrium is reactive: the worst perturbation
grows by a factor of 1.3 over the first step before it decays.

## Raising adult mortality: cycles and chaos

Costantino et al. (1997) set adult mortality at 96% per two weeks in
laboratory cultures and varied the rate at which adults eat pupae,
`c_pa`, taking the LPA model through cycles and chaos. The same design
with the rates above (illustrative; they are not the fitted values of
that experiment):

``` julia
cpas = 0:0.005:1
attractor = map(cpas) do cpa
    c, x = beetle(mu_a = 0.96, cpa = cpa), [100.0, 100, 100]
    for _ in 1:500
        x = advance(c, x, nothing)
    end
    out = Float64[]
    for _ in 1:100
        x = advance(c, x, nothing)
        push!(out, x[1])
    end
    out
end
ρ = map(cpas) do cpa
    c = beetle(mu_a = 0.96, cpa = cpa)
    p = periodic_point(c, [50.0, 40, 30], nothing)
    maximum(abs, eigvals(operator(linearise(c, p.n, nothing))))
end
p1 = scatter(repeat(collect(cpas); inner = 100), reduce(vcat, attractor); ms = 0.8, msw = 0, color = :black,
             legend = false, ylabel = "larvae on the attractor", title = "μa = 0.96")
p2 = plot(cpas, ρ; lw = 2, legend = false, xlabel = "c_pa (adults eating pupae)",
          ylabel = "ρ(J) at the equilibrium")
hline!(p2, [1]; ls = :dash, color = :gray)
plot(p1, p2; layout = (2, 1), size = (700, 550))
```

![](26_resilience_files/figure-commonmark/cell-4-output-1.svg)

The equilibrium is unstable at every `c_pa` (`ρ(J) > 1`,
`resilience < 0`): with this many adults dying, the population never
settles. What it does instead depends on `c_pa`: a two-cycle without
cannibalism of pupae, irregular dynamics at low rates, then an
eight-cycle and, at the highest rates, a three-cycle.

A cycle is analysed in the same way, with the Jacobian of one pass
around it. Its eigenvalues are the cycle’s multipliers (those of
`periodic_point`), and its resilience per step is that of the period map
divided by the period:

``` julia
function cycle(cpa, k)
    c, x = beetle(mu_a = 0.96, cpa = cpa), [100.0, 100, 100]
    for _ in 1:600
        x = advance(c, x, nothing)
    end
    p = periodic_point(c, x, nothing; period = k)
    Jk = linearise(c, p.n, nothing; period = k)
    return (c_pa = cpa, period = k, converged = p.converged, stable = isstable(Jk),
            resilience_per_step = resilience(Jk) / k)
end
[cycle(0.0, 2), cycle(0.5, 8), cycle(1.0, 3)]
```

    3-element Vector{@NamedTuple{c_pa::Float64, period::Int64, converged::Bool, stable::Bool, resilience_per_step::Float64}}:
     (c_pa = 0.0, period = 2, converged = 1, stable = 1, resilience_per_step = 0.09598213936806506)
     (c_pa = 0.5, period = 8, converged = 1, stable = 1, resilience_per_step = 0.11181756035238047)
     (c_pa = 1.0, period = 3, converged = 1, stable = 1, resilience_per_step = 0.07646747962920593)

All three cycles are stable, each settling at a rate comparable to the
stable equilibrium of the long-lived beetle above.

## Two beetles: presses and removals

`MCMState` gives a community the interface that ResilienceCore’s
simulated experiments use. Pressing a species multiplies the mortality
hazard of every one of its stages: a stage that survives a step with
probability `σ` survives with probability `σ^f`. Removing a species
holds its stages at zero. After each perturbation the community returns
to a fixed point.

The two-species model of Edmunds et al. (2003) has a stable coexistence
equilibrium when larvae eat a moderate share of the other species’ eggs
(region A of that paper):

``` julia
two = flour_beetle_competition(merge(TRIBOLIUM_EDMUNDS_2003, (c_el = [0.01 0.018; 0.018 0.01],)))
s = MCMState(two, periodic_point(two, [8.0, 6, 34, 21, 17, 69], nothing).n, nothing)
(species = species_names(two), abundance = round.(state(s); digits = 1), stable = isstable(linearise(s)),
 local_tolerance = local_mortality_tolerances(s), global_tolerance = global_mortality_tolerance(s))
```

    (species = [:confusum, :castaneum], abundance = [49.0, 107.3], stable = true, local_tolerance = [0.1, 0.30000000000000004], global_tolerance = 6.101000000000001)

`local_mortality_tolerances` is the smallest increase in one species’
mortality, in steps of 10%, that loses any species;
`global_mortality_tolerance` the smallest multiplier of every species’
mortality that does. *T. confusum* is lost when its own mortality rises
by 10%, yet no species is lost until the mortality of both rises
6.1-fold: a press on one species alone tips the competition, a press on
both does not. Removing either species releases the other:

``` julia
[(removed = species_names(two)[i], left = round.(state(without_species(s, i)); digits = 1))
 for i in present_species(s)]
```

    2-element Vector{@NamedTuple{removed::Symbol, left::Vector{Float64}}}:
     (removed = :confusum, left = [0.0, 158.3])
     (removed = :castaneum, left = [163.7, 0.0])

## Riparian guilds in a constant flood regime

The five riparian guilds of [vignette
3](03_riparian_guilds.md) compete for free space,
and floods and droughts kill stage-specific fractions of each. The
experiments need a fixed point, so each kind of year is held constant
here, which removes the fluctuations the guilds normally live with:

``` julia
rip = riparian_model(corrected = true)
flood_year = first(e for e in rip.envs if e.flood == 1)
dry_year = first(e for e in rip.envs if e.flood == 0 && e.drought == 0)
wet = MCMState(rip.community, rip.n0, flood_year)
dry = MCMState(rip.community, rip.n0, dry_year)
(flood_every_year = wet, no_floods = dry)
```

    (flood_every_year = MCMState(5 of 5 species present: C, T, W, S, M), no_floods = MCMState(1 of 5 species present: S))

Without floods, sagebrush (`S`) excludes every other guild. With a flood
every year all five persist at a stable equilibrium:

``` julia
guilds = species_names(rip.community)
(abundance = round.(state(wet); sigdigits = 3), largest_multiplier = maximum(abs, eigvals(operator(linearise(wet)))),
 local_tolerance = local_mortality_tolerances(wet))
```

    (abundance = [1.42e6, 650000.0, 2.54e6, 60.0, 2.32e6], largest_multiplier = 0.9516844707836132, local_tolerance = [0.4, 0.9, 1.6, NaN, 1.2000000000000002])

Cottonwood (`C`) is the most sensitive to extra mortality of its own.
Sagebrush cannot be lost (`NaN`): the published model keeps at least 10
individuals in each of its stages. Removing one guild at a time shows
who takes the space:

``` julia
change(t) = round.(state(t) ./ state(wet); digits = 2)
removals = [(removed = guilds[i], abundance_relative_to_before = change(without_species(wet, i)))
            for i in present_species(wet)]
removals
```

    5-element Vector{@NamedTuple{removed::Symbol, abundance_relative_to_before::Vector{Float64}}}:
     (removed = :C, abundance_relative_to_before = [0.0, 1.01, 1.01, 1.0, 1.01])
     (removed = :T, abundance_relative_to_before = [1.01, 0.0, 1.01, 1.0, 1.01])
     (removed = :W, abundance_relative_to_before = [13.06, 13.06, 0.0, 1.0, 13.06])
     (removed = :S, abundance_relative_to_before = [1.0, 1.0, 1.0, 0.0, 1.0])
     (removed = :M, abundance_relative_to_before = [1.0, 1.0, 1.0, 1.0, 0.0])

Removing willow (`W`) multiplies cottonwood, tamarisk and meadow by the
same factor, about 13: every guild’s recruitment is proportional to the
same free space, so freeing space scales them alike. Sagebrush stays at
its floor. Removing any other guild changes the rest by about 1% or
less.

## With EcologicalResilience.jl

`jacobian` and `linearise` are ResilienceCore’s functions, so this
package and others built on ResilienceCore (EcologicalResilience.jl,
FoodWebDynamics.jl) share them. EcologicalResilience also exports its
own `equilibrium` and `subcommunity`; when both packages are loaded,
write `MatrixCommunityModels.equilibrium`. (At present the two cannot be
installed in one environment: MatrixProjectionModels.jl requires
SciMLBase 2 and EcologicalResilience.jl SciMLBase 3.)

## References

Costantino, R. F., Desharnais, R. A., Cushing, J. M. & Dennis, B.
(1997). Chaotic dynamics in an insect population. *Science* 275,
389-391.

Edmunds, J., Cushing, J. M., Costantino, R. F., Henson, S. M., Dennis,
B. & Desharnais, R. A. (2003). Park’s *Tribolium* competition
experiments: a non-equilibrium species coexistence hypothesis. *Journal
of Animal Ecology* 72, 703-712.

Lytle, D. A., Merritt, D. M., Tonkin, J. D., Olden, J. D. & Reynolds, L.
V. (2017). Linking river flow regimes to riparian plant guilds: a
community-wide modeling approach. *Ecological Applications* 27,
1338-1350.
