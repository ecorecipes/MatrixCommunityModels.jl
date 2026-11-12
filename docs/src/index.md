# MatrixCommunityModels.jl

Matrix community models (MCMs; Lytle & Tonkin 2023) are sets of single-species matrix
population models linked by aggregate density dependence. This package builds them by
composition: each species is declared as named transitions between stages, each shared
quantity (space, biomass, aggregate density) is a pool, and a community is the composite
of species along an undirected wiring diagram whose junctions are the pools.

It is built on the [ecorecipes](https://github.com/ecorecipes) packages
MatrixProjectionModels.jl, ProjectionModels.jl and CategoricalPopulationDynamics.jl, and
on Catlab.jl.

## Tutorials

- [Matrix community models by composition](tutorials/01_introduction.md): species, pools,
  composition, projection, and freezing a community into ordinary matrices.
- [Verde River fish](tutorials/02_verde_fish.md): the seven-species model of Rogosch et
  al. (2019), identical to the original R output.
- [Riparian plant guilds](tutorials/03_riparian_guilds.md): the five-guild model of Lytle
  et al. (2017), the original and a corrected variant.
- [Two orchids](tutorials/04_adapt3_cypripedium.md): the community projection of the
  adapt3 R package, reproduced.
- [Cross-species sensitivity](tutorials/05_cross_species_sensitivity.md): interactions
  revealed by perturbing vital rates, and the resulting network.
- [Composing IPMs](tutorials/06_composing_ipms.md): an integral projection model and a
  matrix model sharing a pool.
- [Metacommunity](tutorials/07_metacommunity.md): patches as nested communities, with
  dispersal as a pool.
- [Equilibria, invasion and chance](tutorials/08_equilibria_invasion_chance.md): analytic
  sensitivities, invasion rates in fluctuating environments, demographic stochasticity.
- [Seasons and phenology](tutorials/09_seasons_and_phenology.md): within-year steps, a
  lottery for space, coexistence by timing, and drift.

Each tutorial is a Quarto vignette, also available as
[HTML and PDF in the repository](https://github.com/ecorecipes/MatrixCommunityModels.jl/tree/main/vignettes).

## Quick start

```julia
using MatrixCommunityModels

K = 1000.0
crowding(ctx) = max(0.0, 1 - ctx.pool.space / K)

grass = Species(:grass, [:seedling, :adult],
                [survival(:establish, :seedling => :adult, 0.3),
                 survival(:persist, :adult => :adult, 0.5),
                 fecundity(:seed, :adult => :seedling, ctx -> 6.0 * crowding(ctx))];
                occupancy = [:space => [0.1, 1.0]])
shrub = Species(:shrub, [:seedling, :adult],
                [survival(:establish, :seedling => :adult, 0.1),
                 survival(:persist, :adult => :adult, 0.95),
                 fecundity(:seed, :adult => :seedling, ctx -> 3.0 * crowding(ctx))];
                occupancy = [:space => [0.1, 4.0]])

com = community(grass, shrub)
sol = project(com, [:grass => [50.0, 10.0], :shrub => [5.0, 1.0]], fill(nothing, 200))
sol[:shrub]            # stage-by-time matrix
```
