# MatrixCommunityModels.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ecorecipes.github.io/MatrixCommunityModels.jl/dev/)

Matrix community models (MCMs; [Lytle & Tonkin 2023](https://doi.org/10.1038/s44185-023-00031-5))
are sets of single-species matrix population models linked by aggregate density
dependence. This package builds them by composition. Each species is declared as named
transitions between stages; each shared quantity (space, biomass, aggregate density) is a
pool; a community is the composite of species along an undirected wiring diagram whose
junctions are the pools. A species never refers to another species, only to pool totals,
so the dependency assumption of an MCM is an explicit object rather than code inside a
simulation loop.

## Features

- Built on the ecorecipes stack: a species is a projection net from
  CategoricalPopulationDynamics.jl (named stages and named transitions), given stage
  durations with `expand_stages`, strata such as infection status with `stratify`, and
  *responses* that say how a named transition depends on shared pools and the environment.
- Composition along Catlab undirected wiring diagrams (`compose`), or by pool name
  (`community`). Composites nest. The result is a `PooledProjectionSystem`.
- Projection (`project`) lowers the community to a `CoupledMPMProblem` and solves it with
  MatrixProjectionModels.jl. A species frozen at any state is an ordinary
  `MatrixProjectionModel` (`matrix_model`) or net (`frozen_net`).
- Analysis: `equilibrium`, `jacobian`, `invasion_growth_rate`, `stochastic_invasion_rate`,
  analytic `equilibrium_sensitivity`, the simulation-based cross-species sensitivity of
  Lytle & Tonkin (`abundance_sensitivity`), `InteractionNetwork`, seasonal models
  (`periodic_matrix`), demographic stochasticity (`project(...; rng)`, `ensemble`).
- Published models rebuilt and checked against their original implementations: the Verde
  River fish model (identical output), the Yampa River riparian guild model (within 4e-16,
  plus a corrected variant), and the adapt3 Cypripedium community (within 2e-15).
- A honey bee, *Varroa*, deformed wing virus and bumblebee community, built step by step.

## Installation

The package and its ecorecipes dependencies are not registered. Its `Project.toml` names
their GitHub URLs, so this is enough (Julia 1.11 or later):

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/MatrixCommunityModels.jl")
```

## Quick start

```julia
using MatrixCommunityModels

m = verde_fish_model()                      # Rogosch et al. (2019)
sol = project(m.community, m.n0, m.envs)    # 1964-2017 flow record
sol[:CACL]                                  # biomass by year and stage
```

## Vignettes

| # | Vignette | Description |
|---|----------|-------------|
| 1 | [Introduction](vignettes/01_introduction/01_introduction.md) ([PDF](vignettes/01_introduction/01_introduction.pdf)) | Species, pools, composition, projection |
| 2 | [Verde River fish](vignettes/02_verde_fish/02_verde_fish.md) ([PDF](vignettes/02_verde_fish/02_verde_fish.pdf)) | Seven species linked by total biomass |
| 3 | [Riparian guilds](vignettes/03_riparian_guilds/03_riparian_guilds.md) ([PDF](vignettes/03_riparian_guilds/03_riparian_guilds.pdf)) | Five guilds linked by space; original and corrected |
| 4 | [Two orchids](vignettes/04_adapt3_cypripedium/04_adapt3_cypripedium.md) ([PDF](vignettes/04_adapt3_cypripedium/04_adapt3_cypripedium.pdf)) | An adapt3 community projection reproduced |
| 5 | [Cross-species sensitivity](vignettes/05_cross_species_sensitivity/05_cross_species_sensitivity.md) ([PDF](vignettes/05_cross_species_sensitivity/05_cross_species_sensitivity.pdf)) | Interactions revealed by perturbing vital rates; the interaction network |
| 6 | [Composing IPMs](vignettes/06_composing_ipms/06_composing_ipms.md) ([PDF](vignettes/06_composing_ipms/06_composing_ipms.pdf)) | An integral projection model and a matrix model sharing a pool |
| 7 | [Metacommunity](vignettes/07_metacommunity/07_metacommunity.md) ([PDF](vignettes/07_metacommunity/07_metacommunity.pdf)) | Patches as nested communities, dispersal as a pool |
| 8 | [Equilibria, invasion and chance](vignettes/08_equilibria_invasion_chance/08_equilibria_invasion_chance.md) ([PDF](vignettes/08_equilibria_invasion_chance/08_equilibria_invasion_chance.pdf)) | Analytic sensitivities, invasion rates in fluctuating environments, demographic stochasticity |
| 9 | [Seasons and phenology](vignettes/09_seasons_and_phenology/09_seasons_and_phenology.md) ([PDF](vignettes/09_seasons_and_phenology/09_seasons_and_phenology.pdf)) | Within-year steps, a lottery for space, coexistence by timing, drift |
| 10 | [Bees 1: a honey bee colony](vignettes/10_honeybee_colony/10_honeybee_colony.md) ([PDF](vignettes/10_honeybee_colony/10_honeybee_colony.pdf)) | A colony from modules: life cycle, stage durations, responses, stores |
| 11 | [Bees 2: Varroa and deformed wing virus](vignettes/11_varroa_and_dwv/11_varroa_and_dwv.md) ([PDF](vignettes/11_varroa_and_dwv/11_varroa_and_dwv.pdf)) | Mites as a second species, the virus as a stratum, treatment timing |
| 12 | [Bees 3: a bumblebee population](vignettes/12_bumblebees/12_bumblebees.md) ([PDF](vignettes/12_bumblebees/12_bumblebees.pdf)) | The annual colony cycle from three modules, nest sites |
| 13 | [Bees 4: sharing flowers](vignettes/13_bees_sharing_flowers/13_bees_sharing_flowers.md) ([PDF](vignettes/13_bees_sharing_flowers/13_bees_sharing_flowers.pdf)) | Four forms of competition for nectar and pollen, and virus spillover |

## Status

Early development. The API may change. The Varroa and virus parameters of the bee models
were checked against BEEHAVE; the competition and flower-transmission parameters are
provisional (see those vignettes).

## Data and licences

The package is MIT licensed. `data/verde_fish` is from Rogosch et al. (2019),
<https://doi.org/10.5281/zenodo.1309024> (MIT). `data/riparian/maybell.csv` is from
Lytle et al. (2017), <https://doi.org/10.6084/m9.figshare.4652608> (CC BY 4.0), with
line endings converted. `data/cypripedium` holds matrices estimated with the R packages
lefko3 and adapt3 from the datasets those packages distribute.
