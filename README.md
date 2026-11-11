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

- Declarative species: named survival and fecundity transitions with rates that are
  numbers or functions of the species' own state, the pool totals and the environment.
- Composition along Catlab undirected wiring diagrams (`compose`), or by pool name
  (`community`). Composites nest.
- Projection through any sequence of environments (`project`).
- Freezing a species or community at a state into an ordinary
  `MatrixProjectionModel` or block matrix (`lower`, `community_matrix`).
- Published models rebuilt and checked against their original implementations:
  the Verde River fish model (identical output), the Yampa River riparian guild model
  (within 4e-16, plus a corrected variant), and the adapt3 Cypripedium community
  (within 2e-15).

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

## Status

Early development. Planned: cross-species sensitivity analysis, equilibria and invasion
analysis, integral projection models as species, and metacommunities.

## Data and licences

The package is MIT licensed. `data/verde_fish` is from Rogosch et al. (2019),
<https://doi.org/10.5281/zenodo.1309024> (MIT). `data/riparian/maybell.csv` is from
Lytle et al. (2017), <https://doi.org/10.6084/m9.figshare.4652608>.
