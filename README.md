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
  (`community`). Composites nest. The result is a `PooledProjectionSystem`. The open
  communities and their composition are defined in CategoricalPopulationDynamics.jl and
  re-exported here.
- Integral projection models in the same communities: a stage may carry a continuous
  trait (`ContinuousStage`), with kernels as its named transitions (`kernel_net`), given as
  functions or as the `PKernel` and `FKernel` of IntegralProjectionModels.jl. Such a
  species composes with matrix models through the same pools.
- Projection (`project`) lowers the community to a coupled problem and solves it with the
  stack's own solvers: a community with a continuous stage is an integral projection model,
  lowered to a `CoupledIPMProblem` and solved by IntegralProjectionModels.jl; any other is
  lowered to a `CoupledMPMProblem` and solved by MatrixProjectionModels.jl. Both give the
  same projection, deterministic or with demographic stochasticity. A species frozen at any
  state is an ordinary `MatrixProjectionModel` (`matrix_model`), an `IPMProblem`
  (`ipm_problem`) or a net (`frozen_net`).
- Recruitment into a limited resource within a time step: `lottery` makes the recruits of
  a step share the free space (`min(1, free / demand)` or `free / (free + demand)`), with
  `lottery_demand` adding each species' recruits to a shared demand pool.
- Analysis: `equilibrium`, `jacobian`, `invasion_growth_rate`, `stochastic_invasion_rate`,
  the sensitivity and elasticity of an invasion growth rate to the invader's transitions,
  entries or parameters (`invasion_sensitivity`), analytic `equilibrium_sensitivity`, the
  simulation-based cross-species sensitivity of Lytle & Tonkin (`abundance_sensitivity`),
  `InteractionNetwork`, seasonal models (`periodic_matrix`), demographic stochasticity
  (`project(...; rng)`, `ensemble`).
- Resilience metrics through [ResilienceCore.jl](https://github.com/ecorecipes/ResilienceCore.jl):
  `linearise` gives the Jacobian at a fixed point or around a cycle as a `DiscreteJacobian`,
  `periodic_operator` a rare invader over a cycle of environments as a `PeriodicOperator`, and
  `MCMState` the community interface on which ResilienceCore's press, removal and extinction
  experiments run. `jacobian` and `linearise` are ResilienceCore's functions, extended here.
- Management and decisions: a `ManagedCommunity` gives a community named parameters (a
  `build` function maps the names to vital rates or matrix entries), an action as plain data
  (`remove` a share of a stage after each step, `take` a number, `add` a number, change
  `parameters` or the `environment`) and named outcomes; `project` runs it, deterministic or
  with demographic stochasticity. With
  [DecisionAnalysis.jl](https://github.com/ecorecipes/DecisionAnalysis.jl) loaded (a weak
  dependency, Julia 1.12 or later), an extension gives its adapter verbs methods, so that
  `outcomes(m, actions, draws; replicates, seed)` runs the community over parameter draws,
  actions and replicates with common random numbers, ready for payoff tables, decision rules
  and the value of information.
- Published models rebuilt and checked against their original implementations: the Verde
  River fish model (identical output), the Yampa River riparian guild model (within 4e-16,
  plus a corrected variant), and the adapt3 Cypripedium community (within 2e-15).
- A honey bee, *Varroa*, deformed wing virus and bumblebee community, built step by step.
- A native thistle (the IPM of Rose et al. 2005, whose equilibria are reproduced), musk
  thistle (the matrix model of Jongejans et al. 2008, growth rate 1.748 as published) and
  the flower-head weevil they share.
- Three European trees (Norway spruce, silver fir and beech; `forest_tree`) as IPMs on
  700 mesh cells each, sharing one pool of stand basal area: the model of the R package
  matreex, whose kernels are reproduced to 1e-9 and whose runs to within 1e-3 (1e-5 when
  matreex is run with finer numerical settings).
- A sea fan and its fungal disease (the IPM of Bruno et al. 2011, whose growth rates,
  elasticities and resistance scenario are reproduced): healthy, infected and resistant
  colonies are strata of one life cycle written by size.
- Four perennial plants of the Idaho sagebrush steppe (`sagebrush_plant`,
  `sagebrush_community`): the IPMs of Adler, Ellner & Levine (2010), coupled by
  size-dependent crowding through each species' cover, with the parameters of the R
  package IPMdoit, whose kernels, runs and invasion growth rates are reproduced; each
  species invades the others when rare, and without niche differences none does.
- The length-based multispecies fish model LeMans for 21 species of the North Sea
  (`lemans_community`; Hall et al. 2006, Thorpe et al. 2015): each species in length
  classes, eaten by the others according to predator and prey size through pools with a
  value for each predator length class, with the data and settings of the R package
  LeMaRns, whose runs are reproduced.
- Three coral genera recolonising four fore reefs of Mo'orea (`coral_community`): the
  open-population IPMs of Kayal et al. (2018), with recruits arriving as inflows in
  numbers set by crowding through pools. By default recruitment is read as in the authors'
  script; a reading of its units reconstructed from the published figures
  (`reading = :calibrated`) reproduces the published compositions of the recovering reefs.
- The dominant perennials of five grasslands of the western United States
  (`grassland_community`): the IPMs of Chu & Adler (2015) as run by Tredennick et al.
  (2017), built with the sagebrush steppe machinery; the authors' Kansas runs and the
  published synchrony of growth rates are reproduced.
- Two flour beetles that eat each other's eggs and pupae (`flour_beetle_competition`): the
  competition LPA model of Edmunds et al. (2003) and Desharnais et al. (2005), whose
  equilibria, eigenvalues, bifurcations, coexistence two-cycle and stochastic outcomes are
  reproduced.
- Riparian plants, fish and three aquatic invertebrates under flow regimes designed for
  each (Tonkin et al. 2021; `tonkin_scenario`, `flood_invertebrate`), checked against the
  authors' R code.
- A predator and its prey, each with juveniles and adults (`predator_prey`): the matrix
  model of Barraquand & Gimenez (2019), checked against the authors' difference equations.
- Leopard frogs and the chytrid fungus Bd at four locations (`chytrid_community`; the
  full host-parasite IPM of Wilber et al. 2022): tadpoles, susceptible frogs and infected
  frogs by log Bd load, with a zoospore component joined by pools, weekly steps, seasonal
  temperature and breeding; the authors' Python runs are reproduced to 1e-11.
- The size-spectrum model of the R package mizer for 12 North Sea species and a plankton
  resource (`mizer_community`; Scott, Blanchard & Andersen 2014): mizer's implicit upwind
  step written exactly as a projection matrix of each species, predation through pools
  with a value for each of 1200 predator classes; mizer's runs and steady states are
  reproduced to rounding.
- Eighteen California annual plants with a seed bank (`annual_plant_community`; Godoy,
  Kraft & Levine 2014): seeds and germinants of each species coupled by crowding pools of
  germinants weighted by the interaction coefficients; the niche overlap, fitness ratio and
  predicted outcome of all 153 pairs, checked against the paper's equations run in R.
- Canvasbacks and redheads, two diving ducks of the Prairie Pothole Region
  (`waterfowl_community`; Gibson et al. 2025): yearling and older females and males of each,
  driven year by year by the posterior means of the authors' integrated population model,
  checked against the authors' process equations and stored posterior draws; an extension
  joins fecundity to pools of both species' breeding females (conspecific density and the
  canvasback-to-redhead ratio of nest parasitism).

## Installation

The package and its ecorecipes dependencies are not registered. Its `Project.toml` names
their GitHub URLs, so this is enough (Julia 1.11 or later):

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/MatrixCommunityModels.jl")
```

The decision-analysis extension (see vignette 30) needs Julia 1.12 or later and two more
unregistered packages; without them the rest of the package is unaffected:

```julia
Pkg.add([PackageSpec(url = "https://github.com/ecorecipes/FiniteKernels.jl"),
         PackageSpec(url = "https://github.com/ecorecipes/DecisionAnalysis.jl")])
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
| 6 | [Composing IPMs](vignettes/06_composing_ipms/06_composing_ipms.md) ([PDF](vignettes/06_composing_ipms/06_composing_ipms.pdf)) | An integral projection model and a matrix model sharing a pool; kernels that respond to pools; a seed bank |
| 7 | [Metacommunity](vignettes/07_metacommunity/07_metacommunity.md) ([PDF](vignettes/07_metacommunity/07_metacommunity.pdf)) | Patches as nested communities, dispersal as a pool |
| 8 | [Equilibria, invasion and chance](vignettes/08_equilibria_invasion_chance/08_equilibria_invasion_chance.md) ([PDF](vignettes/08_equilibria_invasion_chance/08_equilibria_invasion_chance.pdf)) | Analytic sensitivities, invasion rates in fluctuating environments, sensitivity of an invasion rate, demographic stochasticity |
| 9 | [Seasons and phenology](vignettes/09_seasons_and_phenology/09_seasons_and_phenology.md) ([PDF](vignettes/09_seasons_and_phenology/09_seasons_and_phenology.pdf)) | Within-year steps, a lottery for space, coexistence by timing, drift |
| 10 | [Bees 1: a honey bee colony](vignettes/10_honeybee_colony/10_honeybee_colony.md) ([PDF](vignettes/10_honeybee_colony/10_honeybee_colony.pdf)) | A colony from modules: life cycle, stage durations, responses, stores |
| 11 | [Bees 2: Varroa and deformed wing virus](vignettes/11_varroa_and_dwv/11_varroa_and_dwv.md) ([PDF](vignettes/11_varroa_and_dwv/11_varroa_and_dwv.pdf)) | Mites as a second species, the virus as a stratum, treatment timing |
| 12 | [Bees 3: a bumblebee population](vignettes/12_bumblebees/12_bumblebees.md) ([PDF](vignettes/12_bumblebees/12_bumblebees.pdf)) | The annual colony cycle from three modules, nest sites |
| 13 | [Bees 4: sharing flowers](vignettes/13_bees_sharing_flowers/13_bees_sharing_flowers.md) ([PDF](vignettes/13_bees_sharing_flowers/13_bees_sharing_flowers.pdf)) | Four forms of competition for nectar and pollen, and virus spillover |
| 14 | [Thistles and a weevil](vignettes/14_thistles_and_weevil/14_thistles_and_weevil.md) ([PDF](vignettes/14_thistles_and_weevil/14_thistles_and_weevil.pdf)) | A published IPM of a native thistle, a published matrix model of a weed and the weevil they share |
| 15 | [A forest on a large mesh](vignettes/15_forest_basal_area/15_forest_basal_area.md) ([PDF](vignettes/15_forest_basal_area/15_forest_basal_area.pdf)) | Three tree species as IPMs on 700 cells each, coupled through stand basal area and checked against the R package matreex |
| 16 | [Sea fans and aspergillosis](vignettes/16_sea_fan_aspergillosis/16_sea_fan_aspergillosis.md) ([PDF](vignettes/16_sea_fan_aspergillosis/16_sea_fan_aspergillosis.pdf)) | A published IPM stratified by infection status: strata, responses by size and pools together |
| 17 | [Four perennials of the sagebrush steppe](vignettes/17_sagebrush_steppe/17_sagebrush_steppe.md) ([PDF](vignettes/17_sagebrush_steppe/17_sagebrush_steppe.pdf)) | Four IPMs coupled by crowding through their covers, checked against the R package IPMdoit; coexistence, invasion when rare and the removal of niche differences |
| 18 | [North Sea fish by length](vignettes/18_north_sea_lemans/18_north_sea_lemans.md) ([PDF](vignettes/18_north_sea_lemans/18_north_sea_lemans.pdf)) | The LeMans model of 21 species in 32 length classes, with predation by predator and prey size through pools indexed by size, checked against the R package LeMaRns; the unfished community and fishing by four fleets |
| 19 | [Corals of Mo'orea](vignettes/19_corals_moorea/19_corals_moorea.md) ([PDF](vignettes/19_corals_moorea/19_corals_moorea.pdf)) | Three open-population coral IPMs whose recruits arrive as inflows set by crowding, checked against the authors' kernel functions and the published simulations of reef recovery |
| 20 | [Five grasslands](vignettes/20_grasslands_synchrony/20_grasslands_synchrony.md) ([PDF](vignettes/20_grasslands_synchrony/20_grasslands_synchrony.pdf)) | Multispecies IPMs of five sites from the sagebrush machinery, the Kansas prairie checked against the authors' code; coexistence, invasion in drawn years and the synchrony of growth rates |
| 21 | [Two flour beetles](vignettes/21_flour_beetles/21_flour_beetles.md) ([PDF](vignettes/21_flour_beetles/21_flour_beetles.pdf)) | Two LPA models coupled by eating each other's eggs and pupae; the published equilibria, bifurcations, coexistence two-cycle and stochastic outcomes reproduced |
| 22 | [Designed flow regimes](vignettes/22_designed_flow_regimes/22_designed_flow_regimes.md) ([PDF](vignettes/22_designed_flow_regimes/22_designed_flow_regimes.pdf)) | Riparian plants, fish and three invertebrates under flows designed for each and the natural regime, checked against the authors' R code; trade-offs and the search for designed flows |
| 23 | [Predator and prey](vignettes/23_predator_and_prey/23_predator_and_prey.md) ([PDF](vignettes/23_predator_and_prey/23_predator_and_prey.pdf)) | Two species with juveniles and adults linked by predation on juveniles, checked against the authors' difference equations; fixed points, oscillations and demographic noise |
| 24 | [Chytrid fungus in leopard frogs](vignettes/24_chytrid_leopard_frogs/24_chytrid_leopard_frogs.md) ([PDF](vignettes/24_chytrid_leopard_frogs/24_chytrid_leopard_frogs.pdf)) | The full host-parasite IPM of Wilber et al. (2022): Bd load as the trait, a zoospore component joined by pools, seasonal temperature and breeding at four locations, checked against the authors' Python code; persistence with demographic stochasticity |
| 25 | [A size-spectrum fish community](vignettes/25_size_spectrum_mizer/25_size_spectrum_mizer.md) ([PDF](vignettes/25_size_spectrum_mizer/25_size_spectrum_mizer.pdf)) | The North Sea model of the R package mizer, its implicit upwind step written as a projection matrix, checked against mizer's runs and steady states; fishing, the size spectrum and a comparison with LeMans |
| 26 | [Resilience metrics](vignettes/26_resilience/26_resilience.md) ([PDF](vignettes/26_resilience/26_resilience.pdf)) | ResilienceCore.jl on matrix communities: a reactive beetle equilibrium, the stability of cycles on the way to chaos, presses and removals in competing beetles and in riparian guilds under a constant flood regime |
| 28 | [Canvasbacks and redheads](vignettes/28_canvasbacks_redheads/28_canvasbacks_redheads.md) ([PDF](vignettes/28_canvasbacks_redheads/28_canvasbacks_redheads.pdf)) | The process model of a two-species integrated population model of diving ducks, run on the posterior-mean yearly rates and checked against the authors' equations, stored draws and transient LTRE; a closed-loop extension through pools of breeding females |
| 29 | [Annual plants with a seed bank](vignettes/29_annual_plants/29_annual_plants.md) ([PDF](vignettes/29_annual_plants/29_annual_plants.pdf)) | Eighteen annual plants of Godoy et al. (2014) as two-stage nets joined by crowding pools; niche and fitness differences of all pairs, predicted and simulated outcomes, and competitive ranks, checked against the paper's equations in R |
| 30 | [Which stage to harvest, and how hard](vignettes/30_harvest_decision/30_harvest_decision.md) ([PDF](vignettes/30_harvest_decision/30_harvest_decision.pdf)) | A managed community and the DecisionAnalysis.jl extension: harvesting the prey of Barraquand & Gimenez (2019) by share or by quota under uncertain vital rates, the decision under several rules, and the value of measuring each rate (illustrative numbers) |
| 31 | [Choosing a river flow regime](vignettes/31_verde_flow_decision/31_verde_flow_decision.md) ([PDF](vignettes/31_verde_flow_decision/31_verde_flow_decision.pdf)) | Environmental flow regime decisions on the Verde River across native fish, riparian cottonwood and aquatic invertebrates under hydrological and carrying-capacity uncertainty: Pareto frontier, stochastic dominance, smooth ambiguity, certainty-equivalent EVPI and PRIM scenario discovery (illustrative uncertainty and rules) |

## Status

Early development. The API may change. The Varroa and virus rules and parameters of the bee
models follow BEEHAVE and were checked against its code. The floral supply rests on
published measurements; the overlap, refuge and virus-on-flowers rates are assumptions
(see vignette 13). In the thistle community the two plant models are published; the
weevil's winter survival is calibrated, and how the species share weevils rests on
assumptions listed in vignette 14. The forest models take their vital rates from the R
package matreex and reproduce its kernels and runs; they leave out its harvesting rules,
disturbances and changing climates, and vignette 15 lists where they and their sources
part.
The sea fan model is published; the size of its
recruits is not given in the paper and is assumed, and infection passed among fans is an
extension (see vignette 16).
The sagebrush steppe models take their parameters and numerical scheme from the R package
IPMdoit and reproduce it; how niche differences are removed is an assumption of this
package, since the papers could not be read (see vignette 17).
The LeMans model takes its equations, settings and North Sea data from the R package
LeMaRns and reproduces its runs; it leaves out LeMaRns's other recruitment and mortality
options and gives each species a single fleet (see vignette 18).
The coral models use the published parameters of Kayal et al. (2018); the units of their
recruitment function are not stated. The default follows the authors' script, which does
not reproduce the published compositions; the calibrated reading reproduces them, but not
the published elasticities to recruitment (see vignette 19).
The grassland models take their parameters from the repository of Tredennick et al.
(2017) and reproduce its Kansas runs; they are not the refitted values printed by Chu &
Adler (2015) (see vignette 20).

The flour beetle parameters are those printed in the papers, which chose them rather than
fitting them; the two-species adult equation of Edmunds et al. (2003) is printed with a
minus sign where a plus is meant (see vignette 21). The invertebrates of the flow-regime
comparison are unstructured logistic populations taken from the authors' code, since
McMullen et al. (2017) could not be read (see vignette 22). The predator-prey parameters are
the simulation values of Barraquand & Gimenez (2019); fitting the integrated community
model is not included (see vignette 23).

The chytrid model takes its equations and fitted values from the authors' Python code for
Wilber et al. (2022) and reproduces its runs; it is the full IPM of their persistence
analysis, not the reduced IPM they fit (see vignette 24).
The mizer model takes its equations, default scheme and North Sea parameters from the R
package mizer 3.4.0 and reproduces its runs; mizer's optional second-order methods are not
written (see vignette 25).
The annual plant model takes its vital rates and interaction coefficients from the authors'
figshare deposit, which do not reproduce the paper's counts of coexisting pairs (10 here, 12
printed), and whose seed production differs from the paper's Table 3 (see vignette 29).

The harvest decision of vignette 30 uses the published predator-prey model, but its prior,
harvest rates, horizon and survival study are illustrative and describe no real population.

The tests (`Pkg.test()`, about seven minutes) compare the models with the authors' code and
published values. Checks that need long runs or ensembles (`test/long`: the five grassland
sites over 2500 years, the forests on 700 cells, the sagebrush runs and exclusion times, the
chytrid runs over 529 weeks, the seeded ensembles of the vignettes) are left out unless asked
for: `MCM_LONG_TESTS=true julia --project -e 'using Pkg; Pkg.test()'` runs them after the
main suite. They also run every week on GitHub Actions (`LongTests.yml`).
The DecisionAnalysis.jl extension is tested in an environment of its own
(`test/decision_analysis`), which `Pkg.test()` runs in a second process on Julia 1.12 or
later; `MCM_DECISION_TESTS=false` leaves it out.

## Data and licences

The package is under mixed licences: the code and most data are MIT, two data
directories, `data/lemans` and `data/mizer`, are GPL-3, and two files of `data/corals` are CC BY 4.0 (see below). `data/verde_fish` is from Rogosch et al. (2019),
<https://doi.org/10.5281/zenodo.1309024> (MIT). `data/riparian/maybell.csv` is from
Lytle et al. (2017), <https://doi.org/10.6084/m9.figshare.4652608> (CC BY 4.0), with
line endings converted. `data/cypripedium` holds matrices estimated with the R packages
lefko3 and adapt3 from the datasets those packages distribute. The fitted vital rates in
`src/forests.jl` are those shipped with the R package matreex 0.4.0
(<https://github.com/gowachin/matreex>, MIT, © 2022 Maxime Jaunatre; its notice is kept in
that file), and `test/reference/forest_*.csv` were produced by running it.
`data/sagebrush` holds parameters copied or extracted from the R package IPMdoit
(<https://github.com/atredennick/IPMdoit>, MIT, © 2015 Andrew Tredennick; its notice is in
`data/sagebrush/LICENSE.txt`), and `test/reference/sagebrush_*.csv` were produced by
running its functions.
`data/lemans` holds the North Sea data sets of the R package LeMaRns 0.1.2
(<https://CRAN.R-project.org/package=LeMaRns>), which is distributed under the **GPL-3**,
not MIT; their provenance and licence are in `data/lemans/LICENSE.txt`.
`test/reference/lemans_*.csv` were produced by running LeMaRns.
`data/corals` holds the parameter estimates of Table S1 of Kayal et al. (2018), typed from
the paper's Supporting Information; the authors' R code carries no licence and is not
included. `data/corals/mcr_lter_observed_colonies.csv` holds the coral colonies observed in 2013
and 2015, binned by size, and `data/corals/mcr_lter_observed_transects.csv` the same colonies
totalled by transect, both derived from MCR LTER data on the Environmental Data Initiative
(knb-lter-mcr.4009.5, <https://doi.org/10.6073/pasta/322a76c6f197bfac60c91d66a3e96fb5>,
**CC BY 4.0**, <https://creativecommons.org/licenses/by/4.0/>). Those two files are under the
same licence; the changes made are listed in `data/corals/LICENSE.txt`.
`test/reference/coral_*.csv` were produced by running the authors' kernel functions, read at run time.
`data/grasslands` holds the fitted parameters of Tredennick et al. (2017), written out from
<https://github.com/atredennick/community_synchrony>, whose DESCRIPTION says "License: MIT"
but which has no LICENSE file; the fits are to the data of Chu & Adler (2015) on Dryad
(CC0). Provenance and notice are in `data/grasslands/LICENSE.txt`.
`test/reference/grassland_*.csv` were produced by running the authors' simulation code.

`test/reference/flow_*.csv` were produced by running the R code of Tonkin et al. (2021)
(<https://github.com/jdtonkin/flow-tradeoffs>, MIT), whose invertebrate parameters are in
`src/flowregimes.jl`. `test/reference/predprey_trajectories.csv` was produced by running the
R code of Barraquand & Gimenez (2019) (<https://github.com/oliviergimenez/predator_prey_icm>,
GPL-3); no code from it is in the package, whose model is written from the paper.

`data/mizer` holds the North Sea parameters of the R package mizer 3.4.0
(<https://CRAN.R-project.org/package=mizer>), also distributed under the **GPL-3**; their
provenance and licence are in `data/mizer/LICENSE.txt`, and `test/reference/mizer_*.csv`
were produced by running mizer. The values in `src/chytrid.jl` are fitted and fixed
parameter values of Wilber et al. (2022), taken from the authors' code
(<https://doi.org/10.5281/zenodo.6407088>), which carries no licence file (its Zenodo
record says "other-open"); none of that code is copied, and
`test/reference/chytrid_*.csv` were produced by running it.
`data/annual_plants` holds the vital rates and interaction coefficients of Godoy, Kraft &
Levine (2014), from the authors' figshare deposit (<https://doi.org/10.6084/m9.figshare.4793488>,
CC0) with line endings converted, and the common competitive effects and responses of the
paper's Appendix S2, Table 2, typed from the paper; see `data/annual_plants/LICENSE.txt`.
`test/reference/annualplants_*.csv` were produced by `test/reference/R/run_annualplants_reference.R`,
which writes the paper's equations in base R.
`data/waterfowl` holds posterior means and data summaries derived from the Dryad dataset of
Gibson et al. (2025) (<https://doi.org/10.5061/dryad.hqbzkh1n9>, CC0); the files here are CC0
too, and the changes are listed in `data/waterfowl/LICENSE.txt`. The model equations and the
LTRE of `src/waterfowl.jl` follow the authors' code (<https://doi.org/10.5281/zenodo.8076163>,
MIT; its notice is in `data/waterfowl/LICENSE.txt`), and `test/reference/waterfowl_*.csv`
were produced by evaluating it, read at run time.
