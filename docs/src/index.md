# MatrixCommunityModels.jl

Matrix community models (MCMs; Lytle & Tonkin 2023) are sets of single-species matrix
population models linked by aggregate density dependence. This package builds them by
composition: each species is a projection net of named transitions, with responses to
shared pools (space, biomass, aggregate density), and a community is the composite of
species along an undirected wiring diagram whose junctions are the pools. Communities are
projected by MatrixProjectionModels.jl, or, when a species has a continuous stage (an
integral projection model), by IntegralProjectionModels.jl.

It is built on the [ecorecipes](https://github.com/ecorecipes) packages
MatrixProjectionModels.jl, IntegralProjectionModels.jl, ProjectionModels.jl and
CategoricalPopulationDynamics.jl, and on Catlab.jl.

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
  matrix model sharing a pool, kernels that respond to pools, and a seed bank.
- [Metacommunity](tutorials/07_metacommunity.md): patches as nested communities, with
  dispersal as a pool.
- [Equilibria, invasion and chance](tutorials/08_equilibria_invasion_chance.md): analytic
  sensitivities, invasion rates in fluctuating environments, what would let an excluded
  species persist (the sensitivity of its invasion rate), demographic stochasticity.
- [Seasons and phenology](tutorials/09_seasons_and_phenology.md): within-year steps, a
  lottery for space, coexistence by timing, and drift.
- Bees: [a honey bee colony](tutorials/10_honeybee_colony.md),
  [Varroa and deformed wing virus](tutorials/11_varroa_and_dwv.md),
  [a bumblebee population](tutorials/12_bumblebees.md) and
  [sharing flowers](tutorials/13_bees_sharing_flowers.md): a community built step by step.
- [Thistles and a weevil](tutorials/14_thistles_and_weevil.md): a published integral
  projection model of a native thistle, a published matrix model of an invasive thistle,
  and the weevil they share.
- [A forest on a large mesh](tutorials/15_forest_basal_area.md): three tree species as
  integral projection models on 700 cells each, coupled through stand basal area and
  checked against the R package matreex.
- [Sea fans and aspergillosis](tutorials/16_sea_fan_aspergillosis.md): a published
  integral projection model stratified by infection status, with healthy, infected and
  resistant colonies as strata of one life cycle.
- [Four perennials of the sagebrush steppe](tutorials/17_sagebrush_steppe.md): the
  integral projection models of Adler, Ellner & Levine (2010), coupled by crowding through
  each species' cover and checked against the R package IPMdoit; coexistence, invasion when
  rare and the removal of niche differences.
- [North Sea fish by length](tutorials/18_north_sea_lemans.md): the length-based
  multispecies model LeMans (Hall et al. 2006) for 21 species, with predation by predator
  and prey size through pools that hold a value for each predator length class, checked
  against the R package LeMaRns; the unfished community and fishing by four fleets.
- [Corals of Mo'orea](tutorials/19_corals_moorea.md): the open-population integral
  projection models of Kayal et al. (2018) for three coral genera recolonising four reefs,
  with recruits arriving as inflows in numbers set by crowding, checked against the
  authors' kernel functions and the published simulations of reef recovery.
- [Five grasslands](tutorials/20_grasslands_synchrony.md): the multispecies integral
  projection models of Chu & Adler (2015) for five sites, built with the sagebrush
  machinery and checked against the authors' Kansas runs; coexistence, invasion in drawn
  years and the synchrony of growth rates (Tredennick et al. 2017).
- [Two flour beetles](tutorials/21_flour_beetles.md): the competition LPA model of Edmunds
  et al. (2003) and Desharnais et al. (2005), two larva-pupa-adult models coupled by eating
  each other's eggs and pupae; the published equilibria, bifurcations, coexistence two-cycle
  and stochastic outcomes.
- [Designed flow regimes](tutorials/22_designed_flow_regimes.md): riparian plants, fish and
  three aquatic invertebrates under flows designed for each and under the natural regime
  (Tonkin et al. 2021), checked against the authors' R code.
- [Predator and prey](tutorials/23_predator_and_prey.md): the two-species, two-stage
  matrix model of Barraquand & Gimenez (2019), checked against the authors' difference
  equations; fixed points, oscillations and demographic noise.
- [Chytrid fungus in leopard frogs](tutorials/24_chytrid_leopard_frogs.md): the full
  host-parasite integral projection model of Wilber et al. (2022), with Bd load as the
  continuous trait, a zoospore component joined by pools, and seasonal temperature and
  breeding at four locations, checked against the authors' Python code; persistence with
  demographic stochasticity.
- [A size-spectrum fish community](tutorials/25_size_spectrum_mizer.md): the North Sea
  model of the R package mizer, its implicit upwind step written exactly as a projection
  matrix, with predation through pools that hold a value for each of 1200 predator
  classes, checked against mizer's runs and steady states; fishing, the size spectrum and
  a comparison with LeMans.
- [Which stage to harvest, and how hard](tutorials/30_harvest_decision.md): a managed
  community (named parameters, an action, named outcomes) and the extension for
  DecisionAnalysis.jl: a harvest of the prey of Barraquand & Gimenez (2019) by share or by
  quota under uncertain vital rates, the decision under several rules, and the value of
  measuring each rate. Every number of the decision problem is illustrative.
- [Choosing a river flow regime across fish, trees and insects under uncertainty](tutorials/31_verde_flow_decision.md):
  environmental flow regimes on the Verde River across native fish, riparian cottonwood
  and aquatic invertebrates under hydrological and carrying-capacity uncertainty
  (`ParetoFront`, `StochasticDominance`, `SmoothAmbiguity`, `ce_evpi`, `prim_boxes`).
- [Managing the North Sea multispecies fishery under recruitment and predation uncertainty](tutorials/32_north_sea_fisheries_decision.md):
  six fleet-effort portfolios on the 21-species North Sea LeMans size-spectrum model
  (`ParetoFront`, `CVaR`, FSD/SSD `StochasticDominance`, `SplineRegression` and
  `GaussianProcessRegression` EVPPI, `prim_boxes` cod-collapse discovery).
- [Managing Varroa, virus spillover and floral competition between honey bees and bumblebees](tutorials/33_pollinator_varroa_decision.md):
  miticide treatment, apiary density caps and wildflower strips in the coupled bee–mite–virus
  community (`InfectionUtility`, `AdditiveUtility`, `ParetoFront`, `survey_or_act`, EVPPI,
  `CredalSet`, `prune_generators`, `SmoothAmbiguity` over all four floral-competition couplings).

Each tutorial is a Quarto vignette, also available as
[HTML and PDF in the repository](https://github.com/ecorecipes/MatrixCommunityModels.jl/tree/main/vignettes).

## Quick start

```julia
using MatrixCommunityModels

K = 1000.0
crowding(ctx) = max(0.0, 1 - ctx.pool.space / K)

# a life cycle is a net of named transitions
grass_life = ValuedProjectionNet([:seedling, :adult],
    :establishing => [(:seedling => :adult) => 0.3],
    :persisting   => [(:adult => :adult) => 0.5],
    :seeding      => [(:adult => :seedling) => 6.0])
shrub_life = ValuedProjectionNet([:seedling, :adult],
    :establishing => [(:seedling => :adult) => 0.1],
    :persisting   => [(:adult => :adult) => 0.95],
    :seeding      => [(:adult => :seedling) => 3.0])

# a response says how a transition depends on a shared pool
grass = species(:grass, respond(grass_life, :seeding => crowding); occupies = [:space => [:seedling => 0.1, :adult => 1.0]])
shrub = species(:shrub, respond(shrub_life, :seeding => crowding); occupies = [:space => [:seedling => 0.1, :adult => 4.0]])

com = community(grass, shrub)
sol = project(com, [:grass => [50.0, 10.0], :shrub => [5.0, 1.0]], fill(nothing, 200))
sol[:shrub]            # stage-by-time matrix
```
