# Composing and projecting communities

Species as open communities, their composition along undirected wiring diagrams, the
state of a community, its projection, and freezing it into ordinary matrices.

## Open communities and their composition

These are defined in CategoricalPopulationDynamics.jl and re-exported by this package.

```@docs
OpenCommunity
species
community
compose(::MatrixCommunityModels.CategoricalPopulationDynamics.UndirectedWiringDiagram, ::AbstractVector{<:OpenCommunity})
subcommunity
system
species_names
stage_names(::OpenCommunity, ::Symbol)
pool_names(::OpenCommunity)
diagram(::OpenCommunity)
state_vector
species_states
pool_totals(::OpenCommunity, ::AbstractVector, ::Any)
advance(::OpenCommunity, ::AbstractVector, ::Any)
```

## Projection and frozen matrices

```@autodocs
Modules = [MatrixCommunityModels]
Pages = ["src/MatrixCommunityModels.jl"]
```

## Recruitment into a limited resource

A lottery in which the recruits of a time step share the free space.

```@autodocs
Modules = [MatrixCommunityModels]
Pages = ["src/recruitment.jl"]
```
