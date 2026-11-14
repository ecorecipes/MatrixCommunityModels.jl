# Integral projection models as projection nets.
#
# A stage may carry a continuous trait (size, weight). It is then represented by mesh cells,
# which are ordinary stages of an ordinary net, so everything else in the package (responses,
# pools, composition, analysis) applies unchanged. The discretisation is the stack's own:
# `ContinuousDomain` from StructuredPopulationCore and `left_kan_extension` from
# CategoricalPopulationDynamics.

"""
    ContinuousStage(name, lower, upper, cells)
    ContinuousStage(name, domain::ContinuousDomain)

A stage whose individuals differ in a continuous trait (size, say) between `lower` and
`upper`. The trait is cut into `cells` mesh cells of equal width. In a net built with
[`kernel_net`](@ref) the stage becomes the stages `name1`, `name2`, ..., one per cell, each
holding the number of individuals in that cell; `name` itself stays usable as their total
(`ctx.n[:plant]`, `stage_total(sol, :thistle, :plant)`).

`meshpoints(stage)` are the trait values at the middle of the cells, `step_size(stage)` is
the cell width and `stage_names(stage)` are the cell names.
"""
struct ContinuousStage
    name::Symbol
    domain::ContinuousDomain{Float64}
    cells::Vector{Symbol}
    index::Dict{Symbol, Int}
    z::Vector{Float64}
end

function ContinuousStage(name::Symbol, domain::ContinuousDomain)
    d = ContinuousDomain(Float64(domain.lower), Float64(domain.upper), domain.n_meshpoints)
    cells = [Symbol(name, i) for i in 1:d.n_meshpoints]
    return ContinuousStage(name, d, cells, Dict(c => i for (i, c) in enumerate(cells)), collect(Float64, meshpoints(d)))
end
ContinuousStage(name::Symbol, lower::Real, upper::Real, cells::Integer) =
    ContinuousStage(name, ContinuousDomain(Float64(lower), Float64(upper), Int(cells)))

stage_names(s::ContinuousStage) = s.cells
meshpoints(s::ContinuousStage) = s.z
step_size(s::ContinuousStage) = step_size(s.domain)
Base.length(s::ContinuousStage) = length(s.cells)
Base.show(io::IO, s::ContinuousStage) =
    print(io, "ContinuousStage(:", s.name, ", ", s.domain.lower, " to ", s.domain.upper, ", ", length(s), " cells)")

# A kernel that is rebuilt every time step from the context.
struct ResponsiveKernel{F}
    f::F
end

"""
    responsive(f)

Mark a kernel of [`kernel_net`](@ref) as depending on the context: `f(ctx)` is called every
time step and returns the kernel to use for that step, in the form the entry would
otherwise take. Use it when the *shape* of a kernel depends on shared pools, on the
species' own state or on the environment, for example growth that slows under crowding:

```julia
:growing => responsive(ctx -> ((z′, z) -> survival(z) * growth(z′, z, ctx.pool.space)))
```

A dependence that only scales a kernel, by a number or by a factor for each trait value,
is cheaper as an ordinary response; see `respond` and [`by_trait`](@ref).
"""
responsive(f) = ResponsiveKernel(f)

const _KernelStage = Union{Symbol, ContinuousStage}
_stage_name(s::Symbol) = s
_stage_name(s::ContinuousStage) = s.name
_stage_cells(s::Symbol) = [s]
_stage_cells(s::ContinuousStage) = s.cells

function _kernel_error(to, from, wanted)
    throw(ArgumentError("the kernel from :$(_stage_name(from)) to :$(_stage_name(to)) must be $wanted"))
end

# The block of the projection matrix for one kernel, rows for the cells of `to`, columns for those of `from`.
# The state holds numbers of individuals per cell, so a density over the destination trait is multiplied
# by the destination cell width.
function _kernel_block(kernel, to::ContinuousStage, from::ContinuousStage, rule::Symbol)
    applicable(kernel, to.z[1], from.z[1]) || _kernel_error(to, from, "a function (z′, z) -> density")
    return left_kan_extension(kernel, to.domain, from.domain; rule) .* (step_size(to) / step_size(from))
end
function _kernel_block(kernel, to::ContinuousStage, from::Symbol, ::Symbol)
    applicable(kernel, to.z[1]) || _kernel_error(to, from, "a function z′ -> density")
    return reshape(Float64[step_size(to) * kernel(z) for z in to.z], :, 1)
end
function _kernel_block(kernel, to::Symbol, from::ContinuousStage, ::Symbol)
    applicable(kernel, from.z[1]) || _kernel_error(to, from, "a function z -> number")
    return reshape(Float64[kernel(z) for z in from.z], 1, :)
end
function _kernel_block(kernel, to::Symbol, from::Symbol, ::Symbol)
    kernel isa Real || _kernel_error(to, from, "a number")
    return fill(Float64(kernel), 1, 1)
end

"""
    kernel_net(stages, transition => kernels, ...; rule = :midpoint)
    kernel_net(stage::ContinuousStage, transition => kernel, ...; rule = :midpoint)

A projection net whose stages may carry a continuous trait: an integral projection model
written, like any other life cycle, as named transitions between stages.

`stages` lists the stages, each a `Symbol` (an ordinary stage) or a
[`ContinuousStage`](@ref). Each transition is `name => [(from => to) => kernel, ...]`, the
form used by `ValuedProjectionNet`, where the kernel is

- between two continuous stages, a function `(z′, z) -> density`: the density of
  individuals with trait `z′` next time step, per individual with trait `z` now;
- from an ordinary stage to a continuous one, a function `z′ -> density`;
- from a continuous stage to an ordinary one, a function `z -> number`;
- between two ordinary stages, a number.

With a single continuous stage a transition can be written `name => kernel`.

Kernels are evaluated at the middle of each cell (`rule`, passed to
`left_kan_extension`, can also be `:trapezoidal` or `:simpson` for kernels between
continuous stages). Individuals that a kernel sends outside the range of a trait are lost,
so the range should be wide enough to hold the population.

A kernel wrapped in [`responsive`](@ref) is rebuilt each time step from the context.

The result is a responsive projection net, so `respond`, `births`, `inflow`, `⊕`,
[`species`](@ref) and the rest apply. [`by_trait`](@ref) writes a response in terms of the
trait and [`trait_weights`](@ref) a pool contribution.

```julia
plant = ContinuousStage(:plant, 0.0, 10.0, 100)
life = kernel_net(plant,
    :growing => (z′, z) -> survival(z) * growth(z′, z),
    :seeding => (z′, z) -> seeds(z) * recruit_size(z′))
with_seed_bank = kernel_net([:seed, plant],
    :dormancy    => [(:seed => :seed) => 0.5],
    :germinating => [(:seed => :plant) => z′ -> 0.1 * recruit_size(z′)],
    :growing     => [(:plant => :plant) => (z′, z) -> survival(z) * growth(z′, z)],
    :seeding     => [(:plant => :seed) => z -> seeds(z)])
```
"""
function kernel_net(stages::AbstractVector, transitions::Pair{Symbol}...; rule::Symbol = :midpoint)
    isempty(stages) && throw(ArgumentError("kernel_net needs at least one stage"))
    all(s -> s isa _KernelStage, stages) ||
        throw(ArgumentError("stages must be Symbols or ContinuousStages"))
    names = Symbol[_stage_name(s) for s in stages]
    allunique(names) || throw(ArgumentError("stage names must be unique, got $names"))
    cells = reduce(vcat, [_stage_cells(s) for s in stages])
    allunique(cells) || throw(ArgumentError("a stage name clashes with the cell names of a continuous stage"))
    position = Dict(c => i for (i, c) in enumerate(cells))
    by_name = Dict{Symbol, Any}(_stage_name(s) => s for s in stages)
    ranges = Dict{Symbol, UnitRange{Int}}()
    start = 1
    for s in stages
        n = length(_stage_cells(s))
        ranges[_stage_name(s)] = start:(start + n - 1)
        start += n
    end
    k = length(cells)

    matrices = Pair{Symbol, Matrix{Float64}}[]
    responses = Pair{Symbol, Any}[]
    for (tname, spec) in transitions
        blocks = if spec isa AbstractVector
            spec
        elseif length(stages) == 1
            [(names[1] => names[1]) => spec]
        else
            throw(ArgumentError("with several stages, give transition :$tname as a list of (from => to) => kernel"))
        end
        M = zeros(Float64, k, k)
        varying = Any[]
        seen = Set{Pair{Symbol, Symbol}}()
        for block in blocks
            (block isa Pair && block.first isa Pair{Symbol, Symbol}) ||
                throw(ArgumentError("transition :$tname: expected (from => to) => kernel, got $block"))
            (from, to), kernel = block.first, block.second
            haskey(by_name, from) && haskey(by_name, to) ||
                throw(ArgumentError("transition :$tname goes from :$from to :$to, but the stages are $names"))
            (from => to) in seen && throw(ArgumentError("transition :$tname lists :$from => :$to twice"))
            push!(seen, from => to)
            rows, cols = ranges[to], ranges[from]
            if kernel isa ResponsiveKernel
                M[rows, cols] .= 1.0
                push!(varying, (rows, cols, by_name[to], by_name[from], kernel.f))
            else
                M[rows, cols] .= _kernel_block(kernel, by_name[to], by_name[from], rule)
            end
        end
        push!(matrices, tname => M)
        isempty(varying) || push!(responses, tname => _kernel_response(varying, position, k, rule))
    end
    groups = Dict{Symbol, Vector{Int}}(s.name => collect(ranges[s.name]) for s in stages if s isa ContinuousStage)
    return ResponsiveProjectionNet(matrix_net(cells, matrices...); responses, groups)
end
kernel_net(stage::ContinuousStage, transitions::Pair{Symbol}...; kwargs...) = kernel_net([stage], transitions...; kwargs...)

# The response of a transition with responsive kernels: its entries have the base value one, and are
# multiplied by this step's kernel values.
function _kernel_response(varying, position::Dict{Symbol, Int}, k::Int, rule::Symbol)
    return function (ctx)
        values = ones(Float64, k, k)
        for (rows, cols, to, from, f) in varying
            values[rows, cols] .= _kernel_block(f(ctx), to, from, rule)
        end
        return (from, to) -> values[position[to], position[from]]
    end
end

"""
    by_trait(stage, f)
    by_trait(stage, factors)

A response value that scales each entry of a transition according to the trait of the
individual it comes from: entries leaving the cell of `stage` at trait value `z` are
multiplied by `f(z)`, and entries leaving other stages are left alone. The factors can
also be given as a vector, one for each cell. Return it from a response:

```julia
respond(life, :seeding => ctx -> by_trait(plant, z -> exp(-0.02 * eggs_on(z, ctx))))
```
"""
by_trait(stage::ContinuousStage, f) = by_trait(stage, Float64[f(z) for z in stage.z])
function by_trait(stage::ContinuousStage, factors::AbstractVector{<:Real})
    length(factors) == length(stage) ||
        throw(DimensionMismatch("by_trait needs one factor for each of the $(length(stage)) cells of :$(stage.name)"))
    index = stage.index
    return function (from, to)
        i = get(index, from, 0)
        return i == 0 ? 1.0 : factors[i]
    end
end

"""
    trait_weights(stage, f)

Pool weights given as a function of the trait: each individual of `stage` with trait `z`
contributes `f(z)`. For the `occupies` argument of [`species`](@ref):

```julia
species(:tree, net; occupies = [:basal_area => trait_weights(tree, d -> pi * (d / 2)^2)])
```
"""
trait_weights(stage::ContinuousStage, f) = Pair{Symbol, Float64}[c => Float64(f(z)) for (c, z) in zip(stage.cells, stage.z)]
