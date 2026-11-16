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
the cell width and `stage_names(stage)` are the cell names. In a response, `ctx.n[stage]`
is the vector of numbers in the cells, in the order of `meshpoints(stage)`.
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

"""
    n[stage::ContinuousStage]

The numbers of individuals in the cells of `stage`, one for each cell in the order of
`meshpoints(stage)`, read from a stage vector such as `ctx.n` in a response. (`ctx.n[:plant]`,
with the stage's name, is their total.)

The cells are found by name, not by position. The result is therefore right whatever other
stages the net holds, and it stays right when the net is refined with `stratify` or
`expand_stages`: each entry is then the total over the stages that descend from that cell.
"""
Base.getindex(n::StageVector, stage::ContinuousStage) = Float64[n[c] for c in stage.cells]

"""
    n[stage::ContinuousStage, stratum]

The numbers of individuals in the cells of `stage` that are in one `stratum` of a net
made with `stratify`, in the order of `meshpoints(stage)`: `ctx.n[plant, :infected]` are
the infected plants by size, where `ctx.n[plant]` are all plants by size. The cells are
found by the names `stratify` gives them (`stratum_stage`, with its default separator).
"""
Base.getindex(n::StageVector, stage::ContinuousStage, stratum::Symbol) =
    Float64[n[stratum_stage(c, stratum)] for c in stage.cells]

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
is cheaper as an ordinary response; see `respond` and [`by_trait`](@ref). A kernel that
depends on one pool alone can be built once at chosen values of the pool and interpolated
between them; see [`tabulated`](@ref).
"""
responsive(f) = ResponsiveKernel(f)

# A kernel that depends on the context only through one pool: built at the grid values of the pool, and
# interpolated between them.
struct TabulatedKernel{G}
    g::G
    pool::Symbol
    grid::Vector{Float64}
    outside::Symbol
end

"""
    tabulated(g; pool, grid, outside = :error)

Like [`responsive`](@ref), a kernel of [`kernel_net`](@ref) whose shape changes from one
time step to the next, but one that depends on the context only through the pool `pool`:
`g(x)` returns the kernel (in the form the entry would otherwise take) when the pool
holds `x`. Rather than being rebuilt every step, the kernel is built once at each value of
`grid` (an increasing list), the first time that value is needed, and kept. Each step the
pool's value `x` falls between two neighbouring grid values `a ≤ x ≤ b`, and the kernel
used is the linear interpolation `(1 - w) K(a) + w K(b)`, `w = (x - a) / (b - a)`, entry
by entry. At a grid value it is the kernel built there. This is how the R package matreex
runs its forests (its grid is the whole basal areas).

A pool value outside the grid, `x < first(grid)` or `x > last(grid)`, is handled according
to `outside`:

- `:error` (the default): throw a `DomainError`, so that the grid can be widened;
- `:clamp`: use the kernel at the nearest end of the grid;
- `:exact`: build the kernel at `x` for that step, as `responsive` would.

A pool value that is not a number (`NaN`) is always an error.

The interpolation error is of the order of the square of the grid spacing times the
curvature of the kernel in `x`; with a growth kernel whose mean shifts with the pool, a
spacing small against the width of the kernel keeps it small. The kernels are kept as the
block of the projection matrix they fill, and only the band of each column that is not
zero is stored, so a growth kernel that moves individuals a few cells takes much less
memory than a dense block (the growth kernels of [`forest_tree`](@ref) on 700 cells take
0.2 to 0.9 MB per grid value, against 3.9 MB). The kept kernels are shared by every copy
of the net and by every species built from it.

```julia
:growing => tabulated(ba -> ((z′, z) -> survival(z, ba) * growth(z′, z, ba)); pool = :basal_area, grid = 0:1:80)
```
"""
function tabulated(g; pool::Symbol, grid, outside::Symbol = :error)
    values = collect(Float64, grid)
    length(values) >= 2 || throw(ArgumentError("a tabulated kernel needs at least two grid values"))
    all(isfinite, values) || throw(ArgumentError("the grid of a tabulated kernel must be finite numbers"))
    all(i -> values[i] < values[i + 1], 1:(length(values) - 1)) ||
        throw(ArgumentError("the grid of a tabulated kernel must be increasing"))
    outside in (:error, :clamp, :exact) ||
        throw(ArgumentError("outside must be :error, :clamp or :exact, got :$outside"))
    return TabulatedKernel(g, pool, values, outside)
end

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
    A = left_kan_extension(kernel, to.domain, from.domain; rule)     # a new matrix, scaled in place
    ratio = step_size(to) / step_size(from)
    ratio == 1 || (A .*= ratio)
    return A
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

# A survival-growth or fecundity kernel of IntegralProjectionModels.jl, with any eviction correction it
# carries. It is discretised there, by the midpoint rule on its own domain, which must be the mesh of the
# continuous stages it joins.
const _SubKernel = Union{IntegralProjectionModels.PKernel, IntegralProjectionModels.FKernel}

function _kernel_block(kernel::_SubKernel, to::ContinuousStage, from::ContinuousStage, rule::Symbol)
    rule === :midpoint || throw(ArgumentError("a $(nameof(typeof(kernel))) of IntegralProjectionModels.jl is " *
        "discretised by the midpoint rule; it cannot be used with rule = :$rule"))
    d = kernel.domain
    same(a) = (Float64(d.lower), Float64(d.upper), d.n_meshpoints) == (a.domain.lower, a.domain.upper, a.domain.n_meshpoints)
    same(to) && same(from) || throw(ArgumentError("the $(nameof(typeof(kernel))) from :$(from.name) to :$(to.name) " *
        "is on the domain $(d.lower) to $(d.upper) with $(d.n_meshpoints) cells, which is not the mesh of both stages"))
    return Matrix{Float64}(IntegralProjectionModels.materialize(kernel))
end
# elsewhere it is an error: such a kernel goes between two continuous stages
_subkernel_error(kernel, to, from) = _kernel_error(to, from,
    "a function or a number: a $(nameof(typeof(kernel))) of IntegralProjectionModels.jl goes between two continuous stages")
_kernel_block(kernel::_SubKernel, to::ContinuousStage, from::Symbol, ::Symbol) = _subkernel_error(kernel, to, from)
_kernel_block(kernel::_SubKernel, to::Symbol, from::ContinuousStage, ::Symbol) = _subkernel_error(kernel, to, from)
_kernel_block(kernel::_SubKernel, to::Symbol, from::Symbol, ::Symbol) = _subkernel_error(kernel, to, from)

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

Between two continuous stages on the same mesh the kernel can also be a `PKernel` or an
`FKernel` of IntegralProjectionModels.jl (with any eviction correction it carries), whose
domain is that mesh; it is discretised there by `materialize`, which uses the midpoint rule,
so `rule` must be `:midpoint`.

With a single continuous stage a transition can be written `name => kernel`.

Kernels are evaluated at the middle of each cell (`rule`, passed to
`left_kan_extension`, can also be `:trapezoidal` or `:simpson` for kernels between
continuous stages). Individuals that a kernel sends outside the range of a trait are lost,
so the range should be wide enough to hold the population.

A kernel wrapped in [`responsive`](@ref) is rebuilt each time step from the context; one
wrapped in [`tabulated`](@ref) is interpolated between kernels built at values of a pool. The
kernels they return may be any of the above, `PKernel` and `FKernel` included.

The continuous stages are kept with the net as its meshes (`meshes(net)`), so a community
that holds the species is solved as a coupled IPM by IntegralProjectionModels.jl (see
[`project`](@ref)).

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
            if kernel isa Union{ResponsiveKernel, TabulatedKernel}
                M[rows, cols] .= 1.0
                push!(varying, (rows, cols, _kernel_source(kernel, by_name[to], by_name[from], rule)))
            else
                M[rows, cols] .= _kernel_block(kernel, by_name[to], by_name[from], rule)
            end
        end
        push!(matrices, tname => M)
        isempty(varying) || push!(responses, tname => _kernel_response(varying, k))
    end
    groups = Dict{Symbol, Vector{Int}}(s.name => collect(ranges[s.name]) for s in stages if s isa ContinuousStage)
    # the continuous stages are recorded as meshes, so that the net can be lowered to a coupled IPM
    meshes = [(name = s.name, domain = s.domain, cells = s.cells) for s in stages if s isa ContinuousStage]
    return ResponsiveProjectionNet(matrix_net(cells, matrices...); responses, groups, meshes)
end
kernel_net(stage::ContinuousStage, transitions::Pair{Symbol}...; kwargs...) = kernel_net([stage], transitions...; kwargs...)

# The response of a transition with responsive or tabulated kernels: its entries have the base value one,
# and are multiplied by this step's kernel values, returned as a matrix laid out like the projection
# matrix. The matrix is lazy (`_KernelValues`): upstream reads it only at the entries of the transition, so
# it holds just this step's blocks, and is one elsewhere.
function _kernel_response(varying, k::Int)
    sources = [(rows, cols, source) for (rows, cols, source) in varying]
    return function (ctx)
        return _KernelValues(k, Tuple((rows, cols, _block(source, ctx)) for (rows, cols, source) in sources))
    end
end

# A kernel rebuilt every step from the context.
struct _Rebuilt{F, T, S}
    f::F
    to::T
    from::S
    rule::Symbol
end
_kernel_source(kernel::ResponsiveKernel, to, from, rule) = _Rebuilt(kernel.f, to, from, rule)
_block(s::_Rebuilt, ctx) = _kernel_block(s.f(ctx), s.to, s.from, s.rule)

# A block kept column by column: in column `j`, the rows `first[j]:last[j]` (an empty range if the column
# is zero), whose values start at `values[offset[j] + first[j]]`. Entries outside are zero.
struct _Envelope
    first::Vector{Int}
    last::Vector{Int}
    offset::Vector{Int}
    values::Vector{Float64}
end
function _Envelope(B::AbstractMatrix{<:Real})
    m, n = size(B)
    lo, hi, offset = ones(Int, n), zeros(Int, n), zeros(Int, n)
    values = Float64[]
    for j in 1:n
        column = view(B, :, j)
        a = findfirst(!iszero, column)
        offset[j] = length(values)
        a === nothing && continue
        b = findlast(!iszero, column)
        lo[j], hi[j] = a, b
        offset[j] = length(values) - a + 1
        append!(values, view(column, a:b))
    end
    return _Envelope(lo, hi, offset, values)
end
@inline Base.getindex(e::_Envelope, i::Int, j::Int) =
    @inbounds(e.first[j] <= i <= e.last[j]) ? @inbounds(e.values[e.offset[j] + i]) : 0.0

# Two kept blocks mixed: `(1 - w) a + w b`.
struct _Mixed
    a::_Envelope
    b::_Envelope
    w::Float64
end
@inline Base.getindex(m::_Mixed, i::Int, j::Int) = (1 - m.w) * m.a[i, j] + m.w * m.b[i, j]

# The kernels of a tabulated kernel between two stages, built on first use, one for each grid value.
struct _KernelTables{K <: TabulatedKernel, T, S}
    kernel::K
    to::T
    from::S
    rule::Symbol
    tables::Vector{Union{Nothing, _Envelope}}
    lock::ReentrantLock
end
_kernel_source(kernel::TabulatedKernel, to, from, rule) =
    _KernelTables(kernel, to, from, rule, Union{Nothing, _Envelope}[nothing for _ in kernel.grid], ReentrantLock())

function _table(t::_KernelTables, i::Int)
    lock(t.lock) do
        table = t.tables[i]
        if table === nothing
            table = _Envelope(_kernel_block(t.kernel.g(t.kernel.grid[i]), t.to, t.from, t.rule))
            t.tables[i] = table
        end
        return table
    end::_Envelope
end

function _block(t::_KernelTables, ctx)
    kernel, grid = t.kernel, t.kernel.grid
    haskey(ctx.pool, kernel.pool) ||
        throw(ArgumentError("a tabulated kernel reads the pool :$(kernel.pool), which the context does not have"))
    x = Float64(getproperty(ctx.pool, kernel.pool))
    isnan(x) && throw(DomainError(x, "the pool :$(kernel.pool) of a tabulated kernel is not a number"))
    if x < first(grid) || x > last(grid)
        kernel.outside === :exact && return _kernel_block(kernel.g(x), t.to, t.from, t.rule)
        kernel.outside === :error && throw(DomainError(x, "the pool :$(kernel.pool) is outside the grid of a " *
            "tabulated kernel, $(first(grid)) to $(last(grid)); widen the grid, or pass `outside = :clamp` or `:exact`"))
        x = clamp(x, first(grid), last(grid))
    end
    i = searchsortedlast(grid, x)                  # grid[i] <= x < grid[i + 1]
    i == length(grid) && return _table(t, i)
    w = (x - grid[i]) / (grid[i + 1] - grid[i])
    w == 0 && return _table(t, i)
    return _Mixed(_table(t, i), _table(t, i + 1), w)
end

# A step's kernel values, laid out like the projection matrix (`k` stages): the entries of each block
# `(rows, cols, block)` are read from it, and every other entry is one.
struct _KernelValues{B <: Tuple} <: AbstractMatrix{Float64}
    k::Int
    blocks::B
end
Base.size(m::_KernelValues) = (m.k, m.k)
Base.getindex(m::_KernelValues, i::Int, j::Int) = _lookup(m.blocks, i, j)
@inline _lookup(::Tuple{}, i::Int, j::Int) = 1.0
@inline function _lookup(blocks::Tuple, i::Int, j::Int)
    rows, cols, block = blocks[1]
    (i in rows && j in cols) || return _lookup(Base.tail(blocks), i, j)
    return Float64(block[i - first(rows) + 1, j - first(cols) + 1])
end

"""
    by_trait(stage, f; stratum = nothing)
    by_trait(stage, factors; stratum = nothing)

A response value that scales each entry of a transition according to the trait of the
individual it comes from: entries leaving the cell of `stage` at trait value `z` are
multiplied by `f(z)`, and entries leaving other stages are left alone. The factors can
also be given as a vector, one for each cell. Return it from a response:

```julia
respond(life, :seeding => ctx -> by_trait(plant, z -> exp(-0.02 * eggs_on(z, ctx))))
```

The response keeps working if the net is refined afterwards with `stratify` or
`expand_stages`: the cells of every stratum take the factor of the cell they descend from.
A response that depends on the stage's own numbers should read them with `ctx.n[stage]`,
which finds the cells by name, rather than by position in `ctx.n`.

A response attached *after* `stratify` sees the cells of the strata, so it names the
stratum the individuals come from: with `stratum`, the factor applies to entries leaving
the cells of `stage` in that stratum (the names `stratum_stage` gives them, with its
default separator, as for [`trait_weights`](@ref)). The stratum they go to is already
fixed by the transition the response is attached to. A life cycle can therefore be
written for one stratum and its fates split after stratifying:

```julia
life = kernel_net(fan, :growing => (y, x) -> survival(x) * growth(y, x))
net = stratify(life, [:healthy, :infected]; coupling = Dict(:growing => [1 0; 1 0]))
net = respond(net,
    :growing_healthy => ctx -> by_trait(fan, x -> 1 - risk(x); stratum = :healthy),
    :growing_healthy_to_infected => ctx -> by_trait(fan, risk; stratum = :healthy))
```

Entries leaving other stages are left alone, but an entry leaving a cell of `stage` that
is not in `stratum` (because the stratum is misspelt, belongs to another transition, or
the net has no strata) is an error, raised when the response is evaluated.
"""
by_trait(stage::ContinuousStage, f; stratum::Union{Nothing, Symbol} = nothing) =
    by_trait(stage, Float64[f(z) for z in stage.z]; stratum)
function by_trait(stage::ContinuousStage, factors::AbstractVector{<:Real}; stratum::Union{Nothing, Symbol} = nothing)
    length(factors) == length(stage) ||
        throw(DimensionMismatch("by_trait needs one factor for each of the $(length(stage)) cells of :$(stage.name)"))
    if stratum === nothing
        index = stage.index
        return function (from, to)
            i = get(index, from, 0)
            return i == 0 ? 1.0 : factors[i]
        end
    end
    cells = Dict(stratum_stage(c, stratum) => i for (i, c) in enumerate(stage.cells))
    return function (from, to)
        i = get(cells, from, 0)
        i == 0 || return factors[i]
        _check_other_stratum(stage, stratum, from)
        return 1.0
    end
end

# An entry leaves `from`, which is not a cell of `stage` in `stratum`: an error if `from` is a cell of `stage`
# all the same, unstratified or in another stratum, since the factor was then meant for it.
function _check_other_stratum(stage::ContinuousStage, stratum::Symbol, from)
    from isa Symbol || return nothing
    name = String(from)
    cell = haskey(stage.index, from) ? from : nothing
    if cell === nothing
        for (i, ch) in pairs(name)
            ch == '_' || continue
            prefix = Symbol(SubString(name, 1, prevind(name, i)))
            haskey(stage.index, prefix) && (cell = prefix; break)
        end
    end
    cell === nothing && return nothing
    throw(ArgumentError("by_trait(:$(stage.name); stratum = :$stratum) met an entry leaving :$from, a cell of " *
        ":$(stage.name) that is not in stratum :$stratum (expected :$(stratum_stage(cell, stratum))); attach it " *
        "to a transition that leaves stratum :$stratum, after `stratify`, or check the stratum's name"))
end

"""
    trait_weights(stage, f; stratum = nothing)

Pool weights given as a function of the trait: each individual of `stage` with trait `z`
contributes `f(z)`. For the `occupies` argument of [`species`](@ref):

```julia
species(:tree, net; occupies = [:basal_area => trait_weights(tree, d -> pi * (d / 2)^2)])
```

In a net made with `stratify` the weights apply to the individuals of every stratum. With
`stratum` they apply to that stratum alone, so that a pool can hold, say, the infected
individuals only:

```julia
species(:host, net; occupies = [:infectious => trait_weights(plant, z -> z^2; stratum = :infected)])
```
"""
trait_weights(stage::ContinuousStage, f; stratum::Union{Nothing, Symbol} = nothing) =
    Pair{Symbol, Float64}[(stratum === nothing ? c : stratum_stage(c, stratum)) => Float64(f(z)) for (c, z) in zip(stage.cells, stage.z)]
