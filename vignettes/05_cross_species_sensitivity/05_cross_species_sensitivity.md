# Revealing interactions: cross-species sensitivity analysis
Simon Frost

- [Interactions as an output](#interactions-as-an-output)
- [One perturbation](#one-perturbation)
- [The full matrix](#the-full-matrix)
- [Which stages matter most](#which-stages-matter-most)
- [The network](#the-network)
- [A note on R²](#a-note-on-r²)
- [Reference](#reference)

``` julia
ENV["GKSwstype"] = "100"
using MatrixCommunityModels
using Plots, Random, Statistics
```

## Interactions as an output

A matrix community model contains no interaction coefficients. Lytle &
Tonkin (2023) recover interactions afterwards: change one vital rate of
one species and stage a little, project the community, and see how the
average abundance of every species and stage responds. Doing this for
every vital rate gives a matrix of net effects, which they summarise as
a network (their Figs. 2 and 3).

This vignette does that for the riparian guild model of the previous
vignette, in its corrected form, perturbing the flood mortality of each
of the 30 guild-stages.

## One perturbation

`abundance_sensitivity` takes a function that builds the model with the
vital rate set to a given value, the values to try, and the
environments. `with_rate` makes a copy of the guild parameters with one
entry changed.

The flow regime is a random draw of years from the record. The same 20
sequences of 300 years are used for every value of the vital rate, so
that differences between runs are due to the vital rate and not to the
luck of the draw. The first 100 years are discarded.

``` julia
base = riparian_model(corrected = true)
rng = Xoshiro(1)
sequences = [[rand(rng, base.envs) for _ in 1:300] for _ in 1:20]

function flood_mortality_sensitivity(guild, stage; envs = r -> sequences[r])
    v0 = RIPARIAN_GUILDS[guild].S[stage]
    values = collect(range(max(0, v0 - 0.05), min(1, v0 + 0.05); step = 0.01))
    build(v) = riparian_model(corrected = true, guilds = with_rate(RIPARIAN_GUILDS, guild, :S, stage, v))
    return abundance_sensitivity(build, values, envs; burnin = 100, replicates = 20), values
end

s, values = flood_mortality_sensitivity(:C, 6)
```

    ((slope = [2.9333329310629065e6, 27630.258619608467, 2139.308508478614, 1671.9234276629882, 4575.162408312951, -16076.338123579006, 1.9980383659485248e6, 53664.25121315984, 14910.945232335136, 6952.269078045299  …  54.56281161324272, 29.519361490532145, 64.73556410076138, 26.68511106397088, 4.546550914646113e6, 143378.03340640862, 2397.3729291399177, 1544.7038042519557, 6738.742026839686, 26108.03894632797], r2 = [0.9723450396997108, 0.970813570996602, 0.9707425978706129, 0.9706395763153273, 0.9706153553918346, 0.9708809728456742, 0.9719809007010939, 0.9708049562950889, 0.9707540011996333, 0.9706585971566825  …  0.6085759407957388, 0.7220297195540397, 0.7870100542927122, 0.8565082104721029, 0.9724336561407138, 0.9706565819473285, 0.9704534039646953, 0.9705647855435483, 0.970405273097389, 0.9694969415992729], means = [2.0739667443201113e6 19346.12579654757 … 4761.161764321525 17967.024797352424; 2.121960545372337e6 19809.105111498346 … 4874.3477177538725 18408.60575491752; … ; 2.2656780072167516e6 21158.246869760467 … 5203.327750668759 19682.464429442836; 2.2834385140301543e6 21323.93705126313 … 5243.659070489611 19837.397321947938]), [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])

The response of mature (stage 6) individuals of each guild to the flood
mortality of mature cottonwood:

``` julia
guilds = [:C, :T, :W, :S, :M]
mature = 6:6:30
plot(values, s.means[:, mature] ./ s.means[1:1, mature]; label = permutedims(String.(guilds)),
     xlabel = "Flood mortality of stage-6 cottonwood", ylabel = "Mean stage-6 abundance (relative)",
     marker = :circle, lw = 2)
```

![](05_cross_species_sensitivity_files/figure-commonmark/cell-4-output-1.svg)

Killing more mature cottonwood in floods reduces cottonwood and lets
every other guild increase: competition, although the model has no
competition coefficient.

``` julia
(slope = round.(s.slope[mature]; digits = 1), r2 = round.(s.r2[mature]; digits = 3))
```

    (slope = [-16076.3, 10136.6, 16485.4, 26.7, 26108.0], r2 = [0.971, 0.969, 0.971, 0.857, 0.969])

## The full matrix

``` julia
labels = [string(g, i) for g in ("HT", "XS", "HS", "DS", "MM") for i in 1:6]
slope = zeros(30, 30)
r2 = zeros(30, 30)
for (a, g) in enumerate(guilds), i in 1:6
    res, _ = flood_mortality_sensitivity(g, i)
    slope[:, (a - 1) * 6 + i] = res.slope
    r2[:, (a - 1) * 6 + i] = res.r2
end
```

Column `j` holds the effect of the flood mortality of guild-stage `j` on
the mean abundance of every guild-stage (rows). To compare stages whose
abundances differ by orders of magnitude, each row is divided by that
stage’s mean abundance in the unperturbed model, giving the proportional
change in abundance per unit of mortality. (Guild labels follow the
paper: HT cottonwood, XS tamarisk, HS willow, DS sagebrush, MM meadow.)

``` julia
baseline = mean(mean_abundance(base.community, base.n0, q; burnin = 100) for q in sequences)
effect = slope ./ baseline
heatmap(labels, labels, clamp.(effect, -3, 3); yflip = true, color = :balance, clims = (-3, 3),
        xrotation = 90, xticks = (0.5:1:29.5, labels), yticks = (0.5:1:29.5, labels), tickfontsize = 5,
        xlabel = "Perturbed flood mortality", ylabel = "Responding abundance", size = (650, 600),
        colorbar_title = "proportional change per unit mortality")
```

![](05_cross_species_sensitivity_files/figure-commonmark/cell-7-output-1.svg)

Red cells are positive effects (the responding stage gains when the
perturbed stage suffers more in floods) and blue cells negative. The
diagonal blocks are within-guild effects. Two features stand out:

``` julia
(positive = count(effect .> 0), negative = count(effect .< 0),
 negative_within_guild = count(effect[i, j] < 0 for i in 1:30, j in 1:30 if (i - 1) ÷ 6 == (j - 1) ÷ 6))
```

    (positive = 706, negative = 194, negative_within_guild = 106)

Most effects between guilds are positive, that is, competitive. The
exceptions involve sagebrush (DS). It declines when willow (HS) suffers
more in floods, a net positive effect of willow on sagebrush that no one
put into the model. And more flood mortality of sagebrush lowers
cottonwood, tamarisk and meadow in about half of their stages. In the
corrected model sagebrush is not held at 10 individuals per stage, as in
the authors’ script, and its numbers swing over orders of magnitude. The
effect of a stage’s mortality on its own abundance is not always
negative:

``` julia
[labels[k] => round(effect[k, k]; digits = 2) for k in (1, 6, 13, 18)]
```

    4-element Vector{Pair{String, Float64}}:
     "HT1" => 0.43
     "HT6" => -7.77
     "HS1" => 0.4
     "HS6" => -3.52

More flood mortality of seedlings (stage 1) *increases* the average
number of seedlings of that guild, because the space freed is
recolonised by new seedlings and recruitment in this model does not
depend on the number of adults.

## Which stages matter most

Counting, for each perturbed stage, how many guild-stages change by more
than half their abundance per unit of mortality gives a simple measure
of influence:

``` julia
influence = vec(sum(abs.(effect) .> 0.5; dims = 1))
bar(labels, influence; xrotation = 90, legend = false, xticks = (0.5:1:29.5, labels), tickfontsize = 6,
    ylabel = "Guild-stages affected (|effect| > 0.5)")
```

![](05_cross_species_sensitivity_files/figure-commonmark/cell-10-output-1.svg)

The older stages of cottonwood (HT5, HT6) and meadow (MM5, MM6), and
stage-5 tamarisk, affect almost the whole community. Lytle & Tonkin
single out mature cottonwood (HT6) as a keystone of this system. Their
analysis used the original model, so the two are not directly
comparable, but the same stage emerges here.

## The network

The matrix can be drawn as a directed network: an arrow from the
perturbed stage to each stage it strongly affects.

``` julia
θ = range(0, 2π; length = 31)[1:30]
x, y = cos.(θ), sin.(θ)
p = plot(; aspect_ratio = 1, axis = false, grid = false, legend = false, size = (600, 600),
         xlims = (-1.3, 1.3), ylims = (-1.3, 1.3))
for j in 1:30, i in 1:30
    (i != j && abs(effect[i, j]) > 1) || continue
    plot!(p, [x[j], x[i]], [y[j], y[i]]; color = effect[i, j] > 0 ? :firebrick : :steelblue,
          alpha = 0.35, lw = min(abs(effect[i, j]), 6) / 2)
end
scatter!(p, x, y; ms = 3 .+ 1.5 .* log10.(baseline .+ 1), color = :white, markerstrokecolor = :black)
annotate!(p, [(1.12 * x[k], 1.12 * y[k], text(labels[k], 7)) for k in 1:30])
p
```

![](05_cross_species_sensitivity_files/figure-commonmark/cell-11-output-1.svg)

Red links are positive effects, blue negative; node size reflects mean
abundance.

`InteractionNetwork` holds the same network and gives the usual summary
measures:

``` julia
net = InteractionNetwork(effect, labels; threshold = 1)
(links = count(net.adjacency), connectance = round(connectance(net); digits = 3),
 link_density = round(link_density(net); digits = 2), reciprocity = round(reciprocity(net); digits = 2))
```

    (links = 124, connectance = 0.143, link_density = 4.13, reciprocity = 0.27)

``` julia
first(links(net), 5)
```

    5-element Vector{@NamedTuple{from::String, to::String, effect::Float64}}:
     (from = "DS1", to = "DS3", effect = -39.83118912586431)
     (from = "DS1", to = "DS4", effect = -37.199751398776144)
     (from = "DS1", to = "DS2", effect = -35.7440324281075)
     (from = "DS1", to = "DS5", effect = -32.16764907481866)
     (from = "DS1", to = "DS6", effect = -27.384900769960044)

``` julia
sort(collect(zip(labels, out_degree(net))); by = last, rev = true)[1:5]
```

    5-element Vector{Tuple{String, Int64}}:
     ("MM6", 29)
     ("HT6", 24)
     ("HS5", 7)
     ("HS3", 6)
     ("HT1", 5)

## A note on R²

The paper keeps only relationships with R² above 0.3. Here most
relationships have an R² close to one; the lowest are those of
sagebrush, whose numbers swing widely:

``` julia
quantile(vec(r2), [0.05, 0.5, 0.95])
```

    3-element Vector{Float64}:
     0.4386682611841727
     0.9975663990581581
     0.9999867120876612

That is a consequence of using the same flow sequences for every value
of the vital rate. If each run draws its own sequence, as in the
original analysis, sampling noise is added to every point and R² falls,
most for the weakest effects:

``` julia
rng2 = Xoshiro(2)
noisy, _ = flood_mortality_sensitivity(:C, 6; envs = r -> [rand(rng2, base.envs) for _ in 1:300])
(common_sequences = round.(s.r2[mature]; digits = 2), independent_sequences = round.(noisy.r2[mature]; digits = 2),
 slope_common = round.(s.slope[mature]; digits = 1), slope_independent = round.(noisy.slope[mature]; digits = 1))
```

    (common_sequences = [0.97, 0.97, 0.97, 0.86, 0.97], independent_sequences = [0.96, 0.57, 0.55, 0.33, 0.56], slope_common = [-16076.3, 10136.6, 16485.4, 26.7, 26108.0], slope_independent = [-16566.7, 8766.6, 26037.1, -509.7, 34222.9])

With common sequences the R² filter is not needed to separate signal
from noise, and far fewer replicates are required; the size of the
effect is the more useful criterion.

## Reference

Lytle, D. A. & Tonkin, J. D. (2023). Matrix community models for ecology
and evolution. *npj Biodiversity* 2, 26.

Tonkin, J. D., Merritt, D. M., Olden, J. D., Reynolds, L. V. & Lytle, D.
A. (2018). Flow regime alteration degrades ecological networks in
riparian ecosystems. *Nature Ecology & Evolution* 2, 86-93.
