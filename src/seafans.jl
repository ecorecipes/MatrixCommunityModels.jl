# Sea fans and aspergillosis: an integral projection model stratified by infection status.
#
# The model is that of Bruno, Ellner, Vu, Kim & Harvell (2011) for the Caribbean sea fan Gorgonia ventalina
# and the fungus Aspergillus sydowii. A colony is described by its size and by whether it is infected. The
# life cycle of a healthy fan is written as one continuous stage, with the kernels of infected fans beside it;
# `stratify` makes the strata, coupling matrices say which stratum each transition leaves and which it enters,
# and the chance of each fate is then attached, by size, to the transitions that leave one stratum.
#
# Size is the cube root of the area of healthy tissue in cm², the time step is a year and numbers are colonies
# per m² of suitable habitat. Equation, table, figure, page and appendix numbers below are the paper's.

"""
    SEA_FAN_PARAMETERS

Default parameters of [`sea_fan`](@ref): the integral projection model of the sea fan
*Gorgonia ventalina* of Bruno, Ellner, Vu, Kim & Harvell (2011), fitted to about 750
colonies followed from 2002 to 2005 on five reefs, two in the Florida Keys and three off
the Yucatán. Size is the cube root of the area of healthy
tissue in cm², and the time step is a year. The rates that changed over the epizootic are
in [`SEA_FAN_PHASES`](@ref). Sources of the values are in the comments of `src/seafans.jl`.
"""
const SEA_FAN_PARAMETERS = (
    size_range = (0.0, 20.0),              # the largest fans measured had a size of about 14 (Figs. 2 and 3); the range holds the population
    height = (1.3, 1.5),                   # height in cm = 1.3 x size^1.5, that is 1.3 x sqrt(area) (Eq. A.3)
    growth = (0.682, 0.945, 0.53),         # healthy fan: size next year is normal, mean a + b x size, with this sd (Eq. 6)
    recovering_growth = (0.27, 0.57),      # infected fan that recovers: mean change in size, and its sd (Eq. 7)
    infected_growth = (-0.33, 0.57),       # infected fan that stays infected: mean change in size, and its sd (Eq. 7)
    tissue_loss = (0.05, -3.2, 0.52),      # share of infected fans that lose much tissue, their mean change in size and its sd (Eq. 7)
    mortality = (-0.55, -0.62),            # logit of the mortality of a healthy fan, a + b x size (p. 130)
    infected_mortality = (1.3, -0.67),     # logit of the mortality of an infected fan (p. 130)
    storm_mortality = (0.093, 2.3, 8.7),   # mortality of large fans in storms, healthy or infected: a x logistic(b x (size - c)) (Eq. B.2)
    large_fan = 5.6,                       # fans above this size run the higher risk of infection (Eq. 5): 175 cm², about 17.5 cm tall
    risk_height = 20.0,                    # cm; at the peak the risk of a large fan rose with its height above this (Appendix D)
    reproductive_height = 20.0,            # cm; shorter fans do not reproduce (p. 130)
    # Recruits for each cm² of healthy tissue on reproductive fans (Eq. 4). The paper prints "2 x 10^4", which cannot be
    # meant: reproductive fans carry some hundreds of cm² of tissue per m², and 0.19 recruits per m² were counted. Read as
    # 2 x 10^-4, a value given to one figure. Calibrated within that rounding: 2.0 x 10^-4 gives growth rates 0.005 to
    # 0.010 below the four that can be read from the paper's Figs. E2 and E3, and 2.4 x 10^-4 gives all four within 0.004.
    fecundity = 2.4e-4,
    # Size of a recruit: mean and sd. NOT IN THE PAPER, which says only that the mortality of recruits predicted from "the
    # estimated recruit size distribution" is 22% (pp. 127 to 128). Assumed: a normal distribution with this sd. Its mean
    # is then set so that 22% of recruits die in their first year. Such a recruit is about 1.6 cm tall.
    recruit_size = (1.17, 0.3))

"""
    SEA_FAN_PHASES

The three phases of the aspergillosis epizootic in [`sea_fan`](@ref), after Bruno et al.
(2011): `pre_epizootic` (before 1994), `peak` (the outbreak, fitted to surveys of 1997 to
1999 in the Florida Keys) and `endemic` (after the peak, from colonies followed in the
Florida Keys and the Yucatán from 2002 to 2005). A phase is the environment of a time step, so an epizootic is
a sequence of phases:

```julia
reef = sea_fan()
project(reef, n0, [fill(SEA_FAN_PHASES.peak, 6); fill(SEA_FAN_PHASES.endemic, 8)])
```

Each phase holds

- `small_fan_risk`, `large_fan_risk`, `risk_per_cm`: the chance that a healthy fan is
  infected within the year. It is `small_fan_risk` up to the size `large_fan` of
  [`SEA_FAN_PARAMETERS`](@ref), and above it `large_fan_risk` plus `risk_per_cm` for each
  cm of height above `risk_height`, which is 20 cm;
- `recovery`: the chance that an infected fan that survives the year is healthy at its end;
- `disease_mortality`: mortality of infected fans of any size, added to the mortality that
  depends on size;
- `recruitment`: a multiplier of the recruits produced, which were far fewer before 2003.

A phase is a `NamedTuple`, so a variant is made with `merge`, for example
`merge(SEA_FAN_PHASES.endemic, (recovery = 0.2,))`.
"""
const SEA_FAN_PHASES = (
    # no infection; recruitment was about 11 times lower than in 2003 (Appendix C)
    pre_epizootic = (small_fan_risk = 0.0, large_fan_risk = 0.0, risk_per_cm = 0.0, recovery = 0.38,
                     disease_mortality = 0.0, recruitment = 1 / 11),
    # Appendix D: ln(risk at 20 cm) = -3.0, slope 0.0083 per cm, ln(disease mortality) = -0.48, recovery 0.10. The risk of
    # small fans was fitted too but is not tabulated; Fig. 8B shows it near 0.01, so the endemic value is used. Recruitment
    # is reduced as before the epizootic, as in the paper's Table E1 and Fig. E2B
    peak = (small_fan_risk = 0.011, large_fan_risk = exp(-3.0), risk_per_cm = 0.0083, recovery = 0.10,
            disease_mortality = exp(-0.48), recruitment = 1 / 11),
    # Eq. 5 and p. 130: risks 0.011 and 0.053, recovery 0.38
    endemic = (small_fan_risk = 0.011, large_fan_risk = 0.053, risk_per_cm = 0.0, recovery = 0.38,
               disease_mortality = 0.0, recruitment = 1.0))

"""
    sea_fan_stage(; parameters = SEA_FAN_PARAMETERS, cells = 100)

The [`ContinuousStage`](@ref) from which [`sea_fan`](@ref) is built: colonies by size,
the cube root of the area of healthy tissue in cm². `meshpoints` of it are the sizes of
the cells. In the stratified net its cell `fan7` becomes `fan7_healthy` and
`fan7_infected`; `ctx.n[stage]` gives the numbers by size whatever the stratum, and
`ctx.n[stage, :infected]` those of the infected fans.
"""
sea_fan_stage(; parameters = SEA_FAN_PARAMETERS, cells::Integer = 100) =
    ContinuousStage(:fan, parameters.size_range[1], parameters.size_range[2], cells)

"""
    sea_fan(; parameters = SEA_FAN_PARAMETERS, cells = 100, resistant = false,
              outside = 1.0, transmission = 0.0, name = :sea_fan)

The sea fan *Gorgonia ventalina* with aspergillosis, as the integral projection model of
Bruno et al. (2011): colonies are classified by size (see [`sea_fan_stage`](@ref)) and by
infection status. Numbers are colonies per m² of suitable habitat and the time step is a
year.

The life cycle is that of a healthy fan, which survives and grows (`growing`) and, if 20
cm or taller, produces recruits (`recruiting`), with the kernels of an infected fan that
survives beside it: one that recovers (`recovering`) and one that stays infected
(`lingering`). `stratify` makes the strata `healthy` and `infected`: a healthy fan that
grows stays healthy or is infected, and `recovering` leads back, so the stages are
`fan1_healthy`, ..., `fan1_infected`, ... and the transitions `growing_healthy`,
`growing_healthy_to_infected`, `recovering_infected_to_healthy`, `lingering_infected` and
`recruiting_healthy`. The chance of each fate at each size is then attached to the
transitions that leave one stratum, with `by_trait(...; stratum = ...)`.

The environment of a time step is the phase of the epizootic, one of
[`SEA_FAN_PHASES`](@ref) or a variant of one; with `nothing` the phase is `endemic`.

- `resistant = true` adds a third stratum, `resistant`, of fans that cannot be infected
  and whose recruits are resistant too: they have the life cycle of a healthy fan
  (`growing_resistant`, `recruiting_resistant`) and are never infected. `resistant = :only`
  gives resistant fans alone, to be composed with the others as a species of their own.
- The species adds its healthy tissue, in cm², to the pool `tissue`, and that of its
  infected fans to the pool `infectious`.
- `outside` and `transmission` say where the infection comes from. The chance of infection
  given by the phase is multiplied by `outside + transmission * infectious / tissue`, read
  from the pools. With the defaults all infection comes from outside the reef and the model
  is the published one, which is linear. `transmission > 0` makes a share of it depend on
  the share of tissue that is infected; the paper has no such term and gives no value.
"""
function sea_fan(; parameters = SEA_FAN_PARAMETERS, cells::Integer = 100, resistant::Union{Bool, Symbol} = false,
                 outside::Real = 1.0, transmission::Real = 0.0, name::Symbol = :sea_fan)
    resistant isa Bool || resistant === :only ||
        throw(ArgumentError("resistant must be false, true or :only, got :$resistant"))
    (outside >= 0 && transmission >= 0) || throw(ArgumentError("outside and transmission cannot be negative"))
    p = parameters
    fan = sea_fan_stage(; parameters, cells)
    height(x) = p.height[1] * x^p.height[2]
    storms(x) = p.storm_mortality[1] * _logistic(p.storm_mortality[2] * (x - p.storm_mortality[3]))
    survival(x) = max(0.0, 1 - _logistic(p.mortality[1] + p.mortality[2] * x) - storms(x))
    infected_survival(x, phase) = max(0.0, 1 - _logistic(p.infected_mortality[1] + p.infected_mortality[2] * x) -
                                           storms(x) - phase.disease_mortality)
    growth(y, x) = _normal_density(y, p.growth[1] + p.growth[2] * x, p.growth[3])
    lost, loss, loss_sd = p.tissue_loss
    # size next year of an infected fan whose size changes by `change` on average, or which loses much tissue
    changing((change, sd)) = (y, x) -> (1 - lost) * _normal_density(y, x + change, sd) + lost * _normal_density(y, x + loss, loss_sd)
    recruits(x) = height(x) >= p.reproductive_height ? p.fecundity * x^3 : 0.0
    recruit_size(y) = _normal_density(y, p.recruit_size[1], p.recruit_size[2])

    phase(ctx) = ctx.env === nothing ? SEA_FAN_PHASES.endemic : ctx.env
    # the life cycle of a healthy fan, a model on its own, and beside it where an infected fan of size x that survives
    # goes if it recovers or if it stays infected
    life = kernel_net(fan,
        :growing => (y, x) -> survival(x) * growth(y, x),
        :recruiting => (y, x) -> recruits(x) * recruit_size(y),
        :recovering => changing(p.recovering_growth),
        :lingering => changing(p.infected_growth))
    life = births(respond(life, :recruiting => ctx -> phase(ctx).recruitment), :recruiting)

    strata = resistant === :only ? [:resistant] : resistant ? [:healthy, :infected, :resistant] : [:healthy, :infected]
    # which stratum each transition leaves and which it enters; coupling matrices are [to, from]. A healthy fan grows
    # the same whether or not it is infected by the next census (p. 128); a resistant fan is never infected
    route(moves...) = Float64[(from => to) in moves for to in strata, from in strata]
    net = stratify(life, strata; coupling = Dict{Symbol, Any}(
        :growing => route(:healthy => :healthy, :healthy => :infected, :resistant => :resistant),
        :recruiting => route(:healthy => :healthy, :resistant => :resistant),
        :recovering => route(:infected => :healthy),
        :lingering => route(:infected => :infected)))

    function risk(x, ph)                                   # chance that a healthy fan of size x is infected within the year
        x <= p.large_fan && return ph.small_fan_risk
        return clamp(ph.large_fan_risk + ph.risk_per_cm * max(0.0, height(x) - p.risk_height), 0.0, 1.0)
    end
    # the inoculum that a healthy fan meets, relative to that of the published model
    function pressure(ctx)
        (transmission > 0 && haskey(ctx.pool, :infectious) && ctx.pool.tissue > 0) || return Float64(outside)
        return outside + transmission * ctx.pool.infectious / ctx.pool.tissue
    end
    infection(x, ctx) = min(1.0, pressure(ctx) * risk(x, phase(ctx)))
    # the chance of each fate at each size, attached to the transitions that leave one stratum
    if :infected in strata
        net = respond(net,
            :growing_healthy => ctx -> by_trait(fan, x -> 1 - infection(x, ctx); stratum = :healthy),
            :growing_healthy_to_infected => ctx -> by_trait(fan, x -> infection(x, ctx); stratum = :healthy),
            :recovering_infected_to_healthy =>
                ctx -> by_trait(fan, x -> infected_survival(x, phase(ctx)) * phase(ctx).recovery; stratum = :infected),
            :lingering_infected =>
                ctx -> by_trait(fan, x -> infected_survival(x, phase(ctx)) * (1 - phase(ctx).recovery); stratum = :infected))
    end

    tissue(x) = x^3                                        # cm² of healthy tissue on a fan of size x
    occupies = Pair{Symbol, Any}[:tissue => trait_weights(fan, tissue)]                  # fans of every stratum
    :infected in strata && push!(occupies, :infectious => trait_weights(fan, tissue; stratum = :infected))
    reads = transmission > 0 && :infected in strata ? [:tissue, :infectious] : Symbol[]
    return species(name, net; occupies, reads)
end
