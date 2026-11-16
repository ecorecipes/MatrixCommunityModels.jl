# A predator and its prey, each with juveniles and adults: the nonlinear matrix model behind the integrated
# community model of Barraquand & Gimenez (2019). Census before breeding, one step a year. Adult predators lower
# the survival of newborn prey (and of their own newborns), juvenile prey raise the fecundity of predators, and
# adult prey lower their own fecundity. The species meet through two pools:
#
#   adult_predators   adult predators, read by the prey
#   juvenile_prey     juvenile prey, read by the predator
#
# Written from the equations and Table 5 of the paper; the authors' code (GPL-3) is not used here.

_logit_survival(a, b, n) = 1 / (1 + exp(-(a + b * n)))

"""
    BARRAQUAND_GIMENEZ_2019

The first parameter set of Barraquand & Gimenez (2019, *Ecol. Model.* 411:108713, Table 5): `alpha` holds
α1 to α8, where newborn predator survival is `logistic(α1 + α2 × adult predators)`, newborn prey survival
`logistic(α3 + α4 × adult predators)`, predator fecundity `exp(α5 + α6 × juvenile prey)` and prey fecundity
`exp(α7 + α8 × adult prey)`. `adult_survival` is the constant adult survival of prey and predator, and `n0` the
initial juvenile and adult prey and juvenile and adult predators.
"""
const BARRAQUAND_GIMENEZ_2019 = (
    alpha = (0.5, -0.01, 0.5, -0.025, 0.0, 0.004, 2.0, -0.005),
    adult_survival = (prey = 0.6, predator = 0.7),
    n0 = (prey = [100.0, 100.0], predator = [20.0, 20.0]))

"""
    BARRAQUAND_GIMENEZ_2019_SET2

The second parameter set of Barraquand & Gimenez (2019, Appendices A and B): as
[`BARRAQUAND_GIMENEZ_2019`](@ref) but `α5 = 0.5`, `α6 = 0.01` and `α7 = 1.5`, a stronger coupling of predator
and prey.
"""
const BARRAQUAND_GIMENEZ_2019_SET2 = merge(BARRAQUAND_GIMENEZ_2019,
    (alpha = (0.5, -0.01, 0.5, -0.025, 0.5, 0.01, 1.5, -0.005),))

"""
    predator_prey(parameters = BARRAQUAND_GIMENEZ_2019)

The two-species, two-stage predator-prey model of Barraquand & Gimenez (2019): a pre-breeding census, so that an
adult's newborns are counted as juveniles a year later only if they survive, and juveniles become adults in one
year. Each species is a net with stages `juvenile` and `adult` and transitions `reproduction` (adults to
juveniles, half the fecundity times newborn survival, marked as births), `maturation` and `adult_survival` (both
at the constant adult survival). The prey occupies the pool `juvenile_prey` and reads `adult_predators`; the
predator occupies `adult_predators` and reads `juvenile_prey`.

Returns `(community, n0)`. Built with births marked and survival probabilities, the community can be projected
with demographic stochasticity (`project(..., rng = ...)`), which gives the binomial survival and Poisson
recruitment of the paper's stochastic model.
"""
function predator_prey(parameters = BARRAQUAND_GIMENEZ_2019)
    α = parameters.alpha
    φ = parameters.adult_survival
    stage(φA) = ValuedProjectionNet([:juvenile, :adult],
        :reproduction => [(:adult => :juvenile) => 0.5],
        :maturation => [(:juvenile => :adult) => φA],
        :adult_survival => [(:adult => :adult) => φA])
    prey = respond(stage(φ.prey),
        :reproduction => ctx -> exp(α[7] + α[8] * ctx.n[:adult]),                        # crowded adults breed less
        :reproduction => ctx -> _logit_survival(α[3], α[4], ctx.pool.adult_predators))    # newborns are eaten
    predator = respond(stage(φ.predator),
        :reproduction => ctx -> exp(α[5] + α[6] * ctx.pool.juvenile_prey),                # food for the young
        :reproduction => ctx -> _logit_survival(α[1], α[2], ctx.n[:adult]))               # crowded newborns die
    com = community(
        species(:prey, births(prey, :reproduction); occupies = [:juvenile_prey => [:juvenile => 1.0]],
                reads = [:adult_predators]),
        species(:predator, births(predator, :reproduction); occupies = [:adult_predators => [:adult => 1.0]],
                reads = [:juvenile_prey]))
    n0 = [:prey => parameters.n0.prey, :predator => parameters.n0.predator]
    return (community = com, n0 = n0)
end
