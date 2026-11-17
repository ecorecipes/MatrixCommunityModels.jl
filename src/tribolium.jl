# Two flour beetles that eat each other's eggs and pupae: the competition LPA model of Edmunds et al. (2003) and
# Desharnais et al. (2005). Each species is the larva-pupa-adult (LPA) net of Dennis et al. (1995), with a time
# step of two weeks. The species meet only through four pools, one per victim stage and species:
#
#   <name>_eggs    the cannibalism hazard on the eggs of species <name>: larvae and adults of both species,
#                  each weighted by its per-capita rate of eating those eggs
#   <name>_pupae   the hazard on its pupae, from the adults of both species
#
# The fraction of eggs that survive to be counted as larvae is exp(-<name>_eggs), and of pupae that emerge as
# adults exp(-<name>_pupae).

"""
    TRIBOLIUM_EDMUNDS_2003

Parameters of the competition LPA model of Edmunds et al. (2003, *J. Anim. Ecol.* 72:703-712), species 1
*Tribolium confusum* (`:confusum`) and species 2 *T. castaneum* (`:castaneum`), in an 8 g habitat (`V = 0.4`,
in units of 20 g of flour). Coefficient matrices are indexed `[victim, eater]`: `c_el[i, j]` is the rate at which
larvae of species `j` eat eggs of species `i`, `c_ea` the same for adults, and `c_pa[i, j]` the rate at which
adults of `j` eat pupae of `i`. The model divides each coefficient by `V`.

The paper varies the interspecific egg-eating by larvae, `c_el[1, 2] = c_el[2, 1]` (its `ceL = cEl`), from
0.016 to 0.032; the value here, 0.01, is the one it lists before doing so. Vignette 21 explores the paper's
range.
"""
const TRIBOLIUM_EDMUNDS_2003 = (
    species = (:confusum, :castaneum),
    V = 0.4,
    b = (10.0, 10.0),            # larval recruits per adult per two weeks without cannibalism
    mu_l = (0.2, 0.2),           # larval mortality per two weeks
    mu_a = (0.02, 0.02),         # adult mortality per two weeks
    c_el = [0.01 0.01; 0.01 0.01],
    c_ea = [0.01 0.01; 0.01 0.01],
    c_pa = [0.01 0.008; 0.005 0.012])   # cpa, cpA; cPa, cPA

"""
    TRIBOLIUM_DESHARNAIS_2005

Parameters of the competition LPA model of Desharnais et al. (2005, *J. Difference Equ. Appl.* 11:311-325,
Table 1), species 1 *Tribolium castaneum* and species 2 *T. confusum*. The coefficients there already include
the habitat size, so `V = 1`. Indexing as in [`TRIBOLIUM_EDMUNDS_2003`](@ref): `c_el[i, j]` is the paper's
`c_el,ij`, eating of species `i`'s eggs by larvae of species `j`. The paper varies `c_el[2, 1]` from 0 to 0.10.
"""
const TRIBOLIUM_DESHARNAIS_2005 = (
    species = (:castaneum, :confusum),
    V = 1.0,
    b = (10.0, 10.0),
    mu_l = (0.2, 0.2),
    mu_a = (0.14, 0.14),
    c_el = [0.012 0.01; 0.062 0.003],
    c_ea = [0.006 0.011; 0.007 0.015],
    c_pa = [0.064 0.017; 0.03 0.01])

"""
    flour_beetle(name; b, mu_l, mu_a, eats = [])

One flour beetle species as the LPA net of Dennis et al. (1995): stages `L` (feeding larvae), `P` (non-feeding
larvae, pupae and callow adults) and `A` (adults), and a time step of two weeks. Transitions are `recruitment`
(adults to larvae, `b` per adult, marked as births), `larval_survival` (`1 - mu_l`), `emergence` (pupae to adults)
and `adult_survival` (`1 - mu_a`).

Recruitment is multiplied by `exp(-pool)` for the pool `<name>_eggs`, and emergence by `exp(-pool)` for the pool
`<name>_pupae`: the shares of eggs and pupae that escape being eaten. The beetle adds to these pools, of its own
species and of others, what it eats: `eats` lists `victim => (larvae_eggs = ..., adults_eggs = ..., adults_pupae
= ...)`, per-capita rates of eating the victim's eggs (by this species' larvae and adults) and pupae (by its
adults). Include `name` itself for cannibalism. Rates are used as given (divide by the habitat size first).
"""
function flour_beetle(name::Symbol; b::Real, mu_l::Real, mu_a::Real, eats = Pair{Symbol, Any}[])
    life = ValuedProjectionNet([:L, :P, :A],
        :recruitment => [(:A => :L) => Float64(b)],
        :larval_survival => [(:L => :P) => 1 - mu_l],
        :emergence => [(:P => :A) => 1.0],
        :adult_survival => [(:A => :A) => 1 - mu_a])
    eggs, pupae = Symbol(name, :_eggs), Symbol(name, :_pupae)
    net = respond(life,
        :recruitment => ctx -> exp(-getproperty(ctx.pool, eggs)),
        :emergence => ctx -> exp(-getproperty(ctx.pool, pupae)))
    net = births(net, :recruitment)
    occupies = Pair{Symbol, Any}[]
    for (victim, c) in eats
        push!(occupies, Symbol(victim, :_eggs) => [:L => Float64(c.larvae_eggs), :A => Float64(c.adults_eggs)])
        push!(occupies, Symbol(victim, :_pupae) => [:A => Float64(c.adults_pupae)])
    end
    own = (eggs, pupae)
    return species(name, net; occupies, reads = [p for p in own if !any(o -> o.first === p, occupies)])
end

"""
    flour_beetle_competition(parameters = TRIBOLIUM_EDMUNDS_2003)

The two-species competition LPA model (Edmunds et al. 2003; Desharnais et al. 2005): two
[`flour_beetle`](@ref)s that eat each other's eggs and pupae as well as their own, composed through the four
pools `<species>_eggs` and `<species>_pupae`. `parameters` is a set such as [`TRIBOLIUM_EDMUNDS_2003`](@ref) or
[`TRIBOLIUM_DESHARNAIS_2005`](@ref); change one with `merge`, for example
`merge(TRIBOLIUM_EDMUNDS_2003, (c_el = [0.01 0.03; 0.03 0.01],))`. The state of each species is `[L, P, A]`.
"""
function flour_beetle_competition(parameters = TRIBOLIUM_EDMUNDS_2003)
    p = parameters
    names = p.species
    beetle(j) = flour_beetle(names[j]; b = p.b[j], mu_l = p.mu_l[j], mu_a = p.mu_a[j],
        eats = [names[i] => (larvae_eggs = p.c_el[i, j] / p.V, adults_eggs = p.c_ea[i, j] / p.V,
                             adults_pupae = p.c_pa[i, j] / p.V) for i in 1:2])
    return community(beetle(1), beetle(2))
end
