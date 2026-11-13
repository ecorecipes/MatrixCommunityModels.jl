# Honey bees, Varroa, deformed wing virus and bumblebees: structure of the composed models
# and the behaviour they were calibrated to (see the bee vignettes for sources).

@testset "bees" begin
    adults(sol) = stage_total(sol, :honeybee, [:hive, :forager])
    year(k) = (365 * (k - 1) + 1):(365 * k + 1)

    @testset "honey bee colony" begin
        colony = honeybee()
        @test species_names(colony) == [:honeybee]
        names = stage_names(colony, :honeybee)
        @test length(names) == 3 + 5 + 12 + 2 + (3 + 7 + 14 + 1) + 2 && names[1] == :egg1 && names[end] == :pollen
        @test length(stage_names(honeybee(drones = false), :honeybee)) == 3 + 5 + 12 + 2 + 2
        @test Set(transition_names(system(colony)[:honeybee])) ⊇ Set([:hatching, :capping, :emerging, :recruiting, :larva_ageing, :honey_kept])
        sol = project(colony, honeybee_state(colony), 1:(365 * 3))
        a = adults(sol)
        @test 30_000 < maximum(a[year(2)]) < 40_000                  # summer peak
        @test 150 < argmax(a[year(2)]) < 250
        @test 7_000 < a[365 * 2 + 1] < 12_000                        # winter cluster
        @test 110_000 < sum(stage_total(sol, :honeybee, :pupa12)[year(2)]) < 160_000   # bees reared in a year
        @test stage_total(sol, :honeybee, :honey)[end] > 30_000      # stores carry the colony through winter
        @test maximum(a[year(3)]) ≈ maximum(a[year(2)]) rtol = 0.01  # a stable annual cycle
        # drones in summer only, a few hundred of them; no brood at all for a few weeks in midwinter
        d = stage_total(sol, :honeybee, :drone)
        @test 200 < maximum(d[year(2)]) < 1500 && d[365 * 2 + 1] < 1
        brood = stage_total(sol, :honeybee, [:egg, :larva, :pupa, :drone_egg, :drone_larva, :drone_pupa])
        @test 10 < count(<(1), brood[year(2)]) < 60
        # a colony too small to function dies, and its stores remain
        small = project(colony, honeybee_state(colony; adults = 500), 1:120)
        @test adults(small)[end] == 0 && stage_total(small, :honeybee, :honey)[end] > 20_000
    end

    @testset "varroa" begin
        c = community(honeybee(), varroa())
        @test species_names(c) == [:honeybee, :varroa]
        @test length(stage_names(c, :varroa)) == 1 + 12 + 15
        sol = project(c, honeybee_state(c; mites = 20), 1:(365 * 2))
        m = abundance(sol, :varroa)
        @test 0.025 < log(m[211] / m[151]) / 60 < 0.04               # about 0.03 per day with drone brood
        @test m[366] < m[300]                                        # decline over winter
        # without drone brood: about 0.02 per day, and of the order of the 12-fold a year of Martin (1998)
        w = community(honeybee(drones = false), varroa(drones = false))
        @test length(stage_names(w, :varroa)) == 13
        mw = abundance(project(w, honeybee_state(w; mites = 20), 1:(365 * 2)), :varroa)
        @test 0.015 < log(mw[211] / mw[151]) / 60 < 0.025
        @test 6 < mw[366] / 20 < 20
        @test m[366] > 3 * mw[366]
        # without virus the mites do not harm the bees in this model
        alone = project(honeybee(), honeybee_state(honeybee()), 1:(365 * 2))
        @test adults(sol) ≈ adults(alone) rtol = 1e-9
        # treatment removes phoretic mites
        t = community(honeybee(), varroa(treatment = d -> 270 <= d < 310 ? 0.115 : 0.0))
        @test abundance(project(t, honeybee_state(t; mites = 20), 1:365), :varroa)[end] < 0.05 * m[366]
    end

    @testset "deformed wing virus" begin
        function run(; treatment = d -> 0.0)
            c = community(honeybee(virus = true), varroa(virus = true, treatment = treatment))
            return c, project(c, honeybee_state(c; mites = 20, viral = 0.1), 1:(365 * 4))
        end
        c, sick = run()
        @test :hive_infected in stage_names(c, :honeybee) && :phoretic_viral in stage_names(c, :varroa)
        @test :hive_carrier in stage_names(c, :honeybee)
        @test Set(transition_names(system(c)[:honeybee])) ⊇ Set([:capping_healthy_to_infected, :hive_life_healthy_to_carrier,
                                                                  :drone_larva_ageing_healthy_to_infected])
        @test Set(transition_names(system(c)[:varroa])) ⊇ Set([:invading_clean_to_viral, :daughters_viral_to_clean,
                                                                :leaving_drone_clean_to_viral])
        @test :viral_phoretic in pool_names(c) && :infected_adults in pool_names(c)
        _, well = run(treatment = d -> 270 <= d < 310 ? 0.115 : 0.0)
        winter = 365 * 3 + 1
        @test adults(sick)[winter] < 5_000 < 8_000 < adults(well)[winter]        # untreated colonies dwindle
        infected(sol) = stage_total(sol, :honeybee, [:hive_infected, :forager_infected])
        @test infected(sick)[winter] / adults(sick)[winter] > 0.3
        @test infected(well)[winter] / adults(well)[winter] < 0.01
        @test stage_total(sick, :varroa, :viral)[winter] / abundance(sick, :varroa)[winter] > 0.95
        # bees infected as adults carry the virus without dying of it
        carriers(sol) = stage_total(sol, :honeybee, [:hive_carrier, :forager_carrier])
        @test carriers(sick)[winter] / adults(sick)[winter] > 0.2
        @test sum(stage_total(sick, :honeybee, [:egg_carrier, :larva_carrier, :pupa_carrier, :drone_carrier])) == 0
        # a cell takes at most four mites, however many there are (a little brood dies after they enter)
        cells, pupae = stage_total(sick, :varroa, :cell12), stage_total(sick, :honeybee, :pupa12)
        @test 3.5 < maximum(cells[pupae .> 1] ./ pupae[pupae .> 1]) < 4.1
        @test maximum(abundance(sick, :varroa)) > 10_000
        # virus without mites goes nowhere
        lone = honeybee(virus = true)
        @test sum(infected(project(lone, honeybee_state(lone), 1:400))) == 0
    end

    @testset "bumblebees" begin
        bb = bumblebee()
        @test length(stage_names(bb, :bumblebee)) == 2 + 28 + 3 + 22 + 1 + 26 + 1 + 30 + 1
        sol = project(bb, honeybee_state(bb; queens = 30), 1:(365 * 2))
        total(stage, r) = sum(stage_total(sol, :bumblebee, stage)[r[1:(end - 1)]])
        founded = total(:founding28, year(1))
        @test 0.08 < founded / 30 < 0.16                             # establishment of about 0.12
        @test 60 < total(:worker_brood22, year(1)) / founded + 8 < 95   # workers reared per colony
        @test 14 < total(:gyne_brood30, year(1)) / founded < 28         # new queens per colony
        @test 18 < total(:male_brood26, year(1)) / founded < 34         # males per colony
        # only hibernating queens survive the winter
        hib = stage_total(sol, :bumblebee, :hibernating)
        @test abundance(sol, :bumblebee)[310] ≈ hib[310] && hib[366] > 30
        # one generation a year: this year's new queens do not found nests
        workers = stage_total(sol, :bumblebee, :worker)
        @test 150 < argmax(workers[1:366]) < 230
        # nest sites cap the number of colonies
        crowded = project(bb, honeybee_state(bb; queens = 3000), 1:365)
        nests = stage_total(crowded, :bumblebee, [:founding, :growing, :reproducing, :ending])
        @test maximum(nests) <= 30 + 1e-9
        @test :infected_worker in stage_names(bumblebee(virus = true), :bumblebee)
    end

    @testset "shared flowers" begin
        function joint(mode; hives, years = 4, virus = false, treatment = d -> 270 <= d < 310 ? 0.115 : 0.0, kw...)
            f = floral_competition(mode; hives = hives, kw...)
            fl = virus ? flower_transmission(hives = hives) : nothing
            c = community(honeybee(virus = virus, forage = f, flowers = fl), varroa(virus = virus, treatment = treatment),
                          bumblebee(virus = virus, forage = f, flowers = fl), f.species...)
            return c, project(c, honeybee_state(c; mites = virus ? 100 : 0, viral = virus ? 0.1 : 0.0, queens = 100),
                              1:(365 * years))
        end
        queens(sol) = stage_total(sol, :bumblebee, :hibernating)[end]
        c, none = joint(:scramble; hives = 0)
        @test Set(pool_names(c)) ⊇ Set([:hb_nectar, :hb_pollen, :bb_nectar, :bb_pollen])
        _, many = joint(:scramble; hives = 10)
        @test queens(many) < queens(none)                            # exploitative competition
        # a refuge that only bumblebees can use protects them from the same hives
        _, refuge = joint(:preemptive; hives = 10, refuge = 0.3)
        @test queens(refuge) > queens(many)
        # two hives per square kilometre are within the default supply: they thrive and bumblebees are unaffected
        _, few = joint(:scramble; hives = 2)
        @test queens(few) ≈ queens(none) rtol = 0.01
        @test maximum(adults(few)[year(4)]) > 30_000
        # a standing crop adds a flowers species and its stock pools
        cs, crop = joint(:standing_crop; hives = 1, years = 2)
        @test :flowers in species_names(cs) && :nectar_stock in pool_names(cs)
        @test all(>=(0), crop.u)
        @test_throws ArgumentError floral_competition(:lottery)

        # virus reaches bumblebees only through flowers visited by infected honey bees
        cv, spill = joint(:scramble; hives = 1, virus = true, treatment = d -> 0.0)
        @test :infectious_visits in pool_names(cv)
        @test maximum(stage_total(spill, :bumblebee, :infected_worker)) > 0
        # about a tenth of bumblebee workers in midsummer, of the order found in field surveys
        day = 365 * 2 + 200
        @test 0.03 < stage_total(spill, :bumblebee, :infected_worker)[day] /
                     stage_total(spill, :bumblebee, [:worker, :infected_worker])[day] < 0.2
        _, clean = joint(:scramble; hives = 0, virus = true)
        @test maximum(stage_total(clean, :bumblebee, :infected_worker)) == 0
    end
end
