## Regenerate the reference for the LeMans model of MatrixCommunityModels.jl (src/lemans.jl).
##
## Runs the length-based multispecies model LeMans (Hall et al. 2006) for the 21 species of the North Sea
## with the R package LeMaRns 0.1.2 (Spence et al. 2020; CRAN, GPL-3), on its own data sets NS_par, NS_tau,
## NS_other, NS_mixed_fish, NS_eta and NS_L50, set up as in the package's vignette:
##
##   LeMansParam(df = NS_par, gdf = NS_mixed_fish, tau = NS_tau, eta = NS_eta, L50 = NS_L50, other = NS_other)
##
## with every other setting at its default (32 length classes, a time step phi_min = 0.1 year, hockey-stick
## recruitment, growth efficiency varying with length, natural mortality "std_RNM").
##
## Usage: Rscript run_lemans_reference.R <output dir> <data dir> [<R library>]
##   writes the reference files lemans_*.csv to <output dir>, and the data files of data/lemans to <data dir>
##   (the North Sea data sets of LeMaRns, written out at full precision).
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
datadir <- args[2]
lib <- if (length(args) >= 3) args[3] else NULL
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(LeMaRns, lib.loc = lib))
stopifnot(as.character(packageVersion("LeMaRns", lib.loc = lib)) == "0.1.2")
num <- function(x) sprintf("%.17g", x)
save <- function(x, name, dir = outdir) write.csv(x, file.path(dir, name), row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------- the data files of data/lemans
nfish <- nrow(NS_par)
sp <- as.character(NS_par$species_names)
gear_of <- as.character(NS_mixed_fish$gear_name)[match(sp, as.character(NS_mixed_fish$catch_species))]
catch_of <- NS_mixed_fish$max_catchability[match(sp, as.character(NS_mixed_fish$catch_species))]
species <- data.frame(species = sp, Linf = num(NS_par$Linf), W_a = num(NS_par$W_a), W_b = num(NS_par$W_b),
                      k = num(NS_par$k), Lmat = num(NS_par$Lmat), a = num(NS_par$a), b = num(NS_par$b),
                      gear = gear_of, max_catchability = num(catch_of), eta = num(NS_eta), L50 = num(NS_L50))
save(species, "north_sea_species.csv", datadir)
diet <- data.frame(predator = rownames(NS_tau))
for (j in seq_len(ncol(NS_tau))) diet[[colnames(NS_tau)[j]]] <- NS_tau[, j]
stopifnot(identical(rownames(NS_tau), sp), identical(colnames(NS_tau), sp))
save(diet, "north_sea_diet.csv", datadir)
save(data.frame(other = num(NS_other)), "north_sea_other_food.csv", datadir)

## ---------------------------------------------------------------- the model
params <- suppressWarnings(LeMansParam(df = NS_par, gdf = NS_mixed_fish, tau = NS_tau, eta = NS_eta,
                                       L50 = NS_L50, other = NS_other))
nsc <- params@nsc
gears <- dimnames(params@Qs)[[3]]
N0 <- get_N0(nsc, nfish, params@mid, params@wgt, params@sc_Linf, intercept = 1e10, slope = -5)
M2_0 <- calc_M2(N0, params@ration, params@wgt, nfish, nsc, params@other, params@sc_Linf, params@suit_M2)

## the quantities of each length class of each species, and the predation mortality at the starting state
classes <- do.call(rbind, lapply(seq_len(nfish), function(i) {
    m <- params@sc_Linf[i]
    Q <- params@Qs[1:m, i, ]
    data.frame(species = sp[i], class = 1:m, lower = num(params@l_bound[1:m]), upper = num(params@u_bound[1:m]),
               mid = num(params@mid[1:m]), weight = num(params@wgt[1:m, i]), leaving = num(params@phi[1:m, i]),
               ration = num(params@ration[1:m, i]), mature = num(params@mature[1:m, i]), M1 = num(params@M1[1:m, i]),
               catchability = num(rowSums(matrix(Q, nrow = m))), N0 = num(N0[1:m, i]), M2_0 = num(M2_0[1:m, i]))
}))
save(classes, "lemans_classes.csv")

## a few suitabilities, the share of a predator's diet that each prey class would make if prey were equally abundant
suit <- do.call(rbind, lapply(c(12, 20, 21), function(i) {
    m <- params@sc_Linf[i]
    do.call(rbind, lapply(c(2, floor(m / 2), m), function(j) {
        s <- params@suit_M2[[i]][j, , ]
        idx <- which(s > 0, arr.ind = TRUE)
        data.frame(predator = sp[i], predator_class = j, prey = sp[idx[, 2]], prey_class = idx[, 1],
                   suitability = num(s[idx]))
    }))
}))
save(suit, "lemans_suitability.csv")

## The steps that are kept from a run: the first 20, and the last of every year. Biomass and spawning stock
## biomass are in grams, after the step; recruits are the numbers added at the start of the step.
kept_steps <- function(steps) steps[steps <= 20 | steps %% 10 == 0]
summarise <- function(label, out) {
    steps <- kept_steps(seq_len(dim(out@N)[3] - 1))
    do.call(rbind, lapply(steps, function(s) {
        N <- out@N[, , s + 1]
        data.frame(run = label, step = s, species = sp, biomass = num(colSums(N * params@wgt)),
                   ssb = num(colSums(N * params@wgt * params@mature)), recruits = num(out@R[s + 1, ]))
    }))
}
final_state <- function(out, label) {
    last <- dim(out@N)[3]
    do.call(rbind, lapply(seq_len(nfish), function(i) {
        m <- params@sc_Linf[i]
        data.frame(run = label, species = sp[i], class = 1:m, N = num(out@N[1:m, i, last]),
                   M2 = num(out@M2[1:m, i, last]))
    }))
}

## 1. unfished, from LeMaRns's starting state, for 100 years (1000 steps)
unfished <- run_LeMans(params, N0 = N0, years = 100, effort = matrix(0, 100, length(gears), dimnames = list(NULL, gears)))
## 2. from there, the four fleets of the North Sea at an effort of 1 each, for 50 years
fleet <- matrix(1, 50, length(gears), dimnames = list(NULL, gears))
fished <- run_LeMans(params, N0 = unfished@N[, , 1001], years = 50, effort = fleet)
## 3. from the unfished state, the otter trawl alone at an effort of 1.5, for 50 years
otter <- matrix(0, 50, length(gears), dimnames = list(NULL, gears)); otter[, "Otter"] <- 1.5
otter_only <- run_LeMans(params, N0 = unfished@N[, , 1001], years = 50, effort = otter)

save(rbind(summarise("unfished", unfished), summarise("fished", fished), summarise("otter", otter_only)),
     "lemans_runs.csv")
save(rbind(final_state(unfished, "unfished"), final_state(fished, "fished"), final_state(otter_only, "otter")),
     "lemans_state.csv")

## the steps at which LeMaRns recruits, over 100 years
recruit_steps <- which(unfished@R[-1, 1] > 0)
save(data.frame(step = recruit_steps), "lemans_recruitment_steps.csv")
cat("gears:", gears, "\n")
cat("steps with recruitment in the first 40 years:", recruit_steps[recruit_steps <= 400], "\n")
