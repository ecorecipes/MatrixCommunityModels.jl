## Regenerate the reference and the parameters for the grassland models of MatrixCommunityModels.jl
## (src/grasslands.jl).
##
## The multispecies integral projection models of Chu & Adler (2015, Ecological Monographs 85:373-392) for the
## dominant perennials of five western North American grasslands, as run by Tredennick, de Mazancourt, Loreau &
## Adler (2017, Ecology 98:971-981) in their repository atredennick/community_synchrony (commit d01b584; its
## DESCRIPTION says "License: MIT"; it has no LICENSE file).
##
## Usage: Rscript run_grassland_reference.R <community_synchrony checkout> <output dir> <data dir>
##   - writes to <data dir> (data/grasslands) the fitted parameters of every site, read from the repository's
##     results/growth_params_list.RDS, surv_params_list.RDS and recruit_parameters.RDS with its own formatting
##     functions (R/import_regression_params_fxns.R), and the mean and variance of log recruit size from its
##     data/<site>/<species>/recSize.csv, as format_recruitment_params computes them;
##   - writes to <output dir> (test/reference) grassland_*.csv: runs of the Kansas community made with the
##     authors' own simulation code (runscripts/run_ipm_source.R, sourced as it stands, with the subroutines of
##     runscripts/run_ipm_simulations.R and R/ipm_kernel_vitalrate_fxns.R), and the authors' saved 2000-year
##     Kansas runs (results/ipm_comp_nocomp_sims.RDS) with the sequence of years that produced them
##     (results/randyr_list.RDS).
## A "mean year" is the authors' own constant environment (their `constant == TRUE`): year effects on growth and
## survival set to zero, the recruitment intercept at its mean (intcpt.mu). Removing interspecific interactions
## is theirs too (`spp_interact == FALSE`: every off-diagonal coefficient set to zero).
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1])
outdir <- args[2]; datadir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(boot))   # inv.logit, used by the authors' survival function
num <- function(x) sprintf("%.17g", x)
save <- function(x, name, dir = outdir) write.csv(x, file.path(dir, name), row.names = FALSE, quote = FALSE)
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }

source(file.path(repo, "R", "import_regression_params_fxns.R"))
source(file.path(repo, "R", "ipm_kernel_vitalrate_fxns.R"))
## the subroutines defined at the top of runscripts/run_ipm_simulations.R (function definitions only)
for (e in parse(file.path(repo, "runscripts", "run_ipm_simulations.R"), keep.source = FALSE))
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        (is.call(e[[3]]) && as.character(e[[3]][[1]]) %in% c("function", "Vectorize")))
        eval(e, globalenv())
stopifnot(exists("make.K.matrix"), exists("WrijG"), exists("sumCover"))

G_all <- readRDS(file.path(repo, "results", "growth_params_list.RDS"))
S_all <- readRDS(file.path(repo, "results", "surv_params_list.RDS"))
R_all <- readRDS(file.path(repo, "results", "recruit_parameters.RDS"))
## mesh of each site, as set in runscripts/run_ipm_simulations.R ("taken from Chu and Adler 2015")
mesh <- list(Arizona = list(c(50, 50), c(170, 40)), Idaho = list(c(50, 75, 50, 50), c(3000, 202, 260, 225)),
             Kansas = list(c(75, 50, 75), c(1656, 823, 2056)), Montana = list(c(75, 50, 5, 50), c(2500, 130, 22, 100)),
             NewMexico = list(c(50, 50), c(600, 1300)))
pars <- function(site) {
    spp <- names(G_all[[site]]); Nyrs <- nrow(G_all[[site]][[1]])
    list(spp = spp, Nyrs = Nyrs,
         G = format_growth_params(site, spp, Nyrs, G_all[[site]]),
         S = format_survival_params(site, spp, Nyrs, S_all[[site]]),
         R = format_recruitment_params(site, spp, Nyrs, R_all[[site]], file.path(repo, "data", site)))
}

## ---------------------------------------------------------------- data/grasslands
species_rows <- list(); year_rows <- list()
for (site in names(mesh)) {
    p <- pars(site); spp <- p$spp
    for (i in seq_along(spp)) {
        row <- data.frame(site = site, species = spp[i], cells = mesh[[site]][[1]][i], largest = mesh[[site]][[2]][i],
                          growth_intercept = num(p$G$intcpt[i]), growth_slope = num(p$G$slope[i]),
                          growth_sigma_a = num(p$G$sigma2.a[i]), growth_sigma_b = num(p$G$sigma2.b[i]),
                          survival_intercept = num(p$S$intcpt[i]), survival_slope = num(p$S$slope[i]),
                          recruitment_intercept = num(p$R$intcpt.mu[i]),
                          recruit_size_mean = num(p$R$sizeMean[i]), recruit_size_variance = num(p$R$sizeVar[i]))
        # neighbour j = 1..4 (blank beyond the site's species): effect of j on i, and the spatial scale of j
        for (j in 1:4) {
            have <- j <= length(spp)
            row[[paste0("growth_crowding_", j)]] <- if (have) num(p$G$nb[i, j]) else ""
            row[[paste0("growth_alpha_", j)]] <- if (have) num(p$G$alpha[i, j]) else ""
            row[[paste0("survival_crowding_", j)]] <- if (have) num(p$S$nb[i, j]) else ""
            row[[paste0("survival_alpha_", j)]] <- if (have) num(p$S$alpha[i, j]) else ""
            row[[paste0("recruitment_crowding_", j)]] <- if (have) num(p$R$dd[i, j]) else ""
        }
        species_rows[[length(species_rows) + 1]] <- row
        year_rows[[length(year_rows) + 1]] <- data.frame(site = site, species = spp[i], year = 1:p$Nyrs,
            growth_intercept = num(p$G$intcpt.yr[, i]), growth_slope = num(p$G$slope.yr[, i]),
            survival_intercept = num(p$S$intcpt.yr[, i]), survival_slope = num(p$S$slope.yr[, i]),
            recruitment_intercept = num(p$R$intcpt.yr[, i]))
    }
}
save(do.call(rbind, species_rows), "species.csv", datadir)
save(do.call(rbind, year_rows), "years.csv", datadir)
say("parameters written")

## ---------------------------------------------------------------- runs of Kansas with the authors' code
site <- "Kansas"
p <- pars(site)
spp_list <- p$spp; Nyrs <- p$Nyrs
randyr <- readRDS(file.path(repo, "results", "randyr_list.RDS"))[[site]]
run <- function(T, constant, interact, years = randyr) {
    Gpars <<- p$G; Spars <<- p$S; Rpars <<- p$R
    if (constant) {
        Rpars$intcpt.yr <<- matrix(Rpars$intcpt.mu, Nyrs, length(spp_list), byrow = TRUE)
        Gpars$intcpt.yr[] <<- 0; Gpars$slope.yr[] <<- 0; Spars$intcpt.yr[] <<- 0; Spars$slope.yr[] <<- 0
    }
    if (!interact) {
        off <- function(m) { d <- diag(m); m[] <- 0; diag(m) <- d; m }
        Rpars$dd <<- off(Rpars$dd); Gpars$nb <<- off(Gpars$nb); Spars$nb <<- off(Spars$nb)
    }
    tlimit <<- T + 1; randyrvec <<- years
    n_spp <<- Nspp <<- length(spp_list); A <<- 10000
    maxSize <<- mesh[[site]][[2]]; bigM <<- mesh[[site]][[1]]; NoOverlap.Inter <<- FALSE
    sink(tempfile()); source(file.path(repo, "runscripts", "run_ipm_source.R")); sink()
    list(cover = covSave, nt = nt, h = h)
}
runs <- list(); states <- list()
keep <- function(label, r) {
    runs[[length(runs) + 1]] <<- data.frame(run = label, year = 0:(nrow(r$cover) - 1), r$cover)
    say(label, paste(signif(100 * r$cover[nrow(r$cover), ], 6), collapse = " "))
    r
}
record <- function(label, r)       # numbers of genets per m2 in each cell (the state of MatrixCommunityModels.jl)
    for (i in seq_along(spp_list))
        states[[length(states) + 1]] <<- data.frame(run = label, species = spp_list[i], cell = seq_along(r$nt[[i]]),
                                                   value = num(r$h[i] * r$nt[[i]]))
keep("mean_30", run(30, TRUE, TRUE))
keep("mean_30_nointer", run(30, TRUE, FALSE))
keep("years_30", run(30, FALSE, TRUE))
keep("years_30_nointer", run(30, FALSE, FALSE))
record("mean_500", keep("mean_500", run(500, TRUE, TRUE)))
record("mean_500_nointer", keep("mean_500_nointer", run(500, TRUE, FALSE)))
out <- do.call(rbind, runs)
names(out) <- c("run", "year", spp_list)
for (s in spp_list) out[[s]] <- num(out[[s]])
save(out, "grassland_kansas_runs.csv")
save(do.call(rbind, states), "grassland_kansas_state.csv")

## ---------------------------------------------------------------- the authors' saved 2000-year runs
saved <- readRDS(file.path(repo, "results", "ipm_comp_nocomp_sims.RDS"))[[site]]
rows <- list()
for (sim in names(saved)) {
    x <- saved[[sim]]
    rows[[length(rows) + 1]] <- data.frame(run = sim, year = 500 + seq_len(nrow(x)), apply(x, 2, num))
}
saved_out <- do.call(rbind, rows); names(saved_out) <- c("run", "year", spp_list)
save(saved_out, "grassland_kansas_saved.csv")
save(data.frame(step = seq_along(randyr), year = randyr), "grassland_kansas_years.csv")
say("done")
