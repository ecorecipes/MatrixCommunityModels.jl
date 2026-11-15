## Regenerate the matreex reference for the forest models of MatrixCommunityModels.jl.
##
## Builds the integral projection models of Norway spruce, silver fir and European beech with the R
## package matreex (0.4.0 used, https://github.com/gowachin/matreex) at the climate of its introductory
## vignette, and writes: values of the vital rates, kernel entries at fixed basal area, runs of one, two
## and three species from young stands with their equilibria, and growth rates of rare invaders.
##
## Usage: Rscript run_forest_reference.R <output dir> <cells> <settings> [<R library holding matreex>]
##   cells: 700 is the mesh of matreex.
##   settings: "default" for matreex's own numerical settings (kernels tabulated at whole basal areas and
##     interpolated between them, 140 quadrature points for each cell a tree can reach by growing less than
##     50 mm, a coarse rule for the 25 cells beyond), or "fine" (kernels every 0.25 m2/ha, 2000 points for
##     each cell within 200 mm, nothing beyond), which removes most of the numerical error of matreex.
##
## The files in test/reference were made with
##   Rscript run_forest_reference.R <dir> 100 fine      (forest_*_100.csv and forest_vital_rates.csv)
##   Rscript run_forest_reference.R <dir> 700 default   (forest_*_700.csv)
## The kernels are kept in <output dir>/ipm_cache_<cells>.rds and reused if that file exists.
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
m <- as.integer(args[2])
fine <- args[3] == "fine"
if (length(args) > 3) .libPaths(c(args[4], .libPaths()))
suppressPackageStartupMessages({library(matreex); library(Matrix)})
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
options(W_matreex_edist = FALSE)
num <- function(x) sprintf("%.17g", x)
short <- function(x) ifelse(is.na(x), "NA", sprintf("%.12g", x))
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
save <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, quote = FALSE)

kinds <- c(spruce = "Picea_abies", fir = "Abies_alba", beech = "Fagus_sylvatica")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
BAgrid <- if (fine) seq(0, 100, by = 0.25) else 0:100
harvest_rate <- 0.006                       # matreex's default: def_harv()
trunk_area <- function(d) pi * (d / 2000)^2 # m2, diameter in mm

## ---------------------------------------------------------------- vital rates, from matreex's own functions
if (fine) {
    covs <- drop(as.matrix(climate))
    rows <- list()
    for (key in names(kinds)) {
        f <- get(data(list = paste0("fit_", kinds[[key]])))
        for (BA in c(0, 20, 45.5)) {
            gr <- matreex:::exp_sizeFun(f$gr$params_m, c(covs, BATOTcomp = BA))
            sv <- matreex:::exp_sizeFun(f$sv$params_m, c(covs, BATOTcomp = BA))
            for (size in c(90, 150.5, 400, 1000)) {
                rows[[length(rows) + 1]] <- data.frame(species = key, rate = "log_increment", a = size, b = BA, value = num(gr(size)))
                rows[[length(rows) + 1]] <- data.frame(species = key, rate = "death", a = size, b = BA,
                                                       value = num(f$sv$family$linkinv(sv(size))))
            }
        }
        rec <- matreex:::exp_recFun(f$rec$params_m, covs)
        for (own in c(0.1, 1, 12.5, 40)) for (others in c(0, 7.25, 30)) {
            rows[[length(rows) + 1]] <- data.frame(species = key, rate = "recruits_per_ha", a = own, b = others,
                                                   value = num(sum(rec(own, others, 1:5, 0.03)) / 0.03))
        }
    }
    save(do.call(rbind, rows), "forest_vital_rates.csv")
}

## ---------------------------------------------------------------- integral projection models
cache <- file.path(outdir, sprintf("ipm_cache_%d.rds", m))
if (file.exists(cache)) {
    ipms <- readRDS(cache)
    say("kernels read from", cache)
} else {
    ipms <- list()
    for (key in names(kinds)) {
        fit <- get(data(list = paste0("fit_", kinds[[key]])))
        mesh_arg <- c(m = m, L = 90, U = get_maxdbh(fit) * 1.1)
        t0 <- Sys.time()
        ipms[[key]] <- if (fine) {
            make_IPM(species = kinds[[key]], climate = climate, fit = fit, clim_lab = "spruce optimum", mesh = mesh_arg,
                     BA = BAgrid, level = 6000, diag_tresh = 200, midbin_tresh = 0)
        } else {
            make_IPM(species = kinds[[key]], climate = climate, fit = fit, clim_lab = "spruce optimum", mesh = mesh_arg, BA = BAgrid)
        }
        say(key, "make_IPM,", m, "cells,", length(BAgrid), "basal areas:", format(difftime(Sys.time(), t0, units = "secs")))
    }
    saveRDS(ipms, cache)
}
delay <- sapply(ipms, function(x) as.integer(x$info["delay"]))
mesh <- lapply(ipms, function(x) x$mesh)          # zeros for the delay classes, then the middles of the cells

## ---------------------------------------------------------------- kernels at fixed basal area
## For five cells of origin: the entries within the cells integrated by Gauss-Legendre quadrature, and the
## sum of the whole column (to = 0), which is the share of the trees that stay on the mesh alive.
rows <- list()
for (key in names(kinds)) for (BA in c(0, 30)) {
    K <- as.matrix(ipms[[key]]$IPM[[as.character(BA)]])
    size <- (delay[[key]] + 1):nrow(K)
    P <- K[size, size]
    band <- as.integer(ipms[[key]]$int["gl_tresh"])
    for (from in c(1, m %/% 4, m %/% 2, 3 * m %/% 4, m)) {
        to <- from:min(m, from + band - 1)
        rows[[length(rows) + 1]] <- data.frame(species = key, basal_area = BA, to = c(to, 0), from = from,
                                               value = num(c(P[to, from], sum(P[, from]))))
    }
}
save(do.call(rbind, rows), sprintf("forest_kernel_%d.csv", m))

## ---------------------------------------------------------------- runs
## a young stand: no seedlings or saplings, trees spread evenly over the cells between 90 and 150 mm
stand <- function(key, basal_area) {
    d <- mesh[[key]]
    present <- d >= 90 & d < 150
    present * basal_area / sum(trunk_area(d[present]))
}
no_harvest <- function(x, species, ...) x * 0
simulate <- function(keys, tlim, harvest = TRUE, settle = TRUE) {
    sp <- lapply(keys, function(key) species(IPM = ipms[[key]], init_pop = def_init_k(stand(key, 1)),
                                             harvest_fun = if (harvest) def_harv else no_harvest))
    t0 <- Sys.time()
    sim <- as.data.frame(sim_deter_forest(forest(species = sp), tlim = tlim, equil_time = if (settle) 50000 else tlim,
                                          equil_dist = 100, equil_diff = 1e-9, SurfEch = 0.03))
    say(paste(keys, collapse = "+"), if (harvest) "" else "without harvest", "simulated to year", max(sim$time[sim$equil]), "in",
        format(difftime(Sys.time(), t0, units = "secs")))
    sim$key <- names(kinds)[match(sim$species, kinds)]
    sim
}
## basal area of each species (m2/ha) every `step` years, and at equilibrium (year -1) if the run settled
basal_areas <- function(sim, run, step = 5, settled = TRUE) {
    ba <- sim[sim$var == "BAsp", ]
    years <- seq(0, max(ba$time[!ba$equil]) - 1, by = step)
    out <- data.frame(run = run, year = c(years, if (settled) -1))
    for (key in names(kinds)) {
        rec <- ba[ba$key == key & !ba$equil, ]
        out[[key]] <- if (nrow(rec) == 0) NA else short(c(rec$value[match(years + 1, rec$time)], if (settled) ba$value[ba$key == key & ba$equil]))
    }
    out
}
## the equilibrium of each species: trees per hectare by class, the delay classes first
equilibrium_state <- function(sim, run) {
    eq <- sim[sim$var == "n" & sim$equil, ]
    eq <- eq[order(match(eq$key, names(kinds)), eq$mesh), ]
    data.frame(run = run, species = eq$key, class = eq$mesh, value = short(eq$value))
}
runs <- list(); states <- list()
for (key in names(kinds)) {                                   # each species alone
    sim <- simulate(key, 500)
    runs[[length(runs) + 1]] <- basal_areas(sim, key)
    if (key == "spruce") states[[length(states) + 1]] <- equilibrium_state(sim, key)
}
sim <- simulate(c("spruce", "fir"), 500)                      # two species
runs[[length(runs) + 1]] <- basal_areas(sim, "pair")
sim <- simulate(c("spruce", "fir", "beech"), 500)             # three species
runs[[length(runs) + 1]] <- basal_areas(sim, "trio")
states[[length(states) + 1]] <- equilibrium_state(sim, "trio")
sim <- simulate("spruce", 3000, harvest = FALSE, settle = FALSE)   # spruce without harvest, which does not settle
runs[[length(runs) + 1]] <- basal_areas(sim, "unharvested", step = 10, settled = FALSE)
save(do.call(rbind, runs), sprintf("forest_runs_%d.csv", m))
save(do.call(rbind, states), sprintf("forest_states_%d.csv", m))

## ---------------------------------------------------------------- a rare invader
## The frozen model of a species that holds `own` m2/ha of a stand of `stand` m2/ha: survival and growth P
## (with the delay classes and the harvest of matreex) and births F, in which recruits are shared among the
## trees in proportion to their basal area, as in Compute_R_0() of the scripts of Guyennon et al. (2023).
if (fine) {
    rows <- list()
    for (key in names(kinds)) for (stand_ba in c(30, 60)) {
        P <- as.matrix(ipms[[key]]$IPM[[as.character(stand_ba)]])
        d <- mesh[[key]]
        P <- P * (1 - harvest_rate * (d > 0))                       # rows of the size classes: felled after growth
        rec <- matreex:::exp_recFun(ipms[[key]]$fit$rec$params_m, drop(as.matrix(climate)))
        own <- 0.1
        recruits <- rec(own, stand_ba - own, d, 0.03) / 0.03        # per hectare, in the first two classes
        Fm <- outer(recruits, trunk_area(d)) / own
        lambda <- max(Re(eigen(P + Fm, only.values = TRUE)$values))
        R0 <- max(Re(eigen(Fm %*% solve(diag(nrow(P)) - P), only.values = TRUE)$values))
        rows[[length(rows) + 1]] <- data.frame(species = key, stand = stand_ba, own = own, lambda = num(lambda), R0 = num(R0))
    }
    save(do.call(rbind, rows), sprintf("forest_invasion_%d.csv", m))
}
say("done")
print(sessionInfo())
