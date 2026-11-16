## Regenerate the reference for the size-spectrum model of MatrixCommunityModels.jl (src/mizer.jl).
##
## Runs the North Sea model of the R package mizer 3.4.0 (Scott, Blanchard & Andersen 2014; CRAN, GPL-3), the
## MizerParams object NS_params that ships with it (12 species, 100 weight classes, a background resource on
## 226 classes, four fleets with knife-edge selectivity), with every setting as stored there: the default
## first-order upwind scheme (`project(..., method = "euler")`), no diffusion, Beverton-Holt reproduction and
## a semichemostat resource. The time step is mizer's default dt = 0.1 year.
##
## Usage: Rscript run_mizer_reference.R <output dir> <data dir> [<R library>]
##   writes the reference files mizer_*.csv to <output dir>, and the data files of data/mizer to <data dir>
##   (the parameters of NS_params, written out at full precision).
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
datadir <- args[2]
lib <- if (length(args) >= 3) args[3] else NULL
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(mizer))
stopifnot(as.character(packageVersion("mizer")) == "3.4.0")
num <- function(x) sprintf("%.17g", x)
save <- function(x, name, dir = outdir) write.csv(x, file.path(dir, name), row.names = FALSE, quote = FALSE)

p <- suppressMessages(validParams(NS_params))
sp <- p@species_params
no_sp <- nrow(sp)
w <- p@w
dw <- p@dw
wf <- p@w_full
dwf <- p@dw_full
stopifnot(all(p@w_min_idx == 1), max(getDiffusion(p)) == 0, p@second_order_w$flux == "upwind",
          !isTRUE(p@second_order_w$bin_average), p@rates_funcs$RDD == "BevertonHoltRDD",
          p@resource_dynamics == "resource_semichemostat", is.null(comment(p@pred_kernel)),
          all(p@ext_encounter == 0), length(p@other_dynamics) == 0, length(p@other_mort) == 0,
          length(p@other_encounter) == 0, all(sp$pred_kernel_type == "lognormal"))

## ---------------------------------------------------------------- the data files of data/mizer
cols <- c("w_min", "w_mat", "w_mat25", "w_max", "w_repro_max", "beta", "sigma", "alpha", "gamma", "q", "h", "n",
          "ks", "p", "k", "m", "z0", "erepro", "R_max", "interaction_resource", "k_vb")
species <- data.frame(species = sp$species)
for (cl in cols) species[[cl]] <- num(sp[[cl]])
save(species, "north_sea_species.csv", datadir)
inter <- data.frame(predator = sp$species)
for (j in seq_len(no_sp)) inter[[sp$species[j]]] <- num(p@interaction[, j])
save(inter, "north_sea_interaction.csv", datadir)
gp <- p@gear_params
save(data.frame(gear = gp$gear, species = gp$species, sel_func = gp$sel_func,
                knife_edge_size = num(gp$knife_edge_size), catchability = num(gp$catchability),
                initial_effort = num(p@initial_effort[as.character(gp$gear)])),
     "north_sea_gears.csv", datadir)
rp <- p@resource_params
save(data.frame(min_w = num(min(w)), max_w = num(max(w)), no_w = length(w), min_w_pp = num(min(wf)),
                kappa = num(rp$kappa), lambda = num(rp$lambda), r_pp = num(rp$r_pp), n = num(rp$n),
                w_pp_cutoff = num(rp$w_pp_cutoff)),
     "north_sea_resource.csv", datadir)
## the stored state of NS_params: numbers per g, by species and weight class, and of the resource
st <- data.frame(component = rep(sp$species, each = length(w)), class = rep(seq_along(w), no_sp),
                 density = num(as.vector(t(p@initial_n))))
st <- rbind(st, data.frame(component = "Resource", class = seq_along(wf), density = num(p@initial_n_pp)))
save(st, "north_sea_initial_state.csv", datadir)

## ---------------------------------------------------------------- derived arrays
long <- function(...) {
    arrays <- list(...)
    out <- data.frame(species = rep(sp$species, each = length(w)), class = rep(seq_along(w), no_sp))
    for (nm in names(arrays)) out[[nm]] <- num(as.vector(t(arrays[[nm]])))
    out
}
save(cbind(data.frame(class = seq_along(wf), w = num(wf), dw = num(dwf)),
           data.frame(rr_pp = num(p@rr_pp), cc_pp = num(p@cc_pp))), "mizer_grid.csv")
phis <- Re(t(mvfft(t(p@ft_pred_kernel_e), inverse = TRUE))) / length(wf)
save(long(search_vol = p@search_vol, intake_max = p@intake_max, metab = p@metab, psi = p@psi, mu_b = p@mu_b,
          maturity = p@maturity, f_mort = getFMort(p, effort = p@initial_effort)), "mizer_arrays.csv")
save(data.frame(species = rep(sp$species, each = length(wf)), offset = rep(seq_along(wf) - 1, no_sp),
                phi = num(as.vector(t(phis)))), "mizer_kernel.csv")

## ---------------------------------------------------------------- rates at the stored state, and one step
r <- getRates(p)
save(long(encounter = r$encounter, feeding_level = r$feeding_level, e_repro = r$e_repro, e_growth = r$e_growth,
          pred_mort = r$pred_mort, mort = r$mort), "mizer_rates.csv")
save(data.frame(species = sp$species, rdi = num(r$rdi), rdd = num(r$rdd)), "mizer_reproduction.csv")
step <- project_simple(p, steps = 1, dt = 0.1, effort = p@initial_effort, method = "euler")
save(data.frame(class = seq_along(wf), resource_mort = num(r$resource_mort), pred_rate_cod = num(r$pred_rate[11, ]),
                n_pp_after_step = num(step$n_pp)), "mizer_resource_rates.csv")
save(long(n_after_step = step$n), "mizer_step.csv")

## ---------------------------------------------------------------- runs
biomass <- function(n) as.vector(n %*% (w * dw))
runs <- NULL
states <- NULL
record <- function(name, sim) {
    times <- as.numeric(dimnames(sim@n)$time)
    for (k in seq_along(times)) {
        runs <<- rbind(runs, data.frame(run = name, year = times[k], species = sp$species,
                                        biomass = num(biomass(sim@n[k, , ]))))
    }
    last <- length(times)
    states <<- rbind(states, data.frame(run = name, component = rep(sp$species, each = length(w)),
                                        class = rep(seq_along(w), no_sp), density = num(as.vector(t(sim@n[last, , ])))))
    states <<- rbind(states, data.frame(run = name, component = "Resource", class = seq_along(wf),
                                        density = num(sim@n_pp[last, ])))
}
## from the stored state with the stored effort (Industrial 0, Pelagic 1, Beam 0.5, Otter 0.5), 20 years
record("stored_effort", project(p, t_max = 20, dt = 0.1, t_save = 1, progress_bar = FALSE))
## without fishing, 30 years
record("unfished", project(p, effort = 0, t_max = 30, dt = 0.1, t_save = 1, progress_bar = FALSE))
## the otter trawl at three times its stored effort, the others as stored, 30 years
record("otter_1.5", project(p, effort = c(Industrial = 0, Pelagic = 1, Beam = 0.5, Otter = 1.5),
                            t_max = 30, dt = 0.1, t_save = 1, progress_bar = FALSE))
save(runs, "mizer_runs.csv")

## the steady state found by running the dynamics (findSteadyState with its default solver "project"), with the
## fishing effort as stored and without fishing
steady_rows <- NULL
for (case in c("stored_effort", "unfished")) {
    eff <- if (case == "unfished") 0 else p@initial_effort
    s <- suppressWarnings(suppressMessages(findSteadyState(p, effort = eff, progress_bar = FALSE, info_level = 0)))
    conv <- attr(s, "convergence")
    steady_rows <- rbind(steady_rows, data.frame(case = case, years = num(conv$years),
                                                 converged = isTRUE(conv$converged),
                                                 distance = num(conv$distance)))
    states <- rbind(states, data.frame(run = paste0("steady_", case), component = rep(sp$species, each = length(w)),
                                       class = rep(seq_along(w), no_sp), density = num(as.vector(t(s@initial_n)))))
    states <- rbind(states, data.frame(run = paste0("steady_", case), component = "Resource", class = seq_along(wf),
                                       density = num(s@initial_n_pp)))
}
save(steady_rows, "mizer_steady.csv")
save(states, "mizer_states.csv")
cat("done\n")
