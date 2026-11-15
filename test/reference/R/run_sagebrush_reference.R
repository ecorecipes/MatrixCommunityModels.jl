## Regenerate the reference for the sagebrush steppe models of MatrixCommunityModels.jl (src/sagebrush.jl).
##
## Runs the multispecies integral projection model of Adler, Ellner & Levine (2010) for four plants of the
## sagebrush steppe in Idaho with the functions and parameters of the R package IPMdoit (Tredennick & Adler,
## https://github.com/atredennick/IPMdoit, commit 94124a4, MIT licence): its make_inits_ipm, crowd_overlap
## (with WrijG and WrijS), make_K_matrix, sum_cover and sum_N, and the vital-rate functions G, S and f of its
## Idaho example (inst/examples/import2ipm_*.r), which read the fitted parameters in inst/extdata.
##
## Usage: Rscript run_sagebrush_reference.R <IPMdoit checkout> <output dir> <data dir>
##   writes the reference files sagebrush_*.csv to <output dir> and the parameter files of data/sagebrush to
##   <data dir> (copies of IPMdoit's Growth_params_*.csv and Surv_params_*.csv, and recruitment.csv and
##   recruitment_years.csv extracted from its recruit_params_GnoG.csv and recSize_*.csv).
##
## Three things are changed from IPMdoit's Idaho example (run_idaho_example.R), which does not run as it stands:
##   - crowd_overlap is given the cover and population vectors, which the example leaves out;
##   - make_K_matrix (unlike make_K_matrix_ss) does not multiply by the mesh width h; that is done here, as
##     the midpoint rule needs (Ellner & Rees 2006);
##   - get_rpa divides the recruits by the species' cover and then f multiplies by the genet's area; the
##     same number is computed here without the division, so that it also holds for a species without
##     genets (an invader). For a species with cover the two agree.
## A "mean year" is added as year 23: the year effects on growth and survival set to zero and the
## recruitment intercept at its posterior mean (intcpt.mu).
args <- commandArgs(trailingOnly = TRUE)
ipmdoit <- normalizePath(args[1])
outdir <- args[2]
datadir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(boot))   # inv.logit, used by IPMdoit's survival function
num <- function(x) sprintf("%.17g", x)
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush(stdout()) }
save <- function(x, name, dir = outdir) write.csv(x, file.path(dir, name), row.names = FALSE, quote = FALSE)

for (f in c("make_inits_ipm.R", "estimate_crowding.R", "get_crowding.R", "util_functions.R", "make_iteration_matrix.R"))
    source(file.path(ipmdoit, "R", f))
## the global settings of IPMdoit's Idaho example
A <- 10000; spp_list <- c("ARTR", "HECO", "POSE", "PSSP"); iter_matrix_dims <- c(75, 75, 50, 50)
max_size <- c(3000, 202, 260, 225); Nyrs <- 22; n_spp <- 4
owd <- setwd(file.path(ipmdoit, "inst", "examples"))
source("import2ipm_Growth.r"); source("import2ipm_Survival.R"); source("import2ipm_Recruitment.r")
setwd(owd)

## ---------------------------------------------------------------- the parameter files of data/sagebrush
ext <- file.path(ipmdoit, "inst", "extdata")
for (s in spp_list) for (kind in c("Growth", "Surv"))
    file.copy(file.path(ext, paste0(kind, "_params_", s, ".csv")), datadir, overwrite = TRUE)
rec <- data.frame(species = spp_list, intercept = num(Rpars$intcpt.mu))
for (j in 1:n_spp) rec[[paste0("crowding_", spp_list[j])]] <- num(Rpars$dd[, j])   # dd[i, j]: effect of j on i
rec$size_mean <- num(Rpars$sizeMean); rec$size_variance <- num(Rpars$sizeVar)
save(rec, "recruitment.csv", datadir)
years <- data.frame(year = 1:Nyrs)
for (j in 1:n_spp) years[[spp_list[j]]] <- num(Rpars$intcpt.yr[, j])
save(years, "recruitment_years.csv", datadir)

## ---------------------------------------------------------------- the model
meanyr <- Nyrs + 1
Gpars$intcpt.yr <- rbind(Gpars$intcpt.yr, 0); Gpars$slope.yr <- rbind(Gpars$slope.yr, 0)
Spars$intcpt.yr <- rbind(Spars$intcpt.yr, 0); Spars$slope.yr <- rbind(Spars$slope.yr, 0)
Rpars$intcpt.yr <- rbind(Rpars$intcpt.yr, Rpars$intcpt.mu)
alphaG <- Gpars$alpha; alphaS <- Spars$alpha       # globals read by WrijG and WrijS
inits <- make_inits_ipm(n_spp, iter_matrix_dims, max_size)

## Without niche differences (an assumption of MatrixCommunityModels.jl, not taken from IPMdoit): a species is
## affected by the cover of every species as by its own. In recruitment, which reads the covers directly, the
## coefficient of every neighbour becomes the species' coefficient on itself. In growth and survival the crowding
## by species j is pi * cover_j / alpha[i, j], on the spatial scale of j, so the coefficient becomes
## nb[i, i] * alpha[i, j] / alpha[i, i]: a unit of cover of any species, placed at random, then crowds as a unit of
## the species' own cover would at its own scale.
no_niche <- function(P, field, scaled) {
    m <- P[[field]]
    for (i in 1:n_spp) m[i, ] <- m[i, i] * (if (scaled) P$alpha[i, ] / P$alpha[i, i] else 1)
    P[[field]] <- m
    P
}
pars <- list(niche = list(G = Gpars, S = Spars, R = Rpars),
             none = list(G = no_niche(Gpars, "nb", TRUE), S = no_niche(Spars, "nb", TRUE), R = no_niche(Rpars, "dd", FALSE)))

rpa_of <- function(cover, yr, Rp)     # get_rpa without the division by cover (see above)
    sapply(1:n_spp, function(i) 100 * exp(Rp$intcpt.yr[yr, i] + sqrt(100 * cover) %*% Rp$dd[i, ]) / A)
crowding <- function(nt) {
    cover <- sum_cover(inits$v, nt, inits$h, A)
    crowd_overlap(A, sum_N(nt, inits$h), NULL, inits$h, alphaG, alphaS, inits$WmatG, inits$WmatS, n_spp,
                  inits$Ctot, inits$Cr, inits$b.r, inits$expv, inits$r.U, inits$v.r, inits$v, cover, nt)
}
kernels <- function(nt, yr, P, only = 1:n_spp) {
    cover <- sum_cover(inits$v, nt, inits$h, A)
    W <- crowding(nt)
    rpa <- rpa_of(cover, yr, P$R)
    lapply(1:n_spp, function(i) if (i %in% only)
        inits$h[i] * make_K_matrix(inits$v[[i]], W$WmatG[[i]], W$WmatS[[i]], P$R, rpa, P$G, P$S, yr, i, inits$h))
}
step <- function(nt, yr, P) {
    cover <- sum_cover(inits$v, nt, inits$h, A)
    K <- kernels(nt, yr, P, only = which(cover > 0))      # a species without cover stays without (as IPMdoit)
    for (i in which(cover > 0)) nt[[i]] <- drop(K[[i]] %*% nt[[i]])
    nt
}
start <- function(present = 1:n_spp) {
    nt <- inits$v
    for (i in 1:n_spp) nt[[i]][] <- if (i %in% present) 0.001 else 0
    nt
}
covers <- function(nt) sum_cover(inits$v, nt, inits$h, A)
lambda <- function(K) max(Mod(eigen(K, only.values = TRUE)$values))
## numbers of genets per m2 in each cell: the state of MatrixCommunityModels.jl
numbers <- function(nt) lapply(1:n_spp, function(i) inits$h[i] * nt[[i]])

run <- function(label, present, T, P, yrs = rep(meanyr, T), every = 1) {
    nt <- start(present)
    rows <- list(data.frame(run = label, year = 0, t(covers(nt))))
    for (t in 1:T) {
        nt <- step(nt, yrs[t], P)
        if (t %% every == 0) rows[[length(rows) + 1]] <- data.frame(run = label, year = t, t(covers(nt)))
    }
    say(label, paste(signif(100 * covers(nt), 6), collapse = " "))
    list(rows = do.call(rbind, rows), nt = nt)
}
out_runs <- list()
keep <- function(r) { out_runs[[length(out_runs) + 1]] <<- r$rows; r$nt }

## ---------------------------------------------------------------- kernels and crowding at one state
states <- list()
record <- function(label, nt)       # a state, as numbers of genets per m2 in each cell
    states[[length(states) + 1]] <<- data.frame(run = label, species = rep(spp_list, iter_matrix_dims),
                                               cell = unlist(lapply(iter_matrix_dims, seq_len)), value = num(unlist(numbers(nt))))
state <- keep(run("community_30", 1:n_spp, 30, pars$niche))
record("community_30", state)
W <- crowding(state)
cw <- list()
for (i in 1:n_spp) for (j in 1:n_spp)
    cw[[length(cw) + 1]] <- data.frame(species = spp_list[i], neighbour = spp_list[j], cell = seq_len(iter_matrix_dims[i]),
                                       growth = num(W$WmatG[[i]][, j]), survival = num(W$WmatS[[i]][, j]))
save(do.call(rbind, cw), "sagebrush_crowding.csv")
kr <- list()
for (label in c("niche", "none")) for (yr in c(meanyr, 3)) {
    K <- kernels(state, yr, pars[[label]])
    for (i in 1:n_spp) {
        pick <- seq(1, iter_matrix_dims[i], by = 7)
        kr[[length(kr) + 1]] <- data.frame(species = spp_list[i], niche = label, year = ifelse(yr == meanyr, 0, yr),
                                           row = rep(pick, times = length(pick)), col = rep(pick, each = length(pick)),
                                           value = num(c(K[[i]][pick, pick])))
        kr[[length(kr) + 1]] <- data.frame(species = spp_list[i], niche = label, year = ifelse(yr == meanyr, 0, yr),
                                           row = 0, col = 0, value = num(lambda(K[[i]])))
    }
}
save(do.call(rbind, kr), "sagebrush_kernels.csv")
say("kernels done")

## ---------------------------------------------------------------- runs in the mean year
eq <- keep(run("community", 1:n_spp, 1500, pars$niche, every = 5))
record("community", eq)
for (i in 1:n_spp) keep(run(paste0("alone_", spp_list[i]), i, 1500, pars$niche, every = 5))
eq_none <- keep(run("community_no_niche", 1:n_spp, 3000, pars$none, every = 10))

## ---------------------------------------------------------------- invasion of the other three, in the mean year
inv <- list()
residents <- list()
for (label in c("niche", "none")) for (i in 1:n_spp) {
    res <- keep(run(paste0("residents_", label, "_without_", spp_list[i]), setdiff(1:n_spp, i), 1500, pars[[label]], every = 10))
    residents[[paste(label, i)]] <- res
    record(paste0("residents_", label, "_without_", spp_list[i]), res)
    K <- kernels(res, meanyr, pars[[label]], only = i)[[i]]
    inv[[length(inv) + 1]] <- data.frame(niche = label, invader = spp_list[i], growth_rate = num(lambda(K)),
                                         t(setNames(num(covers(res)), paste0("cover_", spp_list))))
}
save(do.call(rbind, inv), "sagebrush_invasion.csv")
save(do.call(rbind, states), "sagebrush_state.csv")

## ---------------------------------------------------------------- a fluctuating environment
set.seed(2010)
T <- 600
yrs <- sample(1:Nyrs, T, replace = TRUE)
save(data.frame(step = 1:T, year = yrs), "sagebrush_years.csv")
invisible(keep(run("community_years", 1:n_spp, 300, pars$niche, yrs = yrs)))
## the stochastic growth rate of each invader against the other three, from their mean-year equilibrium
sinv <- list()
burnin <- 100
for (i in 1:n_spp) {
    nt <- residents[[paste("niche", i)]]
    v <- rep(1 / iter_matrix_dims[i], iter_matrix_dims[i])
    total <- 0
    for (t in 1:T) {
        K <- kernels(nt, yrs[t], pars$niche, only = i)[[i]]
        v <- drop(K %*% v); g <- sum(v); v <- v / g
        if (t > burnin) total <- total + log(g)
        nt <- step(nt, yrs[t], pars$niche)
    }
    sinv[[i]] <- data.frame(invader = spp_list[i], burnin = burnin, steps = T, growth_rate = num(exp(total / (T - burnin))))
    say("stochastic invasion", spp_list[i], exp(total / (T - burnin)))
}
save(do.call(rbind, sinv), "sagebrush_stochastic_invasion.csv")

runs <- do.call(rbind, out_runs)
names(runs) <- c("run", "year", spp_list)
for (s in spp_list) runs[[s]] <- num(runs[[s]])
save(runs, "sagebrush_runs.csv")
say("done")
