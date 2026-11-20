## Regenerate the data and the references for the canvasback and redhead model of MatrixCommunityModels.jl
## (src/waterfowl.jl): Gibson, Arnold, Buderman & Koons (2025), Ecological Monographs 95:e1642.
##
## Inputs, read at run time and not copied:
##   <dryad dir>  the files of Dryad doi:10.5061/dryad.hqbzkh1n9 (CC0): estimated_vital_rates.rdata (8000
##                posterior draws of the yearly vital rates and abundances) and model_files.rdata (the data
##                of the integrated population model)
##   <code dir>   the files of Zenodo doi:10.5281/zenodo.8076163 (MIT): IPM_code.R and LTRE_Code.R.R
##
## What it does:
## 1. For each draw, recovers the two constants that the draws do not store but the abundances imply: the
##    breeding success of yearling relative to older females (omega; the stored `bp` has 500000 columns and
##    is not aligned with the 8000 draws) and the ratio of the indirect to the direct harvest hazard of adults
##    (exp(beta.indirect)), both from the step 1961 -> 1962 of the females.
## 2. Writes the posterior means (data/waterfowl/*.csv) that the package ships.
## 3. Evaluates the authors' process equations, the lines of IPM_code.R that define the ducklings (NJ, NJF,
##    NJM, NAF, NAM) and the breeding abundances (N1F, N1M, N2F, N2M), read from the file and evaluated as
##    they stand, at the posterior means, for three single draws, and with the feedback of the package's
##    extension (fecundity responding to the projected female numbers through their departure from the
##    posterior-mean female numbers of the fit, with the slopes of Appendix S1, Table S6). Writes
##    test/reference/waterfowl_projection.csv and the rates of the three draws, waterfowl_draws.csv.
## 4. Runs the authors' transient LTRE (LTRE_Code.R.R, steps 1 to 4, read from the file and evaluated as
##    they stand) on the posterior means instead of the draws, and writes the contributions to
##    test/reference/waterfowl_ltre.csv.
##
## Usage: Rscript run_waterfowl_reference.R <dryad dir> <code dir> <package root>
args <- commandArgs(trailingOnly = TRUE)
dryad <- args[1]; code <- args[2]; root <- args[3]
suppressPackageStartupMessages(library(abind))
datadir <- file.path(root, "data", "waterfowl"); refdir <- file.path(root, "test", "reference")
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)

ev <- new.env(); load(file.path(dryad, "estimated_vital_rates.rdata"), envir = ev)
mf <- new.env(); load(file.path(dryad, "model_files.rdata"), envir = mf)
S <- dim(ev$fecundity)[3]; Y <- 60; years <- 1961:2021
sp <- c("c", "r"); spname <- c("canvasback", "redhead")
fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
rate <- function(name, k) ev[[sub("@", sp[k], name)]]   # years x draws

## the yearly rates of species k, draw j (or the posterior means when j is NULL)
rates_of <- function(k, j = NULL) {
  pick <- function(x) if (is.null(j)) rowMeans(x) else x[, j]
  fec <- if (is.null(j)) rowMeans(ev$fecundity[k, , ]) else ev$fecundity[k, , j]
  data.frame(fecundity = fec,
             sw_hy_f = pick(rate("s.winter_hy_@_f", k)), sw_ahy_f = pick(rate("s.winter_ahy_@_f", k)),
             sw_hy_m = pick(rate("s.winter_hy_@_m", k)), sw_ahy_m = pick(rate("s.winter_ahy_@_m", k)),
             kappa_hy_f = pick(rate("kappa_hy_@_f", k)), kappa_ahy_f = pick(rate("kappa_ahy_@_f", k)),
             kappa_hy_m = pick(rate("kappa_hy_@_m", k)), kappa_ahy_m = pick(rate("kappa_ahy_@_m", k)),
             ss_f = pick(rate("s.summer_@_f", k)), ss_m = pick(rate("s.summer_@_m", k)))
}

## 1. omega and the indirect-harvest ratio of each draw, from the females' step 1961 -> 1962
omega <- indirect <- indirect_m <- matrix(NA, 2, S)
for (k in 1:2) {
  fec <- ev$fecundity[k, 1, ]; sy <- ev$njf[k, , ]; asy <- ev$naf[k, , ]
  swh <- rate("s.winter_hy_@_f", k)[1, ]; swa <- rate("s.winter_ahy_@_f", k)[1, ]
  kh <- rate("kappa_hy_@_f", k)[1, ]; ka <- rate("kappa_ahy_@_f", k)[1, ]; ss <- rate("s.summer_@_f", k)
  omega[k, ] <- (sy[2, ] / (swh * (1 - kh) * sqrt(ss[2, ]) * 0.5 * fec) - asy[1, ]) / sy[1, ]
  indirect[k, ] <- log(asy[2, ] / (sqrt(ss[1, ]) * swa * sqrt(ss[2, ]) * (sy[1, ] + asy[1, ]))) / log(1 - ka)
  # the same ratio from the males, a check that one beta.indirect serves both sexes
  ssm <- rate("s.summer_@_m", k); kam <- rate("kappa_ahy_@_m", k)[1, ]
  indirect_m[k, ] <- log(ev$nam[k, 2, ] / (sqrt(ssm[1, ]) * rate("s.winter_ahy_@_m", k)[1, ] * sqrt(ssm[2, ]) *
                                            (ev$njm[k, 1, ] + ev$nam[k, 1, ]))) / log(1 - kam)
  cat(spname[k], ": omega mean", mean(omega[k, ]), "sd", sd(omega[k, ]), "; exp(beta.indirect) mean",
      mean(indirect[k, ]), "; largest female-male difference", max(abs(indirect[k, ] - indirect_m[k, ])),
      "; mean of bp[, 1:8000] as LTRE_Code.R uses it", mean(ev$bp[k, 1:S]), "\n")
}

## 2. the posterior means that the package ships
vr <- do.call(rbind, lapply(1:2, function(k) cbind(species = spname[k], year = 1961:2020, rates_of(k))))
fmtdf <- function(x) { for (v in names(x)) if (is.numeric(x[[v]]) && v != "year" && v != "draw") x[[v]] <- fmt(x[[v]]); x }
write.csv(fmtdf(vr), file.path(datadir, "vital_rates.csv"), row.names = FALSE, quote = FALSE)
ab <- do.call(rbind, lapply(1:2, function(k) {
  tot <- ev$njf[k, , ] + ev$naf[k, , ] + ev$njm[k, , ] + ev$nam[k, , ]
  data.frame(species = spname[k], year = years,
             sy_f = fmt(rowMeans(ev$njf[k, , ])), asy_f = fmt(rowMeans(ev$naf[k, , ])),
             sy_m = fmt(rowMeans(ev$njm[k, , ])), asy_m = fmt(rowMeans(ev$nam[k, , ])),
             total_q05 = fmt(apply(tot, 1, quantile, 0.05)), total_q95 = fmt(apply(tot, 1, quantile, 0.95)))
}))
write.csv(ab, file.path(datadir, "abundance.csv"), row.names = FALSE, quote = FALSE)
d <- mf$dat
pars <- data.frame(species = spname, omega = fmt(rowMeans(omega)), indirect_harvest = fmt(rowMeans(indirect)),
                   bp_ltre = fmt(rowMeans(ev$bp[, 1:S])),
                   density_mean = fmt(d$mu_BD), density_sd = fmt(d$sd_BD),
                   ratio_mean = fmt(d$mu_HD), ratio_sd = fmt(d$sd_HD))
write.csv(pars, file.path(datadir, "parameters.csv"), row.names = FALSE, quote = FALSE)
cov <- data.frame(year = years, density_canvasback = fmt(d$mean.BD[, 1]), density_redhead = fmt(d$mean.BD[, 2]),
                  ratio = fmt(d$mean.HD),
                  bpop_canvasback = fmt(c(d$y[1, ], NA, NA)), bpop_redhead = fmt(c(d$y[2, ], NA, NA)))
write.csv(cov, file.path(datadir, "covariates.csv"), row.names = FALSE, quote = FALSE)

## 3. the authors' process equations, read from IPM_code.R
ipm <- readLines(file.path(code, "IPM_code.R"))
take <- function(pattern) { i <- grep(pattern, ipm); stopifnot(length(i) == 1); ipm[i] }
births_lines <- sapply(c("^ *NJ\\[k,t\\] <-", "^ *NJF\\[k,t\\] <-", "^ *NJM\\[k,t\\] <-",
                         "^ *NAF\\[k,t\\] <-", "^ *NAM\\[k,t\\] <-"), take)
step_lines <- sapply(paste0("^ *", rep(c("N1F", "N1M", "N2F", "N2M"), 2), "\\[", rep(1:2, each = 4),
                            ",t\\] <- s\\."), take)
births_code <- parse(text = births_lines); step_code <- parse(text = step_lines)

# beta of conspecific female density and of female canvasbacks per female redhead on fecundity, Table S6
beta_density <- c(-0.12, -0.09); beta_ratio <- c(0.15, -0.16)
# the feedback is measured from the posterior-mean breeding females of the fit (2 x 61)
ref_females <- t(sapply(1:2, function(k) rowMeans(ev$njf[k, , ]) + rowMeans(ev$naf[k, , ])))

## run the equations from the 1961 state n1 (2 x 4: sy_f, asy_f, sy_m, asy_m) with rates (list of 2 data
## frames), omega and indirect (2 each); feedback = TRUE scales fecundity by the extension's response
process <- function(n1, rates, om, ind, feedback = FALSE) {
  e <- new.env()
  e$n.species <- 2; e$bp <- om
  for (v in c("N1F", "N2F", "N1M", "N2M", "NJ", "NJF", "NJM", "NAF", "NAM", "theta")) assign(v, matrix(NA, 2, Y), e)
  e$N1F[, 1] <- n1[, 1]; e$N2F[, 1] <- n1[, 2]; e$N1M[, 1] <- n1[, 3]; e$N2M[, 1] <- n1[, 4]
  for (k in 1:2) {
    r <- rates[[k]]; s <- sp[k]
    assign(paste0("s.winter_hy_", s, "_f"), r$sw_hy_f * (1 - r$kappa_hy_f), e)
    assign(paste0("s.winter_hy_", s, "_m"), r$sw_hy_m * (1 - r$kappa_hy_m), e)
    assign(paste0("s.winter_ahy_", s, "_fi"), r$sw_ahy_f * (1 - r$kappa_ahy_f)^ind[k], e)
    assign(paste0("s.winter_ahy_", s, "_mi"), r$sw_ahy_m * (1 - r$kappa_ahy_m)^ind[k], e)
    assign(paste0("s.summer_", s, "_f"), r$ss_f, e); assign(paste0("s.summer_", s, "_m"), r$ss_m, e)
  }
  for (tt in 2:Y) {
    e$t <- tt - 1
    females <- e$N1F[, tt - 1] + e$N2F[, tt - 1]
    for (k in 1:2) {
      m <- 1
      if (feedback) {
        z <- (females[k] - ref_females[k, tt - 1]) / d$sd_BD[k]
        h <- (females[1] / females[2] - ref_females[1, tt - 1] / ref_females[2, tt - 1]) / d$sd_HD
        m <- exp(beta_density[k] * z + beta_ratio[k] * h)
      }
      e$theta[k, tt - 1] <- rates[[k]]$fecundity[tt - 1] * m
      e$k <- k
      for (x in births_code) eval(x, e)
    }
    e$t <- tt
    for (x in step_code) eval(x, e)
  }
  e
}
rows <- function(label, e) do.call(rbind, lapply(1:2, function(k)
  data.frame(label = label, species = spname[k], year = 1961:2020, sy_f = fmt(e$N1F[k, ]),
             asy_f = fmt(e$N2F[k, ]), sy_m = fmt(e$N1M[k, ]), asy_m = fmt(e$N2M[k, ]))))
stored <- function(label, j) do.call(rbind, lapply(1:2, function(k)
  data.frame(label = label, species = spname[k], year = 1961:2020, sy_f = fmt(ev$njf[k, 1:Y, j]),
             asy_f = fmt(ev$naf[k, 1:Y, j]), sy_m = fmt(ev$njm[k, 1:Y, j]), asy_m = fmt(ev$nam[k, 1:Y, j]))))
mean_rates <- lapply(1:2, rates_of)
n1 <- t(sapply(1:2, function(k) c(mean(ev$njf[k, 1, ]), mean(ev$naf[k, 1, ]), mean(ev$njm[k, 1, ]), mean(ev$nam[k, 1, ]))))
out <- rbind(rows("posterior_mean", process(n1, mean_rates, rowMeans(omega), rowMeans(indirect))),
             rows("feedback", process(n1, mean_rates, rowMeans(omega), rowMeans(indirect), feedback = TRUE)))
draws <- c(1, 4000, 8000); dr <- NULL
for (j in draws) {
  rj <- lapply(1:2, rates_of, j = j)
  n1j <- t(sapply(1:2, function(k) c(ev$njf[k, 1, j], ev$naf[k, 1, j], ev$njm[k, 1, j], ev$nam[k, 1, j])))
  e <- process(n1j, rj, omega[, j], indirect[, j])
  err <- max(abs(c(e$N1F / ev$njf[, 1:Y, j], e$N2F / ev$naf[, 1:Y, j], e$N1M / ev$njm[, 1:Y, j], e$N2M / ev$nam[, 1:Y, j]) - 1))
  cat("draw", j, ": largest relative difference of the equations from the stored abundances", err, "\n")
  out <- rbind(out, rows(paste0("draw_", j), e), stored(paste0("stored_", j), j))
  dr <- rbind(dr, do.call(rbind, lapply(1:2, function(k)
    cbind(draw = j, species = spname[k], year = 1961:2020, fmtdf(rj[[k]]),
          omega = fmt(omega[k, j]), indirect_harvest = fmt(indirect[k, j])))))
}
write.csv(out, file.path(refdir, "waterfowl_projection.csv"), row.names = FALSE, quote = FALSE)
write.csv(dr, file.path(refdir, "waterfowl_draws.csv"), row.names = FALSE, quote = FALSE)

## 4. the authors' transient LTRE at the posterior means: two identical copies of the means stand for the
## draws (with one draw R drops a dimension the script indexes)
ltre <- readLines(file.path(code, "LTRE_Code.R.R"))
stopifnot(grepl("^noyears   <- max_year", ltre[27]), grepl("^# Step 2", ltre[115]), grepl("^sd_lamf", ltre[248]),
          grepl("^psr <- .5", ltre[379]), grepl("^\\}", ltre[588]), grepl("^total <- total_f", ltre[711]))
L <- new.env()
two <- function(x) cbind(rowMeans(x), rowMeans(x))
for (v in ls(ev)) if (is.matrix(ev[[v]]) && nrow(ev[[v]]) == Y) assign(v, two(ev[[v]]), L)
arr3 <- function(x) { m <- apply(x, c(1, 2), mean); abind(m, m, along = 3) }
for (v in c("fecundity", "naf", "nam", "njf", "njm")) assign(v, arr3(ev[[v]]), L)
L$bp <- cbind(rowMeans(ev$bp[, 1:S]), rowMeans(ev$bp[, 1:S]))
L$mod <- 2; L$n.chains <- 5; L$samples <- 2; L$max_year <- 60
invisible(capture.output(eval(parse(text = ltre[c(27:248, 379:588, 708:711)]), L)))
names18 <- c("F1", "F2", "Sfj_w", "Sfa_w", "Sfj_s", "Sfa_s", "Sfj_h", "Sfa_h", "nf1", "nf2",
             "Smj_w", "Sma_w", "Smj_s", "Sma_s", "Smj_h", "Sma_h", "nm1", "nm2")
lt <- do.call(rbind, lapply(1:2, function(k) rbind(
  data.frame(quantity = "contribution", parameter = names18, species = spname[k],
             value = fmt(sapply(names18, function(p) L[[paste0("cont_", p)]][1, k]))),
  data.frame(quantity = "sensitivity", parameter = names18, species = spname[k],
             value = fmt(sapply(names18, function(p) L[[paste0("sens_", p)]][1, k]))),
  data.frame(quantity = "total", parameter = "all", species = spname[k], value = fmt(L$total[1, k])),
  data.frame(quantity = "lambda", parameter = c("total", "female"), species = spname[k],
             value = fmt(c(L$mu_lamt[1, k], L$mu_lamf[1, k]))))))
write.csv(lt, file.path(refdir, "waterfowl_ltre.csv"), row.names = FALSE, quote = FALSE)
cat("LTRE at the posterior means: fecundity share", (L$cont_F1[1, ] + L$cont_F2[1, ]) / L$total[1, ],
    "; lambda", L$mu_lamt[1, ], "\n")
