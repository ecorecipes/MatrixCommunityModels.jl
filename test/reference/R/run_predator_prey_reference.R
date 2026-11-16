## Regenerate the reference for the predator-prey model of MatrixCommunityModels.jl (src/predatorprey.jl).
##
## Runs the deterministic difference equations of Barraquand & Gimenez (2019, Appendix A), the pre-breeding
## two-species matrix model, with the authors' own scripts from GitHub oliviergimenez/predator_prey_icm
## (GPL-3), folder difference_equations/predator_prey/pre_breeding. Nothing of those scripts is copied here:
## each is read from the clone at run time, its plotting lines are dropped, its number of years is set, and
## it is evaluated as it stands. The trajectory N (4 x tmax: juvenile prey, adult prey, juvenile predators,
## adult predators) is written out.
##
## Usage: Rscript run_predator_prey_reference.R <clone of predator_prey_icm> <output dir>
##   writes predprey_trajectories.csv to <output dir>, one row per script and year.
args <- commandArgs(trailingOnly = TRUE)
repo <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(Matrix))
pre <- file.path(repo, "difference_equations", "predator_prey", "pre_breeding")

## label = c(script, alpha7 to impose or NA, years)
runs <- list(
  set1_appendix = list(file = "DDmatrix_twospecies_BGpredation_prebreed.R", alpha7 = NA, tmax = 100),
  set1_table5 = list(file = "DDmatrix_twospecies_BGpredation_prebreed.R", alpha7 = 2, tmax = 100),
  set2 = list(file = "DDmatrix_twospecies_BGpredation_prebreed_newparam.R", alpha7 = NA, tmax = 100),
  set1_high = list(file = "increased_fecundity/DDmatrix_twospecies_BGpredation_prebreed.R", alpha7 = NA, tmax = 500),
  set2_high = list(file = "increased_fecundity/DDmatrix_twospecies_BGpredation_prebreed_newparam.R", alpha7 = NA, tmax = 500),
  no_intra = list(file = "DDmatrix_twospecies_BGpredation_prebreed_woutIntraDD.R", alpha7 = NA, tmax = 100))

out <- NULL
for (label in names(runs)) {
  r <- runs[[label]]
  src <- readLines(file.path(pre, r$file))
  src <- src[!grepl("^(pdf|dev\\.off|par|matplot|matlines|library)\\(", src)]
  src <- sub("^tmax *= *[0-9]+", paste("tmax =", r$tmax), src)
  if (!is.na(r$alpha7)) src <- sub("^alpha\\[7\\] *<- *[0-9.]+", paste("alpha[7] <-", r$alpha7), src)
  env <- new.env()
  eval(parse(text = src), envir = env)
  N <- env$N
  a <- env$alpha
  cat(label, ": alpha =", a, "; phiA_V =", env$phiA_V, "; phiA_P =", env$phiA_P,
      "; N at the end =", N[, ncol(N)], "\n")
  out <- rbind(out, data.frame(run = label, t = seq_len(ncol(N)) - 1,
                               prey_juvenile = sprintf("%.17g", N[1, ]), prey_adult = sprintf("%.17g", N[2, ]),
                               predator_juvenile = sprintf("%.17g", N[3, ]), predator_adult = sprintf("%.17g", N[4, ])))
}
write.csv(out, file.path(outdir, "predprey_trajectories.csv"), row.names = FALSE, quote = FALSE)
