## Regenerate the reference for the annual plant model of MatrixCommunityModels.jl (src/annualplants.jl).
##
## The model of Godoy, Kraft & Levine (2014, Ecology Letters 17:836-844), written here directly from the paper's
## eqns. 1-5, in the paper's own form: one number per species, the seeds N in the soil before germination,
##   N[t+1] = (1 - g) s N[t] + g lambda N[t] / (1 + sum_j alpha[i, j] g[j] N[j, t]),
## with the authors' figshare data (doi:10.6084/m9.figshare.4793488, CC0) as shipped in data/annual_plants, and a
## missing coefficient filled with r_i e_j from Appendix S2, Table 2 (data/annual_plants/
## godoy2014_appendix_s2_table2.csv), as in Appendix S4. Base R only.
##
## Usage: Rscript run_annualplants_reference.R <data/annual_plants> <output dir>
##   writes annualplants_pairs.csv (the niche overlap, fitness ratio, invasion growth rates, predicted outcome and,
##   for coexisting pairs, the equilibrium germinants of each of the 153 pairs), annualplants_runs.csv (N by year
##   for three pairs and for all 18 species, 100 years from 100 seeds of each) and annualplants_ranks.csv (the
##   competitive ranks of the paper's sequential simulations, Fig. 2).
args <- commandArgs(trailingOnly = TRUE)
datadir <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

v <- read.csv(file.path(datadir, "species_vital_rates.csv"))
a <- read.csv(file.path(datadir, "alpha_estimates_row_is_target.csv"), check.names = FALSE, row.names = 1)
t2 <- read.csv(file.path(datadir, "godoy2014_appendix_s2_table2.csv"))
sp <- v$species
stopifnot(identical(rownames(a), sp), identical(colnames(a), sp), identical(t2$species, sp))
n <- length(sp)
lambda <- v$lambda; g <- v$g; s <- v$s
A <- as.matrix(a)
missing <- is.na(A)
A[missing] <- outer(t2$response, t2$effect)[missing]          # alpha_ij = r_i e_j
eta <- lambda * g / (1 - (1 - g) * s)                            # eqn. 4

f17 <- function(x) sprintf("%.17g", x)

## the pairs: eqns. 3-5 and the invasion growth rates
pairs <- NULL
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  rho <- sqrt(A[i, j] * A[j, i] / (A[i, i] * A[j, j]))
  ratio <- (eta[j] - 1) / (eta[i] - 1) * sqrt(A[i, j] * A[i, i] / (A[j, j] * A[j, i]))
  inv_i <- (1 - g[i]) * s[i] + g[i] * lambda[i] / (1 + A[i, j] * (eta[j] - 1) / A[j, j])
  inv_j <- (1 - g[j]) * s[j] + g[j] * lambda[j] / (1 + A[j, i] * (eta[i] - 1) / A[i, i])
  outcome <- if (inv_i > 1 && inv_j > 1) "coexistence" else if (inv_i > 1 || inv_j > 1) "exclusion" else "priority_effect"
  winner <- if (outcome == "exclusion") (if (inv_i > 1) sp[i] else sp[j]) else ""
  P <- c(NA, NA)
  if (outcome == "coexistence")                                   # 1 + alpha P = eta for both species
    P <- solve(A[c(i, j), c(i, j)], eta[c(i, j)] - 1)
  pairs <- rbind(pairs, data.frame(i = sp[i], j = sp[j], filled = sum(missing[c(i, j), c(i, j)]),
                                   rho = f17(rho), fitness_ratio = f17(ratio), invasion_i = f17(inv_i),
                                   invasion_j = f17(inv_j), outcome = outcome, winner = winner,
                                   plants_i = f17(P[1]), plants_j = f17(P[2])))
}
write.csv(pairs, file.path(outdir, "annualplants_pairs.csv"), row.names = FALSE, quote = FALSE)
cat("pairs predicted to coexist:", sum(pairs$outcome == "coexistence"), "; priority effects:",
    sum(pairs$outcome == "priority_effect"), "; rho > 1:", sum(as.numeric(pairs$rho) > 1), "\n")

## the paper's model, eqns. 1-2, for the species k
step <- function(N, k) (1 - g[k]) * s[k] * N + g[k] * lambda[k] * N / (1 + A[k, k, drop = FALSE] %*% (g[k] * N))
run <- function(k, years, N0 = 100) {
  N <- matrix(NA, years + 1, length(k)); N[1, ] <- N0
  for (t in 1:years) N[t + 1, ] <- as.vector(step(N[t, ], k))
  N
}

runs <- list(ceme_heco = c("CEME", "HECO"), aghe_agre = c("AGHE", "AGRE"), aghe_erci = c("AGHE", "ERCI"),
             all = sp)
out <- NULL
for (label in names(runs)) {
  k <- match(runs[[label]], sp)
  N <- run(k, 100)
  for (m in seq_along(k))
    out <- rbind(out, data.frame(run = label, species = sp[k[m]], t = 0:100, N = f17(N[, m])))
}
write.csv(out, file.path(outdir, "annualplants_runs.csv"), row.names = FALSE, quote = FALSE)

## competitive ranks (Methods, "Analyzing the relationship between phylogeny and competition"): simulate all
## species from 100 seeds each for 20000 years; the survivors (above 1e-3 of the most abundant) take the next
## rank, tied if several; remove them and repeat
left <- 1:n; rank <- 0; ranks <- rep(NA, n)
while (length(left) > 0) {
  rank <- rank + 1
  N <- run(left, 20000)[20001, ]
  top <- left[N > 1e-3 * max(N)]
  ranks[top] <- rank
  left <- setdiff(left, top)
}
write.csv(data.frame(species = sp, rank = ranks), file.path(outdir, "annualplants_ranks.csv"), row.names = FALSE,
          quote = FALSE)
print(data.frame(species = sp, rank = ranks))
