## Regenerate the adapt3 reference for MatrixCommunityModels.jl.
##
## Runs the R chunks of adapt3's community-projection vignette (comproj.Rmd) up to the
## point where the MPMs, start vectors, density frame and equivalence frame exist, then
## dumps the annual matrices and several project3() runs at full precision.
##
## Usage: Rscript run_adapt3_reference.R <path/to/comproj.Rmd> <output dir>
## Needs R packages lefko3 (6.7.3 used) and adapt3 (2.0.9 used).
suppressPackageStartupMessages({library(lefko3); library(adapt3)})
args <- commandArgs(trailingOnly = TRUE)
rmd <- readLines(args[1])
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Extract R chunks, dropping printed output, and stop before the project3() call.
starts <- grep("^``` r$", rmd); ends <- grep("^```$", rmd)
code <- character(0)
for (s in starts) {
  e <- ends[ends > s][1]
  chunk <- rmd[(s + 1):(e - 1)]
  chunk <- chunk[!grepl("^#>|^rm\\(list", chunk)]
  if (any(grepl("project3\\(", chunk))) break
  code <- c(code, chunk)
}
eval(parse(text = code))

options(digits = 17)
stages <- cypframe_raw$stage
fmt <- function(x) formatC(x, digits = 17, format = "g")

dump_mats <- function(mpm, sp) {
  do.call(rbind, lapply(seq_along(mpm$A), function(k) {
    A <- as.matrix(mpm$A[[k]])
    data.frame(species = sp, matrix = k, year = mpm$labels$year2[k],
               row = as.vector(row(A)), col = as.vector(col(A)), value = fmt(as.vector(A)))
  }))
}
write.csv(rbind(dump_mats(cycamatrix2r, 1), dump_mats(cypamatrix2r, 2)),
          file.path(outdir, "adapt3_matrices.csv"), row.names = FALSE, quote = FALSE)
writeLines(stages, file.path(outdir, "adapt3_stages.txt"))
write.csv(cyp2_dv, file.path(outdir, "adapt3_density_frame.csv"), row.names = FALSE)
write.csv(cyp_eq, file.path(outdir, "adapt3_equivalence_frame.csv"), row.names = FALSE)

dump_proj <- function(proj, tag) {
  comm <- proj$comm_out
  rows <- list()
  for (sp in seq_along(comm)) {
    m <- comm[[sp]][[1]]
    if (is.list(m)) m <- m[[1]]
    m <- as.matrix(m)
    rows[[sp]] <- data.frame(species = sp, time = as.vector(col(m)) - 1,
                             stage = as.vector(row(m)), value = fmt(as.vector(m)))
  }
  write.csv(do.call(rbind, rows), file.path(outdir, paste0("adapt3_", tag, "_stages.csv")),
            row.names = FALSE, quote = FALSE)
  N <- proj$N_out[[1]]
  write.csv(data.frame(time = 0:(ncol(N) - 1), sp1 = fmt(N[1, ]), sp2 = fmt(N[2, ]),
                       agg = fmt(as.vector(proj$agg_density[1, ]))),
            file.path(outdir, paste0("adapt3_", tag, "_N.csv")), row.names = FALSE, quote = FALSE)
}

## (a) the vignette call: cycled annual matrices, floor to integers
dump_proj(project3(mpms = cyp_mpm_list, starts = cyp_start_list, density = cyp_dv_list,
                   equivalence = eq_list, times = 100, integeronly = TRUE), "vignette")
## (b) the same without flooring
dump_proj(project3(mpms = cyp_mpm_list, starts = cyp_start_list, density = cyp_dv_list,
                   equivalence = eq_list, times = 100, integeronly = FALSE), "continuous")
## (c) no equivalence frame: every stage counts 1
dump_proj(project3(mpms = cyp_mpm_list, starts = cyp_start_list, density = cyp_dv_list,
                   times = 100, integeronly = FALSE), "noequiv")
cat("done\n")
