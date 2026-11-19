## Regenerate the reference for the coral models of MatrixCommunityModels.jl (src/corals.jl).
##
## The multispecies, open-population integral projection model of Kayal, Lenihan, Brooks, Holbrook, Schmitt &
## Kendall (2018, Ecology Letters 21:1790-1799) for Acropora, Pocillopora and Porites on four fore reefs of
## Mo'orea (LTER 1, 2, 4 and 5).
##
## Usage: Rscript run_coral_reference.R <Coral-assemblage-IPM.R> <output dir> [<mcr_lter_observed_colonies.csv>]
##   <Coral-assemblage-IPM.R> is the authors' script (https://github.com/MohsenKayal/Coral-assemblage-IPM, commit
##   af547a8; Zenodo doi:10.5281/zenodo.573772). It carries no licence, so it is not copied into the package: this
##   driver reads it at run time and evaluates only its survival, growth, survival-growth kernel and recruit-size
##   functions (S_z*, G_z1z*, P_z1z*, c_z1z*). The script does not run as it stands (it fits models to a data file
##   it does not name, and its recruitment functions call objects that do not exist), so the recruitment of
##   Appendix S3 (equations iii to vi) and the projection loop are written here.
##
## The parameters are those of Table S1 of the paper's Supporting Information, typed below. Recruitment,
## log r = R0 + Rg * gamma, is computed in two readings:
##   "as_written": natural logarithm (as the authors' script, which uses exp), recruits and crowding per 10 m2
##                 (the area of a model iteration), gamma integrated over the mesh as in Appendix S3;
##   "calibrated": base-10 logarithm, R0 for the 25 m2 plots in which recruits were counted (so recruits per
##                 10 m2 are 0.4 * 10^(R0 + Rg * gamma)), and gamma per m2 (a tenth of its value per 10 m2).
##                 This reading is chosen by MatrixCommunityModels.jl because it reproduces the published
##                 compositions (Fig. 3, Fig. S4) and densities (Fig. S2); see the vignette.
## Both cap recruitment at 348 colonies per 10 m2 a year (Appendix S3). Recruits have log10(area + 1) normal with
## mean 0.6 and sd 0.2 (the initial runs). The reef starts empty and the first recruits arrive in the first step.
##
## If the third argument is given (data/corals/mcr_lter_observed_colonies.csv of the package: the colonies observed
## by the MCR LTER in 2013 and 2015, binned onto the mesh), the updated runs of the paper (Fig. 4) are also written to
## coral_observed_runs.csv: from the state observed in 2013, two steps (to 2015), and from that observed in 2015, five
## steps (to 2020), with Acropora recruits of mean 1.2 and sd 0.4 (Appendix S3), in both readings. Acropora was mapped
## over 20 m2 of each reef and the other two genera over 4 m2, as the file's sampled_area_m2 column says.
args <- commandArgs(trailingOnly = TRUE)
authors <- normalizePath(args[1])
outdir <- args[2]
observed_file <- if (length(args) >= 3) args[3] else NA
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
num <- function(x) sprintf("%.17g", x)
save <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------- the authors' functions, read at run time
code <- parse(authors, keep.source = FALSE)
wanted <- c(paste0("S_z", c("Acr", "Poc", "Por")), paste0("G_z1z", c("Acr", "Poc", "Por")),
            paste0("P_z1z", c("Acr", "Poc", "Por")), paste0("c_z1z", c("Acr", "Poc", "Por")))
for (e in code) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% wanted)
        eval(e, globalenv())
}
stopifnot(all(sapply(wanted, exists)))
P_fun <- list(Acr = P_z1zAcr, Poc = P_z1zPoc, Por = P_z1zPor)
c_fun <- list(Acr = c_z1zAcr, Poc = c_z1zPoc, Por = c_z1zPor)

## ---------------------------------------------------------------- Table S1 (Kayal et al. 2018, Supporting Information)
taxa <- c("Acr", "Poc", "Por")
reefs <- c(1, 2, 4, 5)
S1 <- read.table(header = TRUE, text = "
reef taxon recr_int recr_slope surv_int surv_slope grow_int grow_slope grow_sd
1 Acr 0.966  0.040  2.688  0.219 1.542 0.341 0.245
1 Poc 2.688 -0.002  2.400  0.308 1.023 0.748 0.153
1 Por 0.978  0.003  0.758  1.865 1.082 0.637 0.254
2 Acr 0.852  0.040  1.085  0.821 1.168 0.555 0.157
2 Poc 2.475 -0.002  2.234  0.154 0.934 0.765 0.160
2 Por 0.673  0.003 -0.171  1.240 0.644 0.808 0.300
4 Acr 0.707  0.040 -0.350  1.327 1.057 0.617 0.173
4 Poc 1.724 -0.002  0.850  0.699 1.093 0.733 0.163
4 Por 0.139  0.003 -1.417  1.047 1.075 0.634 0.321
5 Acr 0.280  0.040  2.455 -0.451 1.008 0.667 0.199
5 Poc 1.889 -0.002  2.651 -0.564 1.102 0.686 0.184
5 Por 0.669  0.003  0.685  0.759 0.689 0.798 0.273
")
recruit_size <- list(Acr = c(0.6, 0.2), Poc = c(0.6, 0.2), Por = c(0.6, 0.2))
m.par <- function(reef, taxon) {
    r <- S1[S1$reef == reef & S1$taxon == taxon, ]
    c(surv.Intercept = r$surv_int, surv.Slope = r$surv_slope, grow.Intercept = r$grow_int, grow.Slope = r$grow_slope,
      grow.sd = r$grow_sd, recr.Intercept = r$recr_int, recr.Slope = r$recr_slope,
      recr.zmean = recruit_size[[taxon]][1], recr.zsd = recruit_size[[taxon]][2])
}

## ---------------------------------------------------------------- the model (Appendix S3)
L <- 0; U <- 7; m <- 100; h.mesh <- (U - L) / m
meshpts <- L + (1:m) * h.mesh - h.mesh / 2
cap <- 348
kernel <- function(reef, taxon) h.mesh * outer(meshpts, meshpts, P_fun[[taxon]], m.par(reef, taxon))
recruits <- function(reef, taxon, n, reading) {
    # n: densities (per unit log size, per 10 m2) of the three taxa
    p <- m.par(reef, taxon)
    gamma <- if (taxon == "Poc") sum(meshpts * (n$Acr + n$Poc + n$Por)) * h.mesh else sum(n[[taxon]]) * h.mesh
    r <- if (reading == "as_written") exp(p["recr.Intercept"] + p["recr.Slope"] * gamma)
         else 0.4 * 10^(p["recr.Intercept"] + p["recr.Slope"] * gamma / 10)
    min(cap, unname(r))
}
project <- function(reef, reading, T, n = lapply(setNames(taxa, taxa), function(k) rep(0, m))) {
    P <- lapply(setNames(taxa, taxa), function(k) kernel(reef, k))
    cz <- lapply(setNames(taxa, taxa), function(k) c_fun[[k]](meshpts, m.par(reef, k)))
    states <- list(n)
    for (t in 1:T) {
        r <- sapply(taxa, function(k) recruits(reef, k, n, reading))
        n <- lapply(setNames(taxa, taxa), function(k) drop(P[[k]] %*% n[[k]]) + r[[k]] * cz[[k]])
        states[[t + 1]] <- n
    }
    states
}

## ---------------------------------------------------------------- reference files
T <- 40
runs <- list(); cells <- list(); kern <- list()
for (reef in reefs) {
    for (k in taxa) {
        K <- kernel(reef, k)
        pick <- seq(1, m, by = 9)
        kern[[length(kern) + 1]] <- data.frame(reef = reef, taxon = k, row = rep(pick, times = length(pick)),
                                               col = rep(pick, each = length(pick)), value = num(c(K[pick, pick])))
    }
    for (reading in c("as_written", "calibrated")) {
        st <- project(reef, reading, T)
        for (t in 0:T) for (k in taxa) {
            n <- st[[t + 1]][[k]]
            runs[[length(runs) + 1]] <- data.frame(reading = reading, reef = reef, step = t, taxon = k,
                                                   colonies = num(h.mesh * sum(n)), log_size = num(h.mesh * sum(meshpts * n)),
                                                   area = num(h.mesh * sum((10^meshpts - 1) * n)))
        }
        for (t in c(4, 6, 11)) for (k in taxa)
            cells[[length(cells) + 1]] <- data.frame(reading = reading, reef = reef, step = t, taxon = k, cell = 1:m,
                                                     value = num(h.mesh * st[[t + 1]][[k]]))   # colonies per 10 m2 in each cell
        cat(reading, "LTER", reef, "step 4 colonies:", sapply(taxa, function(k) round(h.mesh * sum(st[[5]][[k]]), 2)), "\n")
    }
}
save(do.call(rbind, runs), "coral_runs.csv")
save(do.call(rbind, cells), "coral_state.csv")
save(do.call(rbind, kern), "coral_kernels.csv")

## ---------------------------------------------------------------- the updated runs, from the observed states
if (!is.na(observed_file)) {
    obs <- read.csv(observed_file)
    recruit_size$Acr <- c(1.2, 0.4)
    observed <- function(reef, year) lapply(setNames(taxa, taxa), function(k) {
        n <- rep(0, m)
        rows <- obs[obs$reef == reef & obs$year == year & obs$taxon == k, ]
        n[rows$cell] <- rows$colonies * 10 / rows$sampled_area_m2 / h.mesh   # density per unit log size, per 10 m2
        n
    })
    out <- list()
    for (reef in reefs) for (reading in c("as_written", "calibrated")) for (start in c(2013, 2015)) {
        steps <- if (start == 2013) 2 else 5
        st <- project(reef, reading, steps, observed(reef, start))
        for (t in 0:steps) for (k in taxa) {
            n <- st[[t + 1]][[k]]
            out[[length(out) + 1]] <- data.frame(reading = reading, reef = reef, start = start, step = t, taxon = k,
                                                 colonies = num(h.mesh * sum(n)), log_size = num(h.mesh * sum(meshpts * n)))
        }
    }
    save(do.call(rbind, out), "coral_observed_runs.csv")
}
cat("done\n")
