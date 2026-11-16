## Regenerate the reference for the designed flow regimes of Tonkin et al. (2021), src/flowregimes.jl.
##
## Runs the authors' scripts from GitHub jdtonkin/flow-tradeoffs (MIT): code/dod-3taxa-riparian.R,
## code/dod-3taxa-fish.R, code/dod-3taxa-inverts-scenarios.R and code/dod-3taxa-inverts-findopt.R. Each script
## selects its scenario by commented-out blocks; here the scenario is injected instead (the year types written
## below are those of the blocks for RIP_OPT, FSH_OPT, INV_OPT and NAT), the writing of files and the plots are
## dropped, and the script is otherwise evaluated as it stands. The fish script draws the initial adult biomass
## at random and has no seed; a seed is set here, and the draws are written out so that they can be replayed.
## Two packages the fish script loads are not needed for the outputs used: popbio (for the growth rate of each
## year's matrix, dropped) and mefa (a row repeat, replaced by base R).
##
## Usage: Rscript run_flow_tradeoffs_reference.R <clone of flow-tradeoffs> <output dir>
##   writes flow_riparian.csv, flow_fish_initial.csv, flow_fish_mean.csv, flow_inverts.csv (every 30th day and the
##   first three), flow_findopt.csv and flow_targets.csv to <output dir>. Takes about ten minutes.
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1])
outdir <- normalizePath(args[2], mustWork = FALSE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({library(plyr); library(tidyverse); library(gridExtra); library(RColorBrewer)})
setwd(repo)
num <- function(x) sprintf("%.17g", x)
save <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE, quote = FALSE)
script <- function(name) readLines(file.path(repo, "code", name))
cut_at <- function(src, pattern) src[seq_len(grep(pattern, src)[1] - 1)]
natural <- readRDS("data/all_scenarios_list.rds")$natural.flow
natural200 <- natural[rep(seq_len(nrow(natural)), 5), ][1:200, ]   # the script's mefa:::rep.data.frame(flowdata, 5)
scenarios <- c("RIP_OPT", "FSH_OPT", "INV_OPT", "NAT")

## ---------------------------------------------------------------- riparian plants
riparian_years <- function(fset) {
  if (fset == "RIP_OPT") {
    gap <- 5
    bigflood <- rep(c(rep(0, gap), 1), length.out = 200)
    drought <- rep(c(rep(0, gap - 1), 1, 0), length.out = 200)
    nonflood <- ifelse(bigflood + drought == 1, 0, 1)
  } else if (fset == "FSH_OPT") {
    bigflood <- rep(1, 200); drought <- rep(0, 200); nonflood <- rep(0, 200)
  } else if (fset == "INV_OPT") {
    bigflood <- rep(0, 200); drought <- rep(0, 200); nonflood <- rep(1, 200)
  } else {
    bigflood <- natural200$SP_highflood; drought <- natural200$drought
    nonflood <- ifelse(bigflood + drought == 1, 0, 1)
  }
  list(bigflood = bigflood, drought = drought, nonflood = nonflood)
}
rip <- script("dod-3taxa-riparian.R")
rip <- rip[!grepl("^rm\\(list|^library\\(", rip)]
first <- grep("Simulation choices begin here", rip); last <- grep("Simulation choices end here", rip)
rip <- c(rip[1:first], "fset <- .fset; bigflood <- .years$bigflood; drought <- .years$drought; nonflood <- .years$nonflood",
         rip[last:length(rip)])
rip <- cut_at(rip, "^# OUTPUT ---")
rip_out <- NULL
for (s in scenarios) {
  env <- new.env()
  assign(".fset", s, env); assign(".years", riparian_years(s), env)
  eval(parse(text = rip), envir = env)
  rip_out <- rbind(rip_out, data.frame(scenario = s, year = 1:200, C = num(env$Cnonseedling), T = num(env$Tnonseedling),
                                       W = num(env$Wnonseedling), S = num(env$Snonseedling), M = num(env$Mnonseedling)))
  cat("riparian", s, ": mean cottonwood %K (years 11-200) =", mean(env$Cnonseedling[11:200]), "\n")
}
save(rip_out, "flow_riparian.csv")

## ---------------------------------------------------------------- fish
fish_years <- function(fset) {
  f <- data.frame(SP_highflood = rep(0, 200), SU_highflood = 0, medflood = 0, drought = 0, nonevent = 0)
  if (fset == "RIP_OPT") {
    f$SP_highflood <- rep(c(rep(0, 5), 1), length.out = 200)
    f$drought <- rep(c(rep(0, 4), 1, 0), length.out = 200)
    f$nonevent <- ifelse(f$SP_highflood == 1 | f$drought == 1, 0, 1)
  } else if (fset == "FSH_OPT") {
    f$SP_highflood <- 1
  } else if (fset == "INV_OPT") {
    f$nonevent <- 1
  } else {
    f <- natural200[, c("SP_highflood", "SU_highflood", "medflood", "drought", "nonevent")]
  }
  f
}
fish <- script("dod-3taxa-fish.R")
fish <- fish[!grepl("^rm\\(list|^library\\(|lambda\\.df\\[i,\\] <- ", fish)]
at <- grep("^count <- 200", fish)
fish <- c(fish[1:(at - 1)], "flowdata <- .flow", fish[at:length(fish)])
at <- grep("^for\\(iter in 1:iterations\\)", fish)
fish <- c(fish[1:(at - 1)], "set.seed(.seed)", fish[at:length(fish)])
fish <- cut_at(fish, "^## \\* OUTPUTS")
species <- c("CACL", "GIRO", "LECY", "CAIN", "MIDO", "CYLU", "AMNA")
fish_init <- NULL; fish_mean <- NULL
for (k in seq_along(scenarios)) {
  s <- scenarios[k]
  env <- new.env()
  assign(".flow", fish_years(s), env); assign(".seed", 2020 + k, env)
  eval(parse(text = fish), envir = env)
  rb <- env$replist_biom
  init <- sapply(species, function(sp) rb[[sp]][1, 3, ])
  fish_init <- rbind(fish_init, data.frame(scenario = s, iteration = 1:100, init))
  means <- sapply(species, function(sp) rowMeans(rb[[sp]][, 2, ] + rb[[sp]][, 3, ]))
  fish_mean <- rbind(fish_mean, data.frame(scenario = s, year = 1:200, apply(means, 2, num)))
  native <- 100 * (means[, "CACL"] + means[, "CAIN"] + means[, "GIRO"]) / env$K
  cat("fish", s, ": mean native %K (years 11-200) =", mean(native[11:200]), "\n")
}
save(fish_init, "flow_fish_initial.csv")
save(fish_mean, "flow_fish_mean.csv")

## ---------------------------------------------------------------- invertebrates, the four scenarios
inv <- script("dod-3taxa-inverts-scenarios.R")
inv <- inv[!grepl("^library\\(", inv)]
inv <- sub('^natflow_df <- read.csv\\("export/natflow_mean.csv"\\)',
           "natflow_df <- data.frame(t = c(1:200, 1:165), Q = c(rep(1134, 200), rep(459, 165)))", inv)
inv <- cut_at(inv, "^## exporting")
env <- new.env()
eval(parse(text = inv), envir = env)
fdf <- env$fdf_combined
keep <- fdf$Time %% 30 == 0 | fdf$Time <= 3          # every 30th day, to keep the file small
save(data.frame(taxa = fdf$taxa, scenario = fdf$scenario, Time = fdf$Time, Kpre = num(fdf$Kpre), Nt = num(fdf$Nt))[keep, ],
     "flow_inverts.csv")
inv_target <- fdf %>% filter(Time > 3650, taxa != "Ostracod") %>% group_by(scenario, Time) %>%
  summarise(Nt = mean(Nt), .groups = "drop") %>% group_by(scenario) %>% summarise(mean = mean(Nt))
print(inv_target)

## ---------------------------------------------------------------- invertebrates, the search over flood frequency
opt <- script("dod-3taxa-inverts-findopt.R")
opt <- opt[!grepl("^library\\(", opt)]
at <- grep("summarise\\(ari\\.mean", opt)
for (j in rev(at)) opt <- append(opt, ".found[[length(.found) + 1]] <- means", after = j)
env <- new.env()
assign(".found", list(), env)
eval(parse(text = opt), envir = env)
found <- NULL
for (j in seq_along(env$.found)) {
  m <- env$.found[[j]]
  found <- rbind(found, data.frame(taxa = c("Dragonfly", "Mayfly", "Ostracod")[j], frequency = m$.id,
                                   ari_mean = num(m$ari.mean), geo_mean = num(m$geo.mean)))
}
save(found, "flow_findopt.csv")
print(found)

## ---------------------------------------------------------------- the targets compared in the paper (compile-scenarios.R)
targets <- data.frame(scenario = scenarios,
  plants = sapply(scenarios, function(s) mean(as.numeric(rip_out$C[rip_out$scenario == s])[11:200])),
  fish = sapply(scenarios, function(s) {
    m <- fish_mean[fish_mean$scenario == s, ]
    mean(100 * (as.numeric(m$CACL) + as.numeric(m$CAIN) + as.numeric(m$GIRO))[11:200] / 47660)
  }),
  invertebrates = inv_target$mean[match(scenarios, inv_target$scenario)])
targets$plants <- num(targets$plants); targets$fish <- num(targets$fish); targets$invertebrates <- num(targets$invertebrates)
save(targets, "flow_targets.csv")
print(targets)
