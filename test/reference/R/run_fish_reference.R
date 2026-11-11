## Deterministic driver: patches a copy of fish-model-all-spp.R lines 1-596 (through end of mid loop)
suppressPackageStartupMessages({library(ggplot2); library(plyr); library(tidyverse)})
lambda <- function(A) { ev <- eigen(A)$values; Re(ev[which.max(Re(ev))]) }  # popbio::lambda equivalent
src <- readLines('code/fish-model-all-spp.R')[1:596]
src <- src[!grepl('^library\\(', src)]
variant <- Sys.getenv('VARIANT', 'fixed')
src <- sub("flowdata <- all.scenarios.list\\$natural.flow",
           "flowdata_rds <- all.scenarios.list$natural.flow; flowdata <- read.csv('data/naturalflow_long.csv')", src)
src <- sub("^iterations <- 10", if (variant == 'fixed') "iterations <- 1" else "iterations <- 3", src)
if (variant == 'fixed') {
  ## replace rnbinom(...) draws by the csv value biomXXX3
  for (i in grep('rnbinom', src)) {
    sp <- sub('.*biom([A-Z]{4})2.*', '\\1', src[i-1])
    src[i] <- sub('rnbinom\\(.*\\)\\)', paste0('biom', sp, '3)'), src[i])
  }
} else set.seed(1)
writeLines(src, paste0('patched_', variant, '.R'))
eval(parse(text = src))
cat('first 45 rows equal rds natural.flow: ', isTRUE(all.equal(flowdata[1:45,], flowdata_rds, check.attributes=FALSE)), '\n')
options(digits = 17)
if (variant == 'fixed') {
  out <- do.call(cbind, lapply(sppnames, function(nm) { m <- output.biom.list[[nm]]; colnames(m) <- paste0(nm, '_biom', 1:3); m }))
  out <- data.frame(year = flowdata$year, out, check.names = FALSE)
  write.csv(format(out, digits = 17), 'fish_det_biomass.csv', row.names = FALSE, quote = FALSE)
  outN <- do.call(cbind, lapply(sppnames, function(nm) { m <- output.N.list[[nm]]; colnames(m) <- paste0(nm, '_N', 1:3); m }))
  write.csv(format(data.frame(year = flowdata$year, outN, TotalN = Total.N[,1], check.names=FALSE), digits = 17), 'fish_det_N.csv', row.names = FALSE, quote = FALSE)
  write.csv(format(data.frame(year = flowdata$year, lambda.df), digits = 17), 'fish_det_lambda.csv', row.names = FALSE, quote = FALSE)
  write.csv(format(data.frame(year = flowdata$year, as.data.frame(fec.list)), digits = 17), 'fish_det_fec.csv', row.names = FALSE, quote = FALSE)
  print(head(out, 6)); print(tail(out, 2)); print(head(lambda.df, 4)); print(head(as.data.frame(fec.list), 4))
  cat('final biom vectors after year 54 update:\n'); for (nm in sppnames) print(c(get(paste0('biom', nm))))
} else {
  cat('seeded: stage-3 initial biomass per iteration (CACL..AMNA order of draws):\n')
  print(sapply(sppnames, function(nm) replist[[nm]][1,3,] / get(paste0('den', nm, '3'))))
  print(Total.N[1:5,])
}
