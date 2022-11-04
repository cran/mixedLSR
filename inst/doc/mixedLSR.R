## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mixedLSR)
set.seed(1)

## ----simulate-----------------------------------------------------------------
sim <- simulate_lsr(N = 100, k = 2, p = 30, m = 35)

## ----compute------------------------------------------------------------------
model <- mixed_lsr(sim$x, sim$y, k = 2, alt_iter = 1, anneal_iter = 10, em_iter = 10, verbose = TRUE)

## ----cluster_perf-------------------------------------------------------------
table(sim$true, model$assign)
ari <- mclust::adjustedRandIndex(sim$true, model$assign)
print(paste("ARI:",ari))

## ----visualize----------------------------------------------------------------
plot_lsr(model$a)
plot_lsr(sim$a)

## ----Reproducibility, collapse=TRUE-------------------------------------------
sessionInfo()

