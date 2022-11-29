## setwd("~/research/SoftWare/grfCausalSearch/")
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
library(targets)
## server <- !inherits(try(setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch"),silent = TRUE),"try-error")
# ---------------------------------------------------------------------
# packages
# ---------------------------------------------------------------------
# library("targets");library("future.batchtools");library("tarchetypes");library("future.callr");library("future");library("grf");library("ranger");library("data.table");library("scales");library("riskRegression");library("prodlim");library("survival");library("foreach");library("parallel");library("grid");library("ggplot2");library("gridExtra");library("ggplotify");library("cowplot")
thepackages <- c("targets",
                 "future.batchtools",
                 "tarchetypes",
                 "future.callr",
                 "future",
                 "grf",
                 "ranger",
                 ## "nleqslv"
                 "data.table",
                 "scales",
                 "riskRegression",
                 "prodlim",
                 "survival",
                 "foreach",
                 "parallel",
                 "grid",
                 "ggplot2",
                 "gridExtra",
                 "ggplotify",
                 "cowplot")
targets::tar_option_set(packages = thepackages)
# ---------------------------------------------------------------------
# R functions
# ---------------------------------------------------------------------
for(f in list.files("R",".R$",full.names=TRUE)){source(f)}
for(f in list.files("functions",".R$",full.names=TRUE)){source(f)}
# ---------------------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------------------
source("./setting/simulation_targets.R")

# ---------------------------------------------------------------------
# The target flow
# ---------------------------------------------------------------------

## MCCORES <- 5
## MCCORES <- 50
## MCCORES are set in setting/simulation-targets.R
list(varying_target,
     fixed_target,
     fixed,
     truth_varying,
     truth,
     estimates,
     ate,
     results,
     ranking,
     plotframe)
     ## review_estimates,
     ## review_ate,
     ## review_results)




