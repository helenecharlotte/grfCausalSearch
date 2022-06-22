## setwd("~/research/SoftWare/grfCausalSearch/")
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
server <- !inherits(try(setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch"),silent = TRUE),"try-error")
# ---------------------------------------------------------------------
# packages
# ---------------------------------------------------------------------
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
library(targets)
## library(future)
## library(future.callr)
## library(future.batchtools)
targets::tar_option_set(packages = thepackages)
# ---------------------------------------------------------------------
# R functions
# ---------------------------------------------------------------------
for(f in list.files("R",".R$",full.names=TRUE)){source(f)}
for(f in list.files("functions",".R$",full.names=TRUE)){source(f)}
# ---------------------------------------------------------------------
# Multicore via future
# ---------------------------------------------------------------------
## 
## tweak recources are ignored
## library(future)
## library(future.callr)
## library(batchtools)
## library(future.batchtools)
if (FALSE){
    slurm <- future::tweak(batchtools_slurm,
                           template = 'batchtools.slurm.tmpl')
    future::plan(slurm)
    tar_option_set(
        memory='transient',
        storage='worker',
        resources = tar_resources(
            future = tar_resources_future(
                plan = future::plan(slurm),
                resources = list(
                    template = 'batchtools.slurm.tmpl'
                    ## ,memory = 512,
                    ## ncpus = 1,
                    ## ntasks = 1,
                    ## walltime = 60L
                )
            )
        )
    )
}
## future::plan(batchtools_slurm, template = "slurm.tmpl")
## tar_option_set(
## memory='transient',
## storage='worker',
## resources = tar_resources(
## future = tar_resources_future(
## resources = list(template = 'slurm.tmpl',
## num_cores = 30,
## n_cores = 30),
## )
## )
## )
## }else{
## future::plan(callr)
## }
# ---------------------------------------------------------------------
# Simulation settings
# ---------------------------------------------------------------------
source("./setting/simulation_targets.R")
## source("./setting/dependent_targets.R")
## source("./setting/misspecified_targets.R")
## source("./setting/net_targets.R")

# ---------------------------------------------------------------------
# The target flow
# ---------------------------------------------------------------------

## MCCORES <- 5
MCCORES <- 50
list(tar_target(REPETITIONS, 1:1000),
     varying_target,
     fixed_target,
     fixed,
     truth_varying,
     truth,
     estimates,
     ate,
     results,
     ranking,
     boxplots)




