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
library(future)
library(future.callr)
library(batchtools)
library(future.batchtools)
if (server){
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
}else{
    future::plan(callr)
}
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
list(tar_target(MCCORES, 50),
     tar_target(REPETITIONS, 1:1000),
     varying_target,
     fixed_target,
     fixed,
     truth_varying,
     truth,
     estimates,
     ate,
     results,
     boxplots)
     ## misspecified_varying_target,
     ## misspecified_fixed_target,
     ## misspecified_estimators_target,
     ## misspecified_fixed,
     ## misspecified_truth_varying,
     ## misspecified_truth,
     ## misspecified_estimates,
     ## misspecified_ate,
     ## misspecified_results,
     ## misspecified_boxplots)

## list(tar_target(REPETITIONS, 1:5),

     ## dependent_varying_target,
     ## dependent_fixed_target,
     ## dependent_estimators_target,
     ## dependent_fixed,
     ## dependent_truth_varying,
     ## dependent_truth,
     ## dependent_estimates,
     ## dependent_ate,
     ## dependent_results,
     ## dependent_boxplots,
     ## net_varying_target,
     ## net_fixed_target,
     ## net_estimators_target,
     ## net_fixed,
     ## net_truth_varying,
     ## net_truth,
     ## net_estimates,
     ## net_ate,
     ## net_results,
     ## net_boxplots
     ## )
