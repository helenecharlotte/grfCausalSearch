### run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (13:50) 
## Last-Updated: Jul 11 2022 (10:14) 
##           By: Thomas Alexander Gerds
##     Update #: 66
### Code:
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch"))
library(targets)
library(tarchetypes)
## library(batchtools)
## library(future.batchtools)
# run simulation
## tar_visnetwork(targets_only = TRUE)
## tar_outdated()
tar_make()
## tar_make(reporter = "verbose_positives")
## tar_make(reporter = "summary")
# nohup R CMD BATCH /maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/simulation/run-rao.R &

######################################################################
### run.R ends here
