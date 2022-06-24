### run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (13:50) 
## Last-Updated: Jun 21 2022 (17:10) 
##           By: Thomas Alexander Gerds
##     Update #: 65
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


######################################################################
### run.R ends here
