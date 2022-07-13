### run.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 15 2022 (13:50) 
## Last-Updated: Jul 13 2022 (10:35) 
##           By: Thomas Alexander Gerds
##     Update #: 66
### Code:
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch"))
library(targets)
library(tarchetypes)
## sapply(thepackages,function(x){print(x);do.call("library",list(x))})

## Sys.setenv(DEBUGME = "batchtools")
library(batchtools)
library(future.batchtools)
# run simulation
tar_visnetwork(targets_only = TRUE)
tar_outdated()
x <- tar_manifest()
setDT(x)
x[grepl("INDEPENDENT_ESTIMATE",name)]
## tar_outdated()
tar_make()
tar_make(reporter = "verbose_positives")
tar_make(reporter = "summary")
tar_make_future(workers = 6,reporter = "summary")
tar_make_future(workers = 6,reporter = "verbose_positives")
tar_make(reporter = "summary")

tar_delete(starts_with("SIM_DATA_1_1_1_1.25_0.025_5000"))

# cleaning up
# --------------------------------------------------------------------------
## tar_destroy
tar_delete(starts_with("SIM_DATA"))
tar_delete(starts_with("IND_SIM_DATA"))
tar_delete(starts_with("MISSPECIFIED_SIM_DATA"))
tar_delete(starts_with("DEPENDENT_SIM_DATA"))
tar_delete(starts_with("NET_SIM_DATA"))
tar_prune()


## tar_delete(starts_with("perf"))
## tar_delete(starts_with("estimate"))
## tar_delete(starts_with("sim"))
## tar_delete(starts_with("the"))
## tar_delete(starts_with("truth"))
## tar_delete(starts_with("repetitions"))
## tar_delete(starts_with("fixed"))
## tar_prune()


######################################################################
### run.R ends here
