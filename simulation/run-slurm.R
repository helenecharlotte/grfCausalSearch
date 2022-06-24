### run-slurm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 22 2022 (09:17) 
## Version: 
## Last-Updated: May 27 2022 (16:32) 
##           By: Thomas Alexander Gerds
##     Update #: 15
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch")
library(targets)
library(tarchetypes)
library(batchtools)
library(future.batchtools)
makeClusterFunctionsSlurm("batchtools.slurm.tmpl")
loadRegistry(file.dir = "/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/registry",
             work.dir = "/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/registry",
             writeable = TRUE)
library(future)
library(future.callr)
tar_make_future(workers = 24,reporter = "summary")
## tar_make_future(workers = 24,reporter = "verbose_positives")
## nohup R CMD BATCH /maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/simulation/run-slurm.R &

######################################################################
### run-slurm.R ends here
