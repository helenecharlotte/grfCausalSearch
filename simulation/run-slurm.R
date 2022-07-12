### run-slurm.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 22 2022 (09:17) 
## Version: 
## Last-Updated: Jun 29 2022 (10:39) 
##           By: Thomas Alexander Gerds
##     Update #: 16
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
reg = makeRegistry(file.dir = "/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/registry",
                   work.dir = "/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch/registry")
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
