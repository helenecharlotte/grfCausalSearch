### shower.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 29 2020 (08:16) 
## Version: 
## Last-Updated: Jun 29 2020 (11:24) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


library(ggplot2)
library(data.table)
a <- readRDS("~/research/SoftWare/grfCausalSearch/simulation-results/Sim1c.rds")
## forgot to provide cens input parameter
a[,censpar:=rep(rep(c(0.5,.7,1.5),rep(60,3)),9)]
## 9 * 3 *  6
## n, C, E, 

a[Effect.A2==1.5& n==200,.(n,cens,Effect.A2)]
a[,CensPercent:=round(100*cens)]

cov <- ggplot(a[intervene=="A2"],aes(x=n,y=coverage))+geom_line()+facet_grid(Effect.A2~censpar)
cov+ylim(c(0.85,1))+geom_abline(intercept=.95,slope=0,colour="red")

dev.new()

b <- readRDS("~/research/SoftWare/grfCausalSearch/simulation-results/results-net-weight-formula.rds")

covb <- ggplot(b[intervene=="A2"],aes(x=n,y=coverage))+geom_line()+facet_grid(Effect.A2~censpar)
covb+ylim(c(0.85,1))+geom_abline(intercept=.95,slope=0,colour="red")

######################################################################
### shower.R ends here
