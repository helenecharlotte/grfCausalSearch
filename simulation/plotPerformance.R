### plotPerformance.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 19 2022 (06:05) 
## Version: 
## Last-Updated: Apr 19 2022 (11:24) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(ggplot2)
net <- tar_read(perf_table_net)
crude <- tar_read(perf_table_crude)
crude[,A1_T1 := factor(A1_T1)]
net[,A1_T1 := factor(A1_T1)]

g <- ggplot(crude[intervene == "A1"],aes(x = n, y = bias,color = A1_T1)) + geom_point() + geom_line() + facet_grid(~censored.tau)
g

g <- ggplot(crude[intervene == "A1"],aes(x = n, y = abs.bias,color = A1_T1)) + geom_point() + geom_line() + facet_grid(~censored.tau)
g

g <- ggplot(net[intervene == "A1"&num.tree == 2000],aes(x = n, y = coverage,color = A1_T1)) + geom_point() + geom_line() + facet_grid(~censored.tau)
g

######################################################################
### plotPerformance.R ends here
