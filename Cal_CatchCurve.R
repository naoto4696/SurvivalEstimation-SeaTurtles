#Please install package "HDInterval" and "progress"
library(HDInterval)
library(progress)

#Please delete # if you have "SCL_Radius" form data.
SCL_Radius <- read.csv("~/Desktop/Program/Sample_SR.csv")
source("~/Desktop/Program/SCL_Radius.R")

#Please delete # if you have "Width_Radius" form data.
Width_Radius <- read.csv("~/Desktop/Program/Sample_WR.csv")
source("~/Desktop/Program/Width_Radius.R")

#If you do not have above data, use below parameters as representative values of Japanese loggerhead population
#Pri_M_Lm <- 103.2442
#Pri_S_Lm <- 4.556656
#Pri_M_Gc <- 0.03972297
#Pri_S_Gc  <- 0.004063969

#Calculate parameters using your data
Realdata <- read.csv("~/Desktop/Program/Sample_SCL.csv")
#Realdata <- read.csv("yourdatapass")
source("~/Desktop/Program/SCLdistribution.R")