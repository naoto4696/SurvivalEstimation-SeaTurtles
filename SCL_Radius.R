#simuparameter
realsteps <- 10000000 # number of steps
burnin <- realsteps/100
nsteps <- realsteps + burnin

# Data hundle
SCL <- as.vector(as.ts(SCL_Radius$SCL/10)) # Convert mm to cm
Radius <- as.vector(as.ts(SCL_Radius$Radius/2)) # Convert diameter to radius
SCL <- SCL[is.na(SCL)==F]
Radius <- Radius[is.na(Radius)==F]
Len <- length(SCL)

# Mathematical model
Radius_SCL <- function(C,R) {
  SCL <- C*R
  return(SCL)
}

# Likelihood Function
LikelihoodFunction <- function(C,sigma) {
  LikelihoodFunction <- 1
  Vec <- dnorm(SCL,Radius_SCL(C,Radius),sigma) 
  LikelihoodFunction <- prod(Vec)
  return(LikelihoodFunction)
}  

# Prior distribution
Pi_C <- function(C){
  if(C<0) return(0)
  else return(1)
} 
Pi_sigma <- function(sigma){
  if(sigma<0) return(0)
  else return(1)
} 

# Posterior distribution
joint　<-　function(C,sigma) LikelihoodFunction(C,sigma)*Pi_C(C)*Pi_sigma(sigma)

# MCMC function
Metro<-function(c1,c2) {
  p1 <- c1 + rnorm(1,0,1)
  p2 <- c2 + rnorm(1,0,0.1)
  odds <- joint(p1,p2)/joint(c1,c2)
  pmove <- min(odds,1)
  if(pmove >= runif(1,0,1)) 
    c(p1,p2)
  else
    c(c1,c2)
}

# Init Preparation
for(i in 1:1000000) {
  Init_C <- runif(1,0,5) 
  Init_sigma <- runif(1,0,10)
  if(joint(Init_C,Init_sigma) > 0){
    Init <- c(Init_C,Init_sigma)
    break
  }
}
mcmc_rs <- matrix(0,ncol=3,nrow=nsteps+1)
mcmc_rs[1,1] <- Init[1]
mcmc_rs[1,2] <- Init[2]
print("Init Complete!")

# Prepare progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = nsteps 
)

#MCMC
print("Calculating MCMC...")
for (i in 1:nsteps) {
  dummy <- Metro(mcmc_rs[i,1],mcmc_rs[i,2])
  mcmc_rs[i+1,1] <- dummy[1]
  mcmc_rs[i+1,2] <- dummy[2] 
  pb$tick()
}

# Density
Dens1 <- density(mcmc_rs[-(1:burnin),1])
Dens2 <- density(mcmc_rs[-(1:burnin),2])

# Mean
M_C <- mean(mcmc_rs[-(1:burnin),1])
M_sigma <- mean(mcmc_rs[-(1:burnin),2])

# 95%HDI
HDISCL <- 0.95
LowHDI <- realsteps*(1-HDISCL)/2
HighHDI <- realsteps-LowHDI
LowHDI_C <- hdi(Dens1,ci=HDISCL)[1]
HighHDI_C <- hdi(Dens1,ci=HDISCL)[2]
LowHDI_sigma <- hdi(Dens2,ci=HDISCL)[1]
HighHDI_sigma <- hdi(Dens2,ci=HDISCL)[2]

# Plot
dir.create("~/Desktop/Result_SCL_Radius")
setwd("~/Desktop/Result_SCL_Radius")
print("Plotting graphs...")

Rseq <- seq(0,20,0.001)
tseq <- seq(0,100,0.01)
cplotseq <- seq(burnin,nsteps,by=realsteps/1000)

#fit plot
pdf("Fitplot_SCL_Radius.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(Radius,SCL, type="p", xlab="", ylab="",
     xlim=c(min(Radius),max(Radius)), ylim=c(min(SCL),max(SCL)),
     xaxt="n",yaxt="n", pch=20,
     panel.first=grid(col="grey", lty=2))
for(i in cplotseq){
  par(new=T)
  plot(Rseq,Radius_SCL(mcmc_rs[i,1],Rseq), type="l", xlab="", ylab="",
       ylim=c(min(SCL),max(SCL)), xlim=c(min(Radius), max(Radius)), 
       xaxt="n", yaxt="n", col=rgb(1, 0, 0, alpha=0.003))
}
axis(1, at = pretty(Radius), labels = pretty(Radius))
mtext("Radius(mm)", side = 1, line =3)
axis(2, at = pretty(SCL), labels = pretty(SCL))
mtext("SCL(cm)", side = 2, line = 3)
dev.off()

# Posterior plot
pdf("Posteriorplot_SCL_Radius.pdf", width =7 , height = 3)
par(mfrow=c(1,2))

plot(-100,-100,xlim=c(min(Dens1$x),max(Dens1$x)),ylim=c(min(Dens1$y),max(Dens1$y)),
     xlab="", ylab="", xaxt="n", yaxt="n", main="Regression Coefficient")
lines(Dens1,lwd=1,col=1)
abline(v=M_C,lty=1,col=2,lwd=1)
abline(v=LowHDI_C,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_C,lty=3,col=2,lwd=1.5)
abline(v=Init_C,lty=3,col=3,lwd=1.5)
axis(1, at = pretty(Dens1$x), labels = pretty(Dens1$x))
mtext("C", side = 1, line =3)
axis(2, at = pretty(Dens1$y), labels = pretty(Dens1$y))
mtext("Density", side = 2, line = 3)

plot(-100,-100,xlim=c(min(Dens2$x),max(Dens2$x)),ylim=c(min(Dens2$y),max(Dens2$y)),
     xlab="",ylab="", xaxt="n", yaxt="n", main="Standard Deviation")
lines(Dens2,lwd=1,col=1)
abline(v=M_sigma,lty=1,col=2,lwd=1)
abline(v=LowHDI_sigma,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_sigma,lty=3,col=2,lwd=1.5)
abline(v=Init_sigma,lty=3,col=3,lwd=1.5)
axis(1, at = pretty(Dens2$x), labels = pretty(Dens2$x))
mtext("sigma'", side = 1, line =3)
axis(2, at = pretty(Dens2$y), labels = pretty(Dens2$y))
mtext("Density", side = 2, line = 3)

dev.off()

# Mcmc plot
pdf("mcmcplot_SCL_Radius.pdf", width = 8, height = 6)
mcmcplot <- seq(burnin,nsteps,length.out=10000)
par(mfrow=c(2,1))
plot(mcmc_rs[mcmcplot,1],xlab="Regression Coefficient",ylab="",main="Regression Coefficient",type="l")
plot(mcmc_rs[mcmcplot,2],xlab="sigma",ylab="",main="Sigma",type="l")

dev.off()   

# list
StatisticalResults <- matrix(ncol=4,nrow=2)
colnames(StatisticalResults) <- c("Name","Low_HDI","Mean","High_HDI")
StatisticalResults[,1] <- c("Regression coefficient", "Standard deviation")
StatisticalResults[1,-1] <- round(c(LowHDI_C,M_C,HighHDI_C),digits=5)
StatisticalResults[2,-1] <- round(c(LowHDI_sigma,M_sigma,HighHDI_sigma),digits=5)

write.csv(as.data.frame(StatisticalResults),"StatisticalResults_SCL_Radius.csv")