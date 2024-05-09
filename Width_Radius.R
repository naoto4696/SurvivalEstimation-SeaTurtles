
# Simulation parameters
realsteps <- 10000000 # number of steps
burnin <- realsteps/100
nsteps <- realsteps + burnin

# Data hundle
ID <- as.vector(as.ts(Width_Radius$ID))
Width_dummy <- as.vector(as.ts(Width_Radius$Width)) 
Radius_dummy <- as.vector(as.ts(Width_Radius$Radius)) 
Width <- Width_dummy[is.na(Width_dummy)==F]
Radius <- Radius_dummy[is.na(Width_dummy)==F]
Len <- length(Width)

# Mathematical model
Radius_Width <- function(Hm,Gh,H) {
  Width <- (1-exp(-Gh))*(Hm-H)
  return(Width)
}

# Likelihood Function
LikelihoodFunction <- function(Hm,Gh,sigma) {
  Vec <- dnorm(Width,Radius_Width(Hm,Gh,Radius),sigma) 
  LikelihoodFunction <- prod(Vec)
  if(is.nan(LikelihoodFunction)) return(0)
  else return(LikelihoodFunction)
}  

# Prior distribution
H_max <- max(Radius)
Pi_Hm <- function(Hm){
  if(Hm<H_max) return(0)
  else return(1)
} 
Pi_Gh <- function(Gh){
  if(Gh<0) return(0)
  else return(1)
} 
Pi_sigma <- function(sigma){
  if(sigma<0) return(0)
  else return(1)
} 

# Posterior distribution
joint　<-　function(Hm,Gh,sigma) LikelihoodFunction(Hm,Gh,sigma)*Pi_Hm(Hm)*Pi_Gh(Gh)*Pi_sigma(sigma)

# MCMC function
Metro<-function(c1,c2,c3) {
  p1 <- c1 + rnorm(1,0,1)
  p2 <- c2 + rnorm(1,0,0.005)
  p3 <- c3 + rnorm(1,0,0.1)
  odds <- joint(p1,p2,p3)/joint(c1,c2,c3)
  pmove <- min(odds,1)
  if(pmove >= runif(1,0,1)) 
    c(p1,p2,p3)
  else
    c(c1,c2,c3)
}

# Init Preparation
for(i in 1:1000000) {
  Init_Hm <- runif(1,H_max,50)
  Init_Gh <- runif(1,0,0.5)
  Init_sigma <- runif(1,0,10)
  if(joint(Init_Hm,Init_Gh,Init_sigma) > 0){
      Init <- c(Init_Hm,Init_Gh,Init_sigma)
      break
  }
}

mcmc_rw <- matrix(0,ncol=3,nrow=nsteps+1)
mcmc_rw[1,1] <- Init[1]
mcmc_rw[1,2] <- Init[2]
mcmc_rw[1,3] <- Init[3]
print("Init Complete!")

# Prepare progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = nsteps 
)

#MCMC
print("Calculating MCMC...")
for (i in 1:nsteps) {
  dummy <- Metro(mcmc_rw[i,1],mcmc_rw[i,2],mcmc_rw[i,3])
  mcmc_rw[i+1,1] <- dummy[1]
  mcmc_rw[i+1,2] <- dummy[2] 
  mcmc_rw[i+1,3] <- dummy[3] 
  pb$tick()
}

# Density
Dens1 <- density(mcmc_rw[-(1:burnin),1])
Dens2 <- density(mcmc_rw[-(1:burnin),2])
Dens3 <- density(mcmc_rw[-(1:burnin),3])

# Mean
M_Hm <- mean(mcmc_rw[-(1:burnin),1])
M_Gh <- mean(mcmc_rw[-(1:burnin),2])
M_sigma <- mean(mcmc_rw[-(1:burnin),3])

# 95%HDI
HDIwidth <- 0.95
LowHDI <- realsteps*(1-HDIwidth)/2
HighHDI <- realsteps-LowHDI
LowHDI_Hm <- hdi(Dens1,ci=HDIwidth)[1]
HighHDI_Hm <- hdi(Dens1,ci=HDIwidth)[2]
LowHDI_Gh <- hdi(Dens2,ci=HDIwidth)[1]
HighHDI_Gh <- hdi(Dens2,ci=HDIwidth)[2]
LowHDI_sigma <- hdi(Dens3,ci=HDIwidth)[1]
HighHDI_sigma <- hdi(Dens3,ci=HDIwidth)[2]

# Plot
dir.create("~/Desktop/Result_Width_Radius")
setwd("~/Desktop/Result_Width_Radius")
print("Plotting graphs...")

Rseq <- seq(0,25,0.001)
tseq <- seq(0,60,0.01)
cplotseq <- seq(burnin,nsteps,by=realsteps/1000)

#fit plot
pdf("fitplot_Width_Radius.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(Radius, Width, type="p", xlab="", ylab="", 
     xlim=c(min(Radius),max(Radius)), ylim=c(min(Width),max(Width)), 
     xaxt="n", yaxt="n", pch=20,
     panel.first=grid(col="grey", lty=2))
for(i in cplotseq){
  par(new=T)
  plot(Rseq, Radius_Width(mcmc_rw[i,1],mcmc_rw[i,2],Rseq),type="l", xlab="", ylab=""
       ,xlim=c(min(Radius),max(Radius)), ylim=c(min(Width),max(Width)), 
       xaxt="n", yaxt="n", col=rgb(1, 0, 0, alpha=0.003))
}
axis(1, at = pretty(Radius), labels = pretty(Radius))
mtext("Radius(mm)", side = 1, line =3)
axis(2, at = pretty(Width), labels = pretty(Width))
mtext("Width(mm)", side = 2, line = 3)
dev.off()

# Posterior plot
pdf("Posteriorplot_Radius_Width.pdf", width = 7, height = 6)
par(mfrow=c(2,2))

plot(0,0,xlim=c(min(Dens1$x),max(Dens1$x)),ylim=c(min(Dens1$y),max(Dens1$y)),xlab="",ylab="",main="Maximum Humerus Radius")
lines(Dens1,lwd=1,col=1)
abline(v=M_Hm,lty=1,col=2,lwd=1)
abline(v=LowHDI_Hm,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Hm,lty=3,col=2,lwd=1.5)
abline(v=Init_Hm,lty=3,col=3,lwd=1.5)
axis(1, at = pretty(Dens1$x), labels = pretty(Dens1$x))
mtext("Hm", side = 1, line =3)
axis(2, at = pretty(Dens1$y), labels = pretty(Dens1$y))
mtext("Density", side = 2, line = 3)

plot(0,0,xlim=c(min(Dens2$x),max(Dens2$x)),ylim=c(min(Dens2$y),max(Dens2$y)),
     xaxt="n", yaxt="n", xlab="",ylab="",main="Humerus Growth Coefficient")
lines(Dens2,lwd=1,col=1)
abline(v=M_Gh,lty=1,col=2,lwd=1)
abline(v=LowHDI_Gh,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Gh,lty=3,col=2,lwd=1.5)
abline(v=Init_Gh,lty=3,col=3,lwd=1.5)
axis(1, at = pretty(Dens2$x), labels = pretty(Dens2$x))
mtext("Gh", side = 1, line =3)
axis(2, at = pretty(Dens2$y), labels = pretty(Dens2$y))
mtext("Density", side = 2, line = 3)

plot(0,0,xlim=c(min(Dens3$x),max(Dens3$x)),ylim=c(min(Dens3$y),max(Dens3$y)),
     xaxt="n", yaxt="n", xlab="",ylab="",main="Standard Deviation")
lines(Dens3,lwd=1,col=1)
abline(v=M_sigma,lty=1,col=2,lwd=1)
abline(v=LowHDI_sigma,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_sigma,lty=3,col=2,lwd=1.5)
abline(v=Init_sigma,lty=3,col=3,lwd=3)
axis(1, at = pretty(Dens3$x), labels = pretty(Dens3$x))
mtext("sigma", side = 1, line =3)
axis(2, at = pretty(Dens3$y), labels = pretty(Dens3$y))
mtext("Density", side = 2, line = 3)

dev.off()

# Mcmc plot
pdf("mcmcplot_Width_Radius.pdf", width = 8, height = 12)
mcmcplot <- seq(burnin,nsteps,length.out=10000)
par(mfrow=c(3,1))
plot(mcmc_rw[mcmcplot,1],xlab="steps",ylab="",main="Maximum Humerus Radius",type="l")
plot(mcmc_rw[mcmcplot,2],xlab="steps",ylab="",main="Growth Coefficient",type="l")
plot(mcmc_rw[mcmcplot,3],xlab="steps",ylab="",main="Sigma",type="l")

dev.off()   

# list
StatisticalResults <- matrix(ncol=4,nrow=3)
colnames(StatisticalResults) <- c("Name","Low_HDI","Mean","High_HDI")
StatisticalResults[,1] <- c("Maximum humerus Radius","Growth coefficient","Standard deviation")
StatisticalResults[1,-1] <- round(c(LowHDI_Hm,M_Hm,HighHDI_Hm),digits=5)
StatisticalResults[2,-1] <- round(c(LowHDI_Gh,M_Gh,HighHDI_Gh),digits=5)
StatisticalResults[3,-1] <- round(c(LowHDI_sigma,M_sigma,HighHDI_sigma),digits=5)

write.csv(as.data.frame(StatisticalResults),"StatisticalResults_Width_Radius.csv")

# Preparation

Lm_chain <- mcmc_rs[-(1:burnin),1]*mcmc_rw[-(1:burnin),1]
Pri_M_Lm <- mean(Lm_chain)
Pri_S_Lm <- sqrt(var(Lm_chain))
Pri_M_Gc <- M_Gh
Pri_S_Gc <- (sqrt(var(mcmc_rw[-(1:burnin),2])))
