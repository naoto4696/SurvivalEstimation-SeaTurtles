
# Simulation parameter
count <- 0
realsteps <- 1000000 # number of steps
burnin <- realsteps/100
nsteps <- realsteps + burnin
skipsteps <- 10 #Skip Mcmc steps for move Populaiton coefficient and Survival rate
Nf_th <- 100000  #Threshold of population coefficient

# HMC parameters
elm <- 0.1**(3)
egc <- 0.1**(5)
egv <- 0.1**(4)
ecv <- 0.1**(3)
emr <- 0.1**(5)
enf <- 0.1**(-1)
taulm <- 3
taugc <- 3
taugv <- 3
taucv <- 3
taumr <- 3
taunf <- 3
L <- 100


# Real Data handle
SCL <- as.vector(as.ts(Realdata$SCL))
Number <- as.vector(as.ts(Realdata$Number))
for(i in 1:length(Number)){
  if(Number[i]>0){
    minorder <- i
    break
  }
}
for(i in length(Number):1){
  if(Number[i]>0){
    Maxorder <- i
    break
  }
}
SCL <- SCL[minorder:Maxorder]
Number <- Number[minorder:Maxorder]
Len <- length(SCL)
iter <- max(SCL) - min(SCL) + 1

# Mathematical model
Number_t <- function(Mr,Nf,t) {
  Nt <- Nf*exp(-Mr*t)
}

Berta_curve <- function(Lm,Gc,L) {
  Berta <- (log(Lm/(Lm-L)) )/Gc
} 

Regression_curve <- function(Gv,Cv,t) {
  Rc <- 1/(1+exp(-Gv*(t-Cv)) )
}

Modelfunc <- function(Lm,Gc,Gv,Cv,Mr,Nf,L){
  RL <- Nf*(1-L/Lm)**(Mr/Gc) / (1+exp(Gv*Cv)*(1-L/Lm)**(Gv/Gc)) 
  return(RL)
}

# HMC denominators
HMCfunc <- function(Lm,Gc,Gv,Cv,Mr,Nf) {
  
  theta <- (Nf*(1-SCL/Lm)**(Mr/Gc))/(1+exp(Gv*Cv)*(1-SCL/Lm)**(Gv/Gc)) 
  denom <- 1/(1+exp(Gv*Cv)*(1-SCL/Lm)**(Gv/Gc)) 
  h <- (theta - Number*log(theta)) + (Lm-Pri_M_Lm)**2/(2*Pri_S_Lm**2) +
       (Gc-Pri_M_Gc)**2/(2*Pri_S_Gc**2)
  dhdlm <- ( (Mr/Gc)*(SCL/(Lm*(Lm-SCL))) - 
             (Gv/Gc)*(SCL/Lm**2)*((Lm-SCL)/Lm)**(Gv/Gc-1)*exp(Gv*Cv)*denom )*(theta-Number) +
             (Lm-Pri_M_Lm)/Pri_S_Lm**2
  dhdgc <- ( Gc**(-2)*log((Lm-SCL)/Lm))*(-Mr + Gv*exp(Gv*Cv)*((Lm-SCL)/Lm)**(Gv/Gc)*denom)*(theta-Number) +
            (Gc-Pri_M_Gc)/Pri_S_Gc**2
  dhdgv <- ( -(Cv + log((Lm-SCL)/Lm)/Gc)*exp(Gv*Cv)*((Lm-SCL)/Lm)**(Gv/Gc)*denom )*(theta-Number) 
  dhdcv <- ( -Gv*exp(Gv*Cv)*((Lm-SCL)/Lm)**(Gv/Gc)*denom )*(theta-Number)
  dhdmr <- ( log((Lm-SCL)/Lm)/Gc )*(theta-Number)
  dhdnf <- ( Nf**(-1) )*(theta-Number)
  
  th <- sum(h)
  tdhdlm <- sum(dhdlm)
  tdhdgc <- sum(dhdgc)
  tdhdgv <- sum(dhdgv)
  tdhdcv <- sum(dhdcv)
  tdhdmr <- sum(dhdmr)
  tdhdnf <- sum(dhdnf)
  
  return(c(th,tdhdlm,tdhdgc,tdhdgv,tdhdcv,tdhdmr,tdhdnf))
}

HMCfunc_Mr_Nf <- function(Lm,Gc,Gv,Cv,Mr,Nf) {
  
  theta <- (Nf*(1-SCL/Lm)**(Mr/Gc))/(1+exp(Gv*Cv)*(1-SCL/Lm)**(Gv/Gc)) 
  denom <- 1/(1+exp(Gv*Cv)*(1-SCL/Lm)**(Gv/Gc)) 
  h <- (theta - Number*log(theta)) + (Lm-Pri_M_Lm)**2/(2*Pri_S_Lm**2) + (Gc-Pri_M_Gc)**2/(2*Pri_S_Gc**2)
  dhdmr <- ( log((Lm-SCL)/Lm)/Gc )*(theta-Number)
  dhdnf <- ( Nf**(-1) )*(theta-Number)
  
  th <- sum(h)
  tdhdmr <- sum(dhdmr)
  tdhdnf <- sum(dhdnf)
  
  return(c(th,tdhdmr,tdhdnf))
}

# HMC function
HMC <- function(Lm,Gc,Gv,Cv,Mr,Nf) {
  dLm <- Lm
  dGc <- Gc
  dGv <- Gv
  dCv <- Cv
  dMr <- Mr
  dNf <- Nf
  plm <- rnorm(1,0,taulm)
  pgc <- rnorm(1,0,taugc)
  pgv <- rnorm(1,0,taugv)
  pcv <- rnorm(1,0,taucv)
  pmr <- rnorm(1,0,taumr)
  pnf <- rnorm(1,0,taunf)

  Vect <- HMCfunc(dLm,dGc,dGv,dCv,dMr,dNf)
  HamInit <- Vect[1] + plm**2/2 + pgc**2/2 + pgv**2/2 + pcv**2/2 + pmr**2/2 + pnf**2/2
  
  plm <- plm - elm*Vect[2]/2
  pgc <- pgc - egc*Vect[3]/2
  pgv <- pgv - egv*Vect[4]/2
  pcv <- pcv - ecv*Vect[5]/2
  pmr <- pmr - emr*Vect[6]/2
  pnf <- pnf - enf*Vect[7]/2
  
  for(i in 1:L-1){
    dLm <- dLm + elm*plm
    dGc <- dGc + egc*pgc
    dGv <- dGv + egv*pgv
    dCv <- dCv + ecv*pcv
    dMr <- dMr + emr*pmr
    dNf <- dNf + enf*pnf
    
    if(dLm <= max(SCL)) {
      return(c(Lm,Gc,Gv,Cv,Mr,Nf))
    }
    if(dNf <= sum(Number) | Nf_th <= dNf) {
      return(c(Lm,Gc,Gv,Cv,Mr,Nf))
    }
    
    Vect <- HMCfunc(dLm,dGc,dGv,dCv,dMr,dNf)
    plm <- plm - elm*Vect[2]
    pgc <- pgc - egc*Vect[3]
    pgv <- pgv - egv*Vect[4]
    pcv <- pcv - ecv*Vect[5]
    pmr <- pmr - emr*Vect[6]
    pnf <- pnf - enf*Vect[7]
    
  }
  
  dLm <- dLm + elm*plm
  dGc <- dGc + egc*pgc
  dGv <- dGv + egv*pgv
  dCv <- dCv + ecv*pcv
  dMr <- dMr + emr*pmr
  dNf <- dNf + enf*pnf
  
  if(dLm <= max(SCL)) {
    return(c(Lm,Gc,Gv,Cv,Mr,Nf))
  }
  if(dNf <= sum(Number) | Nf_th <= dNf) {
    return(c(Lm,Gc,Gv,Cv,Mr,Nf))
  }
  
  Vect <- HMCfunc(dLm,dGc,dGv,dCv,dMr,dNf)
  plm <- plm - elm*Vect[2]/2
  pgc <- pgc - egc*Vect[3]/2
  pgv <- pgv - egv*Vect[4]/2
  pcv <- pcv - ecv*Vect[5]/2
  pmr <- pmr - emr*Vect[6]/2
  pnf <- pnf - enf*Vect[7]/2
  
  Vect <- HMCfunc(dLm,dGc,dGv,dCv,dMr,dNf)
  
  HamEnd <- Vect[1] + plm**2/2 + pgc**2/2 + pgv**2/2 + pcv**2/2 + pmr**2/2 + pnf**2/2
  
  odds <- exp(HamInit-HamEnd)
  pmove <- min(odds,1)
  
  if(pmove >= runif(1,0,1)) {
    return(c(dLm,dGc,dGv,dCv,dMr,dNf))
  }
  else {
    return(c(Lm,Gc,Gv,Cv,Mr,Nf))
  }
}

HMC_Mr_Nf <- function(Lm,Gc,Gv,Cv,Mr,Nf) {

  dMr <- Mr
  dNf <- Nf

  pmr <- rnorm(1,0,taumr)
  pnf <- rnorm(1,0,taunf)
  
  Vect <- HMCfunc_Mr_Nf(Lm,Gc,Gv,Cv,dMr,dNf)
  HamInit <- Vect[1] + pmr**2/2 + pnf**2/2
  
  pmr <- pmr - emr*Vect[2]/2
  pnf <- pnf - enf*Vect[3]/2
  
  for(i in 1:L-1){
    
    dMr <- dMr + emr*pmr
    dNf <- dNf + enf*pnf
    
    if(dNf <= sum(Number) | Nf_th <= dNf) {
      return(c(Lm,Gc,Gv,Cv,Mr,Nf))
    }
    
    Vect <- HMCfunc_Mr_Nf(Lm,Gc,Gv,Cv,dMr,dNf)

    pmr <- pmr - emr*Vect[2]
    pnf <- pnf - enf*Vect[3]
    
  }

  dMr <- dMr + emr*pmr
  dNf <- dNf + enf*pnf
  
  if(dNf <= sum(Number) | Nf_th <= dNf) {
    return(c(Lm,Gc,Gv,Cv,Mr,Nf))
  }

  Vect <- HMCfunc_Mr_Nf(Lm,Gc,Gv,Cv,dMr,dNf)
 
  pmr <- pmr - emr*Vect[2]/2
  pnf <- pnf - enf*Vect[3]/2
  
  Vect <- HMCfunc_Mr_Nf(Lm,Gc,Gv,Cv,dMr,dNf)
  
  HamEnd <- Vect[1] + pmr**2/2 + pnf**2/2
  
  odds <- exp(HamInit-HamEnd)
  pmove <- min(odds,1)
  
  if(pmove >= runif(1,0,1)) {
    return(c(Lm,Gc,Gv,Cv,dMr,dNf))
  }
  else {
    return(c(Lm,Gc,Gv,Cv,Mr,Nf))
  }
}

# Init Preparation
Init_Lm <- Pri_M_Lm #known
Init_Gc <- Pri_M_Gc #known
Init_Gv <- 0.5 
Init_Cv <- 25 
Init_Mr <- 0.15
Init_Nf <- 3*sum(Number)

mcmc <- matrix(0,ncol=6,nrow=nsteps+1)
mcmc[1,1] <- Init_Lm
mcmc[1,2] <- Init_Gc
mcmc[1,3] <- Init_Gv
mcmc[1,4] <- Init_Cv
mcmc[1,5] <- Init_Mr
mcmc[1,6] <- Init_Nf

# Prepare progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = nsteps 
)

#MCMC
for (i in 1:nsteps) {
  
  if(i%%skipsteps==0){
  MCMCsample <- HMC(mcmc[i,1],mcmc[i,2],mcmc[i,3],mcmc[i,4],mcmc[i,5],mcmc[i,6])
  mcmc[i+1,1] <- MCMCsample[1]
  mcmc[i+1,2] <- MCMCsample[2]
  mcmc[i+1,3] <- MCMCsample[3]
  mcmc[i+1,4] <- MCMCsample[4]
  mcmc[i+1,5] <- MCMCsample[5]
  mcmc[i+1,6] <- MCMCsample[6]
  }
  else{
    MCMCsample <- HMC_Mr_Nf(mcmc[i,1],mcmc[i,2],mcmc[i,3],mcmc[i,4],mcmc[i,5],mcmc[i,6])
    mcmc[i+1,1] <- MCMCsample[1]
    mcmc[i+1,2] <- MCMCsample[2]
    mcmc[i+1,3] <- MCMCsample[3]
    mcmc[i+1,4] <- MCMCsample[4]
    mcmc[i+1,5] <- MCMCsample[5]
    mcmc[i+1,6] <- MCMCsample[6]
  }
  pb$tick()
}

# Density
Dens1 <- density(mcmc[-(1:burnin),1])
Dens2 <- density(mcmc[-(1:burnin),2])
Dens3 <- density(mcmc[-(1:burnin),3])
Dens4 <- density(mcmc[-(1:burnin),4])
Dens5 <- density(1-mcmc[-(1:burnin),5])
Dens6 <- density(mcmc[-(1:burnin),6])

# Mean
M_Lm <- mean(mcmc[-(1:burnin),1])
M_Gc <- mean(mcmc[-(1:burnin),2])
M_Gv <- mean(mcmc[-(1:burnin),3])
M_Cv <- mean(mcmc[-(1:burnin),4])
M_Mr <- mean(1-mcmc[-(1:burnin),5])
M_Nf <- mean(mcmc[-(1:burnin),6])

# 95%HDI
HDIwidth <- 0.95
LowHDI <- realsteps*(1-HDIwidth)/2
HighHDI <- realsteps-LowHDI
LowHDI_Lm <- hdi(Dens1,ci=HDIwidth)[1]
HighHDI_Lm <- hdi(Dens1,ci=HDIwidth)[2]
LowHDI_Gc <- hdi(Dens2,ci=HDIwidth)[1]
HighHDI_Gc <- hdi(Dens2,ci=HDIwidth)[2]
LowHDI_Gv <- hdi(Dens3,ci=HDIwidth)[1]
HighHDI_Gv <- hdi(Dens3,ci=HDIwidth)[2]
LowHDI_Cv <- hdi(Dens4,ci=HDIwidth)[1]
HighHDI_Cv <- hdi(Dens4,ci=HDIwidth)[2]
LowHDI_Mr <- hdi(Dens5,ci=HDIwidth)[1]
HighHDI_Mr <- hdi(Dens5,ci=HDIwidth)[2]
LowHDI_Nf <- hdi(Dens6,ci=HDIwidth)[1]
HighHDI_Nf <- hdi(Dens6,ci=HDIwidth)[2]

# Plot
dir.create("~/Desktop/Result_SCLdistribution")
setwd("~/Desktop/Result_SCLdistribution")
print("Plotting graphs...")
Lseq <- seq(0,99,1)
tseq <- seq(0,60,0.01)

# Contents plot
pdf("Contentsplot_SCLdistribution.pdf", width = 7, height = 6)
par(mfrow = c(2,2))
cplotseq <- seq(burnin,nsteps,by=realsteps/1000)

for(i in cplotseq){
  plot(tseq,Regression_curve(mcmc[i,3],mcmc[i,4],tseq),type="l",
       xlab="", ylab="", col=rgb(1, 0, 0, alpha=0.003), xlim=c(0,60))
  par(new=T)
}
plot(tseq,Regression_curve(mcmc[realsteps,3],mcmc[realsteps,4],tseq),
     type="l",col=rgb(1, 0, 0, alpha=0.003), xlim=c(0,60),
     xlab="age", ylab="Regression rate", main="Regression curve")

for(i in cplotseq){
  plot(tseq,Number_t(mcmc[i,5],mcmc[i,6],tseq),type="l",
       xlab="",ylab="",log="y",col=rgb(1, 0, 0, alpha=0.003),
       xlim=c(0,60),ylim=c(10,10000))
  par(new=T)
}
plot(tseq,Number_t(mcmc[realsteps,5],mcmc[realsteps,6],tseq),type="l",
     xlab="age",ylab="real_N",log="y",col=rgb(1, 0, 0, alpha=0.003),
     xlim=c(0,60),ylim=c(10,10000),main="Population")

for(i in cplotseq){
    if(mcmc[i,1]>99){
      plot(Lseq,Berta_curve(mcmc[i,1],mcmc[i,2],Lseq),type="l",
           xlab="",ylab="",col=rgb(1, 0, 0, alpha=0.003),
           xlim=c(0,100), ylim=c(0,60))
    par(new=T)
    }
  }
plot(Lseq,Berta_curve(mcmc[realsteps,1],mcmc[realsteps,2],Lseq),type="l",
     xlab="SCL",ylab="age",col=rgb(1, 0, 0, alpha=0.003),
     xlim=c(0,100), ylim=c(0,60), main="Growth Curve")

plot(SCL,Number,type="h",xlab="SCL",ylab="N",
     xlim=c(50,100),ylim=c(0,max(Number)+10),main="SCL distribution")
for(i in cplotseq){
  par(new=T)
  plot(Lseq,Modelfunc(mcmc[i,1],mcmc[i,2],mcmc[i,3],mcmc[i,4],mcmc[i,5],mcmc[i,6],Lseq),type="l",
       col=rgb(1, 0, 0, alpha=0.003),xlab="",ylab="",xlim=c(50,100),ylim=c(0,max(Number)+10))
}
dev.off()

# Fit plot
pdf("Fitplot_SCLdistribution.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(SCL,Number,type="h",xlab="",ylab="",
     xaxt="n", yaxt="n", xlim=c(50,100),ylim=c(0,max(Number)+20))
for(i in cplotseq){
  par(new=T)
  plot(Lseq,Modelfunc(mcmc[i,1],mcmc[i,2],mcmc[i,3],mcmc[i,4],mcmc[i,5],mcmc[i,6],Lseq),type="l",
       col=rgb(1, 0, 0, alpha=0.003),xlab="",ylab="",
       xaxt="n", yaxt="n", xlim=c(50,100), ylim=c(0,max(Number)+20))
}
axis(1, at = pretty(SCL), labels = pretty(SCL))
mtext("SCL(cm)", side = 1, line =3)
axis(2, at = pretty(Number), labels = pretty(Number))
mtext("Observed Number", side = 2, line = 3)
dev.off()

# Posterior plot
pdf("Posteriorplot_SCLdistribution.pdf", width = 7, height = 6)
par(mfrow=c(3,2))

plot(-100,-100,xlim=c(min(Dens1$x),max(Dens1$x)),ylim=c(min(Dens1$y),max(Dens1$y)),xlab="",ylab="Density",main="SCL_Max")
lines(Dens1,lwd=1,col=1)
abline(v=M_Lm,lty=1,col=2,lwd=1)
abline(v=LowHDI_Lm,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Lm,lty=3,col=2,lwd=1.5)
abline(v=Init_Lm,lty=3,col=3,lwd=1.5)

plot(-100,-100,xlim=c(min(Dens2$x),max(Dens2$x)),ylim=c(min(Dens2$y),max(Dens2$y)),xlab="",ylab="Density",main="Growth Coefficient")
lines(Dens2,lwd=1,col=1)
abline(v=M_Gc,lty=1,col=2,lwd=1)
abline(v=LowHDI_Gc,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Gc,lty=3,col=2,lwd=1.5)
abline(v=Init_Gc,lty=3,col=3,lwd=1.5)

plot(-100,-100,xlim=c(min(Dens3$x),max(Dens3$x)),ylim=c(min(Dens3$y),max(Dens3$y)),xlab="",ylab="Density",main="Gain Value")
lines(Dens3,lwd=1,col=1)
abline(v=M_Gv,lty=1,col=2,lwd=1)
abline(v=LowHDI_Gv,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Gv,lty=3,col=2,lwd=1.5)
abline(v=Init_Gv,lty=3,col=3,lwd=1.5)

plot(-100,-100,xlim=c(min(Dens4$x),max(Dens4$x)),ylim=c(min(Dens4$y),max(Dens4$y)),xlab="",ylab="Density",main="Center Value")
lines(Dens4,lwd=1,col=1)
abline(v=M_Cv,lty=1,col=2,lwd=1)
abline(v=LowHDI_Cv,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Cv,lty=3,col=2,lwd=1.5)
abline(v=Init_Cv,lty=3,col=3,lwd=1.5)

plot(-100,-100,xlim=c(min(Dens5$x),max(Dens5$x)),ylim=c(min(Dens5$y),max(Dens5$y)),xlab="",ylab="Density",main="Survival Rate")
lines(Dens5,lwd=1,col=1)
abline(v=M_Mr,lty=1,col=2,lwd=1)
abline(v=LowHDI_Mr,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Mr,lty=3,col=2,lwd=1.5)
abline(v=Init_Mr,lty=3,col=3,lwd=1.5)

plot(-100,-100,xlim=c(min(Dens6$x),max(Dens6$x)),ylim=c(min(Dens6$y),max(Dens6$y)),xlab="",ylab="Density",main="Pupulation Coefficient")
lines(Dens6,lwd=1,col=1)
abline(v=M_Nf,lty=1,col=2,lwd=1)
abline(v=LowHDI_Nf,lty=3,col=2,lwd=1.5)
abline(v=HighHDI_Nf,lty=3,col=2,lwd=1.5)
abline(v=Init_Nf,lty=3,col=3,lwd=3)

dev.off()

# Mcmc plot
pdf("mcmcplot_SCLdistribution.pdf", width = 7, height = 6)
mcmcplot <- seq(burnin,nsteps,by=realsteps/10000)
par(mfrow=c(3,2))
plot(mcmc[mcmcplot,1],xlab="steps",ylab="",main="SCL_max",type="l")
plot(mcmc[mcmcplot,2],xlab="steps",ylab="",main="Growth Coefficient",type="l")
plot(mcmc[mcmcplot,3],xlab="steps",ylab="",main="Gain Value",type="l")
plot(mcmc[mcmcplot,4],xlab="steps",ylab="",main="Center Value",type="l")
plot(1-mcmc[mcmcplot,5],xlab="steps",ylab="",main="Survival Rate",type="l")
plot(mcmc[mcmcplot,6],xlab="steps",ylab="",main="Pupulation Coefficient",type="l")
dev.off()

# list
StatisticalResults <- matrix(ncol=4,nrow=6)
colnames(StatisticalResults) <- c("Name","Low_HDI","Mean","High_HDI")
StatisticalResults[,1] <- c("Maximum carapace length","Growth coefficient","Gain","Threshold","Survival rate","Pupulation coefficient")
StatisticalResults[1,-1] <- round(c(LowHDI_Lm,M_Lm,HighHDI_Lm),digits=3)
StatisticalResults[2,-1] <- round(c(LowHDI_Gc,M_Gc,HighHDI_Gc),digits=4)
StatisticalResults[3,-1] <- round(c(LowHDI_Gv,M_Gv,HighHDI_Gv),digits=3)
StatisticalResults[4,-1] <- round(c(LowHDI_Cv,M_Cv,HighHDI_Cv),digits=2)
StatisticalResults[5,-1] <- round(c(LowHDI_Mr,M_Mr,HighHDI_Mr),digits=3)
StatisticalResults[6,-1] <- round(c(LowHDI_Nf,M_Nf,HighHDI_Nf),digits=-1)

write.csv(as.data.frame(StatisticalResults),"StatisticalResults_SCLdistribution.csv")
