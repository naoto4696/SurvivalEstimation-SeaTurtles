setwd("~/Desktop")

# reading data
SCL_Radius <- read.csv("~/Desktop/Program/SCL_Radius.csv")
Width_Radius <- read.csv("~/Desktop/Program/Width_Radius.csv")
Realdata <- read.csv("~/Desktop/Program/Muroto.csv")

dir.create("~/Desktop/Result_ExistingMethod")
setwd("~/Desktop/Result_ExistingMethod")

# setting
SCL <- as.vector(as.ts(SCL_Radius$SCL/10)) # Convert mm to cm
Radius <- as.vector(as.ts(SCL_Radius$Radius/2)) # Convert diameter to radius
SCL <- SCL[is.na(SCL)==F]
Radius <- Radius[is.na(Radius)==F]

# estimation

lm_result_SR <- lm(SCL~Radius-1) 
Coef_SR <- lm_result_SR$coefficients[1]

C <- Coef_SR

# fit plot

pdf("Fitplot_SCL_Radius_existing.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(Radius,SCL, type="p", xlab="", ylab="",
     xlim=c(min(Radius),max(Radius)), ylim=c(min(SCL),max(SCL)),
     xaxt="n",yaxt="n", pch=20,
     panel.first=grid(col="grey", lty=2))
abline(0, Coef_SR, col="blue", lwd=1.5)
axis(1, at = pretty(Radius), labels = pretty(Radius))
mtext("Radius(mm)", side = 1, line =3)
axis(2, at = pretty(SCL), labels = pretty(SCL))
mtext("SCL(cm)", side = 2, line = 3)
dev.off()

# setting

Width_dummy <- as.vector(as.ts(Width_Radius$Width)) 
Radius_dummy <- as.vector(as.ts(Width_Radius$Radius)) 
Width <- Width_dummy[is.na(Width_dummy)==F]
Radius <- Radius_dummy[is.na(Width_dummy)==F]

# estimation

lm_result_WR <- lm(Width~Radius) 
Intercept_WR <- lm_result_WR$coefficients[1]
Coef_WR <- lm_result_WR$coefficients[2]

Hm <- -Intercept_WR/Coef_WR
Gh <- log(1-Coef_WR)

# fit plot

pdf("fitplot_Width_Radius_existing.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(Radius, Width, type="p", xlab="", ylab="", 
     xlim=c(min(Radius),max(Radius)), ylim=c(min(Width),max(Width)), 
     xaxt="n", yaxt="n", pch=20,
     panel.first=grid(col="grey", lty=2))
abline(Intercept_WR, Coef_WR, col="blue", lwd=1.5)
axis(1, at = pretty(Radius), labels = pretty(Radius))
mtext("Radius(mm)", side = 1, line =3)
axis(2, at = pretty(Width), labels = pretty(Width))
mtext("Width(mm)", side = 2, line = 3)
dev.off()

#setting
SCL <- as.ts(Realdata$SCL)
Number <- as.ts(Realdata$Number)
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
N_max_point <- which.max(Number)

Lm <- C*Hm
Gc <- Gh
age <- (-log((Lm-SCL)/Lm))/Gc

age <- age[N_max_point:length(age)]
Num <- log(Number[N_max_point:length(Number)])

# estimation

lm_result_St <- lm(Num~age)
Intercept_St <- lm_result_St$coefficients[1]
Coef_St <- lm_result_St$coefficients[2]

Sr <- 1+Coef_St

xlongseq <- seq(1,1000,0.01)
tlongseq <- Lm*(1-exp(-Gc*xlongseq))
ylongseq <- exp(Intercept_St+xlongseq*Coef_St)

# fit plot

pdf("Fitplot_SCLdistribution.pdf", width = 7, height = 6)
par(mfrow=c(1,1))
plot(SCL,Number,type="h",xlab="",ylab="",
     xaxt="n", yaxt="n", xlim=c(50,100),ylim=c(0,max(Number)+20))
par(new=TRUE)
plot(tlongseq, ylongseq, col="blue", lwd=1.5, xlab="",ylab="", type="l",
     xaxt="n", yaxt="n", xlim=c(50,100),ylim=c(0,max(Number)+20))
axis(1, at = pretty(SCL), labels = pretty(SCL))
mtext("SCL(cm)", side = 1, line =3)
axis(2, at = pretty(Number), labels = pretty(Number))
mtext("Observed Number", side = 2, line = 3)
dev.off()

# list

StatisticalResults <- matrix(ncol=2,nrow=5)
colnames(StatisticalResults) <- c("Name","Result")
StatisticalResults[,1] <- c("Maximum humerus radius Hm","Regression coefficient C",
                            "Maximum carapace length Lm","Growth coefficient Gc","Survival rate Sr")
StatisticalResults[1,-1] <- round(Hm,digits=5)
StatisticalResults[2,-1] <- round(C,digits=5)
StatisticalResults[3,-1] <- round(Lm,digits=3)
StatisticalResults[4,-1] <- round(Gc,digits=4)
StatisticalResults[5,-1] <- round(Sr,digits=3)

write.csv(as.data.frame(StatisticalResults),"StatisticalResults_ExistingMethod.csv")
