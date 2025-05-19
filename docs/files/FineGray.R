pbc3 <- read.csv("data/pbc3.csv")
pbc3$years <- pbc3$days/365.25
pbc3$log2bili<-log2(pbc3$bili)
pbc3$fstatus <- factor(pbc3$status, 0:2, labels=c("cens", "trans", "death"))

library(survival)

# 1 Aalen-Johansen
print(crsftment<-survfit(Surv(years, fstatus) ~ tment, data = pbc3))
plot(crsftment,col=c(2,2,3,3),lty=c(1,2,1,2))

plot(crsftment[,"death"],col=3,lty=1:2)

lines(crsftment[,"trans"],col=2,lty=1:2)


# 2 Fine-Gray models
# Death
fgdeath <- finegray(Surv(years, fstatus) ~ ., etype="death", data=pbc3)

coxph(Surv(fgstart, fgstop, fgstatus) ~ tment,  
      weight = fgwt, data = fgdeath)

# Trans
fgtrans <- finegray(Surv(years, fstatus) ~ ., etype="trans", data=pbc3)

coxph(Surv(fgstart, fgstop, fgstatus) ~ tment, 
      weight=fgwt, data=fgtrans)


# 3 Fine-Gray models
# Death
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili,  
      weight = fgwt, data = fgdeath)

# Trans
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili, 
      weight=fgwt, data=fgtrans)


# 4
# Non-parametric Aalen-Johansen
crsfsex<-survfit(Surv(years, fstatus) ~ sex, data = pbc3)
plot(crsfsex,col=c(2,2,3,3),lty=c(1,2,1,2))
plot(crsfsex[,"death"],col=3,lty=1:2)
lines(crsfsex[,"trans"],col=2,lty=1:2)


# FG Death
summary(fgfitdeath<-coxph(Surv(fgstart, fgstop, fgstatus) ~ sex, 
                          weight=fgwt, data=fgdeath))
ndata <- data.frame(sex=0:1)
fgsfdeath <- survfit(fgfitdeath, ndata)
plot(fgsfdeath, fun="event", lty=1:2, col=3)
lines(crsfsex[,"death"],col=1,lty=1:2)

# FG Trans
summary(fgfittrans<-coxph(Surv(fgstart, fgstop, fgstatus) ~ sex, 
                          weight=fgwt, data=fgtrans))
fgsftrans <- survfit(fgfittrans, ndata)
plot(fgsftrans, fun="event", lty=1:2, col=2)
lines(crsfsex[,"trans"],col=1,lty=1:2)
