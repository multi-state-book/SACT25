pbc3 <- read.csv("data/pbc3.csv")
pbc3$years <- pbc3$days/365.25
pbc3$fail     <- as.numeric(with(pbc3, status>0))
pbc3$log2bili<-log2(pbc3$bili)
pbc3$fstatus <- factor(pbc3$status, 0:2, 
                       labels=c("cens", "trans", "death"))
library(survival)
library(pseudo)
library(geepack)

## Repetition
# 1.
km <- survfit(Surv(years,status!=0)~tment,data=pbc3)
plot(km, lty=1:2)
plot(km, lty=1:2,fun="event")

# 2. 
riskdiff(km,3)

# 3.
rmstdiff(km,3)

# 4.
summary(coxph(Surv(years,status!=0)~tment,data=pbc3))

# 5.
summary(coxph(Surv(years,status!=0)~tment + alb + log2(bili),data=pbc3))

# 6.
print(crsftment<-survfit(Surv(years, fstatus) ~ tment, data = pbc3))
plot(crsftment,col=c(2,2,3,3),lty=c(1,2,1,2))
plot(crsftment[,"death"],col=3,lty=1:2)
lines(crsftment[,"trans"],col=2,lty=1:2)

# 7.
cidiff(crsftment,etype=1,3)
cidiff(crsftment,etype=2,3)
cidiff(crsftment,etype=3,3)

# 8.
yldiff(crsftment,etype=1,3)
yldiff(crsftment,etype=2,3)
yldiff(crsftment,etype=3,3)

#9.
# Death
fgdeath <- finegray(Surv(years, fstatus) ~ ., etype="death", data=pbc3)
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment,  weight = fgwt, data = fgdeath)

# Trans
fgtrans <- finegray(Surv(years, fstatus) ~ ., etype="trans", data=pbc3)
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment, weight=fgwt, data=fgtrans)


# 10.
# Death
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili,  weight = fgwt, data = fgdeath)

# Trans
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili, weight=fgwt, data=fgtrans)
