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


## Pseudo-values
## Defining two summary functions for summarizing a geese fit 
## (one without `exp(est)` one that does)
posumm<-function(pofit,d=6){
  round(cbind(
    Est   = pofit$beta,
    SD    = sqrt(diag(pofit$vbeta)),
    lo.ci = pofit$beta-1.96*sqrt(diag(pofit$vbeta)),
    up.ci = pofit$beta+1.96*sqrt(diag(pofit$vbeta)),
    Wald  = (pofit$beta/sqrt(diag(pofit$vbeta)))^2,
    PVal  = 2-2*pnorm(abs(pofit$beta/sqrt(diag(pofit$vbeta))))),d)
}

posummExp<-function(pofit,d=6){
  round(cbind(
    est       = pofit$beta,
    SD        = sqrt(diag(pofit$vbeta)),
    exp.est   = exp(pofit$beta),
    exp.lo.ci = exp(pofit$beta-1.96*sqrt(diag(pofit$vbeta))),
    exp.up.ci = exp(pofit$beta+1.96*sqrt(diag(pofit$vbeta))),
    PVal      = 2-2*pnorm(abs(pofit$beta/sqrt(diag(pofit$vbeta))))),d)
}

# 1.
## Add survival POs and failure (event) POs (1-PO) at t = 3 years:
po3 <- pseudosurv(pbc3$years, pbc3$fail,tmax = 3)
pbc3$po3<-as.vector(po3$pseudo)
pbc3$epo3<-as.vector(1-po3$pseudo)

# 2.
## Use *event* POs to model the risk differences.
posumm(geese(epo3 ~ tment, data = pbc3, id = id, 
             mean.link = "identity"))
riskdiff(km,t=3)
posumm(geese(epo3 ~ tment + alb + log2(bili),
             data = subset(pbc3, !is.na(alb)),
             id = id, mean.link = "identity"))

posummExp(geese(epo3 ~ tment, data = pbc3, id = id,
                mean.link = "log"))
posummExp(geese(epo3 ~ tment+alb+log2bili, 
                data = subset(pbc3, !is.na(alb)), 
                id = id, mean.link = "log"))

# 3.
pbc3$rmst3<-pseudomean(pbc3$years, pbc3$fail, tmax = 3)
posumm(geese(rmst3 ~ tment, data = pbc3, 
             id = id, mean.link = "identity"))
rmstdiff(km,t=3)
posumm(geese(rmst3 ~ tment+alb+log2bili, data =subset(pbc3, !is.na(alb)),
             id = id, mean.link = "identity"))

# 4.
## Calculate POs from KM at chosen times and reshape data to long format. 
## Note '1-PO' is created.
potsurv <- pseudosurv(pbc3$years, pbc3$fail,tmax = 1:4)
longpbc3 <- NULL
for(it in 1:length(potsurv$time)){
  longpbc3 <- rbind(longpbc3,
                    cbind(pbc3,
                          pseudo  = 1-potsurv$pseudo[,it],
                          tpseudo = potsurv$time[it],
                          id      = 1:nrow(pbc3)))
}
longpbc3 <- longpbc3[order(longpbc3$id),]

# joint model:
posummExp(geese(pseudo~as.factor(tpseudo)+tment,id=id,
                data=longpbc3, mean.link="cloglog",corstr="independence"))
# Cox model for comparison:
summary(coxph(Surv(years,status!=0)~tment,data=pbc3))

# 5.
# joint model:
posummExp(geese(pseudo~as.factor(tpseudo) + tment + alb + log2(bili), id=id,
                data=subset(longpbc3, !is.na(alb)), 
                mean.link="cloglog",corstr="independence"))
# Cox model for comparison:
summary(coxph(Surv(years,status!=0)~ tment + alb + log2(bili), data=pbc3))


# 6.
## Using pseudoci()
cipo3 <- pseudoci(pbc3$years, pbc3$status,tmax = 3)
pbc3$trans.po3<-as.vector(cipo3$pseudo[[1]])
pbc3$death.po3<-as.vector(cipo3$pseudo[[2]])

# 7.
posumm(geese(trans.po3 ~ tment, data = pbc3, id = id,
             mean.link = "identity"))
cidiff(crsftment,etype=2,3)

posumm(geese(death.po3 ~ tment, data = pbc3, id = id, mean.link = "identity"))
cidiff(crsftment,etype=3,3)

# adjusted
posumm(geese(trans.po3 ~ tment + alb + log2(bili),
             data = subset(pbc3, !is.na(alb)),
             id = id,  mean.link = "identity"))
posumm(geese(death.po3 ~ tment + alb + log2(bili),
             data = subset(pbc3, !is.na(alb)),
             id = id,  mean.link = "identity"))


# 8.
yl3 <- pseudoyl(pbc3$years, pbc3$status,tmax = 3)
pbc3$trans.yl3<-as.vector(yl3$pseudo[[1]])
pbc3$death.yl3<-as.vector(yl3$pseudo[[2]])

posumm(geese(trans.yl3 ~ tment, data = pbc3, id = id, mean.link = "identity"))
yldiff(crsftment,etype=2,3)

posumm(geese(death.yl3 ~ tment, data = pbc3, id = id, mean.link = "identity"))
yldiff(crsftment,etype=3,3)

# 9.

pot <- pseudoci(pbc3$years, pbc3$status,tmax = 1:4)
longpbc3 <- NULL
for(it in 1:length(pot$time)){
  longpbc3 <- rbind(longpbc3,
                    cbind(pbc3,trans.po = pot$pseudo[[1]][,it],
                          pbc3,death.po = pot$pseudo[[2]][,it],
                          tpseudo       = pot$time[it],
                          id            = 1:nrow(pbc3)))
}
longpbc3 <- longpbc3[order(longpbc3$id),]

posummExp(geese(trans.po ~ as.factor(tpseudo)+tment,id=id,data=longpbc3,
                mean.link="cloglog",corstr="independence"))

posummExp(geese(death.po ~ as.factor(tpseudo)+tment,id=id,data=longpbc3,
                mean.link="cloglog",corstr="independence"))

# Fine-Gray models for comparison
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment, weight=fgwt, data=fgtrans)
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment, weight=fgwt, data=fgdeath)

# 10.
posummExp(geese(trans.po ~ as.factor(tpseudo) + tment + alb + log2(bili),id=id,
                data=subset(longpbc3,!is.na(alb)),
                mean.link="cloglog",corstr="independence"))

posummExp(geese(death.po ~ as.factor(tpseudo) + tment + alb + log2(bili),id=id,
                data=subset(longpbc3,!is.na(alb)),
                mean.link="cloglog",corstr="independence"))

# Fine-Gray models for comparison
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili, weight=fgwt, data=fgtrans)
coxph(Surv(fgstart, fgstop, fgstatus) ~ tment+alb+log2bili, weight = fgwt, data = fgdeath)
