pbc3 <- read.csv("data/pbc3.csv")
pbc3$years <- pbc3$days/365.25
library(survival)

# Two-state model ----
## Intensity: ----
### 1. tment ----
summary(coxtment<-coxph(Surv(years,status>0)~tment, data=pbc3, method="breslow"))
print(lr<-survdiff(Surv(years,status!=0) ~ tment, data=pbc3))
coxtment$score
lr$chisq

### 2. bili2 ----
# bili in binary; upper quartile 
pbc3$fbili2<-cut(pbc3$bili,c(0,42.3,454),na.rm=T)
pbc3$bili2<-cut(pbc3$bili,c(0,42.3,454),na.rm=T,labels=F)
summary(fitbili2<-coxph(Surv(years,status>0)~bili2, data=pbc3, method="breslow"))
sfcoxbili2<-survfit(fitbili2, newdata = data.frame(bili2=1:2))
plot(sfcoxbili2, cumhaz = T, col=1:2)

sfbili2<-survfit(Surv(years, status != 0) ~ bili2, data = pbc3)
lines(sfbili2,cumhaz=T,col=3,main="bili")

### 3. tment+bili2 ----
coxph(Surv(years,status>0)~tment+fbili2, data=pbc3, method="breslow")
coxph(Surv(years,status>0)~tment:fbili2+fbili2, data=pbc3, method="breslow")
coxph(Surv(years,status>0)~tment*fbili2, data=pbc3, method="breslow")

### 4. Stratified Cox ----
summary(coxstrata<-coxph(Surv(years,status>0)~tment+strata(fbili2), data=pbc3, method="breslow"))
lr<-survdiff(Surv(years,status!=0) ~ tment+strata(fbili2), data=pbc3)
coxstrata$score
lr$chisq

coxph(Surv(years,status>0)~tment:strata(fbili2), data=pbc3, method="breslow")
coxph(Surv(years,status>0)~tment*strata(fbili2), data=pbc3, method="breslow")

### 5. bili ----
coxph(Surv(years,status>0)~bili, data=pbc3, method="breslow")
coxph(Surv(years,status>0)~bili+I(bili^2), data=pbc3, method="breslow")
coxph(Surv(years,status>0)~I(bili/100)+I((bili/100)^2), data=pbc3, method="breslow")

pbc3$log2bili<-log2(pbc3$bili)

coxph(Surv(years,status>0)~log2bili,data=pbc3, method="breslow")
coxph(Surv(years,status>0)~log2bili+I(log2bili^2), data=pbc3, method="breslow")

### 6. tment+alb+logbili ----
coxph(Surv(years, status != 0)~tment + alb + log2bili,
      data = pbc3, method = "breslow")
coxph(Surv(years, status != 0)~tment + alb + log2bili+ log2bili:tment, 
      data = pbc3, method = "breslow")
coxph(Surv(years, status != 0)~factor(tment) + alb + log2bili:factor(tment), 
      data = pbc3, method = "breslow")


### 7. time-dependent ----
coxph(Surv(years, status != 0) ~ tment + alb + log2bili + tt(tment), 
      data = pbc3, tt = function(x,t, ...) (x)*(t>2), method = "breslow")
coxph(Surv(years, status != 0) ~ tment + alb + log2bili + tt(alb), 
      data = pbc3, tt = function(x,t, ...) (x)*(t>2), method = "breslow")
coxph(Surv(years, status != 0) ~ tment + alb + log2bili + tt(log2bili), 
      data = pbc3, tt = function(x,t, ...) (x)*(t>2), method = "breslow")

## Marginal parameters ----
### 8. Figure 4.5 ----
print(mcox<-coxph(Surv(years,status>0)~tment+alb+log2bili, data=pbc3, method="breslow"))
predcr<-data.frame(tment=0:1,alb=38,log2bili=log2(45))
plot(survfit(mcox,newdata=predcr),lty=1:2,xlab="Years")

# Competing risks ----
## Intensity: ----
### 1. Table 2.13, p. 62 ----
pbc3$fstatus <- factor(pbc3$status, 0:2,
                       labels=c("cens", "trans", "death"))
coxph(Surv(years, fstatus) ~ tment + alb + log2bili + sex + age,
      data = pbc3, method = "breslow", id=id)

coxph(Surv(years, status>0) ~ tment + alb + log2bili + sex + age,
      data = pbc3, method = "breslow", id=id)

### 2. Stratify for sex  ----
coxph(list(Surv(years, fstatus) ~ tment + alb + log2bili + age,
           1:2+1:3~strata(sex)), data = pbc3, method = "breslow", id=id)

### 3. Only sex effect for death  ----
coxph(list(Surv(years, fstatus) ~ tment + alb + log2bili + age,
           1:3~sex), data = pbc3, method = "breslow", id=id)

## Marginal parameters ----
### Figure 4.12 ----
print(fitcr<-coxph(Surv(years,fstatus)~tment+alb+log2bili+sex+age,
             method = "breslow", data = pbc3, id=id))

predcr<-data.frame(tment=0,alb=38,log2bili=log2(45),sex=0,age=40)
plot(survfit(fitcr,newdata=predcr),lty=1:2)

predcr<-data.frame(tment=0,alb=38,log2bili=log2(45),sex=0,age=60)
plot(survfit(fitcr,newdata=predcr),lty=1:2,ylim=c(0,1))