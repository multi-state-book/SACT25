pbc3 <- read.csv("https://multi-state-book.github.io/SACT25/files/pbc3.csv")
#pbc3 <- read.csv("data/pbc3.csv")
pbc3$years <- pbc3$days/365.25
library(survival)

# Two-state model ----
with(pbc3,table(status,tment))
with(pbc3,table(status!=0,tment))
## Intensities ----
### 1. tment ----
print(sftment<-survfit(Surv(years, status != 0) ~ tment, data = pbc3))
plot(sftment,cumhaz=T,lty=1:2,main="tment")
legend("topleft",legend=c("tment=0", "tment=1"),lty=c(1,2))
survdiff(Surv(years,status!=0) ~ tment, data=pbc3)

### 2. bili ----
# bili in quartiles 
pbc3$bili4<-cut(pbc3$bili,quantile(pbc3$bili, na.rm=T), include.lowest=T)
print(sfbili4<-survfit(Surv(years, status != 0) ~ bili4, data = pbc3))
plot(sfbili4,cumhaz=T,col=1:4,main="bili")
survdiff(Surv(years,status!=0) ~ bili4, data=pbc3)

# bili in binary; upper quartile 
pbc3$bili2<-cut(pbc3$bili,c(0,42.3,454),na.rm=T)
print(sfbili2<-survfit(Surv(years, status != 0) ~ bili2, data = pbc3))
plot(sfbili2,cumhaz=T,col=1:2,main="bili")
survdiff(Surv(years,status!=0) ~ bili2, data=pbc3)

### 3. Stratified ----
survdiff(Surv(years,status!=0) ~ tment+strata(bili4), data=pbc3)
survdiff(Surv(years,status!=0) ~ tment+strata(bili2), data=pbc3)
print(sftmentbili2<-survfit(Surv(years, status != 0) ~ tment + bili2, data = pbc3))
plot(sftmentbili2,cumhaz=T,lty=c(1,1,2,2),col=c(1,2,1,2),main="tment by bili")

## Marginal parameters ----
### 4. Kaplan-Meier ----
plot(sftment,lty=1:2,main="KM tment")
legend("bottomleft",legend=c("tment=0", "tment=1"),lty=c(1,2))

plot(sfbili2,col=1:2,,main="KM bili")

plot(sftmentbili2,lty=c(1,1,2,2),col=c(1,2,1,2),main="KM tment by bili")

### 5. RD ----
summary(sftment,t=3)
summary(sftment,t=3)$surv
summary(sftment,t=3)$std.err

riskdiff<-function(km,t){ 
  kmsurv<- summary(km, t)
  risk0 <- 1-kmsurv$surv[1]
  risk1 <- 1-kmsurv$surv[2]
  sd0 <- kmsurv$std.err[1]
  sd1 <- kmsurv$std.err[2]
  rd <-risk1-risk0
  sd <- sqrt(sd0^2 + sd1^2)
  print(paste("Risk0 = ",round(risk0,6)))
  print(paste("Risk1 = ",round(risk1,6)))
  cat("Riskdiff (SD) [95%-CI]\n",rd," (",sd,")", " [", rd-1.96*sd, " to ", rd+1.96*sd,"]", sep="") 
}

riskdiff(sftment,3)

### 6. RMST ----
print(sftment,rmean=3)
summary(sftment,rmean=3)
summary(sftment,rmean=3)$table

rmstdiff<-function(km,t){
  rmst <- summary(km, rmean=t)
  m0 <- rmst$table[1,5]
  m1 <- rmst$table[2,5]
  md <- m1-m0
  sd <- sqrt(rmst$table[1,6]^2 + rmst$table[2,6]^2)
  print(paste("RMST0 = ",round(m0,6)))
  print(paste("RMST1 = ",round(m1,6)))
  cat("RMSTDiff (SD) [95%-CI]\n", md," (",sd,")", " [", md-1.96*sd, " to ", md+1.96*sd,"]", sep="") 
}

rmstdiff(sftment,3)


# Competing risks ----
## Intensities ----
### 1. tment ----
pbc3$fstatus <- factor(pbc3$status, 0:2, labels=c("cens", "trans", "death"))
with(pbc3,table(fstatus,tment))

print(sftrans<-survfit(Surv(years, fstatus=="trans") ~ tment, data = pbc3))
plot(sftrans,cumhaz=T,lty=1:2,main="tment")
print(sfdeath<-survfit(Surv(years, fstatus=="death") ~ tment, data = pbc3))
plot(sfdeath,cumhaz=T,lty=1:2,main="tment")
lines(sftrans,cumhaz=T,lty=1:2)

print(crsftment<-survfit(Surv(years, fstatus) ~ tment, data = pbc3))
print(crsftment,rmean=3)
plot(crsftment,cumhaz=T,col=c(1,1,2,2),lty=c(1,2,1,2),main="tment")

# hmm ...
survdiff(Surv(years,fstatus) ~ tment, data=pbc3)

# Transplantation
survdiff(Surv(years,fstatus=="trans") ~ tment, data=pbc3)
# Death
survdiff(Surv(years,fstatus=="death") ~ tment, data=pbc3)

### 2. bili ----
# bili4
# Transplantation
print(transbili4<-survfit(Surv(years, fstatus=="trans") ~ bili4, data = pbc3))
plot(transbili4,cumhaz=T,col=1:4,main="bili")
survdiff(Surv(years,fstatus=="trans") ~ bili4, data=pbc3)
# Death
print(deathbili4<-survfit(Surv(years, fstatus=="death") ~ bili4, data = pbc3))
plot(deathbili4,cumhaz=T,col=1:4,main="bili")
survdiff(Surv(years,fstatus=="death") ~ bili4, data=pbc3)

# bili2
# Transplantation
print(transbili2<-survfit(Surv(years, fstatus=="trans") ~ bili2, data = pbc3))
plot(transbili2,cumhaz=T,col=1:2,main="bili")
survdiff(Surv(years,fstatus=="trans") ~ bili2, data=pbc3)
# Death
print(deathbili2<-survfit(Surv(years, fstatus=="death") ~ bili2, data = pbc3))
plot(deathbili2,cumhaz=T,col=1:2,main="bili")
survdiff(Surv(years,fstatus=="death") ~ bili2, data=pbc3)

### 3. Stratified ----
# Transplantation
survdiff(Surv(years,fstatus=="trans") ~ tment+strata(bili2), data=pbc3)
# Death
survdiff(Surv(years,fstatus=="death") ~ tment+strata(bili2), data=pbc3)

## Marginal parameters ----
### 4. Aalen-Johansen ----
# tment
print(crsftment)
plot(crsftment,col=c(1,1,2,2),lty=1:2)
plot(crsftment,lty=1:2, noplot = "")

# bili2
print(crsfbili2 <- survfit(Surv(years,fstatus) ~ bili2, data=pbc3),rmean=3)
plot(crsfbili2,col=c(1,1,2,2),lty=1:2)

### 5. RD ----
plot(crsftment,col=c(1,1,2,2),lty=1:2)
summary(crsftment,t=3)
summary(crsftment,t=3)$pstate
summary(crsftment,t=3)$std.err

cidiff<-function(aj,etype=2,t){
  ajt<- summary(aj, t) 
  prob0 <- ajt$pstate[1,etype] 
  prob1 <- ajt$pstate[2,etype]
  sd0 <- ajt$std.err[1,etype] 
  sd1 <- ajt$std.err[2,etype] 
  dd <- prob1-prob0 
  sd <- sqrt(sd0^2 + sd1^2) 
  print(paste("Event type =",aj$state[etype])) 
  print(paste("prob0 =",round(prob0,6)))
  print(paste("prob1 =",round(prob1,6)))
  cat("ProbDiff (SD) [95%-CI]\n",dd," (",sd,")", " [", dd-1.96*sd, " to ", dd+1.96*sd,"]", sep="")
}

cidiff(crsftment,etype=1,3)
riskdiff(sftment,3)
cidiff(crsftment,etype=2,3)
cidiff(crsftment,etype=3,3)


### 6. YL ----
print(crsftment,rmean=3)
summary(crsftment,rmean=3)$table
summary(crsftment,rmean=3)$state

yldiff<-function(aj,etype=1,t){
  yl <- summary(aj,rmean=t)
  m0 <- yl$table[2*etype-1,3]
  m1 <- yl$table[2*etype,3]
  md <- m1-m0
  sd <- sqrt(yl$table[2*etype-1,4]^2 + yl$table[2*etype,4]^2)
  print(paste("Event type =",aj$state[etype]))
  print(paste("mean0 = ",round(m0,6)))
  print(paste("mean1 = ",round(m1,6)))
  cat("MeanDiff (SD) [95%-CI]\n",md," (",sd,")", " [", md-1.96*sd, " to ", md+1.96*sd,"]", sep="")
}

yldiff(crsftment,etype=1,3)
rmstdiff(sftment,3)
yldiff(crsftment,etype=2,3)
yldiff(crsftment,etype=3,3)
