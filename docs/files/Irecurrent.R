library(survival)
rr<-read.csv("rr.csv")
View(rr) # id=172
with(rr,table(status))
with(rr,table(eventno))

#######################################################################
## Intensities
# 1. Nelson-Aalen for recurrent
plot(survfit(Surv(start,stop,status==1)~treat,data=rr),main="recurrent, NAa",
     fun="cumhaz",lty=1:2, xlab="time",ylab="Integrated hazard")

# 2. AG model
coxph(Surv(start,stop,status==1) ~ treat,data=rr)

# 3. AG model given previous event
coxph(Surv(start,stop,status==1) ~ treat+I(eventno>0), data=rr)

# 4. Death
# Nelson-Aalen
plot(survfit(Surv(start,stop,status==2)~treat,data=rr),fun="cumhaz",lty=1:2,
     main="Death: NAa",xlab="time",ylab="Cumulative hazard")
# Cox
coxph(Surv(start,stop,status==2)~treat,data=rr)

# 5. Time to first
# NA time to first event
t2fna<-survfit(Surv(start,stop,status==1)~treat,
               data=subset(rr,eventno==0))
plot(t2fna,lty=1:2,fun="cumhaz",xlab="time",ylab="Integrated hazard",main="T2F")
# Cox model time to first event
coxph(Surv(start,stop,status==1)~treat, data=subset(rr,eventno==0))
# Compare SD to AG model
coxph(Surv(start,stop,status==1) ~ treat, data=rr)

# Time to first composite
# NA
plot(survfit(Surv(start,stop,status>0)~treat,data=subset(rr,eventno==0))
     ,lty=1:2,fun="cumhaz",xlab="time",ylab="Integrated hazard")
# Cox
coxph(Surv(start,stop,status>0) ~treat, data=subset(rr,eventno==0))

### PWP models:
coxph(Surv(start,stop,status>0) ~treat+strata(eventno), data=rr)
coxph(Surv(start,stop,status>0) ~treat*strata(eventno), data=rr)


