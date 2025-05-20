library(mets)
rr<-read.csv("rr.csv")

## Marginal
# 6. Cook-Lawless: Mean no of events
plot(recurrentMarginal(
  Event(start,stop,status)~strata(treat)+cluster(id),
  data=rr,cause=1,death.code=2),col=1,ylim=c(0,0.20),
  main="Cook-Lawless")

# Wrong CMF censor at death
lines(survfit(Surv(start,stop,status==1)~treat,data=rr),
      fun="cumhaz",lty=1:2, col=2,xlab="time",ylab="Mean events")

# 7. Ghosh-Lin
summary(fitGL<-recreg(Event(start,stop,status)~treat+cluster(id),
                      data=rr,cause=1,death.code=2,cens.code=0))

# 8. LWYY
coxph(Surv(start,stop,status==1)~treat+cluster(id),data=rr)

# 9. Composite recurrent event or death
# Mao-Lin
summary(fitGL<-recreg(Event(start,stop,status)~treat+cluster(id),
              data=rr, cause=c(1,2),death.code=2,cens.code=0))

# Cox for death
coxph(Surv(start,stop,status==2)~treat,data=rr)

# 10. Time to first event: Competing risks
# AJ
t2f<-survfit(Surv(start,stop,factor(status))~treat,
             data=subset(rr,eventno==0))
# 1-KM
t2fwrong<-survfit(Surv(start,stop,status==1)~treat,
                  data=subset(rr,eventno==0))
# Plot AJ
plot(t2f[,2],lty=1:2,main="Event: CIF",xlab="time",ylab="Prob")
# add 1-KM
lines(t2fwrong,lty=1:2,fun="event",col=2)

