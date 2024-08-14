# Defining individual to work with
id_owl<-owls_tele$PS0009_OWL09 #Select manually which owl to work with

# name <- "NR87" #identifying which owl you want to work with

dt.plot(id_owl)
title("Sampling Intervals")

SVF <- variogram(id_owl)
level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,12 %#% "hour") # 0-12 hour window
plot(SVF,xlim=xlim,level=level)
title("zoomed in")
plot(SVF,fraction=0.65,level=level)
title("zoomed out")

m.iid <- ctmm(sigma=23 %#% "km^2")
m.ou <- ctmm(sigma=23 %#% "km^2",tau=6 %#% "day")
plot(SVF,CTMM=m.iid,fraction=0.65,level=level,col.CTMM="red")
title("Independent and identically distributed data")
plot(SVF,CTMM=m.ou,fraction=0.65,level=level,col.CTMM="purple")
title("Ornstein-Uhlenbeck movement")

m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
plot(SVF,CTMM=m.ou,level=level,col.CTMM="purple",xlim=xlim)
title("Ornstein-Uhlenbeck movement")
plot(SVF,CTMM=m.ouf,level=level,col.CTMM="blue",xlim=xlim)
title("Ornstein-Uhlenbeck-F movement")
plot(SVF,CTMM=m.ou,fraction=0.65,level=level,col.CTMM="purple")
title("Ornstein-Uhlenbeck movement")
plot(SVF,CTMM=m.ouf,fraction=0.65,level=level,col.CTMM="blue")
title("Ornstein-Uhlenbeck-F movement")


# simulate fake buffalo with the same sampling schedule
willa <- simulate(m.ouf,t=id_owl$t)
plot(willa)
title("simulation")
# now calculate and plot its variogram
SVF2 <- variogram(willa)
plot(SVF2,CTMM=m.ouf,fraction=0.65,level=level,col.CTMM="blue")
title("simulation")

