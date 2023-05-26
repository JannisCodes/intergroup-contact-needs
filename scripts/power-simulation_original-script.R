### Power Analysis ###
### Note:
### Example marginalization
### Becareful with running this. The simulations can only take one CPU core and take forever.

# prepare for re-run:
  cat("\014") # clear console
  rm(list = ls()) # clear workspace
  gc() # run garbage collector
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory
  p_load(simr, data.table, lme4, lmerTest)
  load("df.btw.Rdata")
  
# set random seed for reporducibility
  set.seed(42)

# prepare data for regression (i.e, put into long format and change variable types)
  marginalization = melt(df.btw[!is.na(df.btw$assimilation.post),c("marginalization.pre","marginalization.post",
                                                                   "CtContactNL_c","AvKeyNeedInt_c",
                                                                   "ExternalReference")],
                         id=c("ExternalReference","CtContactNL_c","AvKeyNeedInt_c"))
  marginalization$time = ifelse(grepl(".pre", marginalization$variable), "pre", "post")
  marginalization$variable = gsub(".pre$|.post$", "", marginalization$variable)
  marginalization$time = factor(marginalization$time)
  marginalization$time = factor(marginalization$time,levels(marginalization$time)[c(2,1)])
  marginalization$id = rep(seq(1:20), 2)

# run regresion 
  lm.marg <- lmer(value~
                    AvKeyNeedInt_c*time+
                    CtContactNL_c*time+
                    (1|id), 
                  data=marginalization)
  summary(lm.marg)
  fixef(lm.marg)

# Post-hoc power analysis for single effects (example: key need)
  powerSim(lm.marg, fixed("AvKeyNeed"), nsim=1000)


# Simulate different N.s for all fixed effects
  # extrapolate data:
    sim.N <- extend(lm.marg, along="id", n=500)
  
  # Run and plot simulations
    plt.pwr.contact = powerCurve(sim.N, fixed("CtContactNL_c", "t"), 
                                 along="id", 
                                 breaks = seq(10,100,10), 
                                 nsim=1000)
    plot(plt.pwr.contact)
    
    plt.pwr.time = powerCurve(sim.N, fixed("time"), 
                              along="id", 
                              breaks = seq(10,100,10), 
                              nsim=1000)
    plot(plt.pwr.time)
    
    plt.pwr.Need = powerCurve(sim.N, fixed("AvKeyNeedInt_c", "t"), 
                              along="id", 
                              breaks = seq(10,100,10), 
                              nsim=1000)
    plot(plt.pwr.Need)
    
    plt.pwr.timeCont = powerCurve(sim.N, fixed("time:CtContactNL_c"), 
                                 along="id", 
                                 breaks = seq(10,500,50), 
                                 nsim=1000)
    plot(plt.pwr.timeCont)
    
    plt.pwr.timeNeed = powerCurve(sim.N, fixed("AvKeyNeedInt_c:time"), 
                                  along="id", 
                                  breaks = seq(10,50,5), 
                                  nsim=1000)
    plot(plt.pwr.timeNeed)

# Simulate different number of within measurements
    # extrapolate data:
      sim.within <- extend(lm.marg, within="id", n=12)
      
    # Run and plot simulations
      plt.pwr.t.time = powerCurve(sim.within, fixed("time"), 
                                within="id", 
                                breaks = seq(2,12,1), 
                                nsim=1000)
      plot(plt.pwr.t.time)
      
      plt.pwr.t.contact = powerCurve(sim.within, fixed("CtContactNL_c", "t"), 
                                     within="id", 
                                     breaks = seq(2,12,1), 
                                     nsim=1000)
      plot(plt.pwr.t.contact)
      
      plt.pwr.t.Need = powerCurve(sim.within, fixed("AvKeyNeedInt_c", "t"), 
                                     within="id", 
                                     breaks = seq(2,12,1), 
                                     nsim=1000)
      plot(plt.pwr.t.Need)
      
      plt.pwr.t.timeCont = powerCurve(sim.within, fixed("time:CtContactNL_c"), 
                                within="id", 
                                breaks = seq(2,12,1), 
                                nsim=1000)
      plot(plt.pwr.t.timeCont)
      
      
      plt.pwr.t.timeNeed = powerCurve(sim.within, fixed("AvKeyNeedInt_c:time"), 
                                     within="id", 
                                     breaks = seq(2,12,1), 
                                     nsim=1000)
      plot(plt.pwr.t.timeNeed)

# save simulations, because this takes for freakin ever.      
  save(plt.pwr.contact,plt.pwr.Need,plt.pwr.t.contact,plt.pwr.t.Need,
       plt.pwr.t.time,plt.pwr.t.timeCont,plt.pwr.t.timeNeed,plt.pwr.time,
       plt.pwr.timeCont,plt.pwr.timeNeed, 
       file = "PowerSimPlots.RData")

  #load("PowerSimPlots.RData")
  