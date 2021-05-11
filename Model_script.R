library(lme4)
library(MuMIn)
library(emmeans)
library(RVAideMemoire)
library(DHARMa)
library(ggeffects)
library(bbmle)

model_data<- read.csv("Data/model_data.csv")#load model data

#assess impacts across season type for all four metrics

#assess frich across season (pretty weak but significant model)
m1<-lmer(data = model_data, 
          (FRich)~season + (1|Site))

summary(m1)
r.squaredGLMM(m1)

#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)

#get post hoc info
emm1 = emmeans(m1, specs = pairwise ~ season)
emm1$contrasts

#checking FEve
m1<-lmer(data = model_data, 
         FEve~season + (1|Site))
summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)


#checking FDiv
m1<-lmer(data = model_data, 
         FDiv~season + (1|Site))

summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)

#checking SR
m1<-glmer.nb(data = model_data, 
         SR~season + (1|Site))

summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)



###Assessing fluctuations across time
##Assess SR over time
#negative binomial lead to convergence issues
SR<-glmer.nb(data = trait_table_all, SR~poly(time_numeric,2) + (1|Site))
SR1<-glmer(data = trait_table_all, SR~poly(time_numeric,2) + (1|Site), family = "poisson")
SR3<-glmer(data = trait_table_all, SR~poly(time_numeric,3) + (1|Site), family = "poisson")
SR4<-glmer(data = trait_table_all, SR~poly(time_numeric,4) + (1|Site), family = "poisson")

AICctab(SR, SR1, SR3, weights =T)


summary(SR3)
r.squaredGLMM(SR3)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)


##Assess FEve
#glmer model sucks
m1<-glmer(data = trait_table_all, 
          FEve~poly(time_numeric,2) + (1|Site),family = Gamma(link= "log"))
m1<-lmer(data = trait_table_all, 
          FEve~poly(time_numeric,2) + (1|Site))
summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)


#Assess FDiv
m1<-lmer(data = trait_table_all, 
         FDiv~poly(time_numeric,3) + (1|Site))
summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)




#Assess FRich
m1<-lmer(data = trait_table_all, 
         FRich~poly(time_numeric,2) + (1|Site))
m2<-lmer(data = trait_table_all, 
         FRich~poly(time_numeric,3) + (1|Site))
AIC(m1,m2)

summary(m1)
r.squaredGLMM(m2)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)



library(ggeffects)
plot(ggpredict(m1, terms="time_numeric [all]"))

fdiv_model<-data.frame(ggpredict(m1, terms="time_numeric [all]"))
ggplot(data = NULL) +
  geom_point(data = model_data, mapping = aes(x=time_numeric, y = FDiv)) +
  geom_ribbon(data = fdiv_model, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  geom_line(data = fdiv_model, aes(x= x, y = predicted)) 

plot(ggpredict(m1, terms="time_numeric [all]"))

frich_model<-data.frame(ggpredict(m1, terms="time_numeric [all]"))
ggplot(data = NULL) +
  geom_point(data = model_data, mapping = aes(x=time_numeric, y = FRich)) +
  geom_ribbon(data = frich_model, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  geom_line(data = frich_model, aes(x= x, y = predicted)) 

