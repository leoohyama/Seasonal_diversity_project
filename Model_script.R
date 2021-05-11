library(lme4)
library(MuMIn)
library(emmeans)
library(RVAideMemoire)
library(DHARMa)
library(ggeffects)
library(bbmle)
library(tidyverse)

model_data<- read.csv("Data/model_data.csv")#load model data
model_data$obs<-seq(1:nrow(model_data))



#assess impacts across season type for all four metrics

#assess frich across season (pretty weak but significant model)
m1<-lmer(data = model_data, 
          log(FRich)~season + (1|Site))
hist(model_data$FRich)

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


#checking abundance
m1<-glmer.nb(data = model_data, 
             total~season + (1|Site))
summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)


###Assessing fluctuations across time
##Assess SR over time
#negative binomial lead to convergence issues
SR<-glmer.nb(data = model_data, SR~poly(time_numeric,2) + (1|Site))
SR1<-glmer(data = model_data, SR~poly(time_numeric,2) + (1|Site), family = "poisson")
SR2<-glmer(data = model_data, SR~poly(time_numeric,3) + (1|Site), family = "poisson")
SR3<-glmer(data = model_data, SR~poly(time_numeric,4) + (1|Site), family = "poisson")
nullmod<-glmer(data = model_data, SR~1 + (1|Site), family = "poisson")

AICctab(SR, SR1,SR2, SR3,nullmod, weights =T)


summary(SR1)
r.squaredGLMM(SR1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = SR1, plot = F)
plot(simulationOutput)


#top models had 3 or 4 points of inflection but were overdispersed, including observation level
#random effect led to singulairty in the models. Therefore the best alternative model was the
#2nd degree polynomial model



##Assess FEve
#glmer model sucks
m1<-glmer(data = model_data, 
          FEve~poly(time_numeric,2) + (1|Site),family = Gamma(link= "log"))

Eve1<-lmer(data = model_data, 
         FEve~time_numeric + (1|Site))
Eve2<-lmer(data = model_data, 
          FEve~poly(time_numeric,2) + (1|Site))
Eve3<-lmer(data = model_data, 
         FEve~poly(time_numeric,3) + (1|Site))
Eve4<-lmer(data = model_data, 
         FEve~poly(time_numeric,4) + (1|Site))
Null_mod<-lmer(data = model_data, 
               FEve~1 + (1|Site))

AICctab(Eve1, Eve2,Eve3, Eve4,Null_mod, weights =T)

summary(Eve3)
r.squaredGLMM(Eve3)


#check residuals
simulationOutput <- simulateResiduals(fittedModel = Eve3, plot = F)
plot(simulationOutput)


#Assess FDiv
Div1<-lmer(data = model_data, 
           FDiv~time_numeric + (1|Site))
Div2<-lmer(data = model_data, 
           FDiv~poly(time_numeric,2) + (1|Site))
Div3<-lmer(data = model_data, 
           FDiv~poly(time_numeric,3) + (1|Site))
Div4<-lmer(data = model_data, 
           FDiv~poly(time_numeric,4) + (1|Site))
Null_mod<-lmer(data = model_data, 
           FDiv~1 + (1|Site))

AICctab(Div1, Div2,Div3, Div4,Null_mod, weights =T)


#check residuals
simulationOutput <- simulateResiduals(fittedModel = Div3, plot = F)
plot(simulationOutput)



#Assess Abund
Ab1<-glmer.nb(data = model_data, 
           total~time_numeric + (1|Site))
Ab2<-glmer.nb(data = model_data, 
           total~poly(time_numeric,2) + (1|Site))
Ab3<-glmer.nb(data = model_data, 
           total~poly(time_numeric,3) + (1|Site))
Ab4<-glmer.nb(data = model_data, 
           total~poly(time_numeric,4) + (1|Site))
Null_mod<-glmer.nb(data = model_data, 
               total~1 + (1|Site))

AICctab(Ab1, Ab2,Ab3, Ab4,Null_mod, weights =T)
poly(model_data$time_numeric,4, raw = T)

#check residuals
simulationOutput <- simulateResiduals(fittedModel = Ab4, plot = F)
plot(simulationOutput)


summary(Ab4)
r.squaredGLMM(Ab4)

ab_model<-data.frame(ggpredict(Ab4, terms="time_numeric [all]"))
ggplot(data = NULL) +
  geom_point(data = model_data, mapping = aes(x=time_numeric, y = total)) +
  geom_ribbon(data = ab_model, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  geom_line(data = ab_model, aes(x= x, y = predicted)) 

ply <- contrast(Rich2, "poly")

#Assess FRich

Rich1<-lmer(data = model_data, 
           FRich~time_numeric + (1|Site))
Rich2<-lmer(data = model_data, 
            FRich~poly(time_numeric,2) + (1|Site))
Rich3<-lmer(data = model_data, 
            FRich~poly(time_numeric,3) + (1|Site))
Rich4<-lmer(data = model_data, 
            FRich~poly(time_numeric,4) + (1|Site))
Null_mod<-lmer(data = model_data, 
               FRich~1 + (1|Site))

AICctab(Rich1, Rich2,Rich3, Rich4,Null_mod, weights =T)

summary(Rich4)
r.squaredGLMM(Rich4)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = Rich4, plot = F)
plot(simulationOutput)

#plot model
plot(ggpredict(Rich4, terms="time_numeric [all]"))

Rich4_model<-data.frame(ggpredict(Rich4, terms="time_numeric [all]"))
ggplot(data = NULL) +
  geom_point(data = model_data, mapping = aes(x=time_numeric, y = FRich)) +
  geom_ribbon(data = Rich4_model, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  geom_line(data = Rich4_model, aes(x= x, y = predicted)) 







