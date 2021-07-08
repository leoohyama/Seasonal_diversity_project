library(MuMIn)
library(emmeans)
library(RVAideMemoire)
library(DHARMa)
library(ggeffects)
library(bbmle)
library(tidyverse)
library(glmmTMB)
library(MASS)

model_data<- read.csv("Data/model_data.csv")#load model data
model_data$obs<-seq(1:nrow(model_data))
model_data$SES_fdiv_t<-sqrt(max(model_data$SES_FDiv+1,na.rm = T) - model_data$SES_FDiv) 
model_data<-model_data %>% drop_na(FRich)
hist(model_data$SES_FRich)

#assess impacts across season type for all four metrics



#assess frich across season (pretty weak but significant model)
m1<-glmmTMB(data = model_data, 
          (FRich)~season +  (1|Site) + (1|time_numeric))

m2<-glmmTMB(data = model_data, 
          (SES_FRich)~season + (1|time_numeric)+ (1|Site))

summary(m1)
summary(m2)

r.squaredGLMM(m1)
r.squaredGLMM(m2)

#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
plot(simulationOutput)

#get post hoc info
emm1 = emmeans(m1, specs = pairwise ~ season)
emm1$contrasts




#checking FEve
#assess feve across season (pretty weak but significant model)
m1<-glmmTMB(data = model_data, 
         (FEve)~season + (1|time_numeric)+ (1|Site))

#this model had singulairty
m2<-glmmTMB(data = model_data, 
         (SES_FEve)~season + (1|time_numeric) + (1|Site))

#So we drop the random effect term for site as it is the least concern
m2<-lmer(data = model_data, 
         (SES_FEve)~season + (1|time_numeric))

summary(m1)
r.squaredGLMM(m1)


summary(m2)
r.squaredGLMM(m2)



#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
plot(simulationOutput)




#checking FDiv
#assess FDiv across season 
m1<-glmmTMB(data = model_data, 
         (FDiv)~season + (1|time_numeric) + (1|Site))

m2<-glmmTMB(data = model_data, 
         (SES_FDiv)~season + (1|time_numeric) + (1|Site))


summary(m1)
r.squaredGLMM(m1)

summary(m2)
r.squaredGLMM(m2)



#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel = m2, plot = F)
plot(simulationOutput)


#checking abundance
m1<-glmmTMB(data = model_data, 
             total~season +  (1|Site)+ (1|time_numeric), family = nbinom2)


summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)
emm1 = emmeans(m1, specs = pairwise ~ season)
emm1$contrasts
summary(emm1, type = "response")


#checking SR
m1<-glmmTMB(data = model_data, 
             SR~season + (1|Site), family = nbinom1)

m2<-glmmTMB(data = model_data, 
            SR~season + (1|Site), family = nbinom2)

m3<-glmmTMB(data = model_data, 
            SR~season + (1|Site), family = poisson)

summary(m1)
r.squaredGLMM(m1)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)

#get post hoc info
emm1 = emmeans(m1, specs = pairwise ~ season)
emm1$contrasts
summary(emm1, type = "response")


###Assessing fluctuations across time
##Assess SR over time
#negative binomial lead to convergence issues
SR<-glmmTMB(data = model_data, SR~poly(time_numeric,2) + (1|Site)+ (1|obs), family = poisson)
SR1<-glmmTMB(data = model_data, SR~poly(time_numeric,2) + (1|Site)+ (1|obs), family = poisson)
SR2<-glmmTMB(data = model_data, SR~poly(time_numeric,3) +  (1|Site)+ (1|obs), family = poisson)
SR3<-glmmTMB(data = model_data, SR~poly(time_numeric,4) + (1|Site) + (1|obs), family = poisson)
nullmod<-glmmTMB(data = model_data, SR~1 + (1|Site)+ (1|obs), family = poisson)

AICctab(SR, SR1,SR2, SR3,nullmod, weights =T)


summary(SR2)
r.squaredGLMM(SR2)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = SR2, plot = F)
plot(simulationOutput)




##Assess FEve
#glmer model sucks


Eve1<-glmmTMB(data = model_data, 
         FEve~time_numeric + (1|Site))
Eve2<-glmmTMB(data = model_data, 
          FEve~poly(time_numeric,2) + (1|Site))
Eve3<-glmmTMB(data = model_data, 
         FEve~poly(time_numeric,3) + (1|Site))
Eve4<-glmmTMB(data = model_data, 
         FEve~poly(time_numeric,4) + (1|Site))
Null_mod<-glmmTMB(data = model_data, 
               FEve~1 + (1|Site))

AICctab(Eve1, Eve2,Eve3, Eve4,Null_mod, weights =T)

summary(Eve3)
r.squaredGLMM(Eve3)


#check residuals
simulationOutput <- simulateResiduals(fittedModel = Eve3, plot = F)
plot(simulationOutput)


#Assess FDiv

Div1<-glmmTMB(data = model_data, 
           FDiv~time_numeric + (1|Site))
Div2<-glmmTMB(data = model_data, 
           FDiv~poly(time_numeric,2) + (1|Site))
Div3<-glmmTMB(data = model_data, 
           FDiv~poly(time_numeric,3) + (1|Site))
Div4<-glmmTMB(data = model_data, 
           FDiv~poly(time_numeric,4) + (1|Site))
Null_mod<-glmmTMB(data = model_data, 
           FDiv~1 + (1|Site))

AICctab(Div1, Div2,Div3, Div4,Null_mod, weights =T)


#check residuals
simulationOutput <- simulateResiduals(fittedModel = Div4, plot = F)
plot(simulationOutput)



#Assess Abund
Ab1<-glmmTMB(data = model_data, 
           total~time_numeric + (1|Site), family = poisson)
Ab2<-glmmTMB(data = model_data, 
           total~poly(time_numeric,2) + (1|Site), family = poisson)
Ab3<-glmmTMB(data = model_data, 
           total~poly(time_numeric,3) + (1|Site), family = poisson)
Ab4<-glmmTMB(data = model_data, 
           total~poly(time_numeric,4) + (1|Site) + (1|obs), family = poisson)
Null_mod<-glmmTMB(data = model_data, 
               total~1 + (1|Site), family = poisson)

AICctab(Ab1, Ab2,Ab3, Ab4,Null_mod, weights =T)

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

Rich1<-glmmTMB(data = model_data, 
           FRich~time_numeric + (1|Site))
Rich2<-glmmTMB(data = model_data, 
            FRich~poly(time_numeric,2) + (1|Site))
Rich3<-glmmTMB(data = model_data, 
            FRich~poly(time_numeric,3) + (1|Site))
Rich4<-glmmTMB(data = model_data, 
            FRich~poly(time_numeric,4) + (1|Site))
Null_mod<-glmmTMB(data = model_data, 
               FRich~1 + (1|Site))

AICctab(Rich1, Rich2,Rich3, Rich4,Null_mod, weights =T)

summary(Rich4)
r.squaredGLMM(Rich4)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = Rich3, plot = F)
plot(simulationOutput)

#plot model
plot(ggpredict(Rich4, terms="time_numeric [all]"))

Rich4_model<-data.frame(ggpredict(Rich4, terms="time_numeric [all]"))
ggplot(data = NULL) +
  geom_point(data = model_data, mapping = aes(x=time_numeric, y = FRich)) +
  geom_ribbon(data = Rich4_model, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.3)+
  geom_line(data = Rich4_model, aes(x= x, y = predicted)) 



#Assess Jost



Jost1<-glmer(data = model_data, 
            jost~time_numeric + (1|Site), family = Gamma(link = "log"))

Jost1<-lmer(data = model_data, 
            jost~time_numeric + (1|Site))
Jost2<-lmer(data = model_data, 
            jost~poly(time_numeric,2) + (1|Site))
Jost3<-lmer(data = model_data, 
            jost~poly(time_numeric,3) + (1|Site))
Jost4<-lmer(data = model_data, 
            jost~poly(time_numeric,4) + (1|Site))
Null_mod<-lmer(data = model_data, 
               jost~1 + (1|Site))

AICctab(Rich1, Rich2,Rich3, Rich4,Null_mod, weights =T)

summary(Rich4)
r.squaredGLMM(Rich3)
#check residuals
simulationOutput <- simulateResiduals(fittedModel = Rich4, plot = F)
plot(simulationOutput)



###let's estimate when relative to the average

model_data %>%
  group_by(time_code,Site) %>%
  summarise(diff_FEve = mean(FEve)-FEve)

mean(model_data$total)

model_data %>%
  mutate(mean_SR = mean(SR, na.rm = T),
         mean_total = mean(total, na.rm = T),
         mean_FR = mean(FRich, na.rm = T),
         mean_FD = mean(FDiv, na.rm = T),
         mean_FE = mean(FEve, na.rm = T),
         diff_SR = SR-mean_SR,
         diff_total = total-mean_total,
         diff_Frich = FRich-mean_FR,
         diff_Fdiv = FDiv-mean_FD,
         diff_Feve = FEve-mean_FE) %>%
  group_by(season, time_code) %>%
  summarise(mean_diff_sr = mean(diff_SR, na.rm=T),
            mean_diff_tot = mean(diff_total, na.rm=T),
            mean_diff_frich = mean(diff_Frich, na.rm=T),
            mean_diff_fdiv = mean(diff_Fdiv, na.rm=T),
            mean_diff_feve = mean(diff_Feve, na.rm=T)
            ) %>%
  pivot_longer(cols=c("mean_diff_sr", "mean_diff_tot", "mean_diff_frich",
                        "mean_diff_fdiv", "mean_diff_feve"), names_to = "type", values_to = "values") %>%
  ggplot(.) +
  geom_col(aes(x = time_code, y = values)) +
  facet_wrap(~type, scales = "free_y")



#by season and timecode combining all plots

model_data %>%
  group_by(season, time_code) %>%
  summarise(sum_SR = mean(SR, na.rm = T),
            sum_total = mean(total, na.rm = T),
            sum_FR = mean(FRich, na.rm = T),
            sum_FD = mean(FDiv, na.rm = T),
            sum_FE = mean(FEve, na.rm = T)) 




#SR Frich

model_data$time_numericf<-as.factor(model_data$time_numeric)
m1<-glmmTMB(data = model_data, SR~FRich +(1|time_numericf),family = poisson)
m1<-glm.nb(data = model_data, SR~FRich * time_numericf)

summary(m1)
r.squaredGLMM(m1)

timeplot<-data.frame(ggpredict(m1, terms=c("FRich","time_numericf"), type = "re"))
ggplot(timeplot) +
  geom_line(data = timeplot,aes(x =x, y =  predicted, color = group, group = group)) 
