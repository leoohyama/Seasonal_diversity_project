#Trait analysis Sociometry
library(glmmTMB)
library(MASS)
library(lme4)
library(MuMIn)
library(DHARMa)
library(bbmle)
library(ggeffects)
library(tidyverse)

traits<-read.csv("Data/species_traits_florida.csv")
traits2<-traits
traits<-traits %>% filter(!species == "Solenopsis.invicta")
#make genus column
traits$genus = str_extract(traits$species, "^[\\w]+")
traits$obs<-seq(1:nrow(traits))

m1<-glmer.nb(data = traits, Colony_size~WL + (1|genus))
WL1<-lmer(data = traits, log(Colony_size)~WL + (1|genus))
EL1<-lmer(data = traits, log(Colony_size)~EL + (1|genus))
HW1<-lmer(data = traits, log(Colony_size)~HW + (1|genus))
ML1<-lmer(data = traits, log(Colony_size)~ML + (1|genus))
NULLMOD<-lmer(data = traits, log(Colony_size)~1 + (1|genus))



AICctab(WL1, EL1, HW1,ML1,NULLMOD, weights = T)

summary(EL1)
r.squaredGLMM(EL1)

simulationOutput <- simulateResiduals(fittedModel = m1, plot = F)
plot(simulationOutput)

#plot model
plot(ggpredict(EL1, terms="EL[all]"))

Eyemodel<-data.frame(ggpredict(EL1, terms="EL[all]"), back.transform =FALSE)
ggplot(data = NULL) +
  geom_point(data = traits, mapping = aes(x=EL, y = log(Colony_size))) +
  geom_ribbon(data = Eyemodel, aes(x = x, ymin = log(conf.low), ymax = log(conf.high)), alpha = 0.3)+
  geom_line(data = Eyemodel, aes(x= x, y = log(predicted))) 





m2<-glmer.nb(data = traits, Colony_size~EL + (1|genus))
m2<-glmer.nb(data = traits, Colony_size~HW + (1|genus))
m2<-glmer.nb(data = traits, Colony_size~ML + (1|genus))
summary(m2)
r.squaredGLMM(m2)


m2<-glmer.nb(data = traits, Colony_size~ ML + WL +(1|genus))
