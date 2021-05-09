library(lme4)
library(MuMIn)
library(emmeans)
library(RVAideMemoire)
library(DHARMa)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

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

m1<-glmer(data = trait_table_all, 
          (FEve)~season + (1|Site),family = Gamma(link= "log"))
m1<-glmer(data = trait_table_all, 
          (FDiv)~season + (1|Site), family = Gamma(link= "inverse"))
m1<-glmer.nb(data = trait_table_all, 
             (SR)~season + (1|Site))

hist(trait_table_all$FDiv)
m1<-lmer(data = trait_table_all, FEve~poly(time_numeric,2) + (1|Site))
m1<-lmer(data = trait_table_all, FDiv~poly(time_numeric,2) + (1|Site))

m1.1<-glmer(data = trait_table_all, FEve~poly(time_numeric,2) + (1|Site),family = Gamma(link= "log"))
summary(m1)

AIC(m1,m1.1)
r.squaredGLMM(m1)

##SR
m2<-glmer.nb(data = trait_table_all, SR~poly(time_numeric,2) + (1|Site))
summary(m2)
r.squaredGLMM(m2)


#FDveg
m1.1<-glmer(data = trait_table_all, FRich~poly(time_numeric,2) + (1|Site),family = Gamma(link= "log"))
summary(m1.1)
plot(ggpredict(m1.1, terms="time_numeric [all]"))


library(ggeffects)
plot(ggpredict(m1, terms="time_numeric [all]"))

m1<-glm(data = trait_table_all, (FDiv)~season, family= Gamma(link = "log"))
summary(m1)



plot(trait_table_all$SR,trait_table_all$FEve)
