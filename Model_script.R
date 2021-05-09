library(lme4)
library(MuMIn)

model_data<- read.csv("Data/model_data.csv")
m1<-glmer(data = trait_table_all, (FRich)~season + (1|Site)+(1|obs), family = Gamma(link= "log"))
m1<-glmer(data = trait_table_all, (FEve)~season + (1|Site),family = Gamma(link= "log"))
m1<-glmer(data = trait_table_all, (FDiv)~season + (1|Site), family = Gamma(link= "inverse"))
m1<-glmer.nb(data = trait_table_all, (SR)~season + (1|Site))

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
