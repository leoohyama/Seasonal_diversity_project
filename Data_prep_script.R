library(vegan)
library(tidyverse)
library(zoo)
library(viridis)


#read ecological data
df<-read.csv("~/Desktop/microbeta.csv")
df_c<-df %>%
  filter(Treatment == "C")
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ")
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))

#read morphological trait data
dft<-read.csv("~/Desktop/micro_traits.csv")
#clean up species
dft$Genus.Species<-str_extract(dft$Genus.Species, "^\\w+.\\w+")
#now remove any minors or queens

dft2<-dft %>% filter(!Caste %in% c("minor","queen"))


#get species list to make trait dta.frame
species_list<-df_c %>% dplyr::select(! c(Year, Site, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code"), names_to = "Species", values_to="abund") %>%
  dplyr::select(Species)
species_list<-data.frame(species= unique(species_list$Species))


species_list1<-left_join(species_list, dft2, by = c("species" = "Genus.Species") )
species_list1<-species_list1%>%mutate_at(vars(matches("length")),funs(as.numeric))


#here we calculate the traits and then divide by body size to correct for body size
species_list1$WL<-(species_list1$Weber.s.length_1+species_list1$Weber.s.length_2)/2
species_list1$HW<-((species_list1$Head.width_1+species_list1$Head.width_2)/2)/species_list1$WL
species_list1$EL<-((species_list1$Eye.length_1+species_list1$Eye.length_2)/2)/species_list1$WL
species_list1$ML<-((species_list1$Mandible.length_1+species_list1$Mandible.length_2)/2)/species_list1$WL

species_traits<-species_list1 %>% dplyr::select(species, WL, HW,EL,ML)
setwd("~/Desktop/")
#write.csv(species_traits, "species_traits_florida.csv") write out the trait data
#read in new trait data with colony size
species_traits<-read.csv("species_traits_florida.csv")
species_traits<-species_traits %>% dplyr::select(-Notes) %>%
  mutate(log_colony = log(Colony_size))
dplyr::select(-Colony_size)


#we need species list for everytime code and site
#then we will try calculating FD for each timecode
#this will be how we set up a pipeline to calculate FD for every site for every time_code
#this is used to weight the trait values by abundance

by_sp<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  ungroup() %>%
  filter(time_code == "2018-04-01") %>%
  filter(!total == 0) %>%
  group_by(Species) %>%
  mutate(total_n = sum(total)) %>%
  filter(!duplicated(Species)) %>%
  group_by(Species) %>%
  arrange(Species) %>%
  ungroup()


#now we have the traits
species_list<-by_sp$Species
time1_traits<-by_sp %>% dplyr::select(log_colony,WL,HW,EL,ML)
rownames(time1_traits)<-species_list
#now we have the comm data

time1_comm <- df_c %>% filter(time_code == "2018-04-01") 
site_names<-time1_comm$Site
time1_comm<-time1_comm %>%
  dplyr::select(-c(Year, Site, Month, Treatment,time_code))
time1_comm = time1_comm[,colSums(time1_comm) > 0]
rownames(time1_comm)<-site_names

#rescale trait data
trait_data = scale(time1_traits)
library(vegan)
trait_matrix = vegan::vegdist(trait_data,"euclidean")
branch_lengths<-hclust(trait_matrix,"average")
branch_lengths
library(picante)
trait_tree = ape::as.phylo(branch_lengths)
plot(trait_tree)
str(time1_comm)
FD<-picante::pd(time1_comm, trait_tree)
vegan_tbl<-data.frame(Site = as.vector(rownames(FD)), SR = FD$SR, FD = FD$PD)

#calculate functional richness
trial<-dbFD(time1_traits, time1_comm)
trial$FRic
trial$FEve
trial$FDiv

#program to get all fmetrics

trait_table_all<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site) %>%
  summarise(total = sum(abund))

time_stamp_list<-unique(df_c$time_code)
trait_table_all$FD_veg<-NA
trait_table_all$SR<-NA
trait_table_all$FRich<-NA
trait_table_all$FEve<-NA
trait_table_all$FDiv<-NA


for(i in 1:length(time_stamp_list)){
  
  by_sp<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
    pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
    mutate(genus = str_extract(Species, "^[\\w]+")) %>%
    group_by(time_code, Site,Species) %>%
    summarise(total = sum(abund)) %>%
    left_join(.,species_traits, by = c("Species"="species")) %>%
    ungroup() %>%
    filter(time_code == time_stamp_list[i]) %>%
    filter(!total == 0) %>%
    group_by(Species) %>%
    mutate(total_n = sum(total)) %>%
    filter(!duplicated(Species)) %>%
    group_by(Species) %>%
    arrange(Species) %>%
    ungroup()
  #now we have the traits
  species_list<-by_sp$Species
  time1_traits<-by_sp %>% dplyr::select(log_colony,WL,HW,EL,ML)
  rownames(time1_traits)<-species_list
  #now we have the comm data
  time1_comm <- df_c %>% filter(time_code == time_stamp_list[i]) 
  site_names<-time1_comm$Site
  time1_comm<-time1_comm %>%
    dplyr::select(-c(Year,Site, Month, Treatment,time_code)) 
  rownames(time1_comm)<-site_names
  time1_comm = time1_comm[,colSums(time1_comm) > 0]
  time1_comm = time1_comm[rowSums(time1_comm)> 0, ]
  #remove species removed from comm from traits
  
  
  #rescale trait data
  trait_data = scale(time1_traits)
  library(vegan)
  trait_matrix = vegan::vegdist(trait_data,"euclidean")
  branch_lengths<-hclust(trait_matrix,"average")
  branch_lengths
  library(picante)
  trait_tree = ape::as.phylo(branch_lengths)
  plot(trait_tree)
  str(time1_comm)
  FD<-picante::pd(time1_comm, trait_tree)
  FD
  vegan_tbl<-data.frame(Site = as.vector(rownames(FD)), SR = FD$SR, FD = FD$PD)
  
  for(j in 1:nrow(vegan_tbl)){
    trait_table_all$SR[trait_table_all$Site == paste(vegan_tbl$Site[j]) &
                         trait_table_all$time_code == time_stamp_list[i]]<-vegan_tbl$SR[j]
    trait_table_all$FD_veg[trait_table_all$Site == paste(vegan_tbl$Site[j]) &
                             trait_table_all$time_code == time_stamp_list[i]]<-vegan_tbl$FD[j]
  }
  
  
  #calculate functional richness
  trial<-dbFD(time1_traits, time1_comm)
  func_rich_metrics<-data.frame(Site = rownames(data.frame(trial)), FRic = trial$FRic,
                                FEve = trial$FEve, FDiv = trial$FDiv)
  
  for(k in 1:nrow(func_ruch_metrics)){
    trait_table_all$FRich[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                            trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FRic[k]
    trait_table_all$FEve[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                           trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FEve[k]
    trait_table_all$FDiv[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                           trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FDiv[k]
  }
  
  print(paste(i, "has finished", " "))
}



trait_table_all%>% group_by(time_code) %>%
  summarise(mean_frich = mean(FRich, na.rm = T)) %>%
  ggplot(., aes(x= (time_code), y = mean_frich)) + 
  geom_line()
trait_table_all%>%group_by(time_code) %>%
  summarise(total = n(), mean_frich = mean(FRich, na.rm = T), 
            Frich_std_eror=(sd(FRich, na.rm = T)/sqrt(total)),
            SR = mean(SR, na.rm=T),
            SR_std_eror=(sd(SR, na.rm = T)/sqrt(total))) %>%
  ggplot(.) + 
  geom_line(aes(x= (time_code), y = mean_frich), color = "black") +
  geom_errorbar(aes(x = time_code, ymin = mean_frich-std_eror,
                    ymax = mean_frich+std_eror), width = 0.5)+
  geom_point(aes(x= (time_code), y = mean_frich), color = "black",
             pch = 21, size = 3, fill = "white",stroke = 2.5) +
  labs(x = "Time", y = "Average Functional Richness") +
  theme_bw() 

####facet graph comparinng all metrics
jesuschrist<-trait_table_all %>%
  group_by(time_code) %>%
  summarise(total = n(), mean_FRich = mean(FRich, na.rm = T), 
            std.error.FRich=(sd(FRich, na.rm = T)/sqrt(total)),
            mean_SR = mean(SR, na.rm=T),
            std.error.SR = sd(SR, na.rm=T)/sqrt(total),
            mean_Fdiv = mean(FDiv, na.rm =T),
            std.error.Fdiv = sd(FDiv, na.rm = T)/sqrt(total),
            mean_FEve = mean(FEve, na.rm =T),
            std.error.FEve = sd(FEve, na.rm = T)/sqrt(total)) %>%
  pivot_longer(., cols = c( "mean_SR", "mean_FRich","mean_Fdiv","mean_FEve"
  ), names_to = c("type"),
  values_to= "value") 

jesuschrist$std.error<-NA
for(i in 1:nrow(jesuschrist)){
  if(grepl("FRich", jesuschrist$type[i])){
    jesuschrist$std.error[i]<-jesuschrist$std.error.FRich[i]
  }else{
    if(grepl("SR", jesuschrist$type[i])){jesuschrist$std.error[i]<-jesuschrist$std.error.SR[i]}else{
      if(grepl("Fdiv", jesuschrist$type[i])){
        jesuschrist$std.error[i]<-jesuschrist$std.error.Fdiv[i]
      }else{
        jesuschrist$std.error[i]<-jesuschrist$std.error.FEve[i]
      }
    }
  }
}  

ggplot(data = jesuschrist) + 
  geom_line(aes(x= (time_code), y = value)) +
  geom_errorbar(aes(x = time_code, ymin = value-std.error,
                    ymax = value+std.error), width = 0.5)+
  geom_point(aes(x= (time_code), y = value), 
             pch = 21, size = 3, fill = "white",stroke = 2.5) +
  facet_wrap(~type,scales = "free_y", nrow = 2)



##quick models to assess functional metrics
#June to October Wet Season
month_list = df_c %>% dplyr::select(time_code, Month)
month_list =unique(month_list[,c("time_code", "Month")])


trait_table_all$time_numeric<-NA
trait_table_all$Month<-NA
for(i in 1:nrow(trait_table_all)){
  trait_table_all$Month[i]<-month_list$Month[month_list$time_code == trait_table_all$time_code[i]]
}
trait_table_all$season<-NA
for(i in 1:nrow(trait_table_all)){
  if(trait_table_all$Month[i] %in% c("May","Jun","June", "July", "August", "September")){
    trait_table_all$season[i]<-"WET"
  }else{
    trait_table_all$season[i]<-"DRY"
  }
  
}



for(i in 1:nrow(trait_table_all)){
  if(trait_table_all$time_code[i] == "2018-02-01"){trait_table_all$time_numeric[i]<-1}
  else{if(trait_table_all$time_code[i] == "2018-03-01"){trait_table_all$time_numeric[i]<-2}
    else{if(trait_table_all$time_code[i] == "2018-04-01"){trait_table_all$time_numeric[i]<-3}
      else{if(trait_table_all$time_code[i] == "2018-05-01"){trait_table_all$time_numeric[i]<-4}
        else{if(trait_table_all$time_code[i] == "2018-06-01"){trait_table_all$time_numeric[i]<-5}
          else{if(trait_table_all$time_code[i] == "2018-07-01"){trait_table_all$time_numeric[i]<-6}
            else{if(trait_table_all$time_code[i] == "2018-08-01"){trait_table_all$time_numeric[i]<-7}
              else{if(trait_table_all$time_code[i] == "2018-09-01"){trait_table_all$time_numeric[i]<-8}
                else{if(trait_table_all$time_code[i] == "2018-11-01"){trait_table_all$time_numeric[i]<-9}
                  else{if(trait_table_all$time_code[i] == "2019-02-01"){trait_table_all$time_numeric[i]<-10}
                    else{if(trait_table_all$time_code[i] == "2019-03-01"){trait_table_all$time_numeric[i]<-11}
                      else{if(trait_table_all$time_code[i] == "2019-04-01"){trait_table_all$time_numeric[i]<-12}
                        else{if(trait_table_all$time_code[i] == "2019-05-01"){trait_table_all$time_numeric[i]<-13}
                          else{if(trait_table_all$time_code[i] == "2019-06-01"){trait_table_all$time_numeric[i]<-14}
                            else{print(i)}}}}}}}}}}}}}}
}

hist(trait_table_all$FEve)
obs<-seq(nrow(trait_table_all)) #observation level random effect
library(lme4)
library(MuMIn)

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































#############################
#############################
###treated plots#############################
#############################
#############################
df_t<-df %>%
  filter(Treatment == "T")
df_t$time_code<-paste(df_t$Month,df_t$Year,  sep = " ")
df_t$time_code<-as.Date(as.yearmon(df_t$time_code))


#get species list to make trait dta.frame
species_list<-df_t %>% dplyr::select(! c(Year, Site, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code"), names_to = "Species", values_to="abund") %>%
  dplyr::select(Species)
species_list<-data.frame(species= unique(species_list$Species))


species_list1<-left_join(species_list, dft2, by = c("species" = "Genus.Species") )
species_list1<-species_list1%>%mutate_at(vars(matches("length")),funs(as.numeric))


#here we calculate the traits and then divide by body size to correct for body size
species_list1$WL<-(species_list1$Weber.s.length_1+species_list1$Weber.s.length_2)/2
species_list1$HW<-((species_list1$Head.width_1+species_list1$Head.width_2)/2)/species_list1$WL
species_list1$EL<-((species_list1$Eye.length_1+species_list1$Eye.length_2)/2)/species_list1$WL
species_list1$ML<-((species_list1$Mandible.length_1+species_list1$Mandible.length_2)/2)/species_list1$WL

species_traits<-species_list1 %>% dplyr::select(species, WL, HW,EL,ML)

#we need species list for everytime code and site
#then we will try calculating FD for each timecode
#this will be how we set up a pipeline to calculate FD for every site for every time_code
#this is used to weight the trait values by abundance

by_sp<-df_t %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  ungroup() %>%
  filter(time_code == "2018-04-01") %>%
  filter(!total == 0) %>%
  group_by(Species) %>%
  mutate(total_n = sum(total)) %>%
  filter(!duplicated(Species)) %>%
  group_by(Species) %>%
  arrange(Species) %>%
  ungroup()


#now we have the traits
species_list<-by_sp$Species
time1_traits<-by_sp %>% dplyr::select(WL,HW,EL,ML)
rownames(time1_traits)<-species_list
#now we have the comm data

time1_comm <- df_c %>% filter(time_code == "2018-04-01") 
site_names<-time1_comm$Site
time1_comm<-time1_comm %>%
  dplyr::select(-c(Year, Site, Month, Treatment,time_code))
time1_comm = time1_comm[,colSums(time1_comm) > 0]
rownames(time1_comm)<-site_names

#rescale trait data
trait_data = scale(time1_traits)
library(vegan)
trait_matrix = vegan::vegdist(trait_data,"euclidean")
branch_lengths<-hclust(trait_matrix,"average")
branch_lengths
library(picante)
trait_tree = ape::as.phylo(branch_lengths)
plot(trait_tree)
str(time1_comm)
FD<-picante::pd(time1_comm, trait_tree)
vegan_tbl<-data.frame(Site = as.vector(rownames(FD)), SR = FD$SR, FD = FD$PD)

#calculate functional richness
trial<-dbFD(time1_traits, time1_comm)
trial$FRic
trial$FEve
trial$FDiv

#program to get all fmetrics

trait_table_all_t<-df_t %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site) %>%
  summarise(total = sum(abund))

time_stamp_list<-unique(df_t$time_code)
trait_table_all_t$FD_veg<-NA
trait_table_all_t$SR<-NA
trait_table_all_t$FRich<-NA
trait_table_all_t$FEve<-NA
trait_table_all_t$FDiv<-NA


for(i in 1:length(time_stamp_list)){
  
  by_sp<-df_t %>% dplyr::select(! c(Year, Month, Treatment)) %>%
    pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
    mutate(genus = str_extract(Species, "^[\\w]+")) %>%
    group_by(time_code, Site,Species) %>%
    summarise(total = sum(abund)) %>%
    left_join(.,species_traits, by = c("Species"="species")) %>%
    ungroup() %>%
    filter(time_code == time_stamp_list[i]) %>%
    filter(!total == 0) %>%
    group_by(Species) %>%
    mutate(total_n = sum(total)) %>%
    filter(!duplicated(Species)) %>%
    group_by(Species) %>%
    arrange(Species) %>%
    ungroup()
  #now we have the traits
  species_list<-by_sp$Species
  time1_traits<-by_sp %>% dplyr::select(WL,HW,EL,ML)
  rownames(time1_traits)<-species_list
  #now we have the comm data
  time1_comm <- df_t %>% filter(time_code == time_stamp_list[i]) 
  site_names<-time1_comm$Site
  time1_comm<-time1_comm %>%
    dplyr::select(-c(Year,Site, Month, Treatment,time_code)) 
  rownames(time1_comm)<-site_names
  time1_comm = time1_comm[,colSums(time1_comm) > 0]
  time1_comm = time1_comm[rowSums(time1_comm)> 0, ]
  #remove species removed from comm from traits
  
  
  #rescale trait data
  trait_data = scale(time1_traits)
  library(vegan)
  trait_matrix = vegan::vegdist(trait_data,"euclidean")
  branch_lengths<-hclust(trait_matrix,"average")
  branch_lengths
  library(picante)
  trait_tree = ape::as.phylo(branch_lengths)
  plot(trait_tree)
  str(time1_comm)
  FD<-picante::pd(time1_comm, trait_tree)
  FD
  vegan_tbl<-data.frame(Site = as.vector(rownames(FD)), SR = FD$SR, FD = FD$PD)
  
  for(j in 1:nrow(vegan_tbl)){
    trait_table_all_t$SR[trait_table_all_t$Site == paste(vegan_tbl$Site[j]) &
                           trait_table_all_t$time_code == time_stamp_list[i]]<-vegan_tbl$SR[j]
    trait_table_all_t$FD_veg[trait_table_all_t$Site == paste(vegan_tbl$Site[j]) &
                               trait_table_all_t$time_code == time_stamp_list[i]]<-vegan_tbl$FD[j]
  }
  
  
  #calculate functional richness
  trial<-dbFD(time1_traits, time1_comm)
  func_rich_metrics<-data.frame(Site = rownames(data.frame(trial)), FRic = trial$FRic,
                                FEve = trial$FEve, FDiv = trial$FDiv)
  
  for(k in 1:nrow(func_ruch_metrics)){
    trait_table_all_t$FRich[trait_table_all_t$Site == paste(func_rich_metrics$Site[k]) &
                              trait_table_all_t$time_code == time_stamp_list[i]]<-func_rich_metrics$FRic[k]
    trait_table_all_t$FEve[trait_table_all_t$Site == paste(func_rich_metrics$Site[k]) &
                             trait_table_all_t$time_code == time_stamp_list[i]]<-func_rich_metrics$FEve[k]
    trait_table_all_t$FDiv[trait_table_all_t$Site == paste(func_rich_metrics$Site[k]) &
                             trait_table_all_t$time_code == time_stamp_list[i]]<-func_rich_metrics$FDiv[k]
  }
  
  print(paste(i, "has finished", " "))
}
trait_table_all_t$Treatment<-rep("Removal", nrow(trait_table_all_t))
trait_table_all$Treatment<-rep("Control", nrow(trait_table_all))
trait_table_full<-rbind(trait_table_all_t,trait_table_all)


trait_table_full%>% group_by(time_code, Treatment) %>%
  summarise(mean_frich = mean(SR, na.rm = T)) %>%
  ggplot(., aes(x= (time_code), y = mean_frich, color = Treatment)) + 
  geom_line() +
  geom_point(pch = 21, size = 3, fill = "white",stroke = 2.5) + theme_bw()

trait_table_full%>% 
  ggplot(., aes(x= as.factor(time_code), y = FDiv, color = Treatment)) + 
  geom_boxplot()



###get total abundances by time_code by genus

by_genus<-df_c %>% dplyr::select(! c(Year, Site, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, genus) %>%
  summarise(total = sum(abund)) %>%
  filter(total > 0)


which_genera<-data.frame(generalist = unique(by_genus$genus), keep = NA)

for(i in 1:nrow(which_genera)){
  if(length(unique(by_genus$time_code[by_genus$genus==which_genera[i,1]]))> 7){
    which_genera$keep[i]<-"yes"
  }else{
    which_genera$keep[i]<-"no"
  }
}

which_genera_keep<-which_genera %>% filter(keep == "yes") %>%
  dplyr::select(generalist)

by_genus1<-by_genus %>% filter(genus %in% which_genera_keep$generalist)
ggplot(data = NULL) +
  geom_line(data = by_genus1, aes(x = time_code, y = log(total), color = genus),
            alpha = 1, size = 1.5) 

#' you can look at peaks and vallyes based on genus
#' you can identify maximum overlap time periods
#' how does morphospace change across wet and dry seasons?






average_by_genus<-df_c %>% dplyr::select(! c(Year, Site, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code) %>%
  summarise(total = mean(abund)) 

ggplot(data = NULL) +
  geom_line(data = by_genus, aes(x = time_code, y = log(total), group = genus),
            alpha = 0.4) +
  geom_line(data = average_by_genus, aes(x = time_code, y = log(total)), color = "blue",
            size = 2)





plot_data<-data.frame(matrix(NA, nrow = nrow(callingcard), ncol = 4))


colnames(plot_data)<-c("timestampcode", "Add.beta", "multi_beta", "beta_part")

for(i in 1:nrow(callingcard)){
  data_analysis<-df_c %>% filter(Year == callingcard[i,1] & Month == callingcard[i,2]) %>%
    select(-Year, -Site, -Month, -Treatment)
  alpha = vegan::specnumber(data_analysis)
  spN_allplots = colSums(data_analysis)
  gamma = length(which(spN_allplots > 0))
  plot_data$timestampcode[i]<-paste(callingcard[i,1], callingcard[i,2], sep = "_")
  plot_data$Add.beta[i] = gamma - mean(alpha)
  plot_data$multi_beta[i] = gamma/mean(alpha)
  plot_data$beta_part[i]= 1 - mean(alpha)/gamma
  
}

plot_data$time_period<-seq(1:14)
plot_data$Treatment<-rep("C",14)

#do this with treatment
df_c<-df %>%
  filter(Treatment == "T")


df_c %>% filter()

callingcard<-unique(df_c[,c("Year", "Month")])



plot_data2<-data.frame(matrix(NA, nrow = nrow(callingcard), ncol = 4))
colnames(plot_data2)<-c("timestampcode", "Add.beta", "multi_beta", "beta_part")

for(i in 1:nrow(callingcard)){
  data_analysis<-df_c %>% filter(Year == callingcard[i,1] & Month == callingcard[i,2]) %>%
    select(-Year, -Site, -Month, -Treatment)
  alpha = vegan::specnumber(data_analysis)
  spN_allplots = colSums(data_analysis)
  gamma = length(which(spN_allplots > 0))
  plot_data2$timestampcode[i]<-paste(callingcard[i,1], callingcard[i,2], sep = "_")
  plot_data2$Add.beta[i] = gamma - mean(alpha)
  plot_data2$multi_beta[i] = gamma/mean(alpha)
  plot_data2$beta_part[i]= 1 - mean(alpha)/gamma
  
}

plot_data2$time_period<-seq(1:14)
plot_data2$Treatment<-rep("T",14)
plot_data
plot_data2

#rbind
hello_world<-rbind(plot_data,plot_data2)

ggplot(hello_world) +
  geom_point(aes(x = time_period, y = beta_part, color =Treatment))
