library(vegan)
library(tidyverse)
library(zoo)
library(viridis)
library(picante)
library(vegan)
library(FD)
library(MASS)

#read ecological data
df<-read.csv("Data/microbeta.csv")
df_c<-df %>%
  filter(Treatment == "C") #get only control plots
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ") #convert to date class
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))#convert to date class

#read morphological trait data
dft<-read.csv("Data/micro_traits.csv")
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
#setwd("~/Desktop/")
#write.csv(species_traits, "species_traits_florida.csv") write out the trait data


#read in new trait data with colony size
species_traits<-read.csv("Data/species_traits_florida.csv")


species_traits<-species_traits %>% dplyr::select(-Notes) %>%
  mutate(log_colony = log(Colony_size)) %>% dplyr::select(-Colony_size)


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
trait_table_all$jost<-NA
trait_table_all$mean_FRich_Null<-NA
trait_table_all$mean_FEve_Null<-NA
trait_table_all$mean_FDiv_Null<-NA
trait_table_all$sd_FRich_Null<-NA
trait_table_all$sd_FEve_Null<-NA
trait_table_all$sd_FDiv_Null<-NA

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
  
  #calculate jost diversity
  jost.frame=data.frame(Jost =exp(vegan::diversity(time1_comm, index = "shannon")),
                        Site = rownames(data.frame(exp(vegan::diversity(time1_comm, index = "shannon")))))
  for(m in 1:nrow(jost.frame)){
    trait_table_all$jost[trait_table_all$Site == paste(jost.frame$Site[m]) &
                           trait_table_all$time_code == time_stamp_list[i]]<-jost.frame$Jost[m]
  }
  
  #rescale trait data
  trait_data = scale(time1_traits)
  
  trait_matrix = vegan::vegdist(trait_data,"euclidean")
  branch_lengths<-hclust(trait_matrix,"average")
  branch_lengths
 
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
  
  
  #calculate the ses of functional indices for 1000 null models


  ##################################################################
  #this is where you randomize the community matrix
  ##################################################################
  #first create null model datset to hold all the values
  null_model_results<-data.frame(run = rep(1:1000, each =length(rownames(time1_comm))),
                                 Site = rep(rownames(time1_comm), 1000), FRich.Null = NA,
                                 FEve.Null = NA, FDiv.Null = NA, ID = seq(1:1000),
                                 Time = rep(time_stamp_list[i], 1000))
  for (b in 1:1000) {
    nullmatrix<-time1_comm
    nullmatrix[nullmatrix >0] <- 1
    randommatrixswap<-randomizeMatrix(time1_comm, null.model = c("independentswap"), iterations = 1000)
    randommatrixswap
    #then you calculate calculate functional metrics with randomized comm data for each time
    #period(need to do this 1000 times). 
    trial<-dbFD(time1_traits, randommatrixswap, w.abun =T,calc.FDiv =T)
    func_rich_metrics<-data.frame(Site = rownames(data.frame(trial)), FRic = trial$FRic,
                                  FEve = trial$FEve, FDiv = trial$FDiv, 
                                  Run = rep(paste(b), length(rownames(data.frame(trial))))
                                  )
    
    for(a in 1:nrow(func_rich_metrics)){
      null_model_results$FRich.Null[null_model_results$run == paste(b) &
                                      null_model_results$Site == paste(func_rich_metrics$Site[a])]<-func_rich_metrics$FRic[a]
      null_model_results$FEve.Null[null_model_results$run == paste(b) &
                                      null_model_results$Site == paste(func_rich_metrics$Site[a])]<-func_rich_metrics$FEve[a]
      null_model_results$FDiv.Null[null_model_results$run == paste(b) &
                                      null_model_results$Site == paste(func_rich_metrics$Site[a])]<-func_rich_metrics$FDiv[a]
      null_model_results$time[null_model_results$run == paste(b) &
                                   null_model_results$Site == paste(func_rich_metrics$Site[a])]<-func_rich_metrics$FDiv[a]
    }
    
  }
  null_summary<-null_model_results %>% group_by(Site) %>%
    summarise(mean_FRich.null = mean(FRich.Null, na.rm =T),
              sd_FRich.null = sd(FRich.Null, na.rm =T),
              mean_FEve.null = mean(FEve.Null, na.rm = T),
              sd_FEve.null = sd(FEve.Null, na.rm =T),
              mean_FDiv.null = mean(FDiv.Null, na.rm = T),
              sd_FDiv.null = sd(FDiv.Null, na.rm =T))

  for(z in 1:nrow(null_summary)){
    trait_table_all$mean_FRich_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                      trait_table_all$time_code == time_stamp_list[i]]<-null_summary$mean_FRich.null[z]
    trait_table_all$mean_FEve_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                      trait_table_all$time_code == time_stamp_list[i]]<-null_summary$mean_FEve.null[z]
    trait_table_all$mean_FDiv_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                      trait_table_all$time_code == time_stamp_list[i]]<-null_summary$mean_FDiv.null[z]
    trait_table_all$sd_FRich_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                     trait_table_all$time_code == time_stamp_list[i]]<-null_summary$sd_FRich.null[z]
    trait_table_all$sd_FEve_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                     trait_table_all$time_code == time_stamp_list[i]]<-null_summary$sd_FEve.null[z]
    trait_table_all$sd_FDiv_Null[trait_table_all$Site == paste(null_summary$Site[z]) &
                                   trait_table_all$time_code == time_stamp_list[i]]<-null_summary$sd_FDiv.null[z]
  }

  
  ##########################
  #calculate functional richness
  trial<-dbFD(time1_traits, time1_comm, w.abun =T,calc.FDiv =T)
  func_rich_metrics<-data.frame(Site = rownames(data.frame(trial)), FRic = trial$FRic,
                                FEve = trial$FEve, FDiv = trial$FDiv)
  
  for(k in 1:nrow(func_rich_metrics)){
    trait_table_all$FRich[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                            trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FRic[k]
    trait_table_all$FEve[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                           trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FEve[k]
    trait_table_all$FDiv[trait_table_all$Site == paste(func_rich_metrics$Site[k]) &
                           trait_table_all$time_code == time_stamp_list[i]]<-func_rich_metrics$FDiv[k]
  }
  
  
  ses.mpd.traits = picante::ses.mpd(time1_comm, trait_matrix, null.model = "taxa.labels",
                                    abundance.weighted = TRUE, runs = 500)
  
  print(paste(i, "has finished", " "))
}











trait_table_all %>% mutate(difference = mpd.obs-mpd.rand.mean) %>%
  group_by(time_code) %>%
  summarise(mean_dff = mean(difference, na.rm =T))
min(trait_table_all$mpd.obs.p)
trait_table_all%>% group_by(time_code, Site) %>%
  summarise(mean_metric = mean(mpd.obs.p, na.rm = T)) %>%
  ggplot(., aes(x= (time_code), y = mean_metric, color = Site)) + 
  geom_line()


trait_table_all %>%
  group_by(time_code) %>% summarise(mean = mean(FRich, na.rm =T))
####facet graph comparinng all metrics
facet_data<-trait_table_all %>%
  group_by(time_code) %>%
  summarise(total_plots = n(), 
            mean_FD_veg = mean(FD_veg, na.rm = T),
            std.error.FD_veg=(sd(FD_veg, na.rm = T)/sqrt(total_plots)),
            mean_FRich = mean(FRich, na.rm = T), 
            std.error.FRich=(sd(FRich, na.rm = T)/sqrt(total_plots)),
            mean_SR = mean(SR, na.rm=T),
            std.error.SR = sd(SR, na.rm=T)/sqrt(total_plots),
            mean_Fdiv = mean(FDiv, na.rm =T),
            std.error.Fdiv = sd(FDiv, na.rm = T)/sqrt(total_plots),
            mean_FEve = mean(FEve, na.rm =T),
            std.error.FEve = sd(FEve, na.rm = T)/sqrt(total_plots),
            mean_Jost = mean(jost,na.rm =T),
            std.error.Jost=sd(jost, na.rm = T)/sqrt(total_plots),
            mean_abundance = mean(total),
            std.error.abund = sd(total, na.rm=T)/sqrt(total_plots)) %>%
  pivot_longer(., cols = c("mean_FD_veg", "mean_SR", "mean_FRich","mean_Fdiv","mean_FEve","mean_Jost","mean_abundance"
  ), names_to = c("type"),
  values_to= "value") 

facet_data$std.error<-NA
for(i in 1:nrow(facet_data)){
  if(grepl("FRich", facet_data$type[i])){
    facet_data$std.error[i]<-facet_data$std.error.FRich[i]
  }else{
    if(grepl("SR", facet_data$type[i])){facet_data$std.error[i]<-facet_data$std.error.SR[i]}else{
      if(grepl("Fdiv", facet_data$type[i])){
        facet_data$std.error[i]<-facet_data$std.error.Fdiv[i]
      }else{if(grepl("Jost", facet_data$type[i])){
        facet_data$std.error[i]<-facet_data$std.error.Jost[i]
      }else
      {if(grepl("FEve", facet_data$type[i])){
        facet_data$std.error[i]<-facet_data$std.error.FEve[i]
      }else{if(grepl("FD_veg", facet_data$type[i])){
        facet_data$std.error[i]<-facet_data$std.error.FD_veg[i]
      }else{facet_data$std.error[i]<-facet_data$std.error.abund[i]}
        
      }
      }
        
      }
    }
  }
}  
#getting mean values for null model results:



#this is a key for the geom_rect, displays values for each month
data.frame(date = unique(facet_data$time_code), 
           code = as.numeric(unique(facet_data$time_code)))


facet_data$type<-factor(facet_data$type)
levels(facet_data$type)<-c("Abundance", "Functional Diversity","Functional Divergence", "Functional Evenness", 
                           "Functional Richness", "Jost Diversity", "Species Richness")



ggplot(data = facet_data) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  geom_line(aes(x= (time_code), y = value), size = 1.5) +
  geom_errorbar(aes(x = time_code, ymin = value-std.error,
                    ymax = value+std.error), width = 0.5)+
  geom_point(aes(x= (time_code), y = value), 
             pch = 21, size = 3, fill = "white",stroke = 2.5) +
  facet_wrap(~type,scales = "free_y", nrow = 2) +
  geom_hline(data = subset(facet_data, type == "Abundance"), 
             aes(yintercept = mean(trait_table_all$total)), lty= 2)+
  geom_hline(data = subset(facet_data, type == "Functional Divergence"), 
             aes(yintercept = mean(trait_table_all$FDiv,na.rm=T)), lty= 2) +
  geom_hline(data = subset(facet_data, type == "Functional Evenness"), 
             aes(yintercept = mean(trait_table_all$FEve,na.rm=T)), lty= 2) +
  geom_hline(data = subset(facet_data, type == "Functional Richness"), 
             aes(yintercept = mean(trait_table_all$FRich,na.rm=T)), lty= 2) +
  geom_hline(data = subset(facet_data, type == "Jost Diversity"), 
             aes(yintercept = mean(trait_table_all$jost,na.rm=T)), lty= 2) +
  geom_hline(data = subset(facet_data, type == "Species Richness"), 
             aes(yintercept = mean(trait_table_all$SR,na.rm=T)), lty= 2)  + 
  geom_smooth(data = subset(facet_data, type == "Abundance"), method = "glm.nb", 
              se = TRUE, 
              formula = y ~ poly(x,4), aes(x = time_code, y = value),
              colour = "red") + 
  geom_smooth(data = subset(facet_data, type == "Functional Evenness"), method = "lm", 
              se = TRUE, 
              formula = y ~ poly(x,2), aes(x = time_code, y = value),
              colour = "red") +
  geom_smooth(data = subset(facet_data, type == "Functional Richness"), method = "lm", 
              se = TRUE, 
              formula = y ~ poly(x,3), aes(x = time_code, y = value),
              colour = "red") +
  geom_smooth(data = subset(facet_data, type == "Species Richness"), method = "glm", 
              method.arg = list(family = "poisson"),
              se = TRUE, 
              formula = y ~ poly(x,4), aes(x = time_code, y = value),
              colour = "red") +
  labs(x = "", y =NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5),
        strip.text = element_text(size = 10, face = "bold"))


####facet graph by ses of functional metrics
SES.data<-trait_table_all %>%
  mutate(SES_FRich = (mean(FRich, na.rm = T)- mean(mean_FRich_Null, na.rm = T))/sd_FRich_Null,
         SES_FDiv = (mean(FDiv,na.rm = T)- mean(mean_FDiv_Null,na.rm = T))/sd_FDiv_Null,
         SES_FEve = (mean(FEve,na.rm = T)- mean(mean_FEve_Null,na.rm = T))/sd_FEve_Null) %>% 
  dplyr::select(time_code, Site, SES_FRich, SES_FDiv,SES_FEve)%>%
  group_by(time_code) %>%
  summarise(total_plots = n(), 
            mean_SES_FRich = mean(SES_FRich, na.rm = T),
            std.error.SES_FRich=(sd(SES_FRich, na.rm = T)/sqrt(total_plots)),
            mean_SES_FDiv = mean(SES_FDiv, na.rm = T),
            std.error.SES_FDiv=(sd(SES_FDiv, na.rm = T)/sqrt(total_plots)),
            mean_SES_FEve = mean(SES_FEve, na.rm = T),
            std.error.SES_FEve=(sd(SES_FEve, na.rm = T)/sqrt(total_plots))) %>%
  pivot_longer(.,cols = c("mean_SES_FRich","mean_SES_FDiv", "mean_SES_FEve"), names_to = c("type"),
               values_to= "value") 

SES.data$std.error<-NA
for(i in 1:nrow(SES.data)){
  if(grepl("FRich", SES.data$type[i])){
    SES.data$std.error[i]<-SES.data$std.error.SES_FRich[i]
  }else{
    if(grepl("FDiv", SES.data$type[i])){SES.data$std.error[i]<-SES.data$std.error.SES_FDiv[i]}else{
      SES.data$std.error[i]<-SES.data$std.error.SES_FEve[i]
    }
  }
}  

ggplot(SES.data) +
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  
  geom_line(aes(x= (time_code), y = value), size = 1.5) +
  geom_errorbar(aes(x = time_code, ymin = value-std.error,
                    ymax = value+std.error), width = 0.5)+
  geom_point(aes(x= (time_code), y = value), 
             pch = 21, size = 3, fill = "white",stroke = 2.5) +
  facet_wrap(~type, scales = "free_y", nrow = 3) +
  geom_hline(yintercept = 0, lty = 1) +
  theme_bw()

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
trait_table_all<-trait_table_all %>%
  mutate(SES_FRich = (mean(FRich, na.rm = T)- mean(mean_FRich_Null, na.rm = T))/sd_FRich_Null,
         SES_FDiv = (mean(FDiv,na.rm = T)- mean(mean_FDiv_Null,na.rm = T))/sd_FDiv_Null,
         SES_FEve = (mean(FEve,na.rm = T)- mean(mean_FEve_Null,na.rm = T))/sd_FEve_Null)
write.csv(trait_table_all,"Data/model_data.csv")




















