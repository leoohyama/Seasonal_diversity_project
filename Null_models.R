###This calculate the null models for our data
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

#read  trait data
species_traits<-read.csv("Data/species_traits_florida.csv")

#modify trait data
species_traits<-species_traits %>% dplyr::select(-Notes) %>%
  mutate(log_colony = log(Colony_size)) %>% dplyr::select(-Colony_size)


#program to get all fmetrics

#Set up trait table
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




