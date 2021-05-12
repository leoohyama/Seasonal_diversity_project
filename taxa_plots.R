library(tidyverse)

#Taxaplots


#read ecological data
df<-read.csv("Data/microbeta.csv")
df_c<-df %>%
  filter(Treatment == "C") #get only control plots
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ") #convert to date class
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))#convert to date class


by_sp<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  ungroup() %>%
  filter(!total == 0) %>%
  group_by(time_code,Species) %>%
  mutate(total_n = sum(total)) 