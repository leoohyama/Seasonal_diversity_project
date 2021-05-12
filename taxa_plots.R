library(tidyverse)

#Taxaplots


#read ecological data
df<-read.csv("Data/microbeta.csv")
df_c<-df %>%
  filter(Treatment == "C") #get only control plots
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ") #convert to date class
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))#convert to date class


taxa_plots<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,genus,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  ungroup() %>%
  filter(!total == 0) %>%
  group_by(time_code,Species) %>%
  mutate(total_n = sum(total)) 

taxa_plots %>%
  group_by(time_code, Subfamily, genus) %>%
  summarise(total_f = sum(total)) %>%
  ggplot(.) +
  geom_col(aes(x = time_code, y = log(total_f), fill = genus), position = position_dodge()) +
  facet_wrap(~Subfamily, scales = "free_y")

taxa_plots %>%
  group_by(time_code, Subfamily, genus) %>%
  summarise(total_f = sum(total)) %>%
  ungroup() %>%
  ggplot(.) +
  geom_line(aes(x = time_code, y = log(total_f), color = genus)) +
  geom_point(aes(x = time_code, y = log(total_f), color = genus), pch =21, size = 2, fill = "white") +
  facet_wrap(~Subfamily, scales = "free_y")


# you need to make a graph that shows which genus goes up or down relative to their mean of the year
