library(tidyverse)
library(zoo)
#Taxaplots


#read ecological data
df<-read.csv("Data/microbeta.csv")
species_traits<-read.csv("Data/species_traits_florida.csv")

df_c<-df %>%
  filter(Treatment == "C") #get only control plots
df_c<-df_c[, colSums(df_c != 0) > 0]
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ") #convert to date class
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))#convert to date class


taxa_plots<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,genus,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  ungroup() %>%
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
  geom_line(aes(x = time_code, y = (total_f), color = genus)) +
  geom_point(aes(x = time_code, y = (total_f), color = genus), pch =21, size = 2, fill = "white") +
  facet_wrap(~genus, scales = "free_y")


library(RColorBrewer)
colourCount = length(unique(taxa_plots$Species))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))


taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>%
  ungroup() %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_line(aes(x = time_code, y = (total_f), color = Species)) +
  geom_point(aes(x = time_code, y = (total_f), fill = Species),
             pch =21, size = 2, color = "black") +
  facet_wrap(~genus, scales = "free_y") +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) + theme_dark() +
  theme(axis.text.x = element_blank(),
        legend.position = "top")

taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>% 
  group_by(Species) %>%
  mutate(species_total = sum(total_f)) %>%
  mutate(perc_tot = total_f/species_total) %>%
  ungroup() %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_line(aes(x = time_code, y = (perc_tot), color = Species)) +
  geom_point(aes(x = time_code, y = (perc_tot), fill = Species),
             pch =21, size = 2, color = "black") +
  facet_wrap(~genus) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) + theme_dark() +
  theme(axis.text.x = element_blank(),
        legend.position = "top")

taxa_table<-taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>% 
  group_by(Species) %>%
  mutate(species_total = sum(total_f)) %>%
  filter(!species_total<10) %>%
  mutate(perc_tot = total_f/species_total) %>%
  ungroup() %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_line(aes(x = time_code, y = (perc_tot), color = Species)) +
  geom_point(aes(x = time_code, y = (perc_tot), fill = Species),
             pch =21, size = 2, color = "black") +
  facet_wrap(~genus) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) + theme_dark() +
  theme(axis.text.x = element_blank(),
        legend.position = "top")

taxa_table

#looking at rare species
taxa_table<-taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>% 
  group_by(Species) %>%
  mutate(species_total = sum(total_f)) %>%
  filter(species_total<10) %>%
  mutate(perc_tot = total_f/species_total) %>%
  ungroup() %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_line(aes(x = time_code, y = (perc_tot), color = Species)) +
  geom_point(aes(x = time_code, y = (perc_tot), fill = Species),
             pch =21, size = 2, color = "black") +
  facet_wrap(~genus) +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) + theme_dark() +
  theme(axis.text.x = element_blank(),
        legend.position = "top")

taxa_table


raw_summary<-df_c %>% dplyr::select(! c(Year, Month, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code, Site,genus,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) 
least<-raw_summary %>%
  group_by(Species) %>%
  summarise(total = sum(total)) %>%
  arrange(desc(total))


##delineate between wet and dry


seasonal_summary<-df_c %>% dplyr::select(! c(Year, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site", "Month"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code,Month, Site,genus,Species) %>%
  summarise(total = sum(abund))

seasonal_summary$season<-NA
for(i in 1:nrow(seasonal_summary)){
  if(seasonal_summary$Month[i] %in% c("May","Jun","June", "July", "August", "September")){
    seasonal_summary$season[i]<-"WET"
  }else{
    seasonal_summary$season[i]<-"DRY"
  }
  
}


###seasonal time code species richness
seasonal_summary<-df_c %>% dplyr::select(! c(Year, Treatment)) %>%
  pivot_longer(cols = -c("time_code", "Site", "Month"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code,Month, Site,genus,Species) %>%
  summarise(total = sum(abund)) %>%
  filter(! total == 0) %>%
  group_by(time_code) %>%
  filter(!duplicated(Species)) %>%
  summarise(total_SR = n()) %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=0.03) +
  geom_line(aes(x= (time_code), y = total_SR), size = 1.5) +
  geom_point(aes(x= (time_code), y = total_SR), 
             pch = 21, size = 3, fill = "white",stroke = 2.5) +
  labs(x = NULL, y = "Total species richness") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, face = "bold", vjust = 0.5),
        strip.text = element_text(size = 10, face = "bold"))
  



taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>%
  ungroup() %>%
  ggplot(.) + 
  geom_rect(xmin=17652,
            xmax=17806, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_rect(xmin=18017,
            xmax=Inf, ymin=-Inf, ymax=Inf, fill="lightblue", alpha=1) +
  geom_line(aes(x = time_code, y = (total_f), color = Species, group = Species)) +
  geom_point(aes(x = time_code, y = (total_f), fill = Species),
             pch =21, size = 2, color = "black") +
  facet_wrap(~genus, scales = "free_y") +
  labs(y = "Abundance", x = "Time") +
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_color_manual(values = getPalette(colourCount)) + theme_dark() +
  theme(axis.text.x = element_blank(),
        legend.position = "none" ,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        strip.text = element_text(family = "Helvetica", face= "bold",
                                  size =12))


##getting data prepped for shiny app
shiny_data<-taxa_plots %>%
  group_by(time_code, Subfamily, genus,Species) %>%
  summarise(total_f = sum(total)) %>%
  ungroup() 
str(shiny_data)

saveRDS(shiny_data, "~/Google Drive/Shinyapp/data_ants.rds")
