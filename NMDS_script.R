#NMDS Script
library(vegan)
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
#We have to convert the model_data into nmds format for both traits and abundances
model_data<- read.csv("Data/model_data.csv")#load model data

#let's do abundances first
#read ecological data
df<-read.csv("Data/microbeta.csv")
df_c<-df %>%
  filter(Treatment == "C") #get only control plots
df_c$time_code<-paste(df_c$Month,df_c$Year,  sep = " ") #convert to date class
df_c$time_code<-as.Date(as.yearmon(df_c$time_code))#convert to date class
df_c$season<-NA
for(i in 1:nrow(df_c)){
  if(df_c$Month[i] %in% c("May","Jun","June", "July", "August", "September")){
    df_c$season[i]<-"WET"
  }else{
    df_c$season[i]<-"DRY"
  }
  
}

eco_nmds<-df_c %>% dplyr::select(!c(time_code, Treatment, Month, Site, Year))
eco_nmds<-eco_nmds %>% mutate(total=rowSums(select_if(., is.numeric))) %>%
  filter(!total == 0) 
season<-eco_nmds$season
eco_nmds<-eco_nmds %>% dplyr::select(-c(total, season))
#Abundance sqrt transformed with bray distances 
nmds1=metaMDS(eco_nmds,k=2, trymax=1000, distance = "bray", autotransform = T)


#extracting ellipses for ggplot
ord<-ordiellipse(nmds1, groups = season, draw = "polygon", lty = 1, col = "grey90")
NMDS = data.frame(NMDS1 = nmds1$points[,1], NMDS2 = nmds1$points[,2],group=as.factor(season))

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}




#plotting
scrs <- scores(nmds1, display = 'sites')

scrs <- cbind(as.data.frame(scrs), Season = season)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Season, data = scrs, FUN = mean)

segs <- merge(scrs, setNames(cent, c('Season','oNMDS1','oNMDS2')),
              by = 'Season', sort = FALSE)

ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = Season)) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()    


ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = Season), alpha = 0.8) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=1) +
  geom_segment(data = segs,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders                        # centroids
  geom_point(size = 3) + # sample scores
  geom_point(data = cent, size = 6) + 
  coord_fixed() +
  scale_color_viridis_d(option = "H") +
  theme_bw() +
  theme(legend.text = element_text(size  =12),
        legend.title = element_text(size =12 , face = "bold"))

ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = REALM), alpha = 0.8) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=1) +                      # centroids
  geom_point(size = 3, alpha = 0.5) + # sample scores
  geom_point(data = cent, size = 6) + 
  coord_fixed() +
  scale_color_viridis_d(option = "H") +
  theme_bw() +
  theme(legend.text = element_text(size  =12),
        legend.title = element_text(size =12 , face = "bold"))


#We now need to do an nmds with the weighted mean of each trait for each site

#read in trait data that has already been size corrected
species_traits<-read.csv("Data/species_traits_florida.csv") 
species_traits<-species_traits %>% mutate(log_cs = log(Colony_size)) %>%
  dplyr::select(!c(Colony_size, Notes))


nmds_trait<-df_c %>% dplyr::select(! c(Year, Month, Treatment, season)) %>%
  pivot_longer(cols = -c("time_code", "Site"), names_to = "Species", values_to="abund") %>%
  mutate(genus = str_extract(Species, "^[\\w]+")) %>%
  group_by(time_code,Site,Species) %>%
  summarise(total = sum(abund)) %>%
  left_join(.,species_traits, by = c("Species"="species")) %>%
  filter(!total == 0) %>%
  group_by(time_code, Site) %>%
  summarise(mean_CS = weighted.mean(log_cs, total),
            mean_WL = (weighted.mean(WL, total)),
            mean_HW = (weighted.mean(HW, total)),
            mean_EL = weighted.mean(EL, total),
            mean_ML = (weighted.mean(ML, total)))


nmds_trait$season<-NA
nmds_trait$time_code<- as.character(nmds_trait$time_code) #had to switch to character for this one
for(i in 1:nrow(nmds_trait)){
  if(nmds_trait$time_code[i] %in%  c("2018-05-01", "2019-05-01","2018-06-01",
                                     "2019-06-01", "2018-07-01", "2018-08-01",
                                 "2018-09-01")){
    nmds_trait$season[i]<-"WET"
  }else{
    nmds_trait$season[i]<-"DRY"
  }
  
}

season1<-nmds_trait$season
season2<-nmds_trait$time_code

nmds_trait<-nmds_trait %>% ungroup() %>%
  dplyr::select(-c(time_code, Site, season))

##Trait NMDS
nmds_trait_minus<-nmds_trait %>% dplyr::select(-mean_CS)
nmds_trait_s<-scale(nmds_trait)

nmds1=metaMDS(nmds_trait,k=3, trymax=1000, distance = "mahalanobis", autotransform = F)

#extracting ellipses for ggplot
ord<-ordiellipse(nmds1, groups = season2, draw = "polygon", lty = 1, col = "grey90")
NMDS = data.frame(NMDS1 = nmds1$points[,1], NMDS2 = nmds1$points[,2],group=(season1))

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}




#plotting
scrs <- scores(nmds1, display = 'sites')

scrs <- cbind(as.data.frame(scrs), Season = season1)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Season, data = scrs, FUN = mean)

segs <- merge(scrs, setNames(cent, c('Season','oNMDS1','oNMDS2')),
              by = 'Season', sort = FALSE)

ggplot(scrs, aes(x = NMDS1, y = NMDS2, colour = Season)) + # spiders
  geom_point(data = cent, size = 5) +                         # centroids
  geom_point() +                                              # sample scores
  coord_fixed()    



#do pca
## conduct principal coordinate analysis
pcoatrait<-nmds_trait
b_dist <- vegdist(pcoatrait, method="gower", diag=T)
pcoa1 <- pcoa(b_dist)
nmds_trait$pc1 <- pcoa1$vectors[1:137,1]
nmds_trait$pc2 <- pcoa1$vectors[1:137,2]
nmds_trait$season <- season2
## plot ordination
ggplot(nmds_trait, aes(x=pc1, y=pc2, color=season)) + geom_point() + stat_ellipse() + theme_bw()

?vegdist
