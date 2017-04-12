#exploratory plots
ggplot(PFTC_CO2_2015, aes(site, NEE_lm))+
  geom_boxplot()

ggplot(PFTC_CO2_2015, aes(treatment, NEE_lm, col=treatment))+
  geom_boxplot()+
  facet_grid(. ~ site)

ggplot(PFTC_CO2_2015, aes(origin, NEE_lm))+
  geom_boxplot()

ggplot(PFTC_CO2_2015, aes(treatment, NEE_lm, col=treatment))+
  geom_boxplot()+
  facet_grid(. ~ origin)

# difference between two chambers?
ggplot(PFTC_Reco_final, aes(group, Reco_lm, col=treatment))+
  geom_boxplot()+
  facet_grid(. ~ site)


#===============================================================================================================================
#more eleborate plots
#change order of the treatments
PFTC_Reco_final$treatment<-factor(PFTC_Reco_final$treatment, levels = c("c", "tt0", "otc", "tt1", "tt2", "tt3", "tt4"))

ggplot(Reco_1516, aes(site, Reco_lm))+
  geom_boxplot()+
  facet_grid(. ~ year)+
  theme_bw()+
  ylab("Reco (μmol CO2m–2 s–1)")+
  scale_x_discrete("elevation (m.a.s.l.)", labels = c("elev3000" = "3000", "elev3500" = "3500", "elev3850" = "3850","elev4100" = "4100"))

#plots with same origin in facets
ggplot(PFTC_Reco_final, aes(treatment, Reco_lm))+
  geom_boxplot()+
  facet_grid(. ~ origin)+
  theme_bw()+
  ylab("Reco (μmol CO2m–2 s–1)")+
  xlab("Treatment")
  #scale_x_discrete("treatment", labels = c("c" = "c","tt1" = "tt1", "tt2" = "tt2","tt3" = "tt3", tt4="tt4"))

#overall effect of treatments
ggplot(PFTC_Reco_final, aes(origin, Reco_lm))+
  geom_boxplot()+
  facet_grid(. ~ treatment)+
  theme_bw()+
  ylab("Reco (μmol CO2m–2 s–1)")+
  scale_x_discrete("Treatment")

PFTC_GPP_2016_PAR<- PFTC_GPP_2016 %>%
  filter(PAR.x > 600)
  
ggplot(PFTC_GPP_2016_PAR, aes(origin.x, GPP_lm))+
  geom_boxplot()+
  facet_grid(. ~ treatment)+
  theme_bw()+
  ylab("GPP (μmol CO2m–2 s–1)")+
  scale_x_discrete("Treatment")

ggplot(PFTC_CO2_1516, aes(origin, NEE_lm, col=year))+
  geom_boxplot()+
  facet_grid(. ~ treatment)+
  theme_bw()+
  ylab("Reco (μmol CO2m–2 s–1)")+
  scale_x_discrete("Treatment")
