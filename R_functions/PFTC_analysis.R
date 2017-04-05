
fligner.test(NEE_lm~site, data = PFTC_CO2_1516)
#significant difference in variance! Due to high values at 3850 and less values data at 4100 because of less blocks?
# subset data of different sites
elev3000<-subset(PFTC_CO2data, origin == "elev3000")
elev3500<-subset(PFTC_CO2data, origin == "elev3500")
elev3850<-subset(PFTC_CO2data, origin == "elev3850")
elev4100<-subset(PFTC_CO2data, origin == "elev4100")


model15<- aov(Reco_lm ~ treatment*site, data= PFTC_CO2_1516)
summary(model15)#no interaction between site and treatment >additive
summary.lm(model15)

model16<- aov(Reco_lm ~ treatment*site, data= PFTC_CO2_1516)
summary(model16)# interaction between site and treatment 
summary.lm(model16)

model1516<- aov(Reco_lm ~ treatment*site, data= PFTC_CO2_1516)
summary(model1516) #interaction between site and treatment
summary.lm(model1516)
# most significant differences when combining data of two years