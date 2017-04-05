# PFTC China CO2 flux data analysis
# Load all CO2 flux data
library(plyr)
library(ggplot2)
#install.packages("smatr")
library(smatr)
library(mgcv)

#PFTC_CO2_2016<- read.table("PFTC_CO2flux_all_2016.txt", header = TRUE, sep = "\t", dec = ",") 
PFTC_CO2_2016<- read.table("Data\\CO2\\PFTC_CO2flux_all_new.txt", header = TRUE, sep = "\t", dec = ".") 
str(PFTC_CO2_2016)

# set all numeric parameters to 2 decimals 
is.num <- sapply(PFTC_CO2_2016, is.numeric)
PFTC_CO2_2016[is.num] <- lapply(PFTC_CO2_2016[is.num], round, 2)

# change datetime format and create seperate data and time column from datetime column
PFTC_CO2_2016$datetime<-as.POSIXct(PFTC_CO2_2016$datetime, tz="", format="%d.%m.%Y %H:%M") #format datetime
PFTC_CO2_2016$date<-as.Date(PFTC_CO2_2016$datetime, format="%d.%m.%Y")
PFTC_CO2_2016$time<-format(PFTC_CO2_2016$datetime, format="%H:%M")
PFTC_CO2_2016$year<-format(PFTC_CO2_2016$datetime, format="%Y")
PFTC_CO2_2016$timeframe<-PFTC_CO2_2016$tfinish -PFTC_CO2_2016$tstart

#change order of the treatments
PFTC_CO2_2016$treatment<-factor(PFTC_CO2_2016$treatment, levels = c("c", "tt0", "otc", "tt1", "tt2", "tt3", "tt4"))

# create column for plot origin of transplanted plots 
elevs <- c(3000, 3500, 3850, 4100) 
# tt1 - 1, tt2 + 1, tt3 - 3, tt4 + 3
changeVec <- c(1, -1, 3, -3, 0, 0, 0)
names(changeVec) <- c("tt1", "tt2", "tt3", "tt4", "otc", "c", "tt0")

originVec <- apply(X = cbind(
  as.integer(gsub("elev", "", PFTC_CO2_2016$site, fixed = TRUE)),
  as.character(PFTC_CO2_2016$treatment)),
  FUN = function(inRow, changeVec, elevs) {
    elevIndex <- which(inRow[1] == elevs)
    elevs[elevIndex + changeVec[inRow[2]]]
  }, MARGIN = 1, changeVec = changeVec, elevs = elevs)
PFTC_CO2_2016$origin <- factor(paste("elev", originVec, sep = ""))

#exploratory plots for finding outliers
#ggplot(PFTC_CO2_2016, aes(temp.C, NEE_lm, col=site))+
#  geom_point()+
#  facet_grid(.~ type)

# remove outliers # 2 outliers; inge tt0 _2  elev 3850 (extreme high value), sun elev 3500 otc_7 (+value)
PFTC_CO2_2016<-PFTC_CO2_2016[c(-356, -50), ]

# leave out bad measurements, checked measurement plots individually
#PFTC_CO2_2016<-PFTC_CO2_2016[c(-29, -42, -72, -299, -310, -353, -401, -415), ]

# select only measurements with timeframe of at least 45s
PFTC_CO2_2016<-subset(PFTC_CO2_2016, timeframe>= 45 & rsqd >=0.8)

#seperate photo and res measurements and merge them in new file with new column GPP, 
PFTC_Photo_2016<-subset(PFTC_CO2_2016, type== "photo" )   #& rsqd>=.9 (selecting data with r2>=.9)
PFTC_Reco_2016<-subset(PFTC_CO2_2016, type== "resp" )  #& rsqd>=.9

# calculate GPP from Photo-Reco 
PFTC_Reco_2016$Reco_lm<-PFTC_Reco_2016$NEE_lm*-1
PFTC_Reco_2016$temp.K<-PFTC_Reco_2016$temp.C+273.15
PFTC_GPP_2016<- merge(PFTC_Photo_2016, PFTC_Reco_2016, by=c("site", "block", "treatment", "date", "group"))
PFTC_GPP_2016$GPP_lm<-PFTC_GPP_2016$NEE_lm.x- PFTC_GPP_2016$NEE_lm.y #NEE-Reco

library(dplyr)
# calculate means for each treatment at different elevation (site)
# comparing measurement to their min/max limit, min_limit = mean_Reco*0.5 and max_limit = mean_Reco*2. If exceeding report FALSE 
limits<- PFTC_Reco_2016 %>% 
  group_by(treatment, origin) %>%#site,
  mutate(mn = mean(Reco_lm)) %>%
  mutate(keep = Reco_lm < mn *2 & Reco_lm> mn*0.5)


#rechecked measurements outside of limits, taking out bad measurements and extreme high Reco values compared with other measurement of same plot.
PFTC_Reco_final<-PFTC_Reco_2016[c(-4,-5, -6, -15, -17, -31, -36, -56, -76, -128, -140, -150, -151, -159, -163, -171, -174, -176, -182, -224),]


#Decrease number of columns in df
PFTC_Reco_final<-subset(PFTC_Reco_2016_limit, select=c("year", "date", "time", "timeframe", "site", "block", "treatment", "origin", "temp.K", "PAR", "Reco_lm")) 
                        

#check for different relation with temperature between plots from same origin 
ggplot(PFTC_Reco_final, aes(temp.K, Reco_lm, col=origin))+
  geom_point()+
  geom_smooth(formula = y~x, method = "glm", se = FALSE)


# Calculate flux taking in account biomass at site > use biomass.xlsx. Only possible for controls, because the harvested plots are supposed to be similar to the plots measured. 



#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# NORMALIZING FLUXDATA
#fit Arheniusequation from Loyd&Taylor 1994<- A*exp(E0/T-T0)
#equation from longdoz et al 2000 based on lloyd and taylor 1994  
#Ea<-12970*xtempK.9/(xtempK.9-227.13), Fs =Fs10 * exp {Ea *(Ts - 283.2)/(283.2 * R * Ts)}
require(graphics)

# should do this per origin!! because of relation of Reco with temperature for the different vegetation origins
fit.nls<-nls((Reco_lm~A*exp(-308.56/I(temp.K-227.13))), start=c(A=0 ), data=PFTC_Reco_final)
fit.nls #A=-816.6  residual sum-of-squares: 5972
summary(fit.nls)

#plot residuals vs fitted values
E1 <- resid(fit.nls)   #better: E0 <- rstandard(M0)
F1 <- fitted(fit.nls)
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab = "Residuals",
     cex.lab = 1.5,
     abline(h = 0, lty = 2))  

boxplot(E1 ~ origin, data = PFTC_Reco_final) # heterogeneity?
#Is the ratio of the largest and smalled variance > 4
tapply(E1, INDEX = PFTC_Reco_final$origin, FUN = var)

Recofit<-plot(Reco_lm~temp.K, data = PFTC_Reco_final)
#lines(sort(tempK), fitted(fit.nls)[order(tempK)], col="red") #rough fit
newx<-seq(min(PFTC_Reco_final$temp.K), max(PFTC_Reco_final$temp.K), length=100) #create new x values for fitting line
Reco.pred<-predict(fit.nls, newdata=data.frame(temp.K=newx)) #predict line from newx values 
lines(newx, Reco.pred, col="blue")

#Recalculating flux values taking in account heteroscedasticity 
recalcflux<-PFTC_Reco_final$Reco_lm/fitted(fit.nls)*(coef(fit.nls)*exp(-308.56/I(283.15-227.13)))
PFTC_Reco_final$Reco15<-recalcflux
plot(PFTC_Reco_final$temp.K, PFTC_Reco_final$Reco15)

ggplot(PFTC_Reco_final, aes(site, Reco15, col=factor(treatment)))+
  geom_boxplot()+
  #facet_grid(. ~ site)+
  theme_bw()+
  ylab("Reco (μmol CO2m–2 s–1)")+
  scale_x_discrete("Treatment")



#==============================================================================================================================
#load data PFTC_CO2_2015
library(dplyr)
library(tidyr)
library(raster)

#PFTC_CO2_2016<- read.table("PFTC_CO2flux_all_2016.txt", header = TRUE, sep = "\t", dec = ",") 
PFTC_CO2_2015<- read.table("O:\\PFTC_China\\PFTC_DATA\\CO2_flux\\PFTC_CO2flux_2015.txt", header = TRUE, sep = "\t", dec = ".") 
str(PFTC_CO2_2015)

#remove unimportant columns
PFTC_CO2_2015[13:17]<- NULL

# set all numeric parameters to 2 decimals
is.num <- sapply(PFTC_CO2_2015, is.numeric)
PFTC_CO2_2015[is.num] <- lapply(PFTC_CO2_2015[is.num], round, 2)

#change sign for NEE_lm and NEE_exp
PFTC_CO2_2015$NEE_lm<-PFTC_CO2_2015$NEE_lm*-1
PFTC_CO2_2015$NEE_exp<-PFTC_CO2_2015$NEE_exp*-1

#format datetime column
PFTC_CO2_2015$datetime<-as.POSIXct(PFTC_CO2_2015$datetime, tz="", format="%m/%d/%y %H:%M") #format datetime
PFTC_CO2_2015$datetime<-format(PFTC_CO2_2015$datetime, format="%d.%m.%Y %H:%M")

# create seperate data and time column from datetime column
PFTC_CO2_2015$date<-as.Date(PFTC_CO2_2015$datetime, format="%d.%m.%Y")
PFTC_CO2_2015$year<-format(PFTC_CO2_2015$date, format="%Y")
PFTC_CO2_2015$time<-format(PFTC_CO2_2015$date, format="%H:%M")

# create new column recalculating Reco from NEE
PFTC_CO2_2015$Reco_lm<-PFTC_CO2_2015$NEE_lm*-1

# create column for plot origin of transplanted plots 
elevs <- c(3000, 3500, 3850, 4100) 
# tt1 - 1, tt2 + 1, tt3 - 3, tt4 + 3
changeVec <- c(1, -1, 3, -3, 0, 0, 0)
names(changeVec) <- c("tt1", "tt2", "tt3", "tt4", "otc", "c")

originVec <- apply(X = cbind(
  as.integer(gsub("elev", "", PFTC_CO2_2015$site, fixed = TRUE)),
  as.character(PFTC_CO2_2015$treatment)),
  FUN = function(inRow, changeVec, elevs) {
    elevIndex <- which(inRow[1] == elevs)
    elevs[elevIndex + changeVec[inRow[2]]]
  }, MARGIN = 1, changeVec = changeVec, elevs = elevs)
PFTC_CO2_2015$origin <- factor(paste("elev", originVec, sep = ""))



# combine data of 2015 and 2016
rbind.match.columns <- function(input1, input2) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  
  if (n.input2 > n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}


rbind.match.columns(PFTC_CO2_2015, PFTC_Reco_final)
PFTC_CO2_1516<-rbind.match.columns(PFTC_CO2_2015, PFTC_Reco_2016)
#change order of factor treatment
PFTC_CO2_1516$treatment<-factor(PFTC_CO2_1516$treatment, levels = c("c", "tt0", "otc", "tt1", "tt2", "tt3", "tt4"))
# date of 2016 not in correct format

