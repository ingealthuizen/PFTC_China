#load data PFTC_CO2_2015
library(dplyr)
library(tidyr)
library(raster)

#PFTC_CO2_2016<- read.table("PFTC_CO2flux_all_2016.txt", header = TRUE, sep = "\t", dec = ",") 
PFTC_CO2_2015<- read.table("Data\\CO2\\PFTC_CO2flux_2015.txt", header = TRUE, sep = "\t", dec = ".") 
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
  
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}

rbind.match.columns(PFTC_CO2_2015, PFTC_Reco_final)
PFTC_CO2_1516<-rbind.match.columns(PFTC_CO2_2015, PFTC_CO2_2016)
#change order of factor treatment
PFTC_CO2_1516$treatment<-factor(PFTC_CO2_1516$treatment, levels = c("c", "tt0", "otc", "tt1", "tt2", "tt3", "tt4"))
# date of 2016 not in correct format

#only select data from dark measurements 
Reco_1516<-subset(PFTC_CO2_1516, type== "resp" & rsqd>=.8 ) #& rsqd>=.9& rsqd>=.8
Reco_1516$Reco_lm<-Reco_1516$NEE_lm*-1
Reco_1516$Reco_exp<-Reco_1516$NEE_exp*-1

#subset different years
Reco_15<- subset(Reco_1516, year == "2015")
Reco_16<- subset(Reco_1516, year == "2016")

# excluding tt0 measurements from 2016, because not measured in 2015?
#Reco_1516new<- Reco_1516[Reco_1516$treatment %in% c("c", "otc", "tt1", "tt2", "tt3", "tt4"), ]
# doesn't make much difference

#subset different elevations
Reco_3000<- subset(Reco_1516, origin == "elev3000")
Reco_3500<- subset(Reco_1516, origin == "elev3500")
Reco_3850<- subset(Reco_1516, origin == "elev3850")
Reco_4100<- subset(Reco_1516, origin == "elev4100")

# summary of data
x<-group_by(Reco_16, site & treatment)
# calculate mean, sd and CV per elevation and treatment
y<-summarise(x, mean(Reco_lm), sd(Reco_lm), cv(Reco_lm))

#heteroscedasticity!!


# calculate means for each treatment at different elevation (site)
# comparing measurement to their min/max limit, min_limit = mean_Reco*0.5 and max_limit = mean_Reco*2. If exceeding report FALSE 
x<- Reco_16 %>% 
  group_by(treatment, origin) %>%#site,
  mutate(mn = mean(Reco_lm)) %>%
  mutate(keep = Reco_lm < mn *2 & Reco_lm> mn*0.5)

means<-summarise(x, mean(Reco_lm))
colnames(means)[3] <- "mean_Reco"
summarise(x, var(Reco_lm))


#code not working from mutate because of NA's in means
means %>% 
  spread(key = treatment, value = mean_Reco) %>% # spread Treatments
  mutate( V1 = tt0 - c, V2 = otc - c, V3 = tt1 - c, V4 = tt2 – c, V5 = tt3 – c, V6 = tt4 – c) # %>% #Difference Treatment-Control , V3 = tt1 - c, V4 = tt2 – c, V5 = tt3 – c, V6 = tt4 – c
  #gather(key = T=treatment, value = value, -Control, -...)

means <- join.exp %>% mutate(diff = Expenditure/Inflation)



