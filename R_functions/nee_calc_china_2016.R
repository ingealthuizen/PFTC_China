nee_calc_china_2016 <- function(){
  ##This function is meant to perform the linear and non-linear fitting to nee data.  It will produce a plot of the data, queery the user if they would like to modify the time interval over which the data is fit, then print the paste statement of data values.  Currently, the function is written so that it works for a given filename.  From there, we can extend to a folder/directory.
  
  ## These packages are necessary.
  
  library("proto")
  library("nls2")
  library("minpack.lm")
  
  ####----Insert here the readline() command for querying user for the director of LiCOR files.--------
  ## For querying user to define the working directory in which the LiCOR files are located.
  readline("Please set the working directory to the folder \n that contains the LiCOR files to be analyzed. \n Do so with the upcoming prompt. \n Note that you must choose a(any) file in the \n folder that you want to set as the working directory. \n Please press 'return' to continue.")
  setwd(dirname(file.choose()))
  
  ## Define directory as an object that contains the dir() item names as a vector.
  directory <- dir()
  
  ## For reading the .txt files, replace the ".txt" paste with ".txt"
  photo.names <- directory[grep(paste("[0-9]", ".csv", sep = ""), dir(), ignore.case = TRUE, value = FALSE)]
  # ambient.names <- directory[grep(paste("[0-9]", "a", ".txt", sep = ""), dir(), ignore.case = TRUE, value = FALSE)] ambient will come from the file name
  resp.names <- directory[grep("r", dir(), ignore.case = FALSE, value = FALSE)]
  
  nee.fit <- function(filename){
    
    ## For reading the .txt files, replace read.csv with read.table, and add a skip = 32 parameter to the function.
    input <- read.csv(filename, header = TRUE, sep = ",")[,-1]
    
    ## Quick check for NA's within input.  Due to how the climate variables were combined from different devices, some of the processed data sets have a row of NA's at the bottom (possible elsewhere in the middle too, but none have yet been observed).  The following code will search for rows with NA's and remove them, as the NA's are causing the non-linear fitting algorithms to fail.
    
    input <- input[which(!is.na(input$co2)),]
    
    ## Extract filename information.
    splitname <- strsplit(filename, split = "_")[[1]]
    
    ## Convert date time column.
    input$datetime <- as.POSIXct(input$datetime)

    #  /// define constants for your tent///
    vol = 0.025   # m^3, tent volume
    area = 0.0625   # m^2, tent area
    R = 8.314472   # J/mol K
    
    ## Define vectors to work with

    time <- c(1:length(input$co2)) #s
    co2 <- input$co2 #umol/mol
    h2o <- input$h2o #mmol/mol
    par <- input$par #units?
    press <- input$press #kPa
    temp <- input$temp #C
    
    #  /// Plotting the CO2 vs. Time Curve //

    par(mfrow = c(2,1))
    
    cprime <- co2/(1-(h2o/1000))

    ## Here is where different plots can be added or removed.

    plot(par ~ time, main = "par")
    plot(cprime ~ time, main = filename)

    
    ## Queery user for start time for fitting.  Default is set to 10 in the if() statement
    tstart <- readline("Enter preferred start time for fitting. \n Do not include units. \n Round to nearest integer second. \n Do not use 0. \n  If default of 10s is preferred, press 'return':")
    if(!grepl("^[0-9]+$", tstart)){
      tstart <- 10
    }
    tstart <- as.integer(tstart)
    
    ## Queery user for finish time for fitting.  Default is set to 80 in the if() statement
    tfinish <- readline("Enter preferred finish time for fitting. \n Do not include units. \n Round to nearest integer second. \n  If default of 80s is preferred, press 'return':")
    if(!grepl("^[0-9]+$", tfinish)){
      tfinish <- 80
    }
    tfinish <- as.integer(tfinish)
    
    ## Linear fitting code.
    linear.fit <- lm(cprime[tstart:tfinish]~(time[tstart:tfinish]))

    ## AIC value for linear model for later comparison with non-linear model.
    aic.lm <- AIC(linear.fit)

    # Calculate intercept
    inter<- as.numeric(linear.fit$coeff[1])

    # Calculate slope
    dcdt <- as.numeric(linear.fit$coeff[2])

    # Calculate r-squared (we're not reporting chi-squared significance from the non-linear fit, so this may not be so necessary)
    rsqd <- summary(linear.fit)$r.sq


    # Define these after defining the time interval...
    tav <- mean(temp[tstart:tfinish], na.rm = TRUE)
    pav <- mean(press[tstart:tfinish], na.rm = TRUE)
    wav <- mean(h2o[tstart:tfinish], na.rm = TRUE)
    parav <- mean(par[tstart:tfinish], na.rm = TRUE)

    #ambient
    cprime <- co2/(1-(h2o/1000))
    camb <- as.numeric(strsplit(splitname[6], "amb")[[1]][2])

    ## I think this is just the ideal gas law, but not clear to me why ideal gaw law is necessary for what appears to be a unit conversion.

    camb <- camb*R*(tav+273.15)/(44*pav)    # change to umol/mol or ppm

    # Calculate nee from linear model
    nee_lm <- -(vol*pav*(1000-wav)*dcdt) / (R*area*(tav + 273.15))	# in umol/m2/s
    
    # Make plot of line.
    abline(inter,dcdt, col=6)
    
    ## Non-linear fitting code.
    # Set cnot to the actual first value of cprime used in the fitting time domain.
    cnot = cprime[tstart]
    
    # Define a temporary data frame from which the functional variables come from.
    df = data.frame(cprime, time)
    
    # Define a subset category from the tstart and tfinish variables.
    subsettime <- time > tstart & time < tfinish
    
    # Define boundaries of parameter grid.
    strt <- data.frame(A = c(150, 850), B = c(-1000, 1000))
    
    # Use nls2() to scan through parameter grid, searching for "best" actual starting points.  control variable is set to prevent warnings from ending loop.
    optimize.start <- nls2(cprime ~ (cnot - A)*exp(-time/B) + A, data = df, start=strt, algorithm = "brute-force", control = nls.control(warnOnly = TRUE), subset = subsettime) #(A=375, B=40)
    
    # Run nls() with the optimized starting values from previous nls2().  Control variable is set to prevent warnings from ending loop.  However, they will still be printed at end of run.  When this happens, it is indicative of the fact that the function parameters (A and B) are large (non-physical) for the fitting, yet still produce a fit.  This is worth further investigation.  However, it appears that the nee value produced by the exponential model in such circumstances does not deviate from the linear model by much more than half a percent.  Add a "trace = TRUE" parameter setting to the nls() function to be able to watch the values of A and B change with each iteration.
    tryCatch({
    uptake.fm <- nlsLM(cprime ~ (cnot - A)*exp(-time/B) + A, data = df, start = coef(optimize.start), subset = subsettime, control = nls.control(warnOnly = TRUE))
      
    ##
    sigma <- summary(uptake.fm)$sigma

    ## AIC value for non-linear model for later comparison with linear model.
    aic.nlm <- AIC(uptake.fm)

    Css = summary(uptake.fm)$param[1]  
    tau = summary(uptake.fm)$param[2]
    nee_exp <- ((camb-Css)/(area*tau))*(vol*pav*(1000-wav)/(R*(tav + 273.15))) #equation 4 in Saleska 1999
    
    curve((cnot - Css)*exp(-(x-time[tstart])/tau) + Css, col = 4, add = TRUE)	#equation 3 in Saleska 1999 to plot for visual inspection.
    
    
    ## nee_lm formula produces a negative value of nee_lm, but we are interested in PLANT uptake, so the minus sign must be coerced to positive when saving the values to spreadsheet via the final paste() function call at end of script.
    
    if(length(grep("r.csv", filename, ignore.case = TRUE, value = FALSE)) == 1){
      type <- "resp"
    }else{type <- "photo"}
    
    print(data.frame("tstart" = tstart, "tfinish" = tfinish, "type" = type, "nee_lm" = nee_lm, "nee_exp" = nee_exp, "rsqd" = rsqd, "sigma" = sigma, "aic.lm" = aic.lm, "aic.nlm" = aic.nlm))
    replicate <- readline("Would you like to redo the fitting with \n a different time domain? (y/n)")
    if(replicate == "y"){
      nee.fit(filename)
    } else {
      datetime_info <- as.character(input$datetime[1])
      group <- splitname[1]
      site <- splitname[5]
      treatment <- splitname[7]
      if(type == "resp"){
        block <- strsplit(splitname[length(splitname)], "[a-z].csv")[[1]]
      }else if(type == "photo"){
        block <- strsplit(splitname[length(splitname)], ".csv")[[1]]
      }
      return(c(datetime_info,tstart,tfinish,type,group,site,treatment,block,camb,tav,pav,wav,parav,nee_lm,nee_exp, rsqd,sigma, aic.lm, aic.nlm))
    }
    }, error = function(err) {
      print(paste("Critical Non-Linear Fitting Error: ", err))
      replicate <- readline("Would you like to retry the fitting with \n a different time domain? (y/n) \n A 'no' answer results in no non-linear fitting results being recorded")
      if(replicate == "y"){
        nee.fit(filename)
    } else {
      sigma <- NA
      aic.nlm <- NA
      nee_exp <- NA
      if(length(grep("r.csv", filename, ignore.case = TRUE, value = FALSE)) == 1){
        type <- "resp"
      }else{type <- "photo"}
      datetime_info <- as.character(input$datetime[1])
      group <- splitname[1]
      site <- splitname[5]
      treatment <- splitname[7]
      if(type == "resp"){
        block <- strsplit(splitname[length(splitname)], "[a-z].csv")[[1]]
      }else if(type == "photo"){
        block <- strsplit(splitname[length(splitname)], ".csv")[[1]]
      }
      return(c(datetime_info,tstart,tfinish,type,group,site,treatment,block,camb,tav,pav,wav,parav,nee_lm,nee_exp, rsqd,sigma, aic.lm, aic.nlm))
    }
    })
  }
  
  stats.df <- c()
  if(length(photo.names)>0) {
    for (i in 1:length(photo.names)) {
      stats.df <- rbind(stats.df, nee.fit(photo.names[i]))
      
    }
  }
  
  if (length(resp.names) > 1){
    for (i in 1:length(resp.names)){
      stats.df <- rbind(stats.df, nee.fit(resp.names[i]))
    }
  }

  stats.df <- as.data.frame(stats.df)
  names.vec <- c("datetime","tstart (s)", "tfinish (s)", "type", "group", "site", "treatment", "block", "camb (umol/mol)", "tav (K)", "pav (kPA)", "wav (umol/mol)", "parav (umol/s/m^2)", "nee_lm (umol/s/m^2)", "nee_exp (umol/s/m^2)", "rsqd", "sigma", "aic.lm", "aic.nlm")
  for(i in 1:length(names.vec)){
    names(stats.df)[i] <- names.vec[i]
  }
  
  stats.df
  write.csv(stats.df, file = paste(paste(strsplit(getwd(), "/")[[1]][length(strsplit(getwd(), "/")[[1]])], "summary", sep = " "), ".csv", sep = ""))
  
}

## Updated 11/14/16