##  This is a read me/executable file for running the nee_calc_china_2016.R script for fitting LiCOR CO2 Flux data to both linear and exponential (leaky-fit) models.  For particulars on how the fits are being performed, open and read through the nee_calc_china_2016.R script.  Otherwise, follow the steps listed below for executing.

## FIRST: Execute the nee_calc_china_2016.R file in RStudio so that the function nee_calc_china_2016() exists in your Global Environment.  IMPORTANT!!! THE TENT VOLUME AND GAS EXCHANGE AREA IS FIXED INSIDE THE CODE.  You should check that the volume and area values match those of your tent, and if they do not then correct them inside the code.  See lines 40 and 41.

## SECOND: Run the next line of code below.  You will shortly be promted to set the working directory.  A drop-down/pop-up/toggle box will open.  Use this to navigate to the folder in which your LiCOR .csv files are stored and awaiting analysis.  Select any file in the folder and the folder's directory will be used.

## You will next be asked to define the starting and finishing time intervals, and finally whether or not you would like to re-run the fitting algorithms with different time intervals.  Upon completion, a .csv file will be printed to the working directory, or the same location where your LiCOR files are stored.

## In the current version of this code, once you start analyzing files you can not stop until completed.  That is, partial prints of analysis summaries will not be made.  So it is best to analyze your files in smaller and more manageable chunks.

## Lastly, if there is already a summary spreadsheet that has been printed to your directory (for example, you are running the code a second time through), then remove the summary spreadsheet so that your directory contains only the processed data sets.

nee_calc_china_2016()

## Updated 11/14/16