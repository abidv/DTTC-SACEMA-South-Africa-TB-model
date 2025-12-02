################################################################################
# SA TB Model Version 01 - Read in time-varying data
# Authors: DTTC-SACEMA TB Modelling Group
################################################################################

# This script reads in time-varying parameters related to HIV incidence from 
# and uses an exponential decay model to forecast HIV incidence from 2020 to 2030. 

################################################################################


# Load required libraries
suppressPackageStartupMessages({
	library(readxl)
})


# Define input and output file paths
.args <- if(interactive()){
	c("INPUTS/RSA/RSA_inputdata.xlsx",
		"OUTPUTS/RSA/tv_HIV_ART_data.RDS")
}else{
	commandArgs(trailingOnly = TRUE)
}


# Read time-varying parameter data
timevarpar <- read_excel(.args[1],  
												 sheet = "1_timevarying_parms",
												 range = "A1:E37",       
												 col_types = c("numeric", "numeric", "numeric", "numeric","numeric"), 
												 col_names = T)


# Create exponential decay model to forecast HIV incidence from 2020 to 2030
hivdata = data.frame(t=2010:2019, y=timevarpar$hivinc[timevarpar$calyear>=2010 & timevarpar$calyear<=2019])
hivincmodel = lm(log(y) ~ t, data=hivdata)
timevarpar$hivinc[timevarpar$calyear>=2020 &timevarpar$calyear<=2030] <- exp(predict(hivincmodel, newdata=data.frame(t=2020:2030)))


# Save updated time-varying parameters to RDS
saveRDS(timevarpar, file = tail(.args, 1))
