################################################################################
# SA TB Model Version 01 - Parameter List Generation
# Authors: DTTC-SACEMA TB Modelling Group
################################################################################

# This script creates a list of parameter sets, where each set is sampled from 
# prior distributions, if applicable

################################################################################

# Load required libraries
suppressPackageStartupMessages({
	library(readxl)
})

# Define file paths and number of model iterations
.args <- if(interactive()){
	c("INPUTS/RSA/RSA_inputdata.xlsx",
		"1000",                               
		"OUTPUTS/RSA/initialParametersets.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}

# Set number of model iterations
nrowparammat <- as.numeric(.args[2])

# Read parameter names and prior information from Excel file
PAR <- read_excel(.args[1], 
									sheet = "2_generate_param_list",
									range = "A1:E200",  
									col_types = c("numeric", "text", "numeric", "numeric", "numeric"),
									col_names = T)


# Identify the first empty row, if any, and subset data
first_empty_row <- which(rowSums(is.na(PAR)) == ncol(PAR))[1]
if (is.na(first_empty_row)) {
} else {
  PAR <- PAR[1:(first_empty_row - 1), ]
}


# Initialize empty parameter matrix
ParamMat= data.frame(matrix(ncol = length(PAR$code), nrow = nrowparammat))
colnames(ParamMat)<- PAR$code


# Set seed for reproducibility
set.seed(432784)

# Populate ParamMat with samples from each parameter's prior distribution
for(i in PAR$code){
	local <- PAR[PAR$code==i,]
	if(is.na(local$low)){
		ParamMat[,local$code]= rep(local$best, nrowparammat)
	}  
	else  {
		ParamMat[,local$code]= runif(nrowparammat, local$low,  local$high)
	}
}


# Convert parameter matrix to list of parameter sets
ParamList= vector("list", nrowparammat)
for(i in 1:nrowparammat){
	ParamList[[i]]= c(ParamMat[i,])
}


# Read time vectors for model calibration and projections
years_time_variables <- read_excel(.args[1],  
																	 sheet = "runmodel_years",
																	 range = "A1:B5",       
																	 col_types = c("text", "numeric"), 
																	 col_names = T)

# Define model start and end times for calibration and projections
time_modelstart      <- years_time_variables$value[2]   
time_endcalibration  <- years_time_variables$value[3]  
time_endprojections  <- years_time_variables$value[4]  

timesvect_calibration <- seq(time_modelstart, time_endcalibration, by=1)
timesvect_projections <- seq(time_modelstart, time_endprojections, by=1)

# Save outputs
save(timesvect_calibration, file = "./OUTPUTS/RSA/timesvect_calibration.RData")
save(timesvect_projections, file = "./OUTPUTS/RSA/timesvect_projections.RData")
saveRDS(ParamList, file=tail(.args, 1))
