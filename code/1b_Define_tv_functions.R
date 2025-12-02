################################################################################
# SA TB Model Version 01 - Define time-varying functions
# Authors: SACEMA-DTTC Modelling Group
################################################################################

# This script defines the functions governing the equations for calculating time-varying 
# parameters: HIV incidence, ART initiation rates, ART initiation upon initiating TB treatment, 
# and probability of HIV at birth.

################################################################################


# Define input and output file paths
.args <- if(interactive()){
	c("OUTPUTS/RSA/tv_HIV_ART_data.RDS",
		"OUTPUTS/RSA/tv_functions.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}

# Load time-varying parameters
timevarpar <- readRDS(.args[1])

# Function to return time-varying HIV incidence
hivinc_func <- function(times){
	roundtime = floor(times)
	whichyear= which(timevarpar$calyear==roundtime)
	return(as.numeric(timevarpar[whichyear, "hivinc"]))
}

# Function to return time-varying ART initiation when people are detected for TB
artinit_tb_func <- function(times){
	roundtime = floor(times)
	whichyear= which(timevarpar$calyear==roundtime)
	return(as.numeric(timevarpar[whichyear, "artinit_tb"]))
}

# Function to return time-varying probability of HIV around birth
prob_hiv_func <- function(times){
	roundtime = floor(times)
	whichyear= which(timevarpar$calyear==roundtime)
	return(as.numeric(timevarpar[whichyear, "probhivbirth"]))
}


# Function to return time-varying probability of ART initiation
artgen_general = function(times, L, yearinit, k, mR){
	if(times<2004){return(artinit=0)}
	else{
		artinit = mR + (L-mR)/(1+exp(-k*(times-yearinit)))
		return(artinit)
	}
}

# Save the functions 
save(artgen_general, artinit_tb_func, hivinc_func, prob_hiv_func, file = tail(.args, 1))
