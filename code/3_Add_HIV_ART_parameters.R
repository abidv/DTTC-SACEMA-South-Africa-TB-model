################################################################################
# SA TB Model Version 01 - Extend parameter list to inlcude HIV/ ART
# Authors: DTTC-SACEMA TB Modelling Group
################################################################################

# This script extends the initial parameter list by incorporating the effects of HIV 
# and antiretroviral therapy (ART)

################################################################################

# Load required libraries
suppressPackageStartupMessages({
	library(readxl)
})

# Define file paths and number of model iterations
.args <- if(interactive()){
	c("OUTPUTS/RSA/initialParametersets.RData",
		"INPUTS/RSA/RSA_inputdata.xlsx",
		"OUTPUTS/RSA/Parametersets_added_HIV_ART.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}


# Load initial parameter list and rename as list 2
ParamList <- readRDS(.args[1])
ParamList2 <- ParamList


# Build ParameterList2, extending with ART and HIV parameters
for(i in 1:length(ParamList)){
	
  #read in effect paramaters
	RR1= ParamList[[i]]$hivtb1
	RR2= ParamList[[i]]$hivtb2
	h1RR2pwr= ParamList[[i]]$h1RR2pwr
	h2RR2pwr= ParamList[[i]]$h2RR2pwr
	h3RR2pwr= ParamList[[i]]$h3RR2pwr
	h4RR2pwr= ParamList[[i]]$h4RR2pwr
	
	
	#### Effect of HIV
	
	rapprog_h1 <- ParamList[[i]]$rapprog_h0 * RR1 * (RR2^h1RR2pwr)
	rapprog_h2 <- ParamList[[i]]$rapprog_h0 * RR1 * (RR2^h2RR2pwr) 
	rapprog_h3 <- ParamList[[i]]$rapprog_h0 * RR1 * (RR2^h3RR2pwr)
	rapprog_h4 <- ParamList[[i]]$rapprog_h0 * RR1 * (RR2^h4RR2pwr)
	
	if(rapprog_h1>0.99) {rapprog_h1=0.99}  # Cap values at 0.99
	if(rapprog_h2>0.99) {rapprog_h2=0.99}
	if(rapprog_h3>0.99) {rapprog_h3=0.99}
	if(rapprog_h4>0.99) {rapprog_h4=0.99}
	
	delprog_h1 <- ParamList[[i]]$delprog_h1   
	delprog_h2 <- ParamList[[i]]$delprog_h2   
	delprog_h3 <- ParamList[[i]]$delprog_h3   
	delprog_h4 <- ParamList[[i]]$delprog_h4  
	
	relreinflat_h1 = 1 - ((1-ParamList[[i]]$relreinflat_h0) / (RR1 * (RR2^h1RR2pwr)))
	relreinflat_h2 = 1 - ((1-ParamList[[i]]$relreinflat_h0) / (RR1 * (RR2^h2RR2pwr))) 
	relreinflat_h3 = 1 - ((1-ParamList[[i]]$relreinflat_h0) / (RR1 * (RR2^h3RR2pwr)))
	relreinflat_h4 = 1 - ((1-ParamList[[i]]$relreinflat_h0) / (RR1 * (RR2^h4RR2pwr))) 
	
	relreinfact_h1 = ParamList[[i]]$relreinflat_h1
	relreinfact_h2 = ParamList[[i]]$relreinflat_h2
	relreinfact_h3 = ParamList[[i]]$relreinflat_h3 
	relreinfact_h4 = ParamList[[i]]$relreinflat_h4 
	
	hrelapse_h1 = ParamList[[i]]$hrelapse_h0
	hrelapse_h2 = ParamList[[i]]$hrelapse_h0
	hrelapse_h3 = ParamList[[i]]$hrelapse_h0
	hrelapse_h4 = ParamList[[i]]$hrelapse_h0
	
	relapse_h1 <- ParamList[[i]]$relapse_h1 
	relapse_h2 <- ParamList[[i]]$relapse_h2   
	relapse_h3 <- ParamList[[i]]$relapse_h3  
	relapse_h4 <- ParamList[[i]]$relapse_h4  
	
	recovhl_h1 = ParamList[[i]]$recovhl_h0
	recovhl_h2 = ParamList[[i]]$recovhl_h0
	recovhl_h3 = ParamList[[i]]$recovhl_h0
	recovhl_h4 = ParamList[[i]]$recovhl_h0
	recovhl_h1a = ParamList[[i]]$recovhl_h0
	recovhl_h2a = ParamList[[i]]$recovhl_h0
	recovhl_h3a = ParamList[[i]]$recovhl_h0
	recovhl_h4a = ParamList[[i]]$recovhl_h0
	
	progsymp_h1 <- ParamList[[i]]$progsymp_h0 * (RR1 * (RR2^h1RR2pwr))
	progsymp_h2 <- ParamList[[i]]$progsymp_h0 * (RR1 * (RR2^h2RR2pwr)) 
	progsymp_h3 <- ParamList[[i]]$progsymp_h0 * (RR1 * (RR2^h3RR2pwr))
	progsymp_h4 <- ParamList[[i]]$progsymp_h0 * (RR1 * (RR2^h4RR2pwr))
	
	#### Effect of ART
	
	relbeta_h1a = 1 - ((1-ParamList[[i]]$relbeta_h1) * (1-(ParamList[[i]]$arteff)))
	relbeta_h2a = 1 - ((1-ParamList[[i]]$relbeta_h2) * (1-(ParamList[[i]]$arteff)))
	relbeta_h3a = 1 - ((1-ParamList[[i]]$relbeta_h3) * (1-(ParamList[[i]]$arteff)))
	relbeta_h4a = 1 - ((1-ParamList[[i]]$relbeta_h4) * (1-(ParamList[[i]]$arteff)))
	
	rapprog_h1a <- ParamList[[i]]$rapprog_h0 + (rapprog_h1-ParamList[[i]]$rapprog_h0) * (1-(ParamList[[i]]$arteff))
	rapprog_h2a <- ParamList[[i]]$rapprog_h0 + (rapprog_h2-ParamList[[i]]$rapprog_h0) * (1-(ParamList[[i]]$arteff))
	rapprog_h3a <- ParamList[[i]]$rapprog_h0 + (rapprog_h3-ParamList[[i]]$rapprog_h0) * (1-(ParamList[[i]]$arteff))
	rapprog_h4a <- ParamList[[i]]$rapprog_h0 + (rapprog_h4-ParamList[[i]]$rapprog_h0) * (1-(ParamList[[i]]$arteff))
	
	delprog_h1a <- ParamList[[i]]$delprog_h0 + (delprog_h1-ParamList[[i]]$delprog_h0) * (1-(ParamList[[i]]$arteff))
	delprog_h2a <- ParamList[[i]]$delprog_h0 + (delprog_h2-ParamList[[i]]$delprog_h0) * (1-(ParamList[[i]]$arteff))
	delprog_h3a <- ParamList[[i]]$delprog_h0 + (delprog_h3-ParamList[[i]]$delprog_h0) * (1-(ParamList[[i]]$arteff))
	delprog_h4a <- ParamList[[i]]$delprog_h0 + (delprog_h4-ParamList[[i]]$delprog_h0) * (1-(ParamList[[i]]$arteff))
	
	hrelapse_h1a = ParamList[[i]]$hrelapse_h0
	hrelapse_h2a = ParamList[[i]]$hrelapse_h0
	hrelapse_h3a = ParamList[[i]]$hrelapse_h0
	hrelapse_h4a = ParamList[[i]]$hrelapse_h0
	
	relreinflat_h1a <-  ParamList[[i]]$relreinflat_h0 + (relreinflat_h1 - ParamList[[i]]$relreinflat_h0) * (1-(ParamList[[i]]$arteff))
	relreinflat_h2a <-  ParamList[[i]]$relreinflat_h0 + (relreinflat_h2 - ParamList[[i]]$relreinflat_h0) * (1-(ParamList[[i]]$arteff))
	relreinflat_h3a <-  ParamList[[i]]$relreinflat_h0 + (relreinflat_h3 - ParamList[[i]]$relreinflat_h0) * (1-(ParamList[[i]]$arteff))
	relreinflat_h4a <-  ParamList[[i]]$relreinflat_h0 + (relreinflat_h4 - ParamList[[i]]$relreinflat_h0) * (1-(ParamList[[i]]$arteff))
	
	relreinfact_h1a = relreinflat_h1a
	relreinfact_h2a = relreinflat_h2a
	relreinfact_h3a = relreinflat_h3a 
	relreinfact_h4a = relreinflat_h4a 
	
	progsymp_h1a <-  ParamList[[i]]$progsymp_h0 + (progsymp_h1 - ParamList[[i]]$progsymp_h0) * (1-(ParamList[[i]]$arteff))
	progsymp_h2a <-  ParamList[[i]]$progsymp_h0 + (progsymp_h2 - ParamList[[i]]$progsymp_h0) * (1-(ParamList[[i]]$arteff))
	progsymp_h3a <-  ParamList[[i]]$progsymp_h0 + (progsymp_h3 - ParamList[[i]]$progsymp_h0) * (1-(ParamList[[i]]$arteff))
	progsymp_h4a <-  ParamList[[i]]$progsymp_h0 + (progsymp_h4 - ParamList[[i]]$progsymp_h0) * (1-(ParamList[[i]]$arteff))
	
	
	relapse_h1a <- ParamList[[i]]$relapse_h0 + (relapse_h1-ParamList[[i]]$relapse_h0) * (1-(ParamList[[i]]$arteff))
	relapse_h2a <- ParamList[[i]]$relapse_h0 + (relapse_h2-ParamList[[i]]$relapse_h0) * (1-(ParamList[[i]]$arteff))
	relapse_h3a <- ParamList[[i]]$relapse_h0 + (relapse_h3-ParamList[[i]]$relapse_h0) * (1-(ParamList[[i]]$arteff))
	relapse_h4a <- ParamList[[i]]$relapse_h0 + (relapse_h4-ParamList[[i]]$relapse_h0) * (1-(ParamList[[i]]$arteff))
	
	tbrecov_h1a <- ParamList[[i]]$tbrecov_art  
	tbrecov_h2a <- ParamList[[i]]$tbrecov_art
	tbrecov_h3a <- ParamList[[i]]$tbrecov_art
	tbrecov_h4a <- ParamList[[i]]$tbrecov_art 
	
	# Additional parameters for governing diagnosis 
	accdiagp_h1 = ParamList[[i]]$accdiagp_hx
	accdiagp_h2 = ParamList[[i]]$accdiagp_hx
	accdiagp_h3 = ParamList[[i]]$accdiagp_hx
	accdiagp_h4 = ParamList[[i]]$accdiagp_hx
	accdiagp_h1a = ParamList[[i]]$accdiagp_hx
	accdiagp_h2a = ParamList[[i]]$accdiagp_hx
	accdiagp_h3a = ParamList[[i]]$accdiagp_hx
	accdiagp_h4a = ParamList[[i]]$accdiagp_hx
	
	accdiags_h1 = ParamList[[i]]$accdiags_hx
	accdiags_h2 = ParamList[[i]]$accdiags_hx
	accdiags_h3 = ParamList[[i]]$accdiags_hx
	accdiags_h4 = ParamList[[i]]$accdiags_hx
	accdiags_h1a = ParamList[[i]]$accdiags_hx
	accdiags_h2a = ParamList[[i]]$accdiags_hx
	accdiags_h3a = ParamList[[i]]$accdiags_hx
	accdiags_h4a = ParamList[[i]]$accdiags_hx

	# Append extended parameters to ParamList2
	ParamList2[[i]] <- c(ParamList2[[i]],
											 rapprog_h1= rapprog_h1,  
											 rapprog_h2= rapprog_h2,
											 rapprog_h3= rapprog_h3,
											 rapprog_h4= rapprog_h4,
											 delprog_h1= delprog_h1,  
											 delprog_h2= delprog_h2,
											 delprog_h3= delprog_h3,
											 delprog_h4= delprog_h4,
											 relreinflat_h1= relreinflat_h1,  
											 relreinflat_h2= relreinflat_h2,
											 relreinflat_h3= relreinflat_h3,
											 relreinflat_h4= relreinflat_h4,
											 relreinfact_h1 = relreinfact_h1, 
											 relreinfact_h2 = relreinfact_h2, 
											 relreinfact_h3 = relreinfact_h3, 
											 relreinfact_h4 = relreinfact_h4, 
											 relreinfact_h1a = relreinfact_h1a,
											 relreinfact_h2a = relreinfact_h2a, 
											 relreinfact_h3a = relreinfact_h3a, 
											 relreinfact_h4a = relreinfact_h4a,
											 hrelapse_h1 = hrelapse_h1, 
											 hrelapse_h2 = hrelapse_h2,
											 hrelapse_h3 = hrelapse_h3,
											 hrelapse_h4 = hrelapse_h4,
											 hrelapse_h1a = hrelapse_h1a,
											 hrelapse_h2a = hrelapse_h2a,
											 hrelapse_h3a = hrelapse_h3a,
											 hrelapse_h4a = hrelapse_h4a,
											 recovhl_h1 = recovhl_h1, 
											 recovhl_h2 = recovhl_h2,
											 recovhl_h3 = recovhl_h3,
											 recovhl_h4 = recovhl_h4,
											 recovhl_h1a = recovhl_h1a,
											 recovhl_h2a = recovhl_h2a,
											 recovhl_h3a = recovhl_h3a,
											 recovhl_h4a = recovhl_h4a,
											 progsymp_h1 = progsymp_h1, 
											 progsymp_h2 = progsymp_h2,
											 progsymp_h3 = progsymp_h3,
											 progsymp_h4 = progsymp_h4,
											 relbeta_h1a = relbeta_h1a,
											 relbeta_h2a = relbeta_h2a,
											 relbeta_h3a = relbeta_h3a,
											 relbeta_h4a = relbeta_h4a,
											 rapprog_h1a= rapprog_h1a,  
											 rapprog_h2a= rapprog_h2a,
											 rapprog_h3a= rapprog_h3a,
											 rapprog_h4a= rapprog_h4a,
											 delprog_h1a= delprog_h1a,  
											 delprog_h2a= delprog_h2a,
											 delprog_h3a= delprog_h3a,
											 delprog_h4a= delprog_h4a,
											 relreinflat_h1a= relreinflat_h1a,  
											 relreinflat_h2a= relreinflat_h2a,
											 relreinflat_h3a= relreinflat_h3a,
											 relreinflat_h4a= relreinflat_h4a,
											 progsymp_h1a = progsymp_h1a,
											 progsymp_h2a = progsymp_h2a,
											 progsymp_h3a = progsymp_h3a,
											 progsymp_h4a = progsymp_h4a,
											 relapse_h1 = relapse_h1,
											 relapse_h2 = relapse_h2,
											 relapse_h3 = relapse_h3,
											 relapse_h4 = relapse_h4,
											 relapse_h1a = relapse_h1a,
											 relapse_h2a = relapse_h2a,
											 relapse_h3a = relapse_h3a,
											 relapse_h4a = relapse_h4a,
											 tbrecov_h1a = tbrecov_h1a,
											 tbrecov_h2a = tbrecov_h2a,
											 tbrecov_h3a = tbrecov_h3a,
											 tbrecov_h4a = tbrecov_h4a,
											 accdiagp_h1 = accdiagp_h1,
											 accdiagp_h2 = accdiagp_h2,
											 accdiagp_h3 = accdiagp_h3,
											 accdiagp_h4 = accdiagp_h4,
											 accdiagp_h1a = accdiagp_h1a,
											 accdiagp_h2a = accdiagp_h2a,
											 accdiagp_h3a = accdiagp_h3a,
											 accdiagp_h4a = accdiagp_h4a,
											 accdiags_h1 = accdiags_h1,
											 accdiags_h2 = accdiags_h2,
											 accdiags_h3 = accdiags_h3,
											 accdiags_h4 = accdiags_h4,
											 accdiags_h1a = accdiags_h1a,
											 accdiags_h2a = accdiags_h2a,
											 accdiags_h3a = accdiags_h3a,
											 accdiags_h4a = accdiags_h4a,
											 tbmort_h1a = tbmort_h1a,
											 tbmort_h2a = tbmort_h2a,
											 tbmort_h3a = tbmort_h3a,
											 tbmort_h4a = tbmort_h4a,
											 tbtreatmort_h1a = tbtreatmort_h1a,
											 tbtreatmort_h2a = tbtreatmort_h2a,
											 tbtreatmort_h3a = tbtreatmort_h3a,
											 tbtreatmort_h4a = tbtreatmort_h4a
											 
	)

}

# Save updated parameter list
saveRDS(ParamList2, file=tail(.args, 1))

