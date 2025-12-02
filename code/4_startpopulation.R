################################################################################
# SA TB Model Version 01 - Define the starting population for 1995
# Authors: SACEMA-DTTC Modelling Group
################################################################################

# This script defines the starting population in 1995 for the TB model,
# allowing for variation in the HIV-positive population based on parameter sets.

################################################################################

# Load required libraries
suppressPackageStartupMessages({
	library(readxl)
})

# Define input and output file paths
.args <- if(interactive()){
	c("OUTPUTS/RSA/Parametersets_added_HIV_ART.RData",
		"100000",
		"OUTPUTS/RSA/yinit.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}

ParamList2 <- readRDS(.args[1])
nrowparammat <- .args[2]


# Define empty list for population sizes
yinit= list()

# Sampling to allow for variation in the 1995 starting population for each model iteration
for(i in 1:nrowparammat){
	parm = ParamList2[[i]]   
	parm <- c(parm['prop_h1'], parm['prop_h2'], parm['prop_h3'], parm['prop_h4'],   
						parm['prop_S'], parm['prop_LR'], parm['prop_LD'], parm['prop_IP'], 
						parm['prop_IS'], parm['prop_DP'], parm['prop_DS'], parm['prop_FN'], 
						parm['prop_T'], parm['prop_RH'], parm['prop_R'], 
						parm['pop1995_total'],parm['pop1995_hiv'])
	parm = unlist(parm, use.names=TRUE)
	
	# Normalize HIV population proportions
	hivsum = (parm[['prop_h1']] + parm[['prop_h2']] + parm[['prop_h3']] + parm[['prop_h4']])    
	parm['prop_h1'] = parm[['prop_h1']]/hivsum 
	parm['prop_h2'] = parm[['prop_h2']]/hivsum
	parm['prop_h3'] = parm[['prop_h3']]/hivsum 
	parm['prop_h4'] = parm[['prop_h4']]/hivsum 
	
	
	# Calculate LD proportion as 1 - 'all other proportions'  
	parm['prop_LD'] = 1 - (parm[['prop_S']]+ parm[['prop_LR']]+ parm[['prop_IP']] + parm[['prop_IS']] + parm[['prop_DP']]
												 + parm[['prop_DS']] + parm[['prop_FN']] + parm[['prop_T']] + parm[['prop_RH']] + parm[['prop_R']]) 
	
	# Determine HIV-negative population as total minus positives
	pop1995 = data.frame(Total=parm['pop1995_total'],hivpos=parm['pop1995_hiv'])
	pop1995$hivneg  =pop1995$Total - pop1995$hivpos
	
	# Initialize population structure for the iteration
	yinit[[i]] = c(
		S_h0 = pop1995$hivneg * unname(parm["prop_S"]),       
		LR_h0 = pop1995$hivneg * unname(parm["prop_LR"]),
		LD_h0 = pop1995$hivneg * unname(parm["prop_LD"]),
		IP_h0 = pop1995$hivneg * unname(parm["prop_IP"]),
		IS_h0 = pop1995$hivneg * unname(parm["prop_IS"]),
		DP_h0 = pop1995$hivneg * unname(parm["prop_DP"]),
		DS_h0 = pop1995$hivneg * unname(parm["prop_DS"]),
		FN_h0 = pop1995$hivneg * unname(parm["prop_FN"]),
		T_h0 = pop1995$hivneg * unname(parm["prop_T"]),
		RH_h0 = pop1995$hivneg * unname(parm["prop_RH"]),
		R_h0 = pop1995$hivneg * unname(parm["prop_R"]),
		
		S_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_S"]),
		LR_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_LR"]),
		LD_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_LD"]),
		IP_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_IP"]),
		IS_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_IS"]),
		DP_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_DP"]),
		DS_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_DS"]),
		FN_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_FN"]),
		T_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_T"]),
		RH_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_RH"]),
		R_h1 = pop1995$hivpos * unname(parm["prop_h1"]) * unname(parm["prop_R"]),
		
		S_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_S"]),
		LR_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_LR"]), 
		LD_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_LD"]),
		IP_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_IP"]),
		IS_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_IS"]),
		DP_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_DP"]),
		DS_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_DS"]),
		FN_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_FN"]),
		T_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_T"]),
		RH_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_RH"]),
		R_h2 = pop1995$hivpos * unname(parm["prop_h2"]) * unname(parm["prop_R"]),
		
		S_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_S"]),
		LR_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_LR"]), 
		LD_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_LD"]),
		IP_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_IP"]),
		IS_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_IS"]),
		DP_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_DP"]),
		DS_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_DS"]),
		FN_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_FN"]),
		T_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_T"]),
		RH_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_RH"]),
		R_h3 = pop1995$hivpos * unname(parm["prop_h3"]) * unname(parm["prop_R"]),
		
		S_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_S"]),
		LR_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_LR"]), 
		LD_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_LD"]),
		IP_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_IP"]),
		IS_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_IS"]),
		DP_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_DP"]),
		DS_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_DS"]),
		FN_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_FN"]),
		T_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_T"]),
		RH_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_RH"]),
		R_h4 = pop1995$hivpos * unname(parm["prop_h4"]) * unname(parm["prop_R"]),
		
		S_h1a = 0,
		LR_h1a = 0,
		LD_h1a = 0,
		IP_h1a = 0,
		IS_h1a = 0,
		DP_h1a = 0,
		DS_h1a = 0,
		FN_h1a = 0,
		T_h1a = 0,
		RH_h1a = 0,
		R_h1a = 0,
		
		S_h2a = 0,
		LR_h2a = 0,
		LD_h2a = 0,
		IP_h2a = 0,
		IS_h2a = 0,
		DP_h2a = 0,
		DS_h2a = 0,
		FN_h2a = 0,
		T_h2a = 0,
		RH_h2a = 0,
		R_h2a = 0,
		
		S_h3a = 0,
		LR_h3a = 0,
		LD_h3a = 0,
		IP_h3a = 0,
		IS_h3a = 0,
		DP_h3a = 0,
		DS_h3a = 0,
		FN_h3a = 0,
		T_h3a = 0,
		RH_h3a = 0,
		R_h3a = 0,
		
		S_h4a = 0,
		LR_h4a = 0,
		LD_h4a = 0,
		IP_h4a = 0,
		IS_h4a = 0,
		DP_h4a = 0,
		DS_h4a = 0,
		FN_h4a = 0,
		T_h4a = 0,
		RH_h4a = 0,
		R_h4a = 0,
		
		#These are pseudoparameters "for counters"
		
		CumTBnotif = 0,
		CumFP = 0, 
		CumTBinc = 0,
		CumTBinc_h0 = 0,
		CumTBmort = 0
	)
}


# Save start popaultion sizes 
save(yinit, file = tail(.args, 1))

