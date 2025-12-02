################################################################################################################################
############################### Authors: SACEMA-DTTC Modelling group ###########################################################
################################################################################
# SA TB Model Version 01 - Run model in parallel
# Authors: SACEMA-DTTC Modelling Group
################################################################################

# This script runs the model in parallel across multiple cores using all the previously defined parameters and functions. 
# A function is created to run the model and produce a results list of all key TB indicators and outcome measures.

################################################################################

# Load necessary libraries
suppressPackageStartupMessages({
	library(foreach)
	library(doParallel)
	library(deSolve)
})


# Define input arguments
.args <- if(interactive()){
	c("OUTPUTS/RSA/timesvect_calibration.RData",
		"OUTPUTS/RSA/tv_HIV_ART_data.RDS",
		"OUTPUTS/RSA/tv_functions.RData",
		"OUTPUTS/RSA/cc_functions.RData",
		"OUTPUTS/RSA/TB_model.RData",
		"OUTPUTS/RSA/yinit.RData",
		"OUTPUTS/RSA/Parametersets_added_HIV_ART.Rdata",
		"1000",
		"OUTPUTS/RSA/results_tmp.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}



# Load necessary dependencies

load(.args[1])
timevarpar <- readRDS(.args[2])
load(.args[3])
load(.args[4])
load(.args[5])
load(.args[6])
ParamList2 <- readRDS(.args[7]) 
nrowparammat <- as.numeric(.args[8])

# Create result lis
resultList= list()


# Define a function to run the model
funcMakeResults <- function(model_version, timesvect, yinit_list, param_list){
  
	# Run the model using the ODE solver
  results <- as.data.frame(ode(func=model_version,
                               times=timesvect, 
                               y= yinit_list[[i]],
                               parms=param_list[[i]], 
                               method='rk', hini=0.05))  #hini = adapts the step size
  
  
  # Calculate total population size (N) 
  results$N= results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 + results$RH_h0 + results$R_h0 +
    results$S_h1 + results$LR_h1 + results$LD_h1 + results$IP_h1 + results$IS_h1 + results$DP_h1 + results$DS_h1 + results$FN_h1 + results$T_h1 + results$RH_h1 +  results$R_h1 +
    results$S_h2 + results$LR_h2 + results$LD_h2 + results$IP_h2 + results$IS_h2 + results$DP_h2 + results$DS_h2 + results$FN_h2 + results$T_h2 + results$RH_h2 +  results$R_h2 +
    results$S_h3 + results$LR_h3 + results$LD_h3 + results$IP_h3 + results$IS_h3 + results$DP_h3 + results$DS_h3 + results$FN_h3 + results$T_h3 + results$RH_h3 +  results$R_h3 +
    results$S_h4 + results$LR_h4 + results$LD_h4 + results$IP_h4 + results$IS_h4 + results$DP_h4 + results$DS_h4 + results$FN_h4 + results$T_h4 + results$RH_h4 +  results$R_h4 +
    results$S_h1a + results$LR_h1a + results$LD_h1a + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a + results$T_h1a + results$RH_h1a +  results$R_h1a +
    results$S_h2a + results$LR_h2a + results$LD_h2a + results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$T_h2a + results$RH_h2a +  results$R_h2a +
    results$S_h3a + results$LR_h3a + results$LD_h3a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$T_h3a + results$RH_h3a +  results$R_h3a +
    results$S_h4a + results$LR_h4a + results$LD_h4a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a + results$T_h4a + results$RH_h4a +  results$R_h4a
  
  # Population size (N_h0)  
  results$N_h0= results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 +  results$RH_h0 + results$R_h0
  
  # Population size (N_hx)  
  results$N_hx= results$N - results$N_h0
    
  # Notified TB cases   
  results$notifTB=numeric(length=nrow(results))
  for(k in results$time) {
    results$notifTB[which(results$time==k)] =
      ifelse(k==tail(results$time, n=1),NA,results$CumTBnotif[which(results$time==k+1)]-results$CumTBnotif[which(results$time==k)])
  }

  # False positives 
  results$FP=numeric(length=nrow(results))
  for(k in results$time) {
    results$FP[which(results$time==k)] =
      ifelse(k==tail(results$time, n=1),NA,results$CumFP[which(results$time==k+1)]-results$CumFP[which(results$time==k)])
  }
  
  # Percentage false-positives
  results$p_FP = results$FP / results$notifTB
  
  # Inc TB cases   
  results$IncTB=numeric(length=nrow(results))
  for(k in results$time) {
    results$IncTB[which(results$time==k)] =
      ifelse(k==tail(results$time, n=1),NA,results$CumTBinc[which(results$time==k+1)]-results$CumTBinc[which(results$time==k)])
  }
 
  # Inc TB cases_h0   
  results$IncTB_h0=numeric(length=nrow(results))
  for(k in results$time) {
    results$IncTB_h0[which(results$time==k)] =
      ifelse(k==tail(results$time, n=1),NA,results$CumTBinc_h0[which(results$time==k+1)]-results$CumTBinc_h0[which(results$time==k)])
  } 
  
  # Inc TB cases_hx   
  results$IncTB_hx=results$IncTB-results$IncTB_h0
  
  # Percentage Inc TBHIV   
  results$pIncHIV=results$IncTB_hx/(results$IncTB_hx+results$IncTB_h0)
  
  # TB deaths
  results$TBdeaths=numeric(length=nrow(results))
  for(k in results$time) {
    results$TBdeaths[which(results$time==k)] =
      ifelse(k==tail(results$time, n=1),NA,results$CumTBmort[which(results$time==k+1)]-results$CumTBmort[which(results$time==k)])
  }
  
  results$TBinc100T = results$IncTB / results$N * 100000
  results$TBmort100T = results$TBdeaths / results$N * 100000
  results$TBinc_h0_100T = results$IncTB_h0 / results$N_h0 * 100000
  results$TBinc_hx_100T = results$IncTB_hx / results$N_hx * 100000
  results$hivTBinc_100T = results$IncTB_hx / results$N * 100000
  
  # Size of people in diseased compartments (Prevalence)
  results$TB = results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$IP_h1 + results$IS_h1 + results$DP_h1 + results$DS_h1 + results$FN_h1 + results$IP_h2 + results$IS_h2 + results$DP_h2 + results$DS_h2 + results$FN_h2 +
    results$IP_h3 + results$IS_h3 + results$DP_h3 + results$DS_h3 + results$FN_h3 + results$IP_h4 + results$IS_h4 + results$DP_h4 + results$DS_h4 + results$FN_h4 + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a +
    results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a
  
  # Prevalence of TB disease
  results$p_TB = results$TB / results$N
  
  # Proportion prevalence that are pre-symptomatic
  results$p_presymp = (results$IP_h0 + results$DP_h0 + results$IP_h1 + results$DP_h1 + results$IP_h2 + results$DP_h2 + results$IP_h3 + results$DP_h3 + results$IP_h4 + results$DP_h4 + results$IP_h1a + results$DP_h1a + results$IP_h2a + results$DP_h2a +
                            results$IP_h3a + results$DP_h3a + results$ IP_h4a + results$DP_h4a) / results$TB
  
  # Proportion symptomatic that previously sought care
  results$p_prevcare = (results$FN_h0 + results$FN_h1 + results$FN_h2 + results$FN_h3 + results$FN_h4 + results$FN_h1a + results$FN_h2a +  results$FN_h3a + results$FN_h4a) /
      (results$IS_h0 + results$DS_h0 + results$FN_h0 + results$IS_h1 + results$DS_h1 + results$FN_h1 + results$IS_h2 + results$DS_h2 + results$FN_h2 + results$IS_h3 + results$DS_h3 +
       results$FN_h3 + results$IS_h4 + results$DS_h4 + results$FN_h4 + results$IS_h1a + results$DS_h1a + results$FN_h1a + results$IS_h2a + results$DS_h2a + results$FN_h2a +
      results$IS_h3a + results$DS_h3a + results$FN_h3a  + results$IS_h4a + results$DS_h4a + results$FN_h4a)
  
  
  # HIV prevalence and ART coverage (all ages)
  results$HIV = results$N - (results$S_h0 + results$LR_h0 + results$LD_h0 + results$IP_h0 + results$IS_h0 + results$DP_h0 + results$DS_h0 + results$FN_h0 + results$T_h0 + results$R_h0)
                  
  results$ART =     results$S_h1a + results$LR_h1a + results$LD_h1a + results$IP_h1a + results$IS_h1a + results$DP_h1a + results$DS_h1a + results$FN_h1a + results$T_h1a + results$R_h1a +
    results$S_h2a + results$LR_h2a + results$LD_h2a + results$IP_h2a + results$IS_h2a + results$DP_h2a + results$DS_h2a + results$FN_h2a + results$T_h2a + results$R_h2a +
    results$S_h3a + results$LR_h3a + results$LD_h3a + results$IP_h3a + results$IS_h3a + results$DP_h3a + results$DS_h3a + results$FN_h3a + results$T_h3a + results$R_h3a +
    results$S_h4a + results$LR_h4a + results$LD_h4a + results$IP_h4a + results$IS_h4a + results$DP_h4a + results$DS_h4a + results$FN_h4a + results$T_h4a + results$R_h4a

  results$p_HIV = results$HIV / results$N
  results$p_ART = results$ART / results$HIV
  
  # Mid-year population
  results$Nmid=numeric(length=nrow(results))
  for(k in results$time){
    results$Nmid[which(results$time==k)] = 
      ifelse(k==tail(results$time, n=1),NA, (results$N[which(results$time==k)]+results$N[which(results$time==k+1)])/2) 
  }
  

  # Drop unnecessary columns
  results <- results[,-c(2:105)] 
  return(results)
   
}



# setup parallel backend to use many processors
cores <- detectCores()
cl <- makeCluster(cores - 2)  # Reduce cores to avoid overloading computer
registerDoParallel(cl)

finalMatrix <- foreach(i=1:nrowparammat) %dopar% {        #
  require(deSolve)
  tempMatrix = funcMakeResults(model_version = TB_model, timesvect = timesvect_calibration, yinit_list = yinit, param_list = ParamList2) 
  tempMatrix 
}

stopCluster(cl)

# Save results
save(funcMakeResults, file = "OUTPUTS/RSA/funcMakeResults.RData")
saveRDS(finalMatrix, file=tail(.args, 1))







