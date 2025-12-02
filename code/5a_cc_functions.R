################################################################################################################################
############################### Authors: SACEMA-DTTC Modelling group ###########################################################
################################################################################################################################
################################### SA TB Model Version 01  ####################################################################
################################################################################################################################
###################################### Define the start population 1995   #######################################################
################################################################################################################################


## This script .... 

## INPUTS: 
## OUTPUTS: 


############################################################
###


.args <- if(interactive()){
	c("OUTPUTS/RSA/cc_functions.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}


############################################################
### Function for probability true-positive diagnoses

fun_diagstp_notutt = function(relfoot_phc,pspneg,ltfu_prediag,sputsc,sens_diag,clindiag, ltfu_diag, iltfu_phc, iltfu_other) {

  pdiagstp_notutt =
    relfoot_phc*(1-ltfu_prediag)*(1-pspneg)*(1-sputsc)*sens_diag*(1-ltfu_diag)*(1-iltfu_phc) +  # bacteriologically-confirmed via sputum in PHC
    relfoot_phc*(1-ltfu_prediag)*(1-pspneg)*(1-sputsc)*(1-sens_diag)*clindiag*(1-iltfu_phc) + #bact- but clinically diagnosed in PHC
    relfoot_phc*(1-ltfu_prediag)*(1-pspneg)*(sputsc)*clindiag*(1-iltfu_phc) + # unable to produce sputum, clinically diagnosed in PHC
    relfoot_phc*(1-ltfu_prediag)*pspneg * clindiag *(1-ltfu_diag)*(1-iltfu_phc)+  # sputum negative, other (clinical) diagnosis in PHC
    (1-relfoot_phc)*(1-ltfu_prediag)*(1-pspneg)*(1-sputsc)*sens_diag*(1-ltfu_diag)*(1-iltfu_other) +  # bacteriologically-confirmed via sputum in Other
    (1-relfoot_phc)*(1-ltfu_prediag)*(1-pspneg)*(1-sputsc)*(1-sens_diag)*clindiag*(1-iltfu_other) + #bact- but clinically diagnosed in Other
    (1-relfoot_phc)*(1-ltfu_prediag)*(1-pspneg)*(sputsc)*clindiag*(1-iltfu_other) + # unable to produce sputum, clinically diagnosed in Other
    (1-relfoot_phc)*(1-ltfu_prediag)*pspneg * clindiag *(1-iltfu_other)  # sputum negative, other (clinical) diagnosis in Other

  return(pdiagstp_notutt)
}


############################################################
### Function for for probability false-positive diagnoses

fun_diagsfp_notutt = function(relfoot_phc,ltfu_prediag,ptb,sputsc,spec_diag,pfdiag,ltfu_diag,iltfu_phc, iltfu_other) {

  pdiagsfp_notutt =
    relfoot_phc*(1-ltfu_prediag)*((1-ptb)/ptb)*(1-sputsc)*(1-spec_diag)*(1-ltfu_diag)*(1-iltfu_phc) +  # bacteriologically-confirmed via sputum in PHC
    relfoot_phc*(1-ltfu_prediag)*((1-ptb)/ptb)*(1-sputsc)*spec_diag*pfdiag*(1-iltfu_phc) +
    relfoot_phc*(1-ltfu_prediag)*((1-ptb)/ptb)*sputsc*pfdiag*(1-iltfu_phc) +
    (1-relfoot_phc)*(1-ltfu_prediag)*((1-ptb)/ptb)*(1-sputsc)*(1-spec_diag)*(1-ltfu_diag)*(1-iltfu_other) +  # bacteriologically-confirmed via sputum in Other
    (1-relfoot_phc)*(1-ltfu_prediag)*((1-ptb)/ptb)*(1-sputsc)*spec_diag*pfdiag*(1-iltfu_other) +
    (1-relfoot_phc)*(1-ltfu_prediag)*((1-ptb)/ptb)*sputsc*pfdiag*(1-iltfu_other)
  return(pdiagsfp_notutt)
}

############################################################
save(fun_diagstp_notutt, fun_diagsfp_notutt, file = tail(.args, 1))