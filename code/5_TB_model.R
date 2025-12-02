################################################################################
# SA TB Model Version 01 - Transmission dynamic model equations
# Authors: SACEMA-DTTC Modelling Group
################################################################################

# This script defines the TB model with various compartments and transmission dynamics.

################################################################################

# Define input and output file paths
.args <- if(interactive()){
	c("OUTPUTS/RSA/TB_model.RData")
}else{
	commandArgs(trailingOnly = TRUE)
}


# This function takes the following arguments:
# times:   sequence of times in years that the model runs for
# yinit:   starting values for populations at the first time-point
# pars:    the particular parameter combination for this model run

TB_model <- function(times,yinit,pars){
  with(as.list(c(yinit,pars)), {

    
  #  Calculate Total population
  N <- S_h0 + LR_h0 + LD_h0 + IP_h0 + IS_h0 + DP_h0 + DS_h0 + FN_h0 + T_h0 + RH_h0 + R_h0 +
       S_h1 + LR_h1 + LD_h1 + IP_h1 + IS_h1 + DP_h1 + DS_h1 + FN_h1 + T_h1 + RH_h1 + R_h1 +
       S_h2 + LR_h2 + LD_h2 + IP_h2 + IS_h2 + DP_h2 + DS_h2 + FN_h2 + T_h2 + RH_h2 + R_h2 +
       S_h3 + LR_h3 + LD_h3 + IP_h3 + IS_h3 + DP_h3 + DS_h3 + FN_h3 + T_h3 + RH_h3 + R_h3 +
       S_h4 + LR_h4 + LD_h4 + IP_h4 + IS_h4 + DP_h4 + DS_h4 + FN_h4 + T_h4 + RH_h4 + R_h4 +
       S_h1a + LR_h1a + LD_h1a + IP_h1a + IS_h1a + DP_h1a + DS_h1a + FN_h1a + T_h1a + RH_h1a + R_h1a +
       S_h2a + LR_h2a + LD_h2a + IP_h2a + IS_h2a + DP_h2a + DS_h2a + FN_h2a + T_h2a + RH_h2a + R_h2a +
       S_h3a + LR_h3a + LD_h3a + IP_h3a + IS_h3a + DP_h3a + DS_h3a + FN_h3a + T_h3a + RH_h3a + R_h3a +
       S_h4a + LR_h4a + LD_h4a + IP_h4a + IS_h4a + DP_h4a + DS_h4a + FN_h4a + T_h4a + RH_h4a + R_h4a

      
    #  Define transmission coefficient
    TC = beta*IS_h0/N + 
    	beta*IS_h1/N*relbeta_h1 +
    	beta*IS_h2/N*relbeta_h2 +
    	beta*IS_h3/N*relbeta_h3 +
      beta*IS_h4/N*relbeta_h4 +
      beta*IS_h1a/N*relbeta_h1a +
      beta*IS_h2a/N*relbeta_h2a +
      beta*IS_h3a/N*relbeta_h3a +
      beta*IS_h4a/N*relbeta_h4a +
      beta*IP_h0/N*relbeta_ps +
      beta*IP_h1/N*relbeta_h1*relbeta_ps +
      beta*IP_h2/N*relbeta_h2*relbeta_ps +
      beta*IP_h3/N*relbeta_h3*relbeta_ps +
      beta*IP_h4/N*relbeta_h4*relbeta_ps +
      beta*IP_h1a/N*relbeta_h1a*relbeta_ps +
      beta*IP_h2a/N*relbeta_h2a*relbeta_ps +
      beta*IP_h3a/N*relbeta_h3a*relbeta_ps +
      beta*IP_h4a/N*relbeta_h4a*relbeta_ps +
      beta*DS_h0/N +
      beta*DS_h1/N*relbeta_h1 +
      beta*DS_h2/N*relbeta_h2 +
      beta*DS_h3/N*relbeta_h3 +
      beta*DS_h4/N*relbeta_h4 +
      beta*DS_h1a/N*relbeta_h1a +
    	beta*DS_h2a/N*relbeta_h2a +
      beta*DS_h3a/N*relbeta_h3a +
      beta*DS_h4a/N*relbeta_h4a +
    	beta*DP_h0/N*relbeta_ps +
      beta*DP_h1/N*relbeta_h1*relbeta_ps +
      beta*DP_h2/N*relbeta_h2*relbeta_ps +
      beta*DP_h3/N*relbeta_h3*relbeta_ps +
      beta*DP_h4/N*relbeta_h4*relbeta_ps +
      beta*DP_h1a/N*relbeta_h1a*relbeta_ps +
      beta*DP_h2a/N*relbeta_h2a*relbeta_ps +
      beta*DP_h3a/N*relbeta_h3a*relbeta_ps +
      beta*DP_h4a/N*relbeta_h4a*relbeta_ps +
      beta*FN_h0/N +
      beta*FN_h1/N*relbeta_h1 +
      beta*FN_h2/N*relbeta_h2 +
      beta*FN_h3/N*relbeta_h3 +
      beta*FN_h4/N*relbeta_h4 +
      beta*FN_h1a/N*relbeta_h1a +
      beta*FN_h2a/N*relbeta_h2a +
      beta*FN_h3a/N*relbeta_h3a +
      beta*FN_h4a/N*relbeta_h4a
     

     #  define the time-varying paramaters 
     artinit_h4 <- artgen_general(times, artinit_L_h4, year_art_init_h1, artinit_k, artinit_mr)  # 2008
     artinit_h3 <- artgen_general(times, artinit_L_h3, year_art_init_h2, artinit_k, artinit_mr)  # 2012
     artinit_h2 <- artgen_general(times, artinit_L_h2, year_art_init_h3, artinit_k, artinit_mr)  # 2014
     artinit_h1 <- artgen_general(times, artinit_L_h1, year_art_init_h4, artinit_k, artinit_mr)  # 2018
     
     artinit_tb <- artinit_tb_func(times)

     # hiv incidence as a time-varying parameter
     hivinf <- hivinc_func(times)
     bhivp <- prob_hiv_func(times)
     

     # Rates of leaving diagnostic services as true positive 
     
       roundtime = floor(times)
       if (roundtime<year_xpert_ultra_init) {
         
         p_diagstp_h0 = fun_diagstp_notutt(relfoot_phc= relfoot_phc,
         																		pspneg = 0,
                                            ltfu_prediag = ltfu_prediag_h0,
                                            sputsc = sputsc_symp_h0,
                                            sens_diag = sens_xpert_h0,
                                            clindiag = clindiag_h0,
                                            ltfu_diag = ltfu_diag_ultra_h0,
                                            iltfu_phc = iltfu_phc,
                                            iltfu_other = iltfu_other)
         
         p_diagstp_hx = fun_diagstp_notutt(relfoot_phc=1,
         																	 pspneg = 0,
                                           ltfu_prediag = ltfu_prediag_hx,
                                           sputsc = sputsc_symp_hx,
                                           sens_diag = sens_xpert_hx,
                                           clindiag = clindiag_hx,
                                           ltfu_diag = ltfu_diag_ultra_hx,
                                           iltfu_phc = iltfu_phc,
                                           iltfu_other = iltfu_other)
          
       } else if (roundtime>=year_xpert_ultra_init) {
         
         p_diagstp_h0 = fun_diagstp_notutt(relfoot_phc= relfoot_phc,  
         																	 pspneg = 0,
                                           ltfu_prediag = ltfu_prediag_h0,
                                           sputsc = sputsc_symp_h0,
                                           sens_diag = sens_ultra_h0,  ### Ultra instead of Xpert
                                           clindiag = clindiag_h0,
                                           ltfu_diag = ltfu_diag_ultra_h0,
                                           iltfu_phc = iltfu_phc,
                                           iltfu_other = iltfu_other)
         
         p_diagstp_hx = fun_diagstp_notutt(relfoot_phc=1,
         																	 pspneg = 0,
                                           ltfu_prediag = ltfu_prediag_hx,
                                           sputsc = sputsc_symp_hx,
                                           sens_diag = sens_ultra_hx, ### Ultra instead of Xpert
                                           clindiag = clindiag_hx,
                                           ltfu_diag = ltfu_diag_ultra_hx,
                                           iltfu_phc = iltfu_phc,
                                           iltfu_other = iltfu_other)
        }
         
         diagstp_h0 = (365/diagdur) * p_diagstp_h0
         diagstp_h1 = (365/diagdur) * p_diagstp_hx
         diagstp_h2 = (365/diagdur) * p_diagstp_hx
         diagstp_h3 = (365/diagdur) * p_diagstp_hx
         diagstp_h4 = (365/diagdur) * p_diagstp_hx
         diagstp_h1a = (365/diagdur) * p_diagstp_hx
         diagstp_h2a = (365/diagdur) * p_diagstp_hx
         diagstp_h3a = (365/diagdur) * p_diagstp_hx
         diagstp_h4a = (365/diagdur) * p_diagstp_hx
         
         diagsfn_h0 = (365/diagdur) * (1-p_diagstp_h0)
         diagsfn_h1 = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h2 = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h3 = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h4 = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h1a = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h2a = (365/diagdur) * (1-p_diagstp_hx)
         diagsfn_h3a = (365/diagdur) *(1-p_diagstp_hx)
         diagsfn_h4a = (365/diagdur) * (1-p_diagstp_hx)

         diagptp_h0 = 0
         diagptp_h1 = 0
         diagptp_h2 = 0
         diagptp_h3 = 0
         diagptp_h4 = 0
         diagptp_h1a = 0
         diagptp_h2a = 0
         diagptp_h3a = 0
         diagptp_h4a = 0
         
         diagpfn_h0 = 0
         diagpfn_h1 = 0
         diagpfn_h2 = 0
         diagpfn_h3 = 0
         diagpfn_h4 = 0
         diagpfn_h1a = 0
         diagpfn_h2a = 0
         diagpfn_h3a = 0
         diagpfn_h4a = 0
      
         
        #  False-positive rates
         
         roundtime = floor(times)
         if (roundtime<year_xpert_ultra_init) {
           
           p_diagsfp_h0 = fun_diagsfp_notutt(relfoot_phc=relfoot_phc,
           																	ltfu_prediag = ltfu_prediag_h0,
           																	ptb = ptb_h0,
           																	sputsc = sputsc_symp_h0,
           																	spec_diag = spec_xpert_h0,
           																	pfdiag = pfdiag_h0,
           																	ltfu_diag = ltfu_diag_ultra_h0,
           																	iltfu_phc = iltfu_phc,
           																	iltfu_other = iltfu_other)

           p_diagsfp_hx = fun_diagsfp_notutt(relfoot_phc=1,
                                             ltfu_prediag = ltfu_prediag_hx,
           																	 ptb = ptb_hx,
                                             sputsc = sputsc_symp_hx,
                                             spec_diag = spec_xpert_hx,
           																	 pfdiag = pfdiag_hx,
                                             ltfu_diag = ltfu_diag_ultra_hx,
                                             iltfu_phc = iltfu_phc,
                                             iltfu_other = iltfu_other)

         } else if (roundtime>=year_xpert_ultra_init) {

           p_diagsfp_h0 = fun_diagsfp_notutt(relfoot_phc=relfoot_phc,   
                                             ltfu_prediag = ltfu_prediag_h0,
           																	 ptb = ptb_h0,
                                             sputsc = sputsc_symp_h0,
                                             spec_diag = spec_ultra_h0,  ### Ultra instead of Xpert
           																	 pfdiag = pfdiag_h0,
                                             ltfu_diag = ltfu_diag_ultra_h0,
                                             iltfu_phc = iltfu_phc,
                                             iltfu_other = iltfu_other)

           p_diagsfp_hx = fun_diagsfp_notutt(relfoot_phc=1,
                                             ltfu_prediag = ltfu_prediag_hx,
           																	 ptb = ptb_hx,
                                             sputsc = sputsc_symp_hx,
                                             spec_diag = spec_ultra_hx, ### Ultra instead of Xpert
           																	 pfdiag = pfdiag_hx,
                                             ltfu_diag = ltfu_diag_ultra_hx,
                                             iltfu_phc = iltfu_phc,
                                             iltfu_other = iltfu_other)
         }

         
         diagsfp_h0 = (365/diagdur) * ((1-ptb_h0)/ptb_h0) * p_diagsfp_h0
         diagsfp_h1 = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h2 = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h3 = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h4 = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h1a = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h2a = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h3a = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         diagsfp_h4a = (365/diagdur) * ((1-ptb_hx)/ptb_hx) * p_diagsfp_hx
         
         diagpfp_h0 = 0
         diagpfp_h1 = 0
         diagpfp_h2 = 0
         diagpfp_h3 = 0
         diagpfp_h4 = 0
         diagpfp_h1a = 0
         diagpfp_h2a = 0
         diagpfp_h3a = 0
         diagpfp_h4a = 0
         
            
     ###########################################################################################
     ########### COMPARTMENTS PART ############
     #
     #
     ###  h0 (HIV uninfected)        
     ###  Calculate differentials for the h0 compartments: S LR LD IP IS DP DS FN T RH R
     
     dS_h0  <-  brate*N*(1-bhivp) - ## birth rate
     	TC*S_h0 -                    ## primary infection
     	natmort*S_h0 -               ## natural mortality
     	hivinf*S_h0                  ## HIV infection
     
     dLR_h0 <- TC*S_h0 +            ## primary infection
     	TC*relreinflat_h0*LD_h0 +    ## reinfection (distant latent infection)
     	TC*relreinfact_h0*RH_h0 +      ## reinfection (recovered from active TB)
     	TC*relreinfact_h0*R_h0 -     ## reinfection (recovered from active TB)
     	rapprog_h0*LR_h0 -           ## rapid progression (after recent infection)
     	recdist*LR_h0 -              ## transition from recent to distant infection
     	natmort*LR_h0 -              ## natural mortality
     	hivinf*LR_h0                 ## HIV infection
     
     dLD_h0 <- recdist*LR_h0 -      ## transition from recent to distant infection
     	TC*relreinflat_h0*LD_h0 -    ## reinfection (distant latent infection)
     	delprog_h0*LD_h0 -           ## delayed progression (after distant infection)
     	natmort*LD_h0 -              ## natural mortality
     	hivinf*LD_h0                 ## HIV infection
     
     dIP_h0 <- delprog_h0*LD_h0 +      ## delayed progression (after distant infection)
     	diagptp_h0*iltfu*DP_h0 + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h0*LR_h0 +           ## rapid progression (after recent infection)
     	relapse_h0*R_h0 +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h0*hrelapse_h0*RH_h0 + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h0*DP_h0 -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h0*IP_h0 -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h0*IP_h0 -          ## progression to symptomatic TB
     	accdiagp_h0*IP_h0 -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h0 -              ## natural mortality
     	hivinf*IP_h0                 ## HIV infection
     
     dIS_h0 <- progsymp_h0*IP_h0 -     ## progression to symptomatic TB
     	tbrecov_h0*IS_h0 -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h0*IS_h0 -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h0*IS_h0 -                ## mortality due to untreated TB
     	natmort*IS_h0 -              ## natural mortality
     	hivinf*IS_h0                 ## HIV infection
     
     dDP_h0 <- accdiagp_h0*IP_h0 -    ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h0*DP_h0 -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h0*DP_h0 -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h0 -              ## natural mortality
     	hivinf*DP_h0                 ## HIV infection
     
     dDS_h0 <- accdiags_h0*IS_h0 -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h0*DS_h0 +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h0*FN_h0 -              ## former false-negative people re-entering
     	diagstp_h0*DS_h0 -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h0*DS_h0 -                ## mortality due to untreated TB
     	natmort*DS_h0 -              ## natural mortality
     	hivinf*DS_h0                 ## HIV infection
     
     dFN_h0 <- diagsfn_h0*DS_h0 +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h0*iltfu*DS_h0 - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up
     	fnreltrans*accdiags_h0*FN_h0 -              ## former false-negative people re-entering
     	tbrecov_h0*FN_h0 -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h0*FN_h0 -               ## mortality due to untreated TB
     	natmort*FN_h0 -              ## natural mortality
     	hivinf*FN_h0                 ## HIV infection
     
     dT_h0 <- diagptp_h0*(1-iltfu)*DP_h0 +    ## true-positive diagnosis of pre-symptomatic TB
     	diagstp_h0*(1-iltfu)*DS_h0 - ## true-positive diagnosis of symptomatic TB
     	tfail*T_h0 -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h0 -            ## recovery after TB treatment
     	tbtreatmort_h0*T_h0 -           ## mortality due to TB during treatment
     	natmort*T_h0 -               ## natural mortality
     	hivinf*T_h0                  ## HIV infection
     
     dRH_h0 <- tfail*T_h0 +         ## partial recovery after incomplete treatment
     	tbrecov_h0*FN_h0 +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h0*IP_h0 +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h0*IS_h0 -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h0*RH_h0 -     ## reinfection (recovered from active TB)
     	recovhl_h0*RH_h0 -          ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h0 -              ## natural mortality
     	hivinf*RH_h0                  ## HIV infection
     
     dR_h0 <- recovtreat*T_h0 +        ## recovery after TB treatment
     	recovhl_h0*(1-hrelapse_h0)*RH_h0 -          ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h0*R_h0 -     ## reinfection (revovered from active TB)
     	relapse_h0*R_h0 -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h0 -              ## natural mortality
     	hivinf*R_h0                  ## HIV infection
     
     #  h1 (HIV infected, CD4>500)
     #  Calculate differentials for the h1 compartments: S LR LD IP IS DP DS FN T R
     
     dS_h1  <-  brate*N*(bhivp) -   ## birth rate
     	TC*S_h1 -                    ## primary infection
     	natmort*S_h1 +               ## natural mortality
     	hivinf*S_h0 -                ## HIV infection
     	hivprog_h1h2*S_h1 -          ## Progression from h1 to h2
     	artinit_h1*S_h1              ## ART initiation
     
     dLR_h1 <- TC*S_h1 +            ## primary infection
     	TC*relreinflat_h1*LD_h1 +    ## reinfection (distant latent infection)
     	TC*relreinfact_h0*RH_h0 +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h1*R_h1 -     ## reinfection (revovered from active TB)
     	rapprog_h1*LR_h1 -           ## rapid progression (after recent infection)
     	recdist*LR_h1 -              ## transition from recent to distant infection
     	natmort*LR_h1 +              ## natural mortality
     	hivinf*LR_h0 -               ## HIV infection
     	hivprog_h1h2*LR_h1 -          ## Progression from h1 to h2
     	artinit_h1*LR_h1              ## ART initiation
     
     dLD_h1 <- recdist*LR_h1 -      ## transition from recent to distant infection
     	TC*relreinflat_h1*LD_h1 -    ## reinfection (distant latent infection)
     	delprog_h1*LD_h1 -           ## delayed progression (after distant infection)
     	natmort*LD_h1 +              ## natural mortality
     	hivinf*LD_h0 -               ## HIV infection
     	hivprog_h1h2*LD_h1 -          ## Progression from h1 to h2
     	artinit_h1*LD_h1              ## ART initiation
     
     dIP_h1 <- delprog_h1*LD_h1 +      ## delayed progression (after distant infection)
     	diagptp_h1*iltfu*DP_h1 + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h1*LR_h1 +           ## rapid progression (after recent infection)
     	relapse_h1*R_h1 +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h1*hrelapse_h1*RH_h1 + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h1*DP_h1 -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h1*IP_h1 -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h1*IP_h1 -          ## progression to symptomatic TB
     	accdiagp_h1*IP_h1 -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h1 +              ## natural mortality
     	hivinf*IP_h0  -              ## HIV infection
     	hivprog_h1h2*IP_h1 -          ## Progression from h1 to h2
     	artinit_h1*IP_h1              ## ART initiation
     
     dIS_h1 <- progsymp_h1*IP_h1 -     ## progression to symptomatic TB
     	tbrecov_h1*IS_h1 -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h1*IS_h1 -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h1*IS_h1 -            ## mortality due to untreated TB
     	natmort*IS_h1 +              ## natural mortality
     	hivinf*IS_h0  -              ## HIV infection
     	hivprog_h1h2*IS_h1 -          ## Progression from h1 to h2
     	artinit_h1*IS_h1              ## ART initiation
     
     dDP_h1 <- accdiagp_h1*IP_h1 -    ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h1*DP_h1 -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h1*DP_h1 -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h1 +              ## natural mortality
     	hivinf*DP_h0  -              ## HIV infection
     	hivprog_h1h2*DP_h1 -          ## Progression from h1 to h2
     	artinit_h1*DP_h1              ## ART initiation
     
     dDS_h1 <- accdiags_h1*IS_h1 -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h1*DS_h1 +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h1*FN_h1 -              ## former false-negative people re-entering
     	diagstp_h1*DS_h1 -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h1*DS_h1 -            ## mortality due to untreated TB
     	natmort*DS_h1 +              ## natural mortality
     	hivinf*DS_h0  -              ## HIV infection
     	hivprog_h1h2*DS_h1 -          ## Progression from h1 to h2
     	artinit_h1*DS_h1              ## ART initiation
     
     dFN_h1 <- diagsfn_h1*DS_h1 +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h1*iltfu*DS_h1 - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h1*FN_h1 -              ## former false-negative people re-entering
     	tbrecov_h1*FN_h1 -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h1*FN_h1 -              ## mortality due to untreated TB
     	natmort*FN_h1 +              ## natural mortality
     	hivinf*FN_h0 -              ## HIV infection
     	hivprog_h1h2*FN_h1 -          ## Progression from h1 to h2
     	artinit_h1*FN_h1              ## ART initiation
     
     dT_h1 <- diagptp_h1*(1-iltfu)*(1-artinit_tb)*DP_h1 +       ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	diagstp_h1*(1-iltfu)*(1-artinit_tb)*DS_h1 -           ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	tfail*T_h1 -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h1 -            ## recovery after TB treatment
     	tbtreatmort_h1*T_h1 -        ## mortality due to TB during treatment
     	natmort*T_h1 +               ## natural mortality
     	hivinf*T_h0 -                ## HIV infection
     	hivprog_h1h2*T_h1 -          ## Progression from h1 to h2
     	artinit_h1*T_h1              ## ART initiation
     
     dRH_h1 <- tfail*T_h1 +         ## partial recovery after incomplete treatment
     	tbrecov_h1*FN_h1 +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h1*IP_h1 +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	TC*relreinfact_h1*RH_h1 -     ## reinfection (recovered from active TB)
     	tbrecov_h1*IS_h1 -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	recovhl_h1*RH_h1 -         ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h1 +    ## natural mortality
     	hivinf*RH_h0 -                ## HIV infection
     	hivprog_h1h2*RH_h1 -          ## Progression from h1 to h2
     	artinit_h1*RH_h1              ## ART initiation
     
     dR_h1 <- recovtreat*T_h1 +        ## recovery after TB treatment
     	recovhl_h1*(1-hrelapse_h1)*RH_h1 -  ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h1*R_h1 -     ## reinfection (revovered from active TB)
     	relapse_h1*R_h1 -            ## relapse (endognous reactivation after TB treatment)
     	natmort*R_h1 +    ## natural mortality
     	hivinf*R_h0 -                ## HIV infection
     	hivprog_h1h2*R_h1 -          ## Progression from h1 to h2
     	artinit_h1*R_h1              ## ART initiation
     
     ###  h2 (HIV infected, CD4 350-500)
     ###  Calculate differentials for the h2 compartments: S LR LD IP IS DP DS FN T R
     
     dS_h2  <- TC*S_h2 -            ## primary infection
     	natmort*S_h2 +               ## natural mortality
     	hivprog_h1h2*S_h1 -           ## Progression from h1 to h2
     	hivprog_h2h3*S_h2 -          ## Progression from h2 to h3
     	artinit_h2*S_h2              ## ART initiation
     
     dLR_h2 <- TC*S_h2 +            ## primary infection
     	TC*relreinflat_h2*LD_h2 +    ## reinfection (distant latent infection)
     	TC*relreinfact_h2*RH_h2 +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h2*R_h2 -     ## reinfection (revovered from active TB)
     	rapprog_h2*LR_h2 -           ## rapid progression (after recent infection)
     	recdist*LR_h2 -              ## transition from recent to distant infection
     	natmort*LR_h2 +              ## natural mortality
     	hivprog_h1h2*LR_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*LR_h2 -          ## Progression from h2 to h3
     	artinit_h2*LR_h2              ## ART initiation
     
     dLD_h2 <- recdist*LR_h2 -      ## transition from recent to distant infection
     	TC*relreinflat_h2*LD_h2 -    ## reinfection (distant latent infection)
     	delprog_h2*LD_h2 -           ## delayed progression (after distant infection)
     	natmort*LD_h2 +              ## natural mortality
     	hivprog_h1h2*LD_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*LD_h2 -          ## Progression from h2 to h3
     	artinit_h2*LD_h2              ## ART initiation
     
     dIP_h2 <- delprog_h2*LD_h2 +      ## delayed progression (after distant infection)
     	diagptp_h2*iltfu*DP_h2 + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h2*LR_h2 +           ## rapid progression (after recent infection)
     	relapse_h2*R_h2 +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h2*hrelapse_h2*RH_h2 + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h2*DP_h2 -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h2*IP_h2 -           ## natural recovery from TB (pre-symptomatic to reciovered high-risk)
     	progsymp_h2*IP_h2 -          ## progression to symptomatic TB
     	accdiagp_h2*IP_h2 -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h2 +              ## natural mortality
     	hivprog_h1h2*IP_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*IP_h2 -          ## Progression from h2 to h3
     	artinit_h2*IP_h2              ## ART initiation
     
     dIS_h2 <- progsymp_h2*IP_h2 -     ## progression to symptomatic TB
     	tbrecov_h2*IS_h2 -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h2*IS_h2 -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h2*IS_h2 -            ## mortality due to untreated TB
     	natmort*IS_h2 +              ## natural mortality
     	hivprog_h1h2*IS_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*IS_h2 -          ## Progression from h2 to h3
     	artinit_h2*IS_h2              ## ART initiation
     
     dDP_h2 <- accdiagp_h2*IP_h2 -    ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h2*DP_h2 -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h2*DP_h2 -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h2 +              ## natural mortality
     	hivprog_h1h2*DP_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*DP_h2 -          ## Progression from h2 to h3
     	artinit_h2*DP_h2              ## ART initiation
     
     dDS_h2 <- accdiags_h2*IS_h2 -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h2*DS_h2 +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h2*FN_h2 -              ## former false-negative people re-entering
     	diagstp_h2*DS_h2 -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h2*DS_h2 -            ## mortality due to untreated TB
     	natmort*DS_h2 +              ## natural mortality
     	hivprog_h1h2*DS_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*DS_h2 -          ## Progression from h2 to h3
     	artinit_h2*DS_h2              ## ART initiation
     
     dFN_h2 <- diagsfn_h2*DS_h2 +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h2*iltfu*DS_h2 - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h2*FN_h2 -              ## former false-negative people re-entering
     	tbrecov_h2*FN_h2 -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h2*FN_h2 -               ## mortality due to untreated TB
     	natmort*FN_h2 +              ## natural mortality
     	hivprog_h1h2*FN_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*FN_h2 -          ## Progression from h2 to h3
     	artinit_h2*FN_h2              ## ART initiation
     
     dT_h2 <- diagptp_h2*(1-iltfu)*(1-artinit_tb)*DP_h2 +       ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	diagstp_h2*(1-iltfu)*(1-artinit_tb)*DS_h2 -   ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	tfail*T_h2 -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h2 -            ## recovery after TB treatment
     	tbtreatmort_h2*T_h2 -        ## mortality due to TB during treatment
     	natmort*T_h2 +               ## natural mortality
     	hivprog_h1h2*T_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*T_h2 -          ## Progression from h2 to h3
     	artinit_h2*T_h2              ## ART initiation
     
     dRH_h2 <- tfail*T_h2 +         ## partial recovery after incomplete treatment
     	tbrecov_h2*FN_h2 +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h2*IP_h2 +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h2*IS_h2 -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h2*RH_h2 -     ## reinfection (recovered from active TB)
     	recovhl_h2*RH_h2  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h2 +   ## natural mortality
     	hivprog_h1h2*RH_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*RH_h2 -          ## Progression from h2 to h3
     	artinit_h2*RH_h2              ## ART initiation
     
     dR_h2 <- recovtreat*T_h2  +        ## recovery after TB treatment
     	recovhl_h2*(1-hrelapse_h2)*RH_h2 -  ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h2*R_h2 -     ## reinfection (recovered from active TB)
     	relapse_h2*R_h2 -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h2 +    ## natural mortality
     	hivprog_h1h2*R_h1 -         ## Progression from h1 to h2
     	hivprog_h2h3*R_h2 -          ## Progression from h2 to h3
     	artinit_h2*R_h2              ## ART initiation
     
     ###  h3 (HIV infected, CD4 200-350)
     ###  Calculate differentials for the h3 compartments: S LR LD IP IS DP DS FN T R
     
     dS_h3  <- TC*S_h3 -            ## primary infection
     	natmort*S_h3 +               ## natural mortality
     	hivprog_h2h3*S_h2 -          ## Progression from h2 to h3
     	hivprog_h3h4*S_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*S_h3 - ## Mortality due to HIV at h3
     	artinit_h3*S_h3               ## ART initiation
     
     dLR_h3 <- TC*S_h3 +            ## primary infection
     	TC*relreinflat_h3*LD_h3 +    ## reinfection (distant latent infection)
     	TC*relreinfact_h3*RH_h3 +    ## reinfection (recovered from active TB)
     	TC*relreinfact_h3*R_h3 -     ## reinfection (revovered from active TB)
     	rapprog_h3*LR_h3 -           ## rapid progression (after recent infection)
     	recdist*LR_h3 -              ## transition from recent to distant infection
     	natmort*LR_h3 +              ## natural mortality
     	hivprog_h2h3*LR_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*LR_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*LR_h3 - ## Mortality due to HIV at h3
     	artinit_h3*LR_h3               ## ART initiation
     
     dLD_h3 <- recdist*LR_h3 -      ## transition from recent to distant infection
     	TC*relreinflat_h3*LD_h3 -    ## reinfection (distant latent infection)
     	delprog_h3*LD_h3 -           ## delayed progression (after distant infection)
     	natmort*LD_h3 +              ## natural mortality
     	hivprog_h2h3*LD_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*LD_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*LD_h3 - ## Mortality due to HIV at h3
     	artinit_h3*LD_h3               ## ART initiation
     
     dIP_h3 <- delprog_h3*LD_h3 +      ## delayed progression (after distant infection)
     	diagptp_h3*iltfu*DP_h3 + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h3*LR_h3 +           ## rapid progression (after recent infection)
     	relapse_h3*R_h3 +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h3*hrelapse_h3*RH_h3 + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h3*DP_h3 -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h3*IP_h3 -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h3*IP_h3 -          ## progression to symptomatic TB
     	accdiagp_h3*IP_h3 -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h3 +              ## natural mortality
     	hivprog_h2h3*IP_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*IP_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*IP_h3 - ## Mortality due to HIV at h3
     	artinit_h3*IP_h3               ## ART initiation
     
     dIS_h3 <- progsymp_h3*IP_h3 -     ## progression to symptomatic TB
     	tbrecov_h3*IS_h3 -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h3*IS_h3 -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h3*IS_h3 -            ## mortality due to untreated TB
     	natmort*IS_h3 +              ## natural mortality
     	hivprog_h2h3*IS_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*IS_h3  -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*IS_h3 - ## Mortality due to HIV at h3
     	artinit_h3*IS_h3               ## ART initiation
     
     dDP_h3 <- accdiagp_h3*IP_h3 -     ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h3*DP_h3 -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h3*DP_h3 -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h3 +              ## natural mortality
     	hivprog_h2h3*DP_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*DP_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*DP_h3 - ## Mortality due to HIV at h3
     	artinit_h3*DP_h3               ## ART initiation
     
     dDS_h3 <- accdiags_h3*IS_h3 -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h3*DS_h3 +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h3*FN_h3 -              ## former false-negative people re-entering
     	diagstp_h3*DS_h3 -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h3*DS_h3 -            ## mortality due to untreated TB
     	natmort*DS_h3 +              ## natural mortality
     	hivprog_h2h3*DS_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*DS_h3  -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*DS_h3 - ## Mortality due to HIV at h3
     	artinit_h3*DS_h3               ## ART initiation
     
     dFN_h3 <- diagsfn_h3*DS_h3 +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h3*iltfu*DS_h3 - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h3*FN_h3 -              ## former false-negative people re-entering
     	tbrecov_h3*FN_h3 -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h3*FN_h3 -            ## mortality due to untreated TB
     	natmort*FN_h3 +              ## natural mortality
     	hivprog_h2h3*FN_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*FN_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*FN_h3 - ## Mortality due to HIV at h3
     	artinit_h3*FN_h3               ## ART initiation
     
     dT_h3 <- diagptp_h3*(1-iltfu)*(1-artinit_tb)*DP_h3 +       ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	diagstp_h3*(1-iltfu)*(1-artinit_tb)*DS_h3 -    ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	tfail*T_h3 -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h3 -            ## recovery after TB treatment
     	tbtreatmort_h3*T_h3 -        ## mortality due to TB during treatment
     	natmort*T_h3 +               ## natural mortality
     	hivprog_h2h3*T_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*T_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*T_h3 - ## Mortality due to HIV at h3
     	artinit_h3*T_h3               ## ART initiation
     
     dRH_h3 <- tfail*T_h3 +         ## partial recovery after incomplete treatment
     	tbrecov_h3*FN_h3 +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h3*IP_h3 +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h3*IS_h3 -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h3*RH_h3 -     ## reinfection (recovered from active TB)
     	recovhl_h3*RH_h3  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h3 +    ## natural mortality
     	hivprog_h2h3*RH_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*RH_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*RH_h3 - ## Mortality due to HIV at h3
     	artinit_h3*RH_h3               ## ART initiation
     
     dR_h3 <- recovtreat*T_h3   +        ## recovery after TB treatment
     	recovhl_h3*(1-hrelapse_h3)*RH_h3 - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h3*R_h3 -     ## reinfection (revovered from active TB)
     	relapse_h3*R_h3 -            ## relapse (endognous reactivation after TB treatment)
     	natmort*R_h3 +    ## natural mortality
     	hivprog_h2h3*R_h2 -         ## Progression from h2 to h3
     	hivprog_h3h4*R_h3 -           ## Progression from h3 to h4
     	hivmort_h4*hivmortrat_h3*R_h3 - ## Mortality due to HIV at h3
     	artinit_h3*R_h3               ## ART initiation
     
     ###  h4 (HIV infected, CD4 <200)
     ###  Calculate differentials for the h4 compartments: S LR LD IP IS DP DS FN T R
     
     dS_h4  <- TC*S_h4 -            ## primary infection
     	natmort*S_h4 +               ## natural mortality
     	hivprog_h3h4*S_h3 -          ## Progression from h3 to h4
     	hivmort_h4*S_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*S_h4               ## ART initiation
     
     dLR_h4 <- TC*S_h4 +            ## primary infection
     	TC*relreinflat_h4*LD_h4 +    ## reinfection (distant latent infection)
     	TC*relreinfact_h4*RH_h4 +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h4*R_h4 -     ## reinfection (revovered from active TB)
     	rapprog_h4*LR_h4 -           ## rapid progression (after recent infection)
     	recdist*LR_h4 -              ## transition from recent to distant infection
     	natmort*LR_h4 +              ## natural mortality
     	hivprog_h3h4*LR_h3 -         ## Progression from h3 to h4
     	hivmort_h4*LR_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*LR_h4               ## ART initiation
     
     dLD_h4 <- recdist*LR_h4 -      ## transition from recent to distant infection
     	TC*relreinflat_h4*LD_h4 -    ## reinfection (distant latent infection)
     	delprog_h4*LD_h4 -           ## delayed progression (after distant infection)
     	natmort*LD_h4 +              ## natural mortality
     	hivprog_h3h4*LD_h3 -         ## Progression from h3 to h4
     	hivmort_h4*LD_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*LD_h4               ## ART initiation
     
     dIP_h4 <- delprog_h4*LD_h4 +      ## delayed progression (after distant infection)
     	diagptp_h4*iltfu*DP_h4 + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h4*LR_h4 +           ## rapid progression (after recent infection)
     	relapse_h4*R_h4 +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h4*hrelapse_h4*RH_h4 + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h4*DP_h4 -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h4*IP_h4 -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h4*IP_h4 -          ## progression to symptomatic TB
     	accdiagp_h4*IP_h4 -          ## accessing diagnosis for pre-symptomatic individuals
     	natmort*IP_h4 +              ## natural mortality
     	hivprog_h3h4*IP_h3 -         ## Progression from h3 to h4
     	hivmort_h4*IP_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*IP_h4               ## ART initiation
     
     dIS_h4 <- progsymp_h4*IP_h4 -     ## progression to symptomatic TB
     	tbrecov_h4*IS_h4 -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h4*IS_h4 -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h4*IS_h4 -            ## mortality due to untreated TB
     	natmort*IS_h4 +              ## natural mortality
     	hivprog_h3h4*IS_h3 -         ## Progression from h3 to h4
     	hivmort_h4*IS_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*IS_h4               ## ART initiation
     
     dDP_h4 <- accdiagp_h4*IP_h4 -     ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h4*DP_h4 -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h4*DP_h4 -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h4 +              ## natural mortality
     	hivprog_h3h4*DP_h3 -         ## Progression from h3 to h4
     	hivmort_h4*DP_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*DP_h4               ## ART initiation
     
     dDS_h4 <- accdiags_h4*IS_h4 -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h4*DS_h4 +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h4*FN_h4 -              ## former false-negative people re-entering
     	diagstp_h4*DS_h4 -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h4*DS_h4 -            ## mortality due to untreated TB
     	natmort*DS_h4 +              ## natural mortality
     	hivprog_h3h4*DS_h3 -         ## Progression from h3 to h4
     	hivmort_h4*DS_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*DS_h4               ## ART initiation
     
     dFN_h4 <- diagsfn_h4*DS_h4 +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h4*iltfu*DS_h4 - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h4*FN_h4 -              ## former false-negative people re-entering
     	tbrecov_h4*FN_h4 -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h4*FN_h4 -            ## mortality due to untreated TB
     	natmort*FN_h4 +              ## natural mortality
     	hivprog_h3h4*FN_h3 -         ## Progression from h3 to h4
     	hivmort_h4*FN_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*FN_h4               ## ART initiation
     
     dT_h4 <- diagptp_h4*(1-iltfu)*(1-artinit_tb)*DP_h4 +       ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	diagstp_h4*(1-iltfu)*(1-artinit_tb)*DS_h4 -   ## true-positive diagnosis of pre-symptomatic TB (less those who initate ART)
     	tfail*T_h4 -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h4 -            ## recovery after TB treatment
     	tbtreatmort_h4*T_h4 -        ## mortality due to TB during treatment
     	natmort*T_h4 +               ## natural mortality
     	hivprog_h3h4*T_h3 -          ## Progression from h3 to h4
     	hivmort_h4*T_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*T_h4               ## ART initiation
     
     dRH_h4 <- tfail*T_h4 +         ## partial recovery after incomplete treatment
     	tbrecov_h4*FN_h4 +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h4*IP_h4 +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h4*IS_h4 -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h4*RH_h4 -     ## reinfection (recovered from active TB)
     	recovhl_h4*RH_h4  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h4 +    ## natural mortality
     	hivprog_h3h4*RH_h3 -          ## Progression from h3 to h4
     	hivmort_h4*RH_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*RH_h4               ## ART initiation
     
     dR_h4 <- recovtreat*T_h4   +        ## recovery after TB treatment
     	recovhl_h4*(1-hrelapse_h4)*RH_h4 - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h4*R_h4 -     ## reinfection (recovered from active TB)
     	relapse_h4*R_h4 -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h4 +    ## natural mortality
     	hivprog_h3h4*R_h3 -          ## Progression from h3 to h4
     	hivmort_h4*R_h4 -             ## Mortality due to HIV at h4
     	artinit_h4*R_h4               ## ART initiation
     
     ###  h1a (HIV infected on ART, CD4>500)
     ###  Calculate differentials for the h1a compartments: S LR LD IP IS DP DS FN T R
     
     dS_h1a  <- - TC*S_h1a -        ## primary infection
     	natmort*S_h1a +               ## natural mortality
     	artinit_h1*S_h1             ## ART initiation
     
     dLR_h1a <- TC*S_h1a +            ## primary infection
     	TC*relreinflat_h1a*LD_h1a +    ## reinfection (distant latent infection)
     	TC*relreinfact_h1a*RH_h1a +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h1a*R_h1a -     ## reinfection (revovered from active TB)
     	rapprog_h1a*LR_h1a -           ## rapid progression (after recent infection)
     	recdist*LR_h1a -              ## transition from recent to distant infection
     	natmort*LR_h1a +              ## natural mortality
     	artinit_h1*LR_h1             ## ART initiation
     
     dLD_h1a <- recdist*LR_h1a -      ## transition from recent to distant infection
     	TC*relreinflat_h1a*LD_h1a -    ## reinfection (distant latent infection)
     	delprog_h1a*LD_h1a -           ## delayed progression (after distant infection)
     	natmort*LD_h1a +              ## natural mortality
     	artinit_h1*LD_h1             ## ART initiation
     
     dIP_h1a <- delprog_h1a*LD_h1a +      ## delayed progression (after distant infection)
     	diagptp_h1a*iltfu*DP_h1a + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h1a*LR_h1a +           ## rapid progression (after recent infection)
     	relapse_h1a*R_h1a +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h1a*hrelapse_h1a*RH_h1a + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h1a*DP_h1a -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h1a*IP_h1a -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h1a*IP_h1a -          ## progression to symptomatic TB
     	accdiagp_h1a*IP_h1a -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h1a +              ## natural mortality
     	artinit_h1*IP_h1             ## ART initiation
     
     dIS_h1a <- progsymp_h1a*IP_h1a -     ## progression to symptomatic TB
     	tbrecov_h1a*IS_h1a -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h1a*IS_h1a -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h1a*IS_h1a -            ## mortality due to untreated TB
     	natmort*IS_h1a +              ## natural mortality
     	artinit_h1*IS_h1             ## ART initiation
     
     dDP_h1a <- accdiagp_h1a*IP_h1a -    ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h1a*DP_h1a -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h1a*DP_h1a -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h1a +              ## natural mortality
     	artinit_h1*DP_h1             ## ART initiation
     
     dDS_h1a <- accdiags_h1a*IS_h1a -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h1a*DS_h1a +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h1a*FN_h1a -              ## former false-negative people re-entering
     	diagstp_h1a*DS_h1a -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h1a*DS_h1a -            ## mortality due to untreated TB
     	natmort*DS_h1a +              ## natural mortality
     	artinit_h1*DS_h1             ## ART initiation
     
     dFN_h1a <- diagsfn_h1a*DS_h1a +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h1a*iltfu*DS_h1a - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h1a*FN_h1a -              ## former false-negative people re-entering
     	tbrecov_h1a*FN_h1a -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h1a*FN_h1a -               ## mortality due to untreated TB
     	natmort*FN_h1a +              ## natural mortality
     	artinit_h1*FN_h1             ## ART initiation
     
     dT_h1a <-
     	diagptp_h1*(1-iltfu)*(artinit_tb)*DP_h1 + ## true-positive diagnosis of pre-symptomatic TB and ART initiation
     	diagstp_h1*(1-iltfu)*(artinit_tb)*DS_h1 + ## true-positive diagnosis of symptomatic TB and ART initiation
     	diagptp_h1a*DP_h1a +       ## true-positive diagnosis of pre-symptomatic TB
     	diagstp_h1a*(1-iltfu)*DS_h1a -           ## true-positive diagnosis of symptomatic TB
     	tfail*T_h1a -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h1a -            ## recovery after TB treatment
     	tbtreatmort_h1a*T_h1a -        ## mortality due to TB during treatment
     	natmort*T_h1a +               ## natural mortality
     	artinit_h1*T_h1             ## ART initiation
     
     dRH_h1a <- tfail*T_h1a +         ## partial recovery after incomplete treatment
     	tbrecov_h1a*FN_h1a +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h1a*IP_h1a +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h1a*IS_h1a -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h1a*RH_h1a -     ## reinfection (recovered from active TB)
     	recovhl_h1a*RH_h1a  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h1a +    ## natural mortality
     	artinit_h1*RH_h1             ## ART initiation
     
     dR_h1a <- recovtreat*T_h1a +        ## recovery after TB treatment
     	recovhl_h1a*(1-hrelapse_h1a)*RH_h1a - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h1a*R_h1a -     ## reinfection (recovered from active TB)
     	relapse_h1a*R_h1a -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h1a +    ## natural mortality
     	artinit_h1*R_h1             ## ART initiation
     
     ###  h2a (HIV infected on ART, CD4 350-500)
     ###  Calculate differentials for the h2a compartments: S LR LD IP IS DP DS FN T R
     
     dS_h2a  <- TC*S_h2a -            ## primary infection
     	natmort*S_h2a +               ## natural mortality
     	artinit_h2*S_h2             ## ART initiation
     
     dLR_h2a <- TC*S_h2a +            ## primary infection
     	TC*relreinflat_h2a*LD_h2a +    ## reinfection (distant latent infection)
     	TC*relreinfact_h2a*RH_h2a +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h2a*R_h2a -     ## reinfection (revovered from active TB)
     	rapprog_h2a*LR_h2a -           ## rapid progression (after recent infection)
     	recdist*LR_h2a -              ## transition from recent to distant infection
     	natmort*LR_h2a +              ## natural mortality
     	artinit_h2*LR_h2             ## ART initiation
     
     dLD_h2a <- recdist*LR_h2a -      ## transition from recent to distant infection
     	TC*relreinflat_h2a*LD_h2a -    ## reinfection (distant latent infection)
     	delprog_h2a*LD_h2a -           ## delayed progression (after distant infection)
     	natmort*LD_h2a +              ## natural mortality
     	artinit_h2*LD_h2             ## ART initiation
     
     dIP_h2a <- delprog_h2a*LD_h2a +      ## delayed progression (after distant infection)
     	diagptp_h2a*iltfu*DP_h2a + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h2a*LR_h2a +           ## rapid progression (after recent infection)
     	relapse_h2a*R_h2a +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h2a*hrelapse_h2a*RH_h2a + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h2a*DP_h2a -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h2a*IP_h2a -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h2a*IP_h2a -          ## progression to symptomatic TB
     	accdiagp_h2a*IP_h2a -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h2a +              ## natural mortality
     	artinit_h2*IP_h2             ## ART initiation
     
     dIS_h2a <- progsymp_h2a*IP_h2a -     ## progression to symptomatic TB
     	tbrecov_h2a*IS_h2a -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h2a*IS_h2a -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h2a*IS_h2a -            ## mortality due to untreated TB
     	natmort*IS_h2a +              ## natural mortality
     	artinit_h2*IS_h2             ## ART initiation
     
     dDP_h2a <- accdiagp_h2a*IP_h2a -    ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h2a*DP_h2a -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h2a*DP_h2a -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h2a +              ## natural mortality
     	artinit_h2*DP_h2             ## ART initiation
     
     dDS_h2a <- accdiags_h2a*IS_h2a -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h2a*DS_h2a +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h2a*FN_h2a -              ## former false-negative people re-entering
     	diagstp_h2a*DS_h2a -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h2a*DS_h2a -            ## mortality due to untreated TB
     	natmort*DS_h2a +              ## natural mortality
     	artinit_h2*DS_h2             ## ART initiation
     
     dFN_h2a <- diagsfn_h2a*DS_h2a +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h2a*iltfu*DS_h2a - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h2a*FN_h2a -              ## former false-negative people re-entering
     	tbrecov_h2a*FN_h2a -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h2a*FN_h2a -               ## mortality due to untreated TB
     	natmort*FN_h2a +              ## natural mortality
     	artinit_h2*FN_h2             ## ART initiation
     
     dT_h2a <- diagptp_h2*(1-iltfu)*(artinit_tb)*DP_h2 + ## true-positive diagnosis of pre-symptomatic TB and ART initiation
     	diagstp_h2*(1-iltfu)*(artinit_tb)*DS_h2 + ## true-positive diagnosis of symptomatic TB and ART initiation
     	diagptp_h2a*DP_h2a +       ## true-positive diagnosis of pre-symptomatic TB
     	diagstp_h2a*(1-iltfu)*DS_h2a -           ## true-positive diagnosis of symptomatic TB
     	tfail*T_h2a -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h2a -            ## recovery after TB treatment
     	tbtreatmort_h2a*T_h2a -        ## mortality due to TB during treatment
     	natmort*T_h2a +               ## natural mortality
     	artinit_h2*T_h2             ## ART initiation
     
     dRH_h2a <- tfail*T_h2a +         ## partial recovery after incomplete treatment
     	tbrecov_h2a*FN_h2a +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h2a*IP_h2a +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h2a*IS_h2a -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h2a*RH_h2a -     ## reinfection (recovered from active TB)
     	recovhl_h2a*RH_h2a  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h2a +    ## natural mortality
     	artinit_h2*R_h2             ## ART initiation
     
     dR_h2a <- recovtreat*T_h2a +        ## recovery after TB treatment
     	recovhl_h2a*(1-hrelapse_h2a)*RH_h2a - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h2a*R_h2a -     ## reinfection (recovered from active TB)
     	relapse_h2a*R_h2a -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h2a +    ## natural mortality
     	artinit_h2*R_h2             ## ART initiation
     
     ###  h3a (HIV infected on ART, CD4 200-350)
     ###  Calculate differentials for the h3a compartments: S LR LD IP IS DP DS FN T R
     
     dS_h3a  <- TC*S_h3a -            ## primary infection
     	natmort*S_h3a +               ## natural mortality
     	artinit_h3*S_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*S_h3a  ## Mortality due to HIV at h3 with ART
     
     dLR_h3a <- TC*S_h3a +            ## primary infection
     	TC*relreinflat_h3a*LD_h3a +    ## reinfection (distant latent infection)
     	TC*relreinfact_h3a*RH_h3a +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h3a*R_h3a -     ## reinfection (revovered from active TB)
     	rapprog_h3a*LR_h3a -           ## rapid progression (after recent infection)
     	recdist*LR_h3a -              ## transition from recent to distant infection
     	natmort*LR_h3a +              ## natural mortality
     	artinit_h3*LR_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*LR_h3a  ## Mortality due to HIV at h3 with ART
     
     dLD_h3a <- recdist*LR_h3a -      ## transition from recent to distant infection
     	TC*relreinflat_h3a*LD_h3a -    ## reinfection (distant latent infection)
     	delprog_h3a*LD_h3a -           ## delayed progression (after distant infection)
     	natmort*LD_h3a +              ## natural mortality
     	artinit_h3*LD_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*LD_h3a  ## Mortality due to HIV at h3 with ART
     
     dIP_h3a <- delprog_h3a*LD_h3a +      ## delayed progression (after distant infection)
     	diagptp_h3a*iltfu*DP_h3a + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h3a*LR_h3a +           ## rapid progression (after recent infection)
     	relapse_h3a*R_h3a +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h3a*hrelapse_h3a*RH_h3a + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h3a*DP_h3a -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h3a*IP_h3a -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h3a*IP_h3a -          ## progression to symptomatic TB
     	accdiagp_h3a*IP_h3a -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h3a +              ## natural mortality
     	artinit_h3*IP_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*IP_h3a  ## Mortality due to HIV at h3 with ART
     
     dIS_h3a <- progsymp_h3a*IP_h3a -     ## progression to symptomatic TB
     	tbrecov_h3a*IS_h3a -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h3a*IS_h3a -          ## accessing dagnosis for symptamatic individuals
     	tbmort_h3a*IS_h3a -            ## mortality due to untreated TB
     	natmort*IS_h3a +              ## natural mortality
     	artinit_h3*IS_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*IS_h3a  ## Mortality due to HIV at h3 with ART
     
     dDP_h3a <- accdiagp_h3a*IP_h3a -     ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h3a*DP_h3a -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h3a*DP_h3a -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h3a +              ## natural mortality
     	artinit_h3*DP_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*DP_h3a  ## Mortality due to HIV at h3 with ART
     
     dDS_h3a <- accdiags_h3a*IS_h3a -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h3a*DS_h3a +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h3a*FN_h3a -              ## former false-negative people re-entering
     	diagstp_h3a*DS_h3a -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h3a*DS_h3a -            ## mortality due to untreated TB
     	natmort*DS_h3a +              ## natural mortality
     	artinit_h3*DS_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*DS_h3a  ## Mortality due to HIV at h3 with ART
     
     dFN_h3a <- diagsfn_h3a*DS_h3a +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h3a*iltfu*DS_h3a - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h3a*FN_h3a -              ## former false-negative people re-entering
     	tbrecov_h3a*FN_h3a -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h3a*FN_h3a -            ## mortality due to untreated TB
     	natmort*FN_h3a +              ## natural mortality
     	artinit_h3*FN_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*FN_h3a  ## Mortality due to HIV at h3 with ART
     
     dT_h3a <- diagptp_h3*(1-iltfu)*(artinit_tb)*DP_h3 + ## true-positive diagnosis of pre-symptomatic TB and ART initiation
     	diagstp_h3*(1-iltfu)*(artinit_tb)*DS_h3 + ## true-positive diagnosis of symptomatic TB and ART initiation
     	diagptp_h3a*DP_h3a +       ## true-positive diagnosis of pre-symptomatic TB
     	diagstp_h3a*(1-iltfu)*DS_h3a -           ## true-positive diagnosis of symptomatic TB
     	tfail*T_h3a -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h3a -            ## recovery after TB treatment
     	tbtreatmort_h3a*T_h3a -        ## mortality due to TB during treatment
     	natmort*T_h3a +               ## natural mortality
     	artinit_h3*T_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*T_h3a  ## Mortality due to HIV at h3 with ART
     
     dRH_h3a <- tfail*T_h3a +         ## partial recovery after incomplete treatment
     	tbrecov_h3a*FN_h3a +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h3a*IP_h3a +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h3a*IS_h3a -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h3a*RH_h3a -     ## reinfection (recovered from active TB)
     	recovhl_h3a*RH_h3a  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h3a +    ## natural mortality
     	artinit_h3*RH_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*RH_h3a  ## Mortality due to HIV at h3 with ART
     
     dR_h3a <- recovtreat*T_h3a +        ## recovery after TB treatment
     	recovhl_h3a*(1-hrelapse_h3a)*RH_h3a - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h3a*R_h3a -     ## reinfection (recovered from active TB)
     	relapse_h3a*R_h3a -            ## relapse (endogenous reactivation after TB treatment)
     	natmort*R_h3a +    ## natural mortality
     	artinit_h3*R_h3 -               ## ART initiation
     	hivmort_h4*hivmortrat_h3*hivmortart_ratio_h3*R_h3a  ## Mortality due to HIV at h3 with ART
     
     ###  h4a (HIV infected on ART, CD4 <200)
     ###  Calculate differentials for the h4a compartments: S LR LD IP IS DP DS FN T R
     
     dS_h4a  <- TC*S_h4a -            ## primary infection
     	natmort*S_h4a +               ## natural mortality
     	artinit_h4*S_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*S_h4a ## Mortality due to HIV at h4a
     
     dLR_h4a <- TC*S_h4a +            ## primary infection
     	TC*relreinflat_h4a*LD_h4a +    ## reinfection (distant latent infection)
     	TC*relreinfact_h4a*RH_h4a +     ## reinfection (recovered from active TB)
     	TC*relreinfact_h4a*R_h4a -     ## reinfection (revovered from active TB)
     	rapprog_h4a*LR_h4a -           ## rapid progression (after recent infection)
     	recdist*LR_h4a -              ## transition from recent to distant infection
     	natmort*LR_h4a +              ## natural mortality
     	artinit_h4*LR_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*LR_h4a ## Mortality due to HIV at h4a
     
     dLD_h4a <- recdist*LR_h4a -      ## transition from recent to distant infection
     	TC*relreinflat_h4a*LD_h4a -    ## reinfection (distant latent infection)
     	delprog_h4a*LD_h4a -           ## delayed progression (after distant infection)
     	natmort*LD_h4a +              ## natural mortality
     	artinit_h4*LD_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*LD_h4a ## Mortality due to HIV at h4a
     
     dIP_h4a <- delprog_h4a*LD_h4a +      ## delayed progression (after distant infection)
     	diagptp_h4a*iltfu*DP_h4a + ## true-positive diagnosis of pre-symptomatic TB but initial loss to follow-up
     	rapprog_h4a*LR_h4a +           ## rapid progression (after recent infection)
     	relapse_h4a*R_h4a +            ## relapse (endogenous reactivation after recovery)
     	recovhl_h4a*hrelapse_h4a*RH_h4a + ## relapse (endogenous reactivation after high-risk recovery)
     	diagpfn_h4a*DP_h4a -           ## false-negative diagnosis of pre-symptomatic TB
     	tbrecov_h4a*IP_h4a -           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	progsymp_h4a*IP_h4a -          ## progression to symptomatic TB
     	accdiagp_h4a*IP_h4a -          ## accessing dagnosis for pre-symptamatic individuals
     	natmort*IP_h4a +              ## natural mortality
     	artinit_h4*IP_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*IP_h4a ## Mortality due to HIV at h4a
     
     dIS_h4a <- progsymp_h4a*IP_h4a -     ## progression to symptomatic TB
     	tbrecov_h4a*IS_h4a -           ## natural recovery from TB (symptomatic to pre-symptomatic)
     	accdiags_h4a*IS_h4a -          ## accessing diagnosis for symptomatic individuals
     	tbmort_h4a*IS_h4a -            ## mortality due to untreated TB
     	natmort*IS_h4a +              ## natural mortality
     	artinit_h4*IS_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*IS_h4a ## Mortality due to HIV at h4a
     
     dDP_h4a <- accdiagp_h4a*IP_h4a -     ## accessing dagnosis for pre-symptamatic individuals
     	diagpfn_h4a*DP_h4a -           ## false-negative diagnosis of pre-symptomatic TB
     	diagptp_h4a*DP_h4a -           ## true-positive diagnosis of pre-symptomatic TB
     	natmort*DP_h4a +              ## natural mortality
     	artinit_h4*DP_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*DP_h4a ## Mortality due to HIV at h4a
     
     dDS_h4a <- accdiags_h4a*IS_h4a -     ## accessing dagnosis for symptamatic individuals
     	diagsfn_h4a*DS_h4a +           ## false-negative diagnosis of symptomatic TB
     	fnreltrans*accdiags_h4a*FN_h4a -              ## former false-negative people re-entering
     	diagstp_h4a*DS_h4a -           ## true-positive diagnosis of symptomatic TB
     	tbmort_h4a*DS_h4a -            ## mortality due to untreated TB
     	natmort*DS_h4a +              ## natural mortality
     	artinit_h4*DS_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*DS_h4a ## Mortality due to HIV at h4a
     
     dFN_h4a <- diagsfn_h4a*DS_h4a +      ## false-negative diagnosis of pre-symptomatic TB
     	diagstp_h4a*iltfu*DS_h4a - ## true-positive diagnosis of symptomatic TB but initial loss to follow-up-      ## false-negative diagnosis of pre-symptomatic TB
     	fnreltrans*accdiags_h4a*FN_h4a -              ## former false-negative people re-entering
     	tbrecov_h4a*FN_h4a -           ## natural recovery from TB (FN to recovered high-risk)
     	tbmort_h4a*FN_h4a -            ## mortality due to untreated TB
     	natmort*FN_h4a +              ## natural mortality
     	artinit_h4*FN_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*FN_h4a ## Mortality due to HIV at h4a
     
     dT_h4a <- diagptp_h4*(1-iltfu)*(artinit_tb)*DP_h4 + ## true-positive diagnosis of pre-symptomatic TB and ART initiation
     	diagstp_h4*(1-iltfu)*(artinit_tb)*DS_h4 + ## true-positive diagnosis of symptomatic TB and ART initiation
     	diagptp_h4a*DP_h4a +       ## true-positive diagnosis of pre-symptomatic TB
     	diagstp_h4a*(1-iltfu)*DS_h4a -           ## true-positive diagnosis of symptomatic TB
     	tfail*T_h4a -                 ## treatment failure due to incomplete treatment
     	recovtreat*T_h4a -            ## recovery after TB treatment
     	tbtreatmort_h4a*T_h4a -        ## mortality due to TB during treatment
     	natmort*T_h4a +               ## natural mortality
     	artinit_h4*T_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*T_h4a ## Mortality due to HIV at h4a
     
     dRH_h4a <- tfail*T_h4a +         ## partial recovery after incomplete treatment
     	tbrecov_h4a*FN_h4a +           ## natural recovery from TB (FN to recovered high-risk)
     	tbrecov_h4a*IP_h4a +           ## natural recovery from TB (pre-symptomatic to recovered high-risk)
     	tbrecov_h4a*IS_h4a -           ## natural recovery from TB (symptomatic to recovered high-risk)
     	TC*relreinfact_h4a*RH_h4a -     ## reinfection (recovered from active TB)
     	recovhl_h4a*RH_h4a  -        ## transition from recovery high-risk to recover low-risk or relapse
     	natmort*RH_h4a +    ## natural mortality
     	artinit_h4*RH_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*RH_h4a ## Mortality due to HIV at h4a
     
     dR_h4a <- recovtreat*T_h4a +        ## recovery after TB treatment
     	recovhl_h4a*(1-hrelapse_h4a)*RH_h4a - ## transition from recovery high-risk to recover low-risk
     	TC*relreinfact_h4a*R_h4a -     ## reinfection (revovered from active TB)
     	relapse_h4a*R_h4a -            ## relapse (endognous reactivation after TB treatment)
     	natmort*R_h4a +    ## natural mortality
     	artinit_h4*R_h4 -              ## ART initiation
     	hivmort_h4*hivmortart_ratio_h4*R_h4a ## Mortality due to HIV at h4a
     
     #########################
     ###  CASE DETECTION
     
     ### Cumulative TB cases detected
     
     dCumTBnotif <- diagptp_h0*(1-iltfu)*DP_h0 + diagstp_h0*(1-iltfu)*DS_h0 + # TRUE POSITIVES
     	diagptp_h1*(1-iltfu)*DP_h1 + diagstp_h1*(1-iltfu)*DS_h1 +
     	diagptp_h2*(1-iltfu)*DP_h2 + diagstp_h2*(1-iltfu)*DS_h2 +
     	diagptp_h3*(1-iltfu)*DP_h3 + diagstp_h3*(1-iltfu)*DS_h3 +
     	diagptp_h4*(1-iltfu)*DP_h4 + diagstp_h4*(1-iltfu)*DS_h4 +
     	diagptp_h1a*(1-iltfu)*DP_h1a + diagstp_h1a*(1-iltfu)*DS_h1a +
     	diagptp_h2a*(1-iltfu)*DP_h2a + diagstp_h2a*(1-iltfu)*DS_h2a +
     	diagptp_h3a*(1-iltfu)*DP_h3a + diagstp_h3a*(1-iltfu)*DS_h3a +
     	diagptp_h4a*(1-iltfu)*DP_h4a + diagstp_h4a*(1-iltfu)*DS_h4a +
     	diagpfp_h0*(1-iltfu)*DP_h0 + diagsfp_h0*(1-iltfu)*DS_h0 +            # FALSE POSITIVES
     	diagpfp_h1*(1-iltfu)*DP_h1 + diagsfp_h1*(1-iltfu)*DS_h1 +
     	diagpfp_h2*(1-iltfu)*DP_h2 + diagsfp_h2*(1-iltfu)*DS_h2 +
     	diagpfp_h3*(1-iltfu)*DP_h3 + diagsfp_h3*(1-iltfu)*DS_h3 +
     	diagpfp_h4*(1-iltfu)*DP_h4 + diagsfp_h4*(1-iltfu)*DS_h4 +
     	diagpfp_h1a*(1-iltfu)*DP_h1a + diagsfp_h1a*(1-iltfu)*DS_h1a +
     	diagpfp_h2a*(1-iltfu)*DP_h2a + diagsfp_h2a*(1-iltfu)*DS_h2a +
     	diagpfp_h3a*(1-iltfu)*DP_h3a + diagsfp_h3a*(1-iltfu)*DS_h3a +
     	diagpfp_h4a*(1-iltfu)*DP_h4a + diagsfp_h4a*(1-iltfu)*DS_h4a
     
     dCumFP <- diagpfp_h0*(1-iltfu)*DP_h0 + diagsfp_h0*(1-iltfu)*DS_h0 +            # FALSE POSITIVES
     	diagpfp_h1*(1-iltfu)*DP_h1 + diagsfp_h1*(1-iltfu)*DS_h1 +
     	diagpfp_h2*(1-iltfu)*DP_h2 + diagsfp_h2*(1-iltfu)*DS_h2 +
     	diagpfp_h3*(1-iltfu)*DP_h3 + diagsfp_h3*(1-iltfu)*DS_h3 +
     	diagpfp_h4*(1-iltfu)*DP_h4 + diagsfp_h4*(1-iltfu)*DS_h4 +
     	diagpfp_h1a*(1-iltfu)*DP_h1a + diagsfp_h1a*(1-iltfu)*DS_h1a +
     	diagpfp_h2a*(1-iltfu)*DP_h2a + diagsfp_h2a*(1-iltfu)*DS_h2a +
     	diagpfp_h3a*(1-iltfu)*DP_h3a + diagsfp_h3a*(1-iltfu)*DS_h3a +
     	diagpfp_h4a*(1-iltfu)*DP_h4a + diagsfp_h4a*(1-iltfu)*DS_h4a
     
     dCumTBinc <- delprog_h0*LD_h0 +      ## delayed progression (after distant infection)
     	rapprog_h0*LR_h0 +           ## rapid progression (after recent infection)
     	relapse_h0*R_h0 +            ## relapse (endogenous reactivation after TB treatment)
     	recovhl_h0*hrelapse_h0*RH_h0 +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h1*LD_h1 +      ## delayed progression (after distant infection)
     	rapprog_h1*LR_h1 +           ## rapid progression (after recent infection)
     	relapse_h1*R_h1 +            ## relapse (endogenous reactivation after TB treatment)
     	recovhl_h1*hrelapse_h1*RH_h1 +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h2*LD_h2 +      ## delayed progression (after distant infection)
     	rapprog_h2*LR_h2 +           ## rapid progression (after recent infection)
     	relapse_h2*R_h2 +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h2*hrelapse_h2*RH_h2 +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h3*LD_h3 +      ## delayed progression (after distant infection)
     	rapprog_h3*LR_h3 +           ## rapid progression (after recent infection)
     	relapse_h3*R_h3 +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h3*hrelapse_h3*RH_h3 +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h4*LD_h4 +      ## delayed progression (after distant infection)
     	rapprog_h4*LR_h4 +           ## rapid progression (after recent infection)
     	relapse_h4*R_h4 +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h4*hrelapse_h4*RH_h4 +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h1a*LD_h1a +      ## delayed progression (after distant infection)
     	rapprog_h1a*LR_h1a +           ## rapid progression (after recent infection)
     	relapse_h1a*R_h1a +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h1a*hrelapse_h1a*RH_h1a +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h2a*LD_h2a +      ## delayed progression (after distant infection)
     	rapprog_h2a*LR_h2a +           ## rapid progression (after recent infection)
     	relapse_h2a*R_h2a +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h2a*hrelapse_h2a*RH_h2a +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h3a*LD_h3a +      ## delayed progression (after distant infection)
     	rapprog_h3a*LR_h3a +           ## rapid progression (after recent infection)
     	relapse_h3a*R_h3a +            ## relapse (endognous reactivation after TB treatment)
     	recovhl_h3a*hrelapse_h3a*RH_h3a +## relapse (endogenous reactivation after high-risk recovery)
     	delprog_h4a*LD_h4a +      ## delayed progression (after distant infection)
     	rapprog_h4a*LR_h4a +           ## rapid progression (after recent infection)
     	relapse_h4a*R_h4a +           ## relapse (endogenous reactivation after TB treatment)
     	recovhl_h4a*hrelapse_h4a*RH_h4a ## relapse (endogenous reactivation after high-risk recovery)
     
     dCumTBinc_h0 <-  delprog_h0*LD_h0 +      ## delayed progression (after distant infection)
     	rapprog_h0*LR_h0 +           ## rapid progression (after recent infection)
     	relapse_h0*R_h0  +          ## relapse (endogenous reactivation after TB treatment)
     	recovhl_h0*hrelapse_h0*RH_h0 ## relapse (endogenous reactivation after high-risk recovery)
     
     dCumTBmort <- tbmort_h0*(IS_h0+FN_h0+DS_h0) + tbmort_h1*(IS_h1+FN_h1+DS_h1) + tbmort_h2*(IS_h2+FN_h2+DS_h2) + 
     	tbmort_h3*(IS_h3+FN_h3+DS_h3) + tbmort_h4*(IS_h4+FN_h4+DS_h4) + tbmort_h1a*(IS_h1a+FN_h1a+DS_h1a) + tbmort_h2a*(IS_h2a+FN_h2a+DS_h2a) + ### TB HIV mortality halve attributed to TB
     	tbmort_h3a*(IS_h3a+FN_h3a+DS_h3a) + tbmort_h4a*(IS_h4a+FN_h4a+DS_h4a) +
     	tbtreatmort_h0*T_h0 + tbtreatmort_h1*T_h1 + tbtreatmort_h2*T_h2 + tbtreatmort_h3*T_h3 + tbtreatmort_h4*T_h4 + tbtreatmort_h1a*T_h1a + 
     	tbtreatmort_h2a*T_h2a + tbtreatmort_h3a*T_h3a + tbtreatmort_h4a*T_h4a
     
     
      
     
     return(list(c(dS_h0, dLR_h0, dLD_h0, dIP_h0, dIS_h0, dDP_h0, dDS_h0, dFN_h0, dT_h0, dRH_h0, dR_h0,
     							dS_h1, dLR_h1, dLD_h1, dIP_h1, dIS_h1, dDP_h1, dDS_h1, dFN_h1, dT_h1, dRH_h1,  dR_h1,
     							dS_h2, dLR_h2, dLD_h2, dIP_h2, dIS_h2, dDP_h2, dDS_h2, dFN_h2, dT_h2, dRH_h2,  dR_h2,
     							dS_h3, dLR_h3, dLD_h3, dIP_h3, dIS_h3, dDP_h3, dDS_h3, dFN_h3, dT_h3, dRH_h3,  dR_h3,
     							dS_h4, dLR_h4, dLD_h4, dIP_h4, dIS_h4, dDP_h4, dDS_h4, dFN_h4, dT_h4, dRH_h4,  dR_h4,
     							dS_h1a, dLR_h1a, dLD_h1a, dIP_h1a, dIS_h1a, dDP_h1a, dDS_h1a, dFN_h1a, dT_h1a, dRH_h1a,  dR_h1a,
     							dS_h2a, dLR_h2a, dLD_h2a, dIP_h2a, dIS_h2a, dDP_h2a, dDS_h2a, dFN_h2a, dT_h2a, dRH_h2a,  dR_h2a,
     							dS_h3a, dLR_h3a, dLD_h3a, dIP_h3a, dIS_h3a, dDP_h3a, dDS_h3a, dFN_h3a, dT_h3a, dRH_h3a,  dR_h3a,
     							dS_h4a, dLR_h4a, dLD_h4a, dIP_h4a, dIS_h4a, dDP_h4a, dDS_h4a, dFN_h4a, dT_h4a, dRH_h4a,  dR_h4a,
     							dCumTBnotif, dCumFP, dCumTBinc, dCumTBinc_h0, dCumTBmort
     )))})
}


# Save TB model as a function
save(TB_model, file = tail(.args, 1))
