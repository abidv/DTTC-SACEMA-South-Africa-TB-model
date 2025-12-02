## DTTC-SACEMA South Africa Tuberculosis Model 

This repository provides the R code for a deterministic, compartmental transmission-dynamic model of the tuberculosis (TB) epidemic in South Africa. The model was developed through a collaboration between the Desmond Tutu TB Centre (DTTC) and the South African Centre for Epidemiological Modelling and Analysis (SACEMA), and builds on an earlier version developed between 2017 and 2020.
The model structure and underlying assumptions are described in detail in the accompanying publication. In brief, the model incorporates established TB transmission dynamics alongside additional features to represent transitions and losses along the TB care cascade. It distinguishes between subclinical and clinical TB disease states and further differentiates individuals based on their care-seeking status. HIV infection, disease progression, and antiretroviral therapy are also explicitly modelled.

## Related publication

**Impact of prolonged tuberculosis healthcare service disruptions following the COVID-19 pandemic in South Africa – a mathematical modelling study***

Authors:  Abigail K. de Villiers, Sue-Ann Meehan, Muhammad Osman, Rory Dunbar, James A. Seddon, Cari van Schalkwyk, Gerald J. Maarman, Karen Du Preez, Anneke C. Hesseling, Florian M. Marx, for the DTTC-SACEMA TB Modelling Group

[![DOI](https://)


## License 

Impact of prolonged tuberculosis healthcare service disruptions following the COVID-19 pandemic in South Africa – a mathematical modelling study © 2025 by Abigail K. de Villiers, Florian M. Marx is licensed under CC BY-NC-SA 4.0. ![License: CC BY-NC-SA 4.0](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)  
Unless otherwise noted, all content associated with this work is licensed under CC BY-NC-SA 4.0.   
To view a copy of this license, visit [Creative Commons BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). 



## Directory structure

```
.
├── code
├── makefile
├── INPUTS
├── OUTPUTS
```

- *code* : model R code
- *makefile* : (optional) script to run the model using the makefile system 
- *INPUTS* : excel sheet containing all input data 
- *OUTPUTS* :  functions and data generated during the analyses


## Required packages 

The following R packages are essential for running this model:
- readxl: For reading data from Excel files.
- tictoc: For timing code execution, allowing you to monitor runtime.
- foreach: Provides a simple syntax for parallel loop execution.
- doParallel: Supports parallel backend registration for foreach loops.
- deSolve: Enables solving differential equations, which is central to this transmission model.
- ggplot2: For generating high-quality data visualizations.
- stringr: For working with and manipulating string data.
- plyr: Simplifies data manipulation tasks, especially for the data frames.


## Analysis order


#### 1. Read in Time-Varying Parameters and Define Time-Varying Functions
Reading in necessary parameter values and defining functions governing the equations for calculating time-varying parameters: HIV incidence, ART initiation rates, ART initiation upon initiating TB treatment, and probability of HIV at birth.

- **Scripts:** `./code/1a_Read_tv_data.R` & `./code/1b_Define_tv_functions.R`
- **Dependencies:**
  - For `tv_HIV_ART_data.RDS`: `./code/1a_Read_tv_data.R` & `./INPUTS/RSA/RSA_inputdata.xlsx` (Sheet: `1_time_varying_parms`)
  - For `tv_functions.RData`: `./code/1b_Define_tv_functions.R` & `./OUTPUTS/RSA/tv_HIV_ART_data.RDS`
- **Outputs:**
  - `./OUTPUTS/RSA/tv_HIV_ART_data.RDS`
  - `./OUTPUTS/RSA/tv_functions.RData`


#### 2. Generate Parameter List
This step generates a sampled list of parameters used in the model. Users can specify the number of parameter sets they want to sample (eg. nrowparammat = 1000). 
Users can define the years for which they want to run the TB model either in this script or read them in from `./INPUTS/RSA/RSA_inputdata.xlsx` (sheet: "runmodel_years")

- **Scripts:** `./code/2_Generate_parameter_List.R`
- **Dependencies:** 
  - `./INPUTS/RSA/RSA_inputdata.xlsx`  (Sheet: `2_generate_param_list` and Sheet: "runmodel_years")
- **Outputs:**
  - `./OUTPUTS/RSA/initialParametersets.RData`: Contains the generated list of initial parameter sets.


#### 3. Add Parameters Governing HIV and ART Status
This step adds parameters that govern HIV and ART status to the parameter list generated above. Parameters values are determined based on HIV and ART status affects the natural history of TB.

- **Scripts:** `./code/3_Add_HIV_ART_parameters.R`
- **Dependencies:** 
  - `./OUTPUTS/RSA/initialParametersets.RData`
  - `./INPUTS/RSA/RSA_inputdata.xlsx` (Sheet: `2_generate_param_list`)
- **Outputs:**
  - `./OUTPUTS/RSA/Parametersets_added_HIV_ART.RData`: Contains the updated parameter sets that now include parameters related to HIV and ART status.


#### 4. Defining Model Start Population
This step defines the initial population sizes for the model for the specified start year e.g. 1995. 

- **Scripts:** `./code/4_startpopulation.R`
- **Dependencies:** 
  - `./OUTPUTS/RSA/Parametersets_added_HIV_ART.RData`
- **Outputs:**
  - `./OUTPUTS/RSA/yinit.RData`: Contains the initialized population data for the model.


#### 5. Transmission Dynamic Model Equations & Care Cascade Function
This step involves defining the model functions which included the differential equations governing transmission dynamics as well as functions defining the TB care cascade.

- **Scripts:** 
  - `./code/5_TB_model.R`
  - `./code/5a_cc_functions.R`
- **Outputs:**
  - `./OUTPUTS/RSA/TB_model.RData`: Contains the model differential equations related to TB dynamics.
  - `./OUTPUTS/RSA/cc_functions.RData`: Contains functions used in the care cascade model.


#### 6. Run Model "Offline" - Run Model in Parallel on Local Machine
This step runs the model in parallel across multiple cores using all the previously defined parameters and functions. 
A function is created to run the model and produce a results list of all key TB indicators and outcome measures. 


- **Scripts:** `./code/6_Run_model_parallel.R`

- **Dependencies:** 
  - **tv_targets:** 
    - `./OUTPUTS/RSA/tv_HIV_ART_data.RDS`
    - `./OUTPUTS/RSA/tv_functions.RData`
  - **param_targets:** 
    - `./OUTPUTS/RSA/initialParametersets.RData`
    - `./OUTPUTS/RSA/Parametersets_added_HIV_ART.RData`
    - `./OUTPUTS/RSA/yinit.RData`
  - **fun_targets:** 
    - `./OUTPUTS/RSA/TB_model.RData`
    - `./OUTPUTS/RSA/cc_functions.RData`

- **Outputs:**
  - `./OUTPUTS/RSA/results.RData`: Contains the results generated by running the model.




## Additional included files 

- `./OUTPUTS/RSA Paramlist.f.RData`: Contains the final 1,000 parameter sets selected after model calibration 

-  `./OUTPUTS/RSA /yinit.f.RData`:  Contains the final 1,000 initial population sizes for the calibrated model. 

These files can be used to run the final calibrated baseline TB model. 

