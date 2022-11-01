# Multi-zone_Hospital_conc_AJE
This repository reproduces the research contained in the manuscript "A Mathematical Model for Assessing Transient Airborne Infection Risk in Hospital Ward"; Alexander J. Edwards, Lee Benson, Zeyu Guo, Martin Lopez-Garcia, Catherine J. Noakes, Daniel Peckham, Marco-Felipe King.

# Software
This code is written using Python in Spyder 4.1.4. Users will also require CONTAM 3.4.0.3 to reproduce the airflow simulations.

# Code

The script entitled 'FixedInf_9Zone_Conc_AJE.py' will reproduce the results used in Case study 1. This script uses pre-defined functions which can be found in the folder 'Function_scripts_AJE'. In particular, 'Function_scripts_AJE/Vent_Set_A_AJE.py' is a script which defines the set-up of Ventilation Setting A (Originally described Noakes et. al. [2009]). It then solves governing equations, which are set out in the script file ''Function_scripts_AJE/SE_Conc_Eqn_AJE.py'. 


The script entitled 'InfHCW_12zone_Conc_AJE.py' will reproduce the results used in Case study 2. To run this script you must define the file path for the file containing the exported CONTAM results at the beginning of the code (Line 84 in the code). This script uses pre-defined functions which can be found in the folder 'Function_scripts_AJE'. In particular, 'Function_scripts_AJE/Contam_flows_12zone_AJE.py' defines the ventialtion matrix set-up. This uses the script 'Contam_export_scripts_AJE/boundary_flow_Contam_12zone_AJE' which has previously calculated and defined the inter-zonal flow values using exported CONTAM airflow simulations. The outputs from this script are produced by calling a pre-defined output script which can be found in 'Function_scripts_AJE/output_AJE.py'.

The scripts contained in 'Contam_export_scripts_AJE' provide the functions which convert exported airflow results from CONTAM simulations into the inter-zonal flow values to be used within the transmission model. The exported results used in this case can be seen in 'Contam_flow_results.csv'. Note that the script 'Contam_export_scripts_AJE/boundary_flow_Contam_12zone_AJE' indexes according to the particular set-up of the results found in 'Contam_flow_results.csv'. Should you use an alternative results format, this may need to be altered. 

The scripts contained within the folder 'Function_scripts_AJE' contains functions which are used within the the two main scripts; 'FixedInf_9Zone_Conc_AJE.py' and 'InfHCW_12zone_Conc_AJE.py'. 'Function_scripts_AJE/Contam_flows_12zone_AJE.py' defines the ventialtion matrix set-up using inter-zonal flow values which have been exported from airflow simulations in CONTAM. 'Function_scripts_AJE/output_AJE.py' is predefined to produce particular outputs used in Case Study 2. 'Function_scripts_AJE/Vent_Set_A_AJE.py' is a script which defines the set-up of Ventilation Setting A (Originally described Noakes et. al. [2009]). 'Function_scripts_AJE/SE_Conc_Eqn_AJE.py' defines the governing equations for the transmission model.

The scripts contained within the folders 'Contam_export_scripts_AJE' and 'Function_scripts_AJE' must be callable file paths when running 'FixedInf_9Zone_Conc_AJE.py' and 'InfHCW_12zone_Conc_AJE.py', to reproduce the research in Case Study 1 and Case Study 2, respectively.
