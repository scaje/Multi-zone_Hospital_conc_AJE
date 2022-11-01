# -*- coding: utf-8 -*-
"""
@author: scaje
 This code aims to solve for concentration of pathogen in the air and 
the SE epidemic model for a n zone hopsital ward (layout defined in geometry 
                                                  matrix)

The code calls on functions previously defined for the concentration and 
epidemic model. 

This code solves for the setup for a 12 zone hospital ward with zones defined as
Zone 1 = Bay 1 - Volume= 98.35m^3  
Zone 2 = Bay 2  -  Volume= 98.35m^3 
Zone 3 = Single Room 1  -  Volume=   28.57m^3
Zone 4 = Single Room 2  -  Volume= 28.57m^3   
Zone 5 = Nurse Room 1  -  Volume= 36.35m^3 
Zone 6 = Corridor 1  -  Volume= 31.74m^3  
Zone 7 = Corridor 2  -  Volume= 31.73m^3 
Zone 8 = Corridor 3  -  Volume= 31.74m^3  
Zone 9 = Sluice 1 ( a contamintaion room) Volume= 47.24m^3  
Zone 10 = Nurse Room 2 (staff room) Volume= 50.46m^3  
Zone 11 = Clinic room 1 Volume= 43.42m^3  
Zone 12 = Doctors office 1 Volume= 46.94m^3  


Each Bay has 4 fixed patients, the corridor has 4 suscpetibles (C1=1,C2=2,C3=1) and the nurse
room has 5 HCWs in total, 4 suscpetible and 1 infector, nobody is present in
the sluice, 3 nurses in staff room, 2 poeple in consulting room and 1 in doctors
office.

Ventilation setting for a 12zone ward is imported and represents

 

This code illustrates the following scenario
Infector HCW in nurse station = 30mins    ZONE 5
Infector HCW in bay 1 = 1 35mins     ZONE 1
Infector HCW in Bay 2 = 2 35mins      ZONE 2
Infector HCW in SR1 = 10mins     ZONE 3
Infector HCW in SR2= 10mins     ZONE 4
Infector HCW at nurse station between schedule = 30mins      ZONE 5
Then repeat round
Infector HCW in bay 1 = 1 35mins     ZONE 1
Infector HCW in Bay 2 = 2 35mins      ZONE 2
Infector HCW in SR1 = 10mins     ZONE 3
Infector HCW in SR2= 10mins     ZONE 4
Total = 4hrs
Any interactionin the corridor is ignored and assumed to be sufficiently
 small to exclude from the simulation.


User must define infected individuals and coresponding zone, and any time
periods or zones which the infectors move to, below.

Created 13/04/2022 AJE
"""


import matplotlib.pyplot as plt
import random, math
import numpy as np
from scipy.integrate import odeint #for solving odes
import matplotlib.colors as mcolors #colour name package
from matplotlib.pyplot import cm #colour map package
from pywaffle import Waffle #visual package for visuallising icons
import time #for live run time of code
import Contam_flows_12zone_AJE as VentMatrix #imports setup for 6 zone ward ventilaton setting from another file where it is already defined
from Contam_flows_12zone_AJE import VentilationMatrix #import function which defines ventilation matrix
from Contam_flows_12zone_AJE import InvVentilationMatrix #imports function which defines inverse ventilation matrix
from SE_Conc_Eqn_AJE import odes #imports predefined SE ode functions for transient concentration
from SE_Conc_Eqn_AJE import steadyodes ##imports predefined SE ode functions for steady concentration
from output_AJE import output_SE_Ct #This imports the plotting ode for all possible outputs for multizonal transient concentreation SE model
from boundary_flow_Contam_12zone_AJE import boundary_flow_contam #this import the function which changes the boundary flow values based on output from contam simulation


start_time = time.time() #time code started running
############################################################################
#Initial values from other simulations for reference
########################################################
###########################Initial values###################################


## DEFINE FILEPATH FOR EXPORTED CONTAM AIRFLOW SIMULATION RESULTS
filepath = r"C:\Users\scaje\OneDrive - University of Leeds\UNIV. OF LEEDS\PhD PROJECT\Ward Transmission\Code\Multi-zone_Hospital_Conc_AJE\Contam_export_scripts_AJE\Contam_flow_results_AJE.csv"




#outbreak parameters

#quanta rate = quanta/min . person (as 0.5 quanta per min)
q=0.5
#pulmonary rate = volume/min ( as 0.01volume/min)
p=0.01



##############################################################################
##############################################################################
############################## ZONAL SETUP ###################################
##############################################################################

######################Run ventilation setting############################ 
n=12
VentMatrix

###############################################################

#Pulmonary rate in each zone
p_zonal = np.zeros(n)
#for when  is the same in each room
for i in range(n):
    p_zonal[i] = p
#If volumes are different
#p_zonal[0]=
#p_zonal[1]=
#p_zonal[2]=
#   .
#   .
#   .
#p_zonal[n]=
print("Pulmonary rate p_zonal = " + str(p_zonal))
############################################################################

#Zonal quanta
q_zonal = np.zeros(n)
#for when  is the same in each room
for i in range(n):
    q_zonal[i] = q
#If volumes are different
#q_zonal[0]=
#q_zonal[1]=
#q_zonal[2]=
#   .
#   .
#   .
#q_zonal[n]=
print("quanta q_zonal = " + str(q_zonal))
############################################################################


 
#Zonal infections
#Zonal volumes little v
I_zonal = np.zeros(n)
#for when is the same in each room
#for i in range(n):
#    I_zonal[i] = I0
#If volumes are different
I_zonal[0]=0
I_zonal[1]=0
I_zonal[2]=0
I_zonal[3]=0
I_zonal[4]=1
I_zonal[5]=0

#   .
#   .
#   .
#I_zonal[n]= 
print("infections I_zonal = " + str(I_zonal))
##########################################################################

#############################################################################
#############################################################################
######################### DEFINE  Transient ODES #############################
#############################################################################


#################################################################
# declare the time vector in which to solve ODEs
#A grid of time points (in minutes)
#14400 time steps used, this is split up proportional to the length of each time period
#defined below - This is amount of seconds in a 4hr = 4x60x60 simulation
#t = np.linspace(0,240,14400) is total time period and vector length.



###########################################################################
################ transient infector #############################
##############################################################################
#Thinking about a transient infector i.e only present for a certain amount of time

#time periods of scenario
t1 = np.linspace(0,30,1800) # for infector present in nurse station - Zone 5

t2 = np.linspace(30,65,2100) # infector starts drug round in bay 1- zone 1

t3 = np.linspace(65,100,2100) # infector continues drug round bay 2- zone 2

t4 = np.linspace(100,110,600) # infector continues drug round single room 1- zone 3

t5 = np.linspace(110,120,600) # infector continues drug round single room 2 - zone 4

t6 = np.linspace(120,150,1800) # for infector present in nurse station - Zone 5

t7 = np.linspace(150,185,2100) # infector starts drug round in bay 1- zone 1

t8 = np.linspace(185,220,2100) # infector continues drug round bay 2- zone 2

t9 = np.linspace(220,230,600) # infector continues drug round single room 1- zone 3

t10 = np.linspace(230,240,600) # infector continues drug round single room 2 - zone 4

#NOTE: - time stepping for each second in a 4 hour simulation 
#that is 4 x 60 x 60 = 14,400 time steps, or 3600 time steps per hour



#infection vectors to correspon with scenario
I_zonal_t1 = I_zonal # for infector present in nurse station - Zone 5

I_zonal_t2 = np.zeros(n)
I_zonal_t2[0] = 1    # infector starts drug round bay 1- zone 1

I_zonal_t3 = np.zeros(n)
I_zonal_t3[1] = 1    # carries out drug round bay 2- zone 2

I_zonal_t4 = np.zeros(n)
I_zonal_t4[2] = 1    # carries out drug round single room 1 - zone 3

I_zonal_t5 = np.zeros(n)
I_zonal_t5[3] = 1    # carries out drug round single room 2- zone 4

I_zonal_t6 = np.zeros(n)
I_zonal_t6[4] = 1   # infector returns to nurse station - zone 5

I_zonal_t7 = np.zeros(n)
I_zonal_t7[0] = 1    # infector starts drug round bay 1- zone 1

I_zonal_t8 = np.zeros(n)
I_zonal_t8[1] = 1    # carries out drug round bay 2- zone 2

I_zonal_t9 = np.zeros(n)
I_zonal_t9[2] = 1    # carries out drug round single room 1 - zone 3

I_zonal_t10 = np.zeros(n)
I_zonal_t10[3] = 1    # carries out drug round single room 2- zone 4

#combining for solution loop

t=[t1,t2,t3,t4,t5,t6,t7,t8,t9,t10]
I_zonal_t=[I_zonal_t1,I_zonal_t2,I_zonal_t3,I_zonal_t4,I_zonal_t5,I_zonal_t6,I_zonal_t7,I_zonal_t8,I_zonal_t9,I_zonal_t10]
###########################################################################
############################################################################
##############################################################################
#################################################################

#Loops below calculate the solution for transient infector over any specified time periods

###########################################
##############  Transient #################
###########################################

#looping solutions over different time periods for transient model
C0 = np.zeros(n) #inital concentration
E0 = np.zeros(n) #inital exposed
S0 = VentMatrix.K_zonal - I_zonal_t1 - E0 #inital suceptibles
Ct = np.empty((0,n))
St = np.empty((0,n))
Et = np.empty((0,n))
#combining intial conditions
X0 = np.hstack( (C0, S0, E0) )
print(X0)

for i in range(len(t)):
    

    V_zonal = VentilationMatrix(n, t[i], VentMatrix.Q_zonal, VentMatrix.geometry, filepath)
    print("V_zonal = " + str(V_zonal)) #print bounday flow to check its updating each step

    
    #solving
    x = odeint(odes, X0, t[i], args=(n, V_zonal, I_zonal_t[i], q_zonal, p_zonal, VentMatrix.v_zonal)) #args=()
    
    #re-defining initial conditions
    C0 = x[:, 0:n]
    S0 = x[:, n:2*n]
    E0 = x[:, 2*n:3*n]
    
    
    X0 = np.hstack( (C0[-1,:], S0[-1,:], E0[-1,:]) )
    print(X0)
    
    
    #storing results in a vector
    Ct = np.vstack((Ct, C0))
    St = np.vstack((St,S0))
    Et = np.vstack((Et,E0))
    
    

    #End
#########################################
############ Steady State ###############
#########################################

E0star = np.zeros(n)
S0star = VentMatrix.K_zonal - I_zonal_t1 - E0star

Cstar = np.empty((0,n))
Ststar = np.empty((0,n))
Etstar = np.empty((0,n))

#initial condition
X0star = np.hstack((S0star, E0star))
print(X0star)

for j in range(len(t)):


    V_zonal = VentilationMatrix(n, t[j], VentMatrix.Q_zonal, VentMatrix.geometry, filepath)
    print("V_zonal = " + str(V_zonal)) #print bounday flow to check its updating each step
    V_zonal_inv = InvVentilationMatrix(V_zonal)
    print( "V_zonal_inv = " + str(V_zonal_inv))
    

    #solving steady state system with steady concentration
    x = odeint(steadyodes, X0star, t[j], args=(n, V_zonal_inv, I_zonal_t[j], q_zonal, p_zonal, VentMatrix.v_zonal))
    
    #rdfining initial conditions with stored solution
    S0star = x[:, 0:n]
    E0star = x[:, n:]
    
    
    #redefine initial conditions
    X0star = np.hstack((S0star[-1,:], E0star[-1,:]))
    print(X0star)
    
    #store values in a vector
    Ststar = np.vstack((Ststar, S0star))
    Etstar = np.vstack((Etstar, E0star))


    #defining the concentration value for each zone at each time period (cols represent zones, rows represent time periods)
    Cstar_t = np.matmul(V_zonal_inv, I_zonal_t[j]) * q_zonal
    Cstar_t = np.tile(Cstar_t, (len(t[j]), 1))
    Cstar = np.vstack((Cstar,Cstar_t))
    

#end

###########################################################################
###########################################################################
####################### population values###################################


#transinet version
St_pop = np.sum(St, axis=1) #axis=1 does rows, axis=0 does columns
Et_pop = np.sum(Et, axis=1)
S0_pop = np.sum(S0)
E0_pop = np.sum(E0)
I0_pop = np.sum(I_zonal)

#steady state
Ststar_pop = np.sum(Ststar, axis=1) #axis=1 does rows, axis=0 does columns
Etstar_pop = np.sum(Etstar, axis=1)
S0star_pop = np.sum(S0star)
E0star_pop = np.sum(E0star)






############################################################################
############################################################################
######################### Plotting #########################################
############################################################################
############################################################################

#define t for plotting 
#plotting in hours 
#t_hours = t/60

t_hours=np.empty((0,0))
for i in range(10):#note range runs for number of time periods defined
    t_hours = np.append(t_hours,t[i]/60) #to make all times plottable in hours not minutes
    



############################################################################
############################################################################

#uses a predefined function to plot all of the required outputs for this model
output_SE_Ct(5, t_hours, Ct, Cstar, St, Et, Ststar, Etstar, St_pop, Et_pop, Ststar_pop, Etstar_pop, I0_pop, start_time)



boundary_flow_test = boundary_flow_contam( r"C:\Users\scaje\OneDrive - University of Leeds\UNIV. OF LEEDS\PhD PROJECT\Ward Transmission\St_James_ward_info\CONTAM_Model_Zeyu\Contam_data\Windows_closed.csv", n, t[3], VentMatrix.geometry)
print("boundary flow test" + str(boundary_flow_test))