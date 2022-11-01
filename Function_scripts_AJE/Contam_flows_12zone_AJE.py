# -*- coding: utf-8 -*-
"""

@author: scaje

This code sets out the setup for a 12 zone hospital ward  with 
the geometry taken as a subset of a UK adult Respiratory ward with 
the zones defined as
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

total ward volume=573.41


Each Bay has 4 fixed patients, the corridor has 4 suscpetibles (C1=1,C2=2,C3=1) and the nurse
room has 5 HCWs in total (4 suscpetible and 1 infector), nobody is present in
the sluice, 3 nurses in staff room, 2 poeple in consulting room and 1 in doctors
office to simulate expected average occupancy.


This code is designed to be imported into code which solves the governing equations.

This file is set up to deal with imported boundary flows from the results of a contam simulation.
This fuction is called within the ventilation matrix function.

Created 08/03/2022 AJE 

"""

import numpy as np
import scipy as sp
from boundary_flow_Contam_12zone_AJE import boundary_flow_contam #imports boundary flow function which uses cintam results to define them 


##############################################################################
##############################################################################
############################## ZONAL SETUP ###################################
##############################################################################
##############################################################################
#number of zones n
n=12

##############################################################################
############################## geometry matrix SETUP #########################
##############################################################################
""" the aim of this nxn matrix, geometry(nxn), is to characterise the geometry
of the zonal set-up and so if zone i is connected to zone j then entry 
geometry[i,j]=1, if zone i is not connected to zone j then geometry[1,j]=0."""
 
#defined in such away that input should be [i,j] where i<j
geometry=np.zeros((n,n))
geometry[0,5]=1
geometry[1,7]=1
geometry[2,4]=1
geometry[3,4]=1
geometry[4,6]=1
geometry[5,6]=1
geometry[6,7]=1
geometry[6,8]=1
geometry[6,9]=1
geometry[7,10]=1
geometry[7,11]=1




for i in range(n):
    for j in range(n):    
        geometry[j,i] = geometry[i,j]

print("geometry matric geometry" +str(geometry))

##############################################################################
#Zonal volumes little v
v_zonal = np.zeros(n)
#for when volume is the same in each room
#for i in range(n):
#    v_zonal[i] = V
#If volumes are different
v_zonal[0]=98.35
v_zonal[1]=98.35
v_zonal[2]=28.57
v_zonal[3]=28.57
v_zonal[4]=36.35
v_zonal[5]=31.74
v_zonal[6]=31.74
v_zonal[7]=31.74
v_zonal[8]=47.24
v_zonal[9]=50.46
v_zonal[10]=43.42
v_zonal[11]=46.94

#   .
#   .
#   .
#v_zonal[n]= 
print("volume v_zonal = " + str(v_zonal))
##########################################################################

#number of people in each zone is not the same due to scenario

K_zonal = np.zeros(n)
K_zonal[0]=4
K_zonal[1]=4
K_zonal[2]=1
K_zonal[3]=1
K_zonal[4]=5
K_zonal[5]=1
K_zonal[6]=2
K_zonal[7]=1
K_zonal[8]=0
K_zonal[9]=3
K_zonal[10]=2
K_zonal[11]=1

##############################################################################
###########################################################################

#zonal ventialtion rate m^3/min
Q_zonal = np.zeros(n)

#NOTE: Q_zonal[i] here has been calculated as a proportion of the volume size
#in order to lead to specific ACH rates - the ones used in this study have been
#pre-calculated and including in the commenting for each zone below.


Q_zonal[0]=4.91 #for 3ach= 4.91 #for1.5ach = 2.45 for 0.5ach = 0.81 #for 6ACH=9.82
Q_zonal[1]=4.91 #for 3ach=4.91 #for1.5ach = 2.45 for 0.5ach = 0.81 #for 6ACH=9.82
Q_zonal[2]=1.42 #for 3ach=1.42 #for1.5ach = 0.71 for 0.5ach = 0.23 #for 6ACH=2.84
Q_zonal[3]=1.42 #for 3ach=1.42 #for1.5ach = 0.71 for 0.5ach = 0.23 #for 6ACH=2.84
Q_zonal[4]=1.81 #for 3ach=1.81 #for1.5ach = 0.9 for 0.5ach = 0.3 #for 6ACH=3.62
Q_zonal[5]=1.59 #for 3ach=1.59 #for1.5ach = 0.795 for 0.5ach =0.265  #for 6ACH=3.18
Q_zonal[6]=1.59 #for 3ach=1.59 #for1.5ach = 0.795 for 0.5ach =0.265  #for 6ACH=3.18
Q_zonal[7]=1.59#for 3ach=1.59 #for1.5ach = 0.795 for 0.5ach =0.265  #for 6ACH=3.18
Q_zonal[8]=2.36 #for 3ach=2.36 #for1.5ach = 1.18 for 0.5ach = 0.4 #for 6ACH=4.72
Q_zonal[9]=2.52 #for 3ach=2.52 #for1.5ach = 1.26 for 0.5ach = 0.42 #for 6ACH=5.04
Q_zonal[10]=2.17 #for 3ach=2.17 #for1.5ach = 1.08 for 0.5ach = 0.36 #for 6ACH=4.34
Q_zonal[11]=2.32 #for 3ach=2.32 #for1.5ach = 1.16 for 0.5ach = 0.39 #for 6ACH=4.64

#add more for more zones, currently for 3 zones
print("ventialtion Q_zonal =" + str(Q_zonal))
###########################################################################







###########################################################################
#The function below uses the above set-up, alongside inported boundary flow values
#to define and caluclate the ventialtion matrix, required to solve the governing equations
#The inverse ventialtion matrix is also calculated below

def VentilationMatrix(n, t, Q_zonal, geometry, filepath):
    
    
    #define boundary flow matrix from boundary flow contam func (uses contam flows)
    boundary_flow = boundary_flow_contam(filepath, n, t, geometry)
    
    
    #VENTILATION MATRIX
    V_zonal = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i==j:
                V_zonal[i,i] = Q_zonal[i]
            
                for k in range(n):
                    if geometry[i,k] > 0:
                        
                        V_zonal[i,i] = V_zonal[i,i] + boundary_flow[i,k]
                    else:
                        V_zonal[i,i] = V_zonal[i,i]
                    
                    
            else:
                if geometry[i,j]>0:
                        
                    V_zonal[i,j] = - boundary_flow[j,i]
                else:
                    V_zonal[i,j] = 0
                
    print("Ventilation Matrix V = " + str(V_zonal))   
    
    
    return V_zonal




def InvVentilationMatrix(V_zonal):
    
    
    #calculate inverse of ventiation matrix for steady state calculation
    V_zonal_inv = sp.linalg.inv(V_zonal)
    print("Inverse Ventilation Matrix V = " + str(V_zonal_inv))
    
    return V_zonal_inv
###########################################################################