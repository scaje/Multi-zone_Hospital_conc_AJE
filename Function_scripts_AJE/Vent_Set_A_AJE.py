# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 14:40:50 2022

@author: scaje

This code involved the setup of the ventilation matrix for ventialtion 
setting A from Noakes et al [2009]. 

The values shown here represent those consistent with ventialton setting A
and can be modified accordingly.


Below represents homogenous mixing  on a 9 zone hopital ward with 
global mixing rate of 9m^3/min, a ward ventilation rate of 3ACH and an 
infectious person in zone 1 with 3 people in each zone. 


Created 10/02/2022 AJE 

"""

import numpy as np
import scipy as sp

##############################################################################
##############################################################################
############################## ZONAL SETUP ###################################
##############################################################################
##############################################################################
#number of zones n
n=9

##############################################################################
############################## geometry matrix SETUP #########################
##############################################################################
""" the aim of this nxn matrix, geometry(nxn), is to characterise the geometry
 and so if zone i is connected to zone j then entry geometry[i,j]=1,
 if zone i is not connected to zone j then geometry[1,j]=0."""
 

geometry=np.zeros((n,n))
geometry[0,1]=1
geometry[1,2]=1
geometry[2,5]=1
geometry[4,5]=1
geometry[3,4]=1
geometry[5,8]=1
geometry[7,8]=1
geometry[6,7]=1

for i in range(n):
    for j in range(n):    
        geometry[j,i] = geometry[i,j]

print("geometry matric geometry" +str(geometry))

##############################################################################
##############################################################################
###########################################################################
###############################################################

V=60
#Zonal volumes little v
v_zonal = np.zeros(n)
#for when volume is the same in each room
for i in range(n):
    v_zonal[i] = V
#If volumes are different
#v_zonal[0]=
#v_zonal[1]=
#v_zonal[2]=
#   .
#   .
#   .
#v_zonal[n]= 
print("volume v_zonal = " + str(v_zonal))
############################################################################
############################################################################
#zonal ventialtion rate
Q_zonal = np.zeros(n)
#for when ventilation rate is same in each zone
for i in range(n):
    Q_zonal[i] = 3 # ventilation rate is the same in each zone 

#if ventilation rates are different
#Q_zonal[0]=6
#Q_zonal[1]=0
#Q_zonal[2]=3
#Q_zonal[3]=6
#Q_zonal[4]=0
#Q_zonal[5]=3
#Q_zonal[6]=6
#Q_zonal[7]=0
#Q_zonal[8]=3
#add more for more zones, currently for 3 zones
print("ventialtion Q_zonal =" + str(Q_zonal))
###########################################################################

# boundary_flow defines flow from zone i to j
boundary_flow = np.zeros((n,n))
#Diagonal terms
for i in range(n):
    boundary_flow[i,i] = 0

 
 
#For setting A from noakes 2009
boundary_flow[0,1]=9
boundary_flow[1,0]=9
boundary_flow[1,2]=9
boundary_flow[2,1]=9
boundary_flow[2,5]=9
boundary_flow[5,2]=9
boundary_flow[3,4]=9
boundary_flow[4,3]=9
boundary_flow[4,5]=9
boundary_flow[5,4]=9
boundary_flow[5,8]=9
boundary_flow[8,5]=9
boundary_flow[6,7]=9
boundary_flow[7,6]=9          
boundary_flow[7,8]=9
boundary_flow[8,7]=9
#add more for more zones

print("flow from zone i to zone k boundary_flow =" + str(boundary_flow))
###########################################################################

###########################################################################
###########################################################################
#This defined the ventialtion matrix which is needed for the governing equations.
#This uses the parameters defined above. The inverse is also caluclated and used
#within the governing equations.


def VentilationMatrix(n, t, Q_zonal, geometry, boundary_flow):
    
    
    
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
