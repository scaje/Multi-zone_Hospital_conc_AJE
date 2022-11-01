# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 12:19:56 2021

@author: scaje
This code aims to solve for concentration of pathogen in the air and 
the SE epidemic model for a n zone hopsital ward (layout defined in geometry 
                                                  matrix)

The code calls on functions previously defined for the concentration and 
epidemic model. 

The ventialtion setting is imported also. 

This code shows ventialtion setting A being run for specific parameters but this
is to be used as an example and can be altered to any parameters or ventialtion
settings for any number of zones. 
"""


import matplotlib.pyplot as plt
import random, math
import numpy as np
from scipy.integrate import odeint #for solving odes
import matplotlib.colors as mcolors #colour name package
from matplotlib.pyplot import cm #colour map package
from pywaffle import Waffle #visual package for visuallising icons
from matplotlib.ticker import StrMethodFormatter #to set decimal places in axis
import time #for live run time of code
import Vent_Set_A_AJE as VentMatrix #imports setup for 6 zone ward ventilaton setting from another file where it is already defined
from Vent_Set_A_AJE import VentilationMatrix #import function which defines ventilation matrix
from Vent_Set_A_AJE import InvVentilationMatrix #imports function which defines inverse ventilation matrix
from SE_Conc_Eqn_AJE import odes #imports predefined SE ode functions for transient concentration
from SE_Conc_Eqn_AJE import steadyodes ##imports predefined SE ode functions for steady concentration
""" To run the code, one needs to define the initial parameters below and then 
define the setup matrices for each specific zone. A time period and timestep needs to be defined
and then the transient concentrations can be solved and the graphs can be produced. 

before running the code, itial conditions for the epidemic model also need to be defined.
These can be found after the ode definition function.
"""

start_time = time.time() #time code started running
############################################################################
#Initial values from other simulations for reference
########################################################
###########################Initial values###################################




#quanta rate = quanta/min . person (as 0.5 quanta per min)
q=0.5
#pulmonary rate = volume/min ( as 0.01volume/min)
p=0.01
#Ventilaation rate = volume/ min (taken as 27m^3/min)
Q=3 #3m^3 in each zone 
#volume of indoor space V m^3
V=60
#desired number of poeple in each of the zones
K=3






##############################################################################
##############################################################################
############################## ZONAL SETUP ###################################
##############################################################################

######################Run ventilation setting############################ 
n=9


##########################################################################


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
I_zonal = np.zeros(n)
#for when is the same in each room
#for i in range(n):
#    I_zonal[i] = I0
#If volumes are different
I_zonal[0]=1
I_zonal[1]=0
I_zonal[2]=0
I_zonal[3]=0
I_zonal[4]=0
I_zonal[5]=0
I_zonal[6]=0
I_zonal[7]=0
I_zonal[8]=0
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
t = np.linspace(0,240,14400)
#################################################################

C0 = np.zeros(n) #inital concentration
E0 = np.zeros(n) #inital exposed
S0 = K - I_zonal - E0 #inital suceptibles


 

#combining intial conditions
X0 = np.hstack( (C0, S0, E0) )
print(X0)
#Solving the ODEs
V_zonal = VentilationMatrix(n, t, VentMatrix.Q_zonal, VentMatrix.geometry, VentMatrix.boundary_flow)

x = odeint(odes, X0, t, args=(n, V_zonal, I_zonal, q_zonal, p_zonal, VentMatrix.v_zonal))

Ct = x[:, 0:n]
St = x[:, n:2*n]
Et = x[:, 2*n:3*n]


print(Ct.shape)
print(St.shape)
print(Et.shape)



###########################################################################
###########################################################################
####################### Steady State ######################################


#Define Steady State Concentration 
#Cstar = np.matmul(np.matmul(V_zonal_inv , q_zonal),I_zonal)

V_zonal = VentilationMatrix(n, t, VentMatrix.Q_zonal, VentMatrix.geometry, VentMatrix.boundary_flow)
print("V_zonal = " + str(V_zonal)) #print bounday flow to check its updating each step
V_zonal_inv = InvVentilationMatrix(V_zonal)
print( "V_zonal_inv = " + str(V_zonal_inv))

Cstar = np.matmul(V_zonal_inv, I_zonal) * q_zonal


E0star = np.zeros(n)
S0star = K - I_zonal - E0star

X0star = np.hstack((S0star, E0star))
print(X0star)
#solving steady odes



x = odeint(steadyodes, X0star, t, args=(n, V_zonal_inv, I_zonal, q_zonal, p_zonal, VentMatrix.v_zonal))

Ststar = x[:, 0:n]
Etstar = x[:, n:]
    

###########################################################################
###########################################################################
####################### population values###################################

#Total population
N = K*n #K people in each zone
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
t_hours = t/60

#Define Colours array for plotting 
#cm is imported colour maps from matplotlib.pyplot library
#tab20 is selected color map palette
# n is number of different colours needed
#This has to be applied before each use
#colour = iter(cm.tab20(np.linspace(0, 1, n)))



######################### Concentration #########################################
############################################################################
###### Plotting transient solution to concentration 
for i in range(n):
    plt.plot(t_hours,Ct[:,i], label='Zone %s' %(i+1))
plt.title("Concentration of Pathogen C(t)")
plt.xlabel("Time [hr]")
plt.ylabel("Concentration [$qm^{-3}$]")
#plt.legend(title='Concentration')
plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Concentration')
plt.show()


############################################################################
###### Plotting transient concentration with steady state
#plt.figure(dpi=1000)#set dots per inch for better quality images

colour = iter(cm.tab20(np.linspace(0, 1, n))) #defining colour map for new loop
for i in range(n):
    c = next(colour)#choosing next random colour for plotting
    plt.plot(t_hours, Ct[:,i],color=c, label = 'Zone %s' %(i+1))
    plt.hlines(Cstar[i],0,t_hours[-1],color=c, linestyle='dashed' , label = 'Zone %s - Steady State' %(i+1) )#, label = 'Zone %s - Steady State' %(i+1) 
plt.title("Concentration of Pathogen C(t)")
plt.xlabel("Time [hr]")
plt.ylabel("Concentration [$qm^{-3}$]")
plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Concentration')
plt.show()








###########################################################################
######################### population plotting #############################
###########################################################################
#normal population
plt.plot(t_hours,St_pop, label='Susceptible')
plt.plot(t_hours,Et_pop, label = 'Exposed')
plt.title("Epidemic model of total population")
plt.xlabel("Time/Hours")
plt.ylabel("Number of people")
plt.legend()
plt.show()

#full epidemic model plot suscepible vs exposed
plt.plot(t_hours,St_pop,'tab:blue', label='Susceptible')
plt.plot(t_hours, Ststar_pop, 'tab:blue', linestyle='dashed' , label='Susceptible - Steady C*')
plt.plot(t_hours,Et_pop, 'tab:red', label = 'Exposed')
plt.plot(t_hours, Etstar_pop, 'tab:red', linestyle='dashed' , label='Exposed - Steady C*')
plt.title("Epidemic model of total population (steady vs transient)")
plt.xlabel("Time/Hours")
plt.ylabel("Number of people")
plt.legend()
plt.show()

#full population susceptible plot
plt.plot(t_hours,St_pop,'tab:blue', label='Susceptible')
plt.plot(t_hours, Ststar_pop, 'tab:blue', linestyle='dashed' , label='Susceptible - Steady C*')
plt.title("Susceptible model of total population (steady vs transient)")
plt.xlabel("Time/Hours")
plt.ylabel("Number of people")
plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Suscpetible')
#plt.legend()
plt.show()

#full population exposed plot
#plt.figure(dpi=1000)#set dots per inch for better quality images

plt.plot(t_hours,Et_pop,'tab:red', label = 'Transient C(t)')
plt.plot(t_hours, Etstar_pop, 'tab:red', linestyle='dashed' , label='Steady-state C*')
#plt.title("Exposed model of total population")
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # to set the decimal places used in yaxis
plt.xlabel("Time [hr]")
plt.ylabel("Number of Exposed E(t)")
plt.legend(fontsize=15)#loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Exposed')
#plt.legend()
plt.show()


print("--- Run Time = %s seconds ---" % (time.time() - start_time))

