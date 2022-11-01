# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 10:30:59 2022

@author: scaje
 
This code defines a function which is designed to provide all of the plotting 
output for the Susceptible-Exposed (SE) Transient concentration model and 
Transient vs. steady state bersion of this model. 
Thisis designed to work with an number of n multizonal model within an allocated
time period. The ouput is as follows:
    - Transient concentration over time
    -transient vs steady state concentration
    -popn level eposure and susceptible trans vs steady state

    
    Created 14/03/2022 AJE
"""
import matplotlib.pyplot as plt
import random, math
import numpy as np
import matplotlib.colors as mcolors #colour name package
from matplotlib.pyplot import cm #colour map package
from matplotlib.ticker import StrMethodFormatter #to set decimal places in axis
import time #for live run time of code

#Define Colours array for plotting 
#cm is imported colour maps from matplotlib.pyplot library
#tab20 is selected color map palette
# n is number of different colours needed
#This has to be applied before each use
#colour = iter(cm.tab20(np.linspace(0, 1, n)))

def output_SE_Ct(n, t, Ct, Cstar, St, Et, Ststar, Etstar, St_pop, Et_pop, Ststar_pop, Etstar_pop, I0_pop, start_time):
    
    ######################### Concentration #########################################
    ############################################################################
    ###### Plotting transient solution to concentration 
    for i in range(n):
        plt.plot(t, Ct[:,i], label='Zone %s' %(i+1))
    plt.title("Concentration of Pathogen C(t)")
    plt.xlabel("Time [hr]")
    plt.ylabel("Concentration [$qm^{-3}$]")
    #plt.legend(title='Concentration')
    plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Concentration')
    plt.show()
    
    
    ############################################################################
    ###### Plotting transient concentration with steady state
    #plt.figure(dpi=2000)#set dots per inch for better quality images
    
    colour = iter(cm.tab10(np.linspace(0, 1, n))) #defining colour map for new loop
    #for more than 10 items use tab20 or tab10 when less than 10
    for i in range(n):
        c = next(colour)#choosing next random colour for plotting
        plt.plot(t, Ct[:,i],color=c, label = 'Zone %s' %(i+1))
        plt.plot(t, Cstar[:,i],color=c, linestyle='dashed')#, label = 'Zone %s - Steady State' %(i+1) )
    #plt.title("Concentration of Pathogen C(t)")
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # to set the decimal places used in yaxis
    plt.xlabel("Time [hr]")
    plt.ylabel("Concentration [$qm^{-3}$]")
    plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Concentration')
    plt.show()
    
    
    ###########################################################################
    ######################### population plotting #############################
    ###########################################################################
    #normal population
    plt.plot(t,St_pop, label='Susceptible')
    plt.plot(t,Et_pop, label = 'Exposed')
    plt.title("Epidemic model of total population")
    plt.xlabel("Time [hr]")
    plt.ylabel("People")
    plt.legend()
    plt.show()
    
    #full epidemic model plot suscepible vs exposed
    plt.plot(t,St_pop,'tab:blue', label='Susceptible')
    plt.plot(t, Ststar_pop, 'tab:blue', linestyle='dashed' , label='Susceptible - Steady C*')
    plt.plot(t,Et_pop, 'tab:red', label = 'Exposed')
    plt.plot(t, Etstar_pop, 'tab:red', linestyle='dashed' , label='Exposed - Steady C*')
    plt.title("Epidemic model of total population (steady vs transient)")
    plt.xlabel("Time [hr]")
    plt.ylabel("People")
    plt.legend()
    plt.show()
    
    #full population susceptible plot
    plt.plot(t,St_pop,'tab:blue', label='Susceptible')
    plt.plot(t, Ststar_pop, 'tab:blue', linestyle='dashed' , label='Susceptible - Steady C*')
    plt.title("Susceptible model of total population (steady vs transient)")
    plt.xlabel("Time [hr]")
    plt.ylabel("People")
    plt.legend(loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Suscpetible')
    #plt.legend()
    plt.show()
    
    #full population exposed plot
    #plt.figure(dpi=2000)#set dots per inch for better quality images
    
    plt.plot(t,Et_pop,'tab:red', label = 'Transient C(t)')
    plt.plot(t, Etstar_pop, 'tab:red', linestyle='dashed' , label='Steady C*')
    #plt.title("Number of Exposures E(t) in the Population")
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # to set the decimal places used in yaxis
    plt.xlabel("Time [hr]")
    plt.ylabel("Number of Exposed E(t)")
    plt.legend(fontsize=15)#loc='center left',prop={'size': 8}, bbox_to_anchor=(1, 0.5), title='Exposed')
    #plt.legend()
    plt.show()
    
    
    print("--- Run Time = %s seconds ---" % (time.time() - start_time))
