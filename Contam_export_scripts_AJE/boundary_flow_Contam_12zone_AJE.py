# -*- coding: utf-8 -*-
"""
Created on Mon May  9 15:58:09 2022

@author: scaje
"""
from read_csv_AJE import ReadCSV
import numpy as np



def boundary_flow_contam(filepath,n,t,geometry):
    flows = ReadCSV(filepath)
    print(flows)
    
    
    # boundary_flow defines flow from zone i to j
    boundary_flow = np.zeros((n,n))
    #Diagonal terms are zero
    for i in range(n):
        boundary_flow[i,i] = 0
    
    #other entries in boundary flow matrix
    #For homogenous mixing
    #for i in range(n-1):#n-1 since considering i+1 in indices below 
     #       boundary_flow[i,i+1]=9
     #       boundary_flow[i+1,i]=9
     
    
     
    # t is a vector of vectors e.g t[i][0] calls the first entry of the ith vector in t
    # then since the format of the excel data includes 2 entries for each flow value we take steps in 2 rather than 1
    #the second index then selects the colum which refers to the boundary

    boundary_flow[0,5] = flows[int((t[0])*2),6] + flows[int((t[0])*2+1),6]
    boundary_flow[1,7] = flows[int((t[0])*2),8] + flows[int((t[0])*2+1),8]
    boundary_flow[2,4] = flows[int((t[0])*2),4] + flows[int((t[0])*2+1),4]
    boundary_flow[3,4] = flows[int((t[0])*2),5] + flows[int((t[0])*2+1),5]
    boundary_flow[4,6] = flows[int((t[0])*2),7] + flows[int((t[0])*2+1),7]
    boundary_flow[5,6] = flows[int((t[0])*2),10] + flows[int((t[0])*2+1),10]
    boundary_flow[6,7] = flows[int((t[0])*2),11] + flows[int((t[0])*2+1),11]    
    boundary_flow[6,8] = -(flows[int((t[0])*2+1),13] +flows[int((t[0])*2),13])
    boundary_flow[6,9] = -(flows[int((t[0])*2+1),14] + flows[int((t[0])*2),14])
    boundary_flow[7,10] = -(flows[int((t[0])*2+1),15] +flows[int((t[0])*2),15])
    boundary_flow[7,11] = -(flows[int((t[0])*2+1),16] +flows[int((t[0])*2),16])
    
    #The following loop adjusts for the positive flow direction which is set in...
    #...contam to avoid negative flows i.e if the boundary flow value, i to j is...
    #... negative, then this flow should be for flow j to i
    for i in range(n):
        for j in range(n):
            if geometry[i,j] > 0:
                if boundary_flow[i,j] < 0:
                    boundary_flow[j,i] = - boundary_flow[i,j]
                    boundary_flow[i,j] = 0
                

    
    print("flow from zone i to zone k boundary_flow =" + str(boundary_flow))
    ###########################################################################
    
    
    return boundary_flow