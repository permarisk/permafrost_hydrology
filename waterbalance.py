#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 09:04:24 2022

@author: mlanger
"""

import numpy as np

def SaturationVapourPressure(T):
    #saturation vapour pressure [Pa] - August-Roche-Magnus equation
    e_sat = 100 * 6.1094 * np.exp(17.625*T/(T+243.04))
    return e_sat #[Pa]
       
def SaturationVapourPressureDerivative(T):
    #derivative of the saturation vapour pressure [Pa/K] - August-Roche-Magnus equation
    de_sat_dT = 2.61701e6 * np.exp(17.625*T/(T+243.04)) * (T+243.04)**(-2)
    return de_sat_dT #[Pa/K]

def WaterVapourAerodynamicResistance(z,Uz,hc):
    z_w = z #wind speed measurement hight [m]
    k = 0.41 #von Karman constant [-]
    
    d = 2./3. * hc #displacment hight estmated based on vegetation hight [m]
    z_0m = 0.123 * hc #roughness length momentum [m]
    z_0v = 0.1 * z_0m #roughness length surface [m]
    
    r_av = np.log((z_w-d)/z_0m) * np.log((z_w-d)/z_0v) * (k**2 * Uz)**(-1) 
    
    return r_av
    
def PenmanMonteith(e0, ea, dedT, R_net, G, z, Uz, hc, r_s):
    
    gamma = 0.066 #psychrometric constant [kPa °C -1] (at sea level  101.325 kPa!)
    rho_air = 1.2250 #density of air at seal level ( 101.325 kPa) and 15°C [kg/m³]
    capacity_air = 1.013e-3 #heat capacity of air [MJ kg -1 °C -1]
        
    r_av = WaterVapourAerodynamicResistance(z, Uz, hc)
    
    LET1 = dedT * (R_net - G) + (86400 * rho_air * capacity_air * (e0 - ea)) / r_av
    LET = LET1 / (dedT + gamma * (1. + r_s/r_av) )
    
    return LET



#==============================================================================


#load meteorological observations
file = '/home/mlanger/Dokumente/Lehre/PermafrostHydrology_2022/Data/Boike_2013/datasets/Samoylov_2002-2011_meteorology.tab'
import pandas as pd
dataset = pd.read_csv(file,delimiter="\t")
dataset['DateTime'] = pd.to_datetime(dataset['Date/Time'])

dataset = dataset.set_index(['DateTime'])
dataset = dataset.loc['2011-1-1':'2011-8-25'] 

df_mean = dataset.groupby(by=pd.Grouper(freq='D')).mean()
df_min = dataset.groupby(by=pd.Grouper(freq='D')).min()
df_max = dataset.groupby(by=pd.Grouper(freq='D')).max()
 
meteorology_mean = df_mean   
meteorology_min = df_min   
meteorology_max = df_max   

#Note net radiation and groud heat flux must be proved as [MJ m^2 / day] 

#as a initial guess we can approximate the ground heat flux to be 2% of the net radiation 
G_av = Rnet_av*0.02 
T_av = meteorology_mean['T2 [°C]'] 
RH_av = meteorology_mean['RH [%]'] / 100.





