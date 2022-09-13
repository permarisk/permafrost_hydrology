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

def ThermalConductivity(Water, Mineral, Organic):
    
    #thermal conductivty estimated according to 
    ka = 0.025;       #air [Hillel(1982)]
    kw = 0.57;        #water [Hillel(1982)]
    ko = 0.25;        #organic [Hillel(1982)]
    km = 3.8;         #mineral [Hillel(1982)]
    ki = 2.2;         #ice [Hillel(1982)]

    air = 1.0 - Water - Mineral - Organic
    TC = (Water * kw**0.5 + Mineral * km**0.5 + Organic * ko**0.5 + air * ka**0.5)**2.0
    return TC
    
def ThawDepth(Tav,Theta_w, Theta_m, Theta_o):
    #estmates thaw depth based on Stefans equation
    L_sl = 334.0e3 #J/kg laten heat of fusion
    rho_w = 1.0e3 #[kg/m³] denstiy of water
    dt = 86400 #[s] seconds per day
    
    J = np.cumsum(Tav*(Tav>0.0)) * dt #[Ks] integrated thawing degree days 
    TC_t = ThermalConductivity(Theta_w, Theta_m, Theta_o) #[W/mK] thermal conductivity
        
    d = ((2. * TC_t)/(rho_w * L_sl * Theta_w) * np.abs(J))**0.5 #[m] thaw depth
    
    T_f = 0. #[°C]
    
    G = TC_t * (Tav - T_f) / d * dt / 1e6 #[MJ/day] ground heat flux per day
    G[np.isinf(G)] = 0.0 #set inf to zero if d == 0
    G[np.isnan(G)] = 0.0 #set nan to zero if d and Tav == 0
    
    return d, G

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

#load wind speed measurments
file = '/home/mlanger/Dokumente/Lehre/PermafrostHydrology_2022/Data/Boike_2013/datasets/Samoylov_2002-2011_wind_speed.tab'
import pandas as pd
dataset = pd.read_csv(file,delimiter="\t")
dataset['DateTime'] = pd.to_datetime(dataset['Date/Time'])

dataset = dataset.set_index(['DateTime'])
dataset = dataset.loc['2011-1-1':'2011-8-25'] 

windspeed_mean = dataset.groupby(by=pd.Grouper(freq='D')).mean()


#variables
Uz_av = windspeed_mean['ff [m/s]']
Rnet = np.array([0,200,600])#[W/m²s]

#Note net radiation and groud heat flux must be proved as [MJ m^2 / day] 
dt = 86400
Rnet_av = meteorology_mean['NET [W/m**2]']
Rnet_av = Rnet_av * dt / 1e6

G_av = Rnet_av*0.02 
T_av = meteorology_mean['T2 [°C]'] 
RH_av = meteorology_mean['RH [%]'] / 100.
ea = RH_av * SaturationVapourPressure(T_av)/1000.
dedT = SaturationVapourPressureDerivative(T_av)/1000.
T_min = meteorology_min['T2 [°C]']
T_max = meteorology_max['T2 [°C]']
e0 = 0.5 * (SaturationVapourPressure(T_min) + SaturationVapourPressure(T_max))/1000.

Theta_w = 0.6
Theta_m = 0.2
Theta_o = 0.05

[thaw_depth, G] = ThawDepth(T_av,Theta_w, Theta_m, Theta_o)
thaw_depth.plot()

#input parameter for Penman Monteith    
z = 3.0 #[m] wind speed mesurmeent hight above ground
hc = 0.10 #[m] crop (vegetation) hight above ground
rs = 200 #[-] surface resitance to evapotranspiration

#Penman Monteith equation 
LET = PenmanMonteith(e0, ea, dedT, Rnet_av, G_av, z, Uz_av, hc, rs) #[MJ/m²day]
L = 2.45 #[MJ kg -1]
ET = LET/L #[mm/day]

#Note that percipitation was measured as mm/h which must be intergrated over 24h 
P = meteorology_mean['Precip [mm/h]']*24.


#thaw_depth.plot()
ET.cumsum().plot()
(P.cumsum()+100.).plot()
