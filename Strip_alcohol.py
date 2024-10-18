import sys
import os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Add the directory containing the package to the system path
sys.path.append('/Users/samdenhartog/Documents/Documenten - MacBook Pro van Sam/mep')


# Correct import path
from biochemical_engineering_relations.mass_transfer.kLa import *




def Henry(T, compound='Ethanol'): # https://doi.org/10.1016/j.atmosenv.2006.06.024
    if compound.lower() == 'ethanol':
        ln_Kh = -(15.87) + 6274.0/T
        Kh = np.exp(ln_Kh) * 1000   # mol/(m3*atm)
    return Kh

# Column properties
H_L = 1 # m
D_L = 0.19 # m
r_L = 0.19/2
A = np.pi * r_L**2 # m2
V_L = A*H_L # m3

# Experimental conditions
T_C = 20 # degrees C
T_K = T_C + 273 # K
P_tot = 101325 # Pa (1 atm = 101325 Pa)
Q_G = 0.003 # m3/s = 180 L/min
u_G = Q_G/A # m/s
C_max = 21 # g/L article states 20.7 g/L https://www.sciencedirect.com/science/article/pii/S1385894718307903

eps_G = 0.4 # gas holdup (estimation by Rik)
V_G = (eps_G*V_L)/(1-eps_G)

# Use Empirical relation: Deckwer1, Deckwer2, Jackson_Shen, Heijnen_VanRiet
k_H = Henry(T_K)/101325 # mol/(m3*Pa)

kLa = Empirical.Deckwer1(u_G) # s-1

# Constants
R = 8.314 # J/(mol*K)
MW_eth = 46.068 # g/mol

F_G = Q_G*P_tot/(R*T_K) # mol/s
C_begin = C_max*1000/MW_eth # mol/m3
print('Start')

def alcohol_conc(t, var):
    F = np.empty(2)

    C_alc = var[0]
    gam = var[1]
    #C_gas = 1/k_H * C_alc
    #r = -kLa * (C_alc-(P_vap*101325/(R*T_K)))#mol/m3/s
    #r = -kLa * (C_alc-(P_vap/k_H))#mol/m3/s
    CL_star = k_H * gam * P_tot
    r = kLa*(CL_star-C_alc)
    F[0] = r # mol/(m3*s) --> dC_L/dt
    F[1] = (-F_G*gam - r*V_L)*R*T_K/(V_G*P_tot) # mol_eth/(mol_gas*s) --> dy_eth/dt
    
    return F
    

total_time = 8 * 3600 # s


print(kLa*3600)
gam_0 = 0
sol = solve_ivp(alcohol_conc, (0, total_time), y0=[C_begin, gam_0], method='LSODA')

fig, ax1 = plt.subplots()

color='tab:red'
ax1.plot(sol.t/3600, sol.y[0]*MW_eth/1000, color=color)
ax1.set_ylabel('Concentration (g/L)', color=color)
ax1.set_xlabel('Time (h)')
ax1.tick_params(axis='y', labelcolor=color)

color='tab:blue'
ax2 = ax1.twinx()
ax2.plot(sol.t/3600, sol.y[1], color='tab:blue')
ax2.set_ylabel('Fraction of ethanol in offgas', color=color)
ax2.tick_params(axis='y', labelcolor=color)
plt.show()