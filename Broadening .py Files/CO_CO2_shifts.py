# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating delta_CO2
of CO transitions in the HITRAN database as described in:
Tan et al, "H$_2$, He, and CO$_2$ line-broadening coefficients, and 
temperature-dependence exponents for the HITRAN database. Part II: 
CO2, N2O, CO, SO2, OH, OCS, H2CO, HCN, PH3, H2S and GeH4", 
Astrophysical Journal Supplement Series, 2022
Sample input and output files are included.
Note that the program can be easily modified for any other file formats.
'''

from astropy.io import ascii
from astropy.table import Table, Column
import numpy as np
import sys
import os

#--------------read CO HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for CO2-shifts of CO:')

savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close

total = ascii.read(table, Reader=ascii.FixedWidthNoHeader,
                         names=('hitpar','Br','J','v_f','v_i'),
                         col_starts=(0,  117,118, 79, 96),
                         col_ends = (159,117,120, 86, 99),
                         )

# unused columns in next steps
Parline= np.array(total['hitpar'])

# column used in next steps
Branch= np.array(total['Br'])
J= np.array(total['J'])
v1_f= np.array(total['v_f'])
v1_i= np.array(total['v_i'])

#-----------------calcuating |m| for CO lines----------------------------------
# *Note that m stands for |m| which is related to the lower J rotational quantum number as follows:
# P branch: m = -J" (However in this work we are using |m| so for P branches this is just J")
# Q branch: m = J"
# R branch: m = J" + 1

ms = []
inx = []
for i in range(len(Branch)):
    if Branch[i]=='R':
        ms.append(J[i]+1) 
        inx.append(-1)
    elif Branch[i]=='P':
        ms.append(J[i])
        inx.append(1)
    elif Branch[i]=='Q':
        ms.append(J[i])
        inx.append(0)
    else:
        pass

#-----------------calcuating multipliers for CO lines----------------------------------
a_CO2 = [0.5] #VP

for i in range(len(Branch)):
    multipliers_CO2 = a_CO2*(v1_f-v1_i)
    
#-------------Function for generating shift values for CO broadened by CO2 -----------------------------
def dCO2(x, y, z): 
    alph1_rot = 1.25396
    alph2_rot = -2.05688
    alph3_rot = 0.803285
    beta2_rot = 0.001053
    beta3_rot = 0.002796
    alph1_vib = 0.01503
    alph2_vib = 0.02691
    alph3_vib = -0.04405
    beta2_vib = 0.02746
    beta3_vib = 0.008576
    ddCO2 = (y*(alph1_rot+alph2_rot*np.exp(-x*beta2_rot)+alph3_rot*np.exp(-x*beta3_rot))+
            (z*(alph1_vib+alph2_vib*np.exp(-x*beta2_vib)+alph3_vib*np.exp(-x*beta3_vib))))
    return ddCO2# x in this calculation stands for |m|, y stands for inx values, and z are the multiplier values

#--------------Fill empty lists with calculated broadening-------------------------------

CO2_shifts = []
err_CO2 = []
ref_CO2 = []

for i in range(len(Branch)):
    x = ms[i]
    y = inx[i]
    z = multipliers_CO2[i]
    CO2_shifts.append(float(dCO2(x, y, z)))# CO2 pressure-induced line shifts
    err_CO2.append("3")                     # CO2 shifts uncertainty code
    ref_CO2.append("1345")# CO2 shifts data references: Described in Tan et al. 2022 For the CO-CO2 system, the measured data from Hashemi et al. 2016 http://dx.doi.org/10.1016/j.jms.2016.02.014 is used to extrapolate the broadening for all the transitions
    
#------------create new HITRAN data file with CO2 shifts for CO--------
with open(savepath,'w') as out:
    for i in range(len(ms)):
        out.write("%160s, %9.6f, %3s, %3s \n" %( Parline[i], CO2_shifts[i], err_CO2[i], ref_CO2[i]))
    else:
        print('end for calculation: output "160.par + delta_CO2" ')