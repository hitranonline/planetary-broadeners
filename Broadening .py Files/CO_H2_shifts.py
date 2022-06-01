# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating delta_H2
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

readpath = input('input HITRAN 160 .par file to do the calculation for H2-shifts of CO:')

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
a_H2 = [0.345] #VP

for i in range(len(Branch)):
    multipliers_H2 = a_H2*(v1_f-v1_i)
    
#-------------Function for generating shift values for CO broadened by H2 -----------------------------
def dH2(x, y, z):
    alph1_rot = 0.06963
    alph2_rot = -0.243263
    alph3_rot = 0.173377
    beta2_rot = 0.002443
    beta3_rot = 0.00350517
    alph1_vib = -0.00628
    alph2_vib = -0.00223
    alph3_vib = 0.001072
    beta2_vib = 1.15326
    beta3_vib = 0.18625
    ddH2 = (y*(alph1_rot+alph2_rot*np.exp(-x*beta2_rot)+alph3_rot*np.exp(-x*beta3_rot))+
            (z*(alph1_vib+alph2_vib*np.exp(-x*beta2_vib)+alph3_vib*np.exp(-x*beta3_vib))))
    return ddH2# x in this calculation stands for |m|, y stands for inx values, and z are the multiplier values

#--------------Fill empty lists with calculated broadening-------------------------------

H2_shifts = []
err_H2 = []
ref_H2 = []

for i in range(len(Branch)):
    x = ms[i]
    y = inx[i]
    z = multipliers_H2[i]
    H2_shifts.append(float(dH2(x, y, z)))# H2 pressure-induced line shifts
    err_H2.append("3")                   # He shifts uncertainty code
    ref_H2.append("1345")# He shifts data references: Described in Tan et al. 2022 CO-H2 broadening were obtained by fitting the Padé approximation on data from Malathy Devi et al. 2004 https://dx.doi.org/10.1016/j.jms.2004.05.006 and Sung and Varanasi 2004 https://dx.doi.org/10.1016/S0022-4073(03)00202-4
    
#------------create new HITRAN data file with H2 shifts for CO--------
with open(savepath,'w') as out:
    for i in range(len(ms)):
        out.write("%160s, %9.6f, %3s, %3s \n" %( Parline[i], H2_shifts[i], err_H2[i], ref_H2[i]))
    else:
        print('end for calculation: output "160.par + delta_H2" ')