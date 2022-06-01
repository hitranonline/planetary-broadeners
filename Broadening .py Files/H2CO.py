# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2 
and n_H2 of H2CO transitions in the HITRAN database as described in:
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

#--------------read H2CO HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He- and H2-broadening and temperature dependence of H2CO:')

savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close

total = ascii.read(table, Reader = ascii.FixedWidthNoHeader,
                         names = ('parline', 'J_', 'Ka_', 'air'),
                         col_starts = (0, 113, 116, 35),
                         col_ends = (159, 114, 117, 40),
                         ) # This work assumes the HITRAN .par format is the input data

# columns unused in next steps
Parline = np.array(total['parline'])

# columns used in next steps
J = np.array(total['J_'])
Ka = np.array(total['Ka_'])
Air_broadening = np.array(total['air'])

#-----------------calcuating J+0.2Ka for H2CO lines----------------------------------
JKa = list()
for i in range(len(J)):
    xx = J[i] + 0.2 * Ka[i]
    JKa.append(xx)
    if JKa[i]==0:
        JKa[i]=1
    else:
        pass

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = -24.09414
    a1 = 32.4839
    a2 = 2.97868
    a3 = 0.47408
    b1 = 4.07669
    b2 = 31.84113
    b3 = -3.37705
    b4 = 0.18356

    ggHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))*Air_broadening[i]
    return ggHe # This currently populates He broadening (also x is J+0.2Ka)
                # To get He/Air broadening values, remove the double () and *Air_broadening[i]

def err_gHe(x):
    if x<15:
        err = 5
    elif x>=15 and x<30:
        err = 4
    elif x>=30:
        err = 3
    return err    

#-----------------define function for gH2---------------------------------------
def gH2(x):
    a0 = 27.529045  
    a1 = -103.93252
    a2 = 26.695497  
    a3 = 1.630053   
    b1 = -80.069841 
    b2 = 23.497867  
    b3 = 1.010394   
    b4 = 0.005558   
    
    ggH2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))*Air_broadening[i]
    return ggH2 # This currently populates H2 broadening (also x is J+0.2Ka)
                # To get H2/Air broadening values, remove the double () and *Air_broadening[i]

def err_gH2(x):
    if x<16:
        err = 5
    elif x>=16 and x<30:
        err = 4
    elif x>=30:
        err = 3
    return err
  
 #--------------Fill empty lists with calculated broadening-------------------------------
ref_air = []

gamma_He = []
n_He = []
ref_n_He = []
err_n_He = []
err_He = []
ref_He = []

gamma_H2 = []
n_H2 = []
ref_n_H2 = []
err_n_H2 = []
err_H2 = []
ref_H2 = []

#--The reference numbers below correspond to "global reference IDs" in the HITRAN database. The mapping is also provided here in the code."
for i in range(len(J)):
    jka = JKa[i]
    gamma_He.append(float(gHe(jka))) # He broadening
    err_He.append(str(err_gHe(jka))) # He uncertainty code
    ref_He.append("1427") # He Broadening Data References: Tan et al. 2022
    n_He.append("0.75") # He Temperature Dependence
    ref_n_He.append("1436") # He Temperature Dependence Reference: Due to a lack of available measurements a default value of 0.75 for He-temperature dependence values have been assigned.
    err_n_He.append("3") # He Temperature Dependence uncertainty code
    gamma_H2.append(float(gH2(jka))) # H2 broadening
    err_H2.append(str(err_gH2(jka))) # H2 uncertainty code
    ref_H2.append("1427") # H2 Broadening Data References: Tan et al. 2022
    n_H2.append("0.75") # H2 Temperature Dependence
    ref_n_H2.append("1436") # H2 Temperature Dependence Reference: Due to a lack of available measurements a default value of 0.75 for H2-temperature dependence values have been assigned.
    err_n_H2.append("3") # H2 Temperature Dependence uncertainty code
    ref_air.append("825") # Air Broadening Data Reference: Jacquemart et al. 2010 https://doi.org/10.1016/j.jqsrt.2010.02.004
        

#------------create new HITRAN data file with H2 and He broadening and temperature dependence for H2CO--------
with open(savepath,'w') as out:
    for i in range(len(J)):
        out.write("%160s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s \n" %(Parline[i], ref_air[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i]))
    else:
        print('end for calculation: output "160.par + ref_air + gamma_He + n_He + gamma_H2 + n_H2" ')
