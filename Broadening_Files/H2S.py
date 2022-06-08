# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2
and n_H2 of H2S transitions in the HITRAN database as described in:
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

#--------------read H2S HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He- and H2-broadening and temperature dependence of H2S:')
savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close


total = ascii.read(table, Reader = ascii.FixedWidthNoHeader,
                         names = ('parline', 'J_', 'Ka_'),
                         col_starts = (0, 113, 116),
                         col_ends = (159, 114, 117),
                         )# This work assumes the HITRAN .par format is the input data

# unused columns in next steps
Parline = np.array(total['parline'])

# columns used in next steps
J = np.array(total['J_'])
Ka = np.array(total['Ka_'])

#-----------------calcuating J+0.2Ka of H2S lines for He----------------------------------
JKa = list()
for i in range(len(J)):
    xx = J[i] + 0.2 * Ka[i]
    JKa.append(xx)
    
#-----------------calcuating J+0.2Ka of H2S lines for H2----------------------------------
JKa_H2 = list()
for i in range(len(J)):
    xxx = J[i] + 0.2 * Ka[i]
    JKa_H2.append(xxx)
    if JKa_H2[i]<=1.2:
        JKa_H2[i]=2
    else:
        pass

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = 18.04211    
    a1 = 13.10827 
    a2 = -2.96011    
    a3 = 0.70801 
    b1 = 405.14936   
    b2 = -0.36953    
    b3 = -4.27884   
    b4 = 1.77897   

    if x>1 and x<30:
        ggHe = (a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4)
    elif x <=1:
        ggHe = (a0+a1*1+a2*1**2+a3*1**3)/(1+b1*1+b2*1**2+b3*1**3+b4*1**4) 
    elif x>=30:
        ggHe = (a0+a1*30+a2*30**2+a3*30**3)/(1+b1*30+b2*30**2+b3*30**3+b4*30**4)
    else:
        print ("error in calculation")    
    return ggHe # x in these calculations stands for J+0.2Ka

def err_gHe(x):
    if x<18 :
        err = 4
    elif x>=18 and x<30:
        err = 3
    else:
        err =2    
    return err    
    
#-----------------define function for gH2---------------------------------------

def gH2(x):
    a0 = 0.01908 
    a1 = 1.25017 
    a2 = -1.52728
    a3 = 0.93939 
    b1 = -1.89026
    b2 = -4.80047
    b3 = 6.22612 
    b4 = 0.81255 
    ggH2 = (a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4)
    return ggH2

def err_gH2(x):
    if x<12:
        err = 5
    elif x>=12 and x<30:
        err = 4
    elif x>=30:
        err = 3
    return err
    
#--------------Fill empty lists with calculated broadening-------------------------------

gamma_He = []
err_He = []
ref_He = []
n_He = []
err_n_He = []
ref_n_He = []

gamma_H2 = []
err_H2 = []
ref_H2 = []
n_H2 = []
err_n_H2 = []
ref_n_H2 = []

#--The reference numbers below correspond to "global reference IDs" in the HITRAN database. The mapping is also provided here in the code."
for i in range(len(J)):
    jka = JKa[i]
    jka_h2 = JKa_H2[i]
    gamma_He.append(float(gHe(jka)))  # He broadening
    err_He.append(str(err_gHe(jka)))  # He uncertainty code
    ref_He.append("1427") # He Broadening Data Reference: Tan et al. 2022
    n_He.append("0.46")# He Temperature Dependence
    err_n_He.append("4")# He Temperature Dependence uncertainty code
    ref_n_He.append("1514")# He Temperature Dependence reference: Tan et al. 2022 He-H2S temperature dependence values were calculated using the first equation under the Results section in Flatin et al. 1994 https://dx.doi.org/10.1006/jmsp.1994.1086 by using their broadening values.
    gamma_H2.append(float(gH2(jka_h2)))  # H2 broadening
    err_H2.append(str(err_gH2(jka)))  # H2 uncertainty code
    ref_H2.append("1427") # H2 Broadening Data Reference: Tan et al. 2022
    n_H2.append("0.70")# H2 Temperature Dependence
    err_n_H2.append("4")# H2 Temperature Dependence uncertainty code
    ref_n_H2.append("1514")# H2 Temperature Dependence reference: Tan et al. 2022 H2-H2S temperature dependence values were calculated using the first equation under the Results section in Flatin et al. 1994 https://dx.doi.org/10.1006/jmsp.1994.1086 by using their broadening values.

#------------create new HITRAN format with H2 and He broadening and temperature dependence for H2S--------
with open(savepath,'w') as out:
    for i in range(len(J)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s \n" %( Parline[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2" ')