# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2,
n_H2, gamma_CO2 and n_CO2 of CO2 transitions in the HITRAN database as described in:
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

#--------------read CO2 HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He-, H2- and CO2-broadening and temperature dependence of CO2:')

savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close


total = ascii.read(table, Reader = ascii.FixedWidthNoHeader,
                         names = ('parline', 'branch', 'J'),
                         col_starts = (0, 117, 118),
                         col_ends = (159, 117, 120),
                         ) # This work assumes the HITRAN .par format is the input data

# unused columns in next steps
Parline = np.array(total['parline'])

# column used in next steps
Branch = np.array(total['branch'])
J = np.array(total['J'])

#-----------------calcuating |m| for CO2 lines----------------------------------
# *Note that m stands for |m| which is related to the lower J rotational quantum number as follows:
# P branch: m = -J" (However in this work we are using |m| so for P branches this is just J")
# Q branch: m = J"
# R branch: m = J" + 1

m = []
for i in range(len(Branch)):
    if Branch[i]=='R':
        m.append(J[i]+1)
    elif Branch[i]=='P': 
        m.append(J[i])
    elif Branch[i]=='Q':
        m.append(J[i])
    else:
        pass

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = 0.07206
    a1 = -0.02269
    a2 = 0.10172
    a3 = 0.01168
    b1 = -0.3246
    b2 = 1.43332
    b3 = 0.21907
    b4 = 8.94019E-5

    ggHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggHe # x in this calculation stands for |m|

def err_gHe(x):
    if x<40:
        err = 5
    elif x>=40:
        err = 4
    return err    
    
#-----------------define function for gH2---------------------------------------
def gH2(x):
    a0 = 0.30051
    a1 = 1.99925
    a2 = -0.02836
    a3 = 6.34937E-4
    b1 = 14.15000
    b2 = 0.02731
    b3 = -9.28600E-4
    b4 = 6.25400E-5 

    ggH2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggH2 

def err_gH2(x):
    if x<40:
        err = 5
    elif x>=40:
        err = 4
    return err    

#-----------------define function for gCO2---------------------------------------
def gCO2(x):
    a0 = 1.312E-1
    a1 = 1.320E-2
    a2 = -3.851E-4
    a3 = 4.312E-6
    b1 = 1.396E-1
    b2 = -3.00E-3
    b3 = 2.635E-5
    b4 = 1.954E-7
    
    ggCO2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggCO2 

def err_gCO2(x):
    if x<40:
        err = 5
    elif x>=40:
        err = 4
    return err  

#-----------------define function for nCO2---------------------------------------
def nCO2(x):
    a0 = 7.926E-1 
    a1 = -5.339E-2
    a2 = 5.805E-5
    a3 = 6.916E-5
    b1 = -4.258E-2
    b2 = -2.530E-3
    b3 = 1.644E-4
    b4 = -1.619E-7
    
    nnCO2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return nnCO2 

def err_nCO2(x):
    if x<90:
        err = 5
    elif x>=90:
        err = 4
    return err  
  
#--------------Fill empty lists with calculated broadening-------------------------------

gamma_He = []
n_He = []
err_n_He = []
ref_n_He = []
err_He = []
ref_He = []

gamma_H2 = []
n_H2 = []
err_n_H2 = []
ref_n_H2 = []
err_H2 = []
ref_H2 = []

gamma_CO2 = []
n_CO2 = []
err_n_CO2 = []
ref_n_CO2 = []
err_CO2 = []
ref_CO2 = []

#--The reference numbers below correspond to "global reference IDs" in the HITRAN database. The mapping is also provided here in the code."
for i in range(len(m)):
    xx = m[i]
    gamma_He.append(float(gHe(xx)))   # He broadening
    n_He.append('0.3000')             # He temperature dependence
    err_n_He.append('4')              # He temperature dependence uncertainty code
    ref_n_He.append('1501')# He temperature dependence Data Reference: Nakamichi et al. 2006 https://doi.org/10.1039/B511772K
    err_He.append(str(err_gHe(xx)))   # He uncertainty code
    ref_He.append("1511")  # He Broadening Data Reference: Tan et al. 2022 Padé fit to data from Nakamichi et al. 2006 https://doi.org/10.1039/B511772K
    gamma_H2.append(float(gH2(xx)))   # H2 broadening
    n_H2.append('0.5800')             # H2 temperature dependence
    err_n_H2.append('4')              # H2 temperature dependence uncertainty code
    ref_n_H2.append('1499')# H2 temperature dependence Data Reference: Hanson and Whitty 2014 https://doi.org/10.2172/1222583
    err_H2.append(str(err_gH2(xx)))   # H2 uncertainty code
    ref_H2.append("1509")  # H2 Broadening Data Reference: Tan et al. 2022 average value of H2/air; H2 data from Padmanabhan et al. 2014 https://doi.org/10.1016/j.jqsrt.2013.07.016 
    gamma_CO2.append(float(gCO2(xx))) # CO2 broadening
    n_CO2.append(float(nCO2(xx)))     # CO2 temperature dependence
    err_n_CO2.append(str(err_nCO2(xx)))# CO2 temperature dependence uncertainty code
    ref_n_CO2.append('1273')# CO2 temperature dependence Data Reference: Hashemi et al. 2020 https://doi.org/10.1016/j.jqsrt.2020.107283
    err_CO2.append(str(err_gCO2(xx))) # CO2 uncertainty code
    ref_CO2.append("1359")# CO2 Broadening Data Reference: Tan et al. 2022 Padé fit to data from Hashemi et al. 2013 https://dx.doi.org/10.1139/cjp-2013-0051 and Predoi-Cross et al. 2007 https://doi.org/10.1016/j.jms.2007.07.004

#------------create new HITRAN data file with He, H2 and CO2 broadening and temperature dependence for CO2--------
with open(savepath,'w') as out:
    for i in range(len(m)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s \n" %( Parline[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i],
        gamma_CO2[i], err_CO2[i], ref_CO2[i], n_CO2[i], err_n_CO2[i], ref_n_CO2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2 + gamma_CO2 + n_CO2" ')