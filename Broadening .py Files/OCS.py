# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2
n_H2, gamma_CO2 and n_CO2 of OCS transitions in the HITRAN database as described in:
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

#--------------read OCS HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He and H2-broadening and temperature dependence of OCS:')

savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close

total = ascii.read(table, Reader = ascii.FixedWidthNoHeader,
                         names = ('parline', 'branch', 'J'),
                         col_starts = (0, 117, 118),
                         col_ends = (159, 117, 120),
                         )# This work assumes the HITRAN .par format is the input data
                         
# unused columns in next steps
Parline = np.array(total['parline'])

# columns used in next steps
Branch = np.array(total['branch'])
J = np.array(total['J'])

#-----------------calcuating |m| of OCS lines for H2----------------------------------
# *Note that m in this work stands for |m| which is related to the lower J rotational quantum number as follows:
# P branch: m = -J" (However in this work we are using |m| so for P branches this is just J")
# Q branch: m = J"
# R branch: m = J" + 1

m_H2 = []
for i in range(len(Branch)):
    if Branch[i]=='R':
        m_H2.append(J[i]+1)
    elif Branch[i]=='P':
        m_H2.append(J[i])
    elif Branch[i]=='Q':
        m_H2.append(J[i])
    if m_H2[i]<=1:
        m_H2[i]=2
    elif m_H2[i]==61:
        m_H2[i]=57
    else:
        pass
       
#-----------------calcuating |m| of OCS lines for He----------------------------------
# *Note that m in this work stands for |m| which is related to the lower J rotational quantum number as follows:
# P branch: m = -J" (However in this work we are using |m| so for P branches this is just J")
# Q branch: m = J"
# R branch: m = J" + 1

m_He = []
for i in range(len(Branch)):
    if Branch[i]=='R':
        m_He.append(J[i]+1)
    elif Branch[i]=='P':
        m_He.append(J[i])
    elif Branch[i]=='Q':
        m_He.append(J[i])
    else:
        pass

#-----------------define function for gH2---------------------------------------

def gH2(x):
    a0 = -8.02672
    a1 = 4.87015
    a2 = 2.44905
    a3 = -0.04140
    b1 = -9.36773
    b2 = 25.58158
    b3 = -0.34727
    b4 = -0.00113
    
    ggH2 = (a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4)
    return ggH2 # x in this calculation stands for |m|

def err_gH2(x):
    if x<12:
        err = 5
    elif x>=12 and x<30:
        err = 4
    elif x>=30:
        err = 3
    return err
    
#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = -4.48798
    a1 = 6.50867
    a2 = 5.60066
    a3 = 1.36104
    b1 = 3.86063
    b2 = 87.3008
    b3 = 15.66005
    b4 = 0.03454

    ggHe = (a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4)
    return ggHe

def err_gHe(x):
    if x<70:
        err = 5
    elif x>=70 and x<90:
        err = 4
    elif x>=90:
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
for i in range(len(m_He)):
    xx = m_H2[i]
    x = m_He[i]
    gamma_He.append(float(gHe(x))) # He broadening
    err_He.append(str(err_gHe(x))) # He uncertainty code
    ref_He.append("1427") # He Broadening Data References: Tan et al. 2022
    n_He.append("0.75") # He Temperature Dependence
    err_n_He.append("3") # He Temperature Dependence uncertainty code
    ref_n_He.append("951") # He Temperature Dependence reference: OCS-He temperature dependence values for all transitions set to 0.75 due to lack of data
    gamma_H2.append(float(gH2(xx))) # H2 broadening
    err_H2.append(str(err_gH2(x))) # H2 uncertainty code
    ref_H2.append("1512") # H2 Broadening Data Reference: Tan et al. 2022 Padé Approximation fit to data from Broquier et al. 1986 https://doi.org/10.1063/1.450421
    n_H2.append("0.75") # H2 Temperature Dependence
    err_n_H2.append("3") # H2 Temperature Dependence uncertainty code
    ref_n_H2.append("992") # H2 Temperature Dependence reference: Default value of 0.75 for OCS-H2 temperature dependence exponents
        
#------------create new HITRAN format with H2 and He broadening and temperature dependence for OCS--------
with open(savepath,'w') as out:
    for i in range(len(m_He)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s \n" %(Parline[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2"')