# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He and n_He of 
N2O transitions in the HITRAN database as described in:
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

#--------------read N2O HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He-broadening and temperature dependence of N2O:')

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

#-----------------calcuating |m| for N2O lines----------------------------------
# m stands for |m| which is related to the lower J rotational quantum number as follows:
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
    if m[i]>40:
        m[i]=40
    else:
        pass

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = 29.92585
    a1 = 275.28681
    a2 = -21.0512
    a3 = 0.78324
    b1 = 4411.70782
    b2 = -370.09121
    b3 = 15.49987
    b4 = -0.03189

    ggHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggHe # x in this calculation stands for |m|

def err_gHe(x):
    if x<38:
        err = 5
    elif x>=38:
        err = 4
    return err    
    
#--------------Fill empty lists with calculated broadening-------------------------------

gamma_He = []
err_He = []
ref_He = []
n_He = []
ref_n_He = []
err_n_He = []

#--The reference numbers below correspond to "global reference IDs" in the HITRAN database. The mapping is also provided here in the code."
for i in range(len(m)):
    xx = m[i]
    gamma_He.append(float(gHe(xx))) # He broadening
    err_He.append(str(err_gHe(xx))) # He uncertainty code
    ref_He.append("1504") # He Broadening Data References: The He broadening data from Nakayama et al. 2007 https://doi.org/10.1016/j.chemphys.2007.03.001 and from Tasinato et al. 2010 https://doi.org/10.1063/1.3386385 were used to fit the Padé approximant in Tan et al. 2022
    n_He.append("0.30")   # He Temperature Dependence
    ref_n_He.append("1515")# He Temperature Dependence Reference: As stated in Tan et al. 2022, due to the lack of He-temperature dependence data for N2O, the He-temperature dependence value from Nakamichi et al. https://doi.org/10.1039/b511772k for CO2 lines is used.
    err_n_He.append("3")  # He Temperature Dependence Uncertainty Code

#------------create new HITRAN data file with He broadening for N2O--------
with open(savepath,'w') as out:
    for i in range(len(m)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s \n" %( Parline[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He" ')