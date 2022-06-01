# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2, n_H2, 
gamma_CO2 and n_CO2 of CO transitions in the HITRAN database as described in:
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

readpath = input('input HITRAN 160 .par file to do the calculation for He-, H2- and CO2-broadening and temperature dependence of CO:')

savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close

total = ascii.read(table, Reader=ascii.FixedWidthNoHeader,
                         names=('hitpar','Br','J_'),
                         col_starts=(0,117,118),
                         col_ends = (159,117,120),
                         )

# unused columns in next steps
Parline= np.array(total['hitpar'])

# column used in next steps
Branch= np.array(total['Br'])
J= np.array(total['J_'])

#-----------------calcuating |m| for CO lines----------------------------------
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

#-----------------define function for gH2---------------------------------------
def gH2(x):
    a0 = 0.08228
    a1 = -0.07411
    a2 = 0.10795
    a3 = 0.00211
    b1 = -1.0
    b2 = 1.53458
    b3 = 0.03054
    b4 = 6.9468E-5

    ggH2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggH2# x in this calculation stands for |m|

def errgH2(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3

#-----------------define function for nH2---------------------------------------
def nH2(x):
    a0 = 0.64438
    a1 = 0.49261
    a2 = -0.0748
    a3 = 0.0032
    b1 = 0.69861
    b2 = -0.09569
    b3 = 0.003
    b4 = 5.7852E-5

    nnH2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return nnH2

def errnH2(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = 0.0809
    a1 = 0.3641
    a2 = -0.04025
    a3 = 0.00178
    b1 = 8.1769
    b2 = -0.9105
    b3 = 0.0397
    b4 = 2.556E-6

    ggHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggHe

def errgHe(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3
        
#-----------------define function for nHe---------------------------------------
def nHe(x):
    a0 = 0.5393
    a1 = 0.1286
    a2 = -0.0129
    a3 = 0.00175
    b1 = 0.3146
    b2 = -0.0417
    b3 = 0.00403
    b4 = -6.589E-6

    nnHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return nnHe

def errnHe(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3


#-----------------define function for gCO2---------------------------------------
def gCO2(x):
    a0 = 0.12106
    a1 = 0.05433
    a2 = -0.00851
    a3 = 6.90673E-4
    b1 = 0.63012
    b2 = -0.07902
    b3 = 0.006
    b4 = 1.703E-4

    ggCO2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggCO2

def errgCO2(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3

#-----------------define function for nCO2---------------------------------------
def nCO2(x):
    a0 = 0.70343
    a1 = -0.10857
    a2 = 0.00407
    a3 = 1.112E-4
    b1 = -0.14755
    b2 = 0.00528
    b3 = 1.3829E-4
    b4 = 1.4546E-6

    nnCO2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return nnCO2

def errnCO2(x):
    if x<=101:
        return 5
    elif x>101 and x<=121:
        return 4
    else:
        return 3
        
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
    n_He.append(float(nHe(xx)))       # He temperature dependence
    err_n_He.append(str(errnHe(xx)))  # He temperature dependence uncertainty code
    ref_n_He.append('1345')# He temperature dependence Data Reference: Described in Tan et al. 2022 For CO-He the data from Predoi-Cross et al. 2016 https://doi.org/10.1016/j.jqsrt.2016.08.007, Sinclair et al. 1998 https://doi.org/10.1006/jmsp.1998.7628, Luo et al. 2001 https://doi.org/10.1063/1.1383049, Thibault et al. 1992 http://dx.doi.org/10.1063/1.463865 were used
    err_He.append(str(errgHe(xx)))    # He uncertainty code
    ref_He.append("1345")  # He Broadening Data References: Described in Tan et al. 2022 For CO-He the data from Predoi-Cross et al. 2016 https://doi.org/10.1016/j.jqsrt.2016.08.007, Sinclair et al. 1998 https://doi.org/10.1006/jmsp.1998.7628, Luo et al. 2001 https://doi.org/10.1063/1.1383049, Thibault et al. 1992 http://dx.doi.org/10.1063/1.463865 were used
    gamma_H2.append(float(gH2(xx)))   # H2 broadening
    n_H2.append(float(nH2(xx)))       # H2 temperature dependence
    err_n_H2.append(str(errnH2(xx)))  # H2 temperature dependence uncertainty code
    ref_n_H2.append('1345')# H2 temperature dependence Data Reference: Described in Tan et al. 2022 CO-H2 broadening were obtained by fitting the Padé approximation on data from Malathy Devi et al. 2004 https://dx.doi.org/10.1016/j.jms.2004.05.006 and Sung and Varanasi 2004 https://dx.doi.org/10.1016/S0022-4073(03)00202-4
    err_H2.append(str(errgH2(xx)))    # H2 uncertainty code
    ref_H2.append("1345")  # H2 Broadening Data References: Described in Tan et al. 2022 CO-H2 broadening were obtained by fitting the Padé approximation on data from Malathy Devi et al. 2004 https://dx.doi.org/10.1016/j.jms.2004.05.006 and Sung and Varanasi 2004 https://dx.doi.org/10.1016/S0022-4073(03)00202-4
    gamma_CO2.append(float(gCO2(xx))) # CO2 broadening
    n_CO2.append(float(nCO2(xx)))     # CO2 temperature dependence
    err_n_CO2.append(str(errnCO2(xx)))# CO2 temperature dependence uncertainty code
    ref_n_CO2.append('1345')# CO2 temperature dependence Data Reference: Described in Tan et al. 2022 For the CO-CO2 system, the measured data from Hashemi et al. 2016 http://dx.doi.org/10.1016/j.jms.2016.02.014 is used to extrapolate the broadening for all the transitions
    err_CO2.append(str(errgCO2(xx)))  # CO2 uncertainty code
    ref_CO2.append("1345")# CO2 Broadening Data Reference: Described in Tan et al. 2022 For the CO-CO2 system, the measured data from Hashemi et al. 2016 http://dx.doi.org/10.1016/j.jms.2016.02.014 is used to extrapolate the broadening for all the transitions
    
#------------create new HITRAN data file with He, H2 and CO2 broadening and temperature dependence for CO--------
with open(savepath,'w') as out:
    for i in range(len(m)):
        out.write("%160s, %8.4f, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s \n" %( Parline[i],
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i],
        gamma_CO2[i], err_CO2[i], ref_CO2[i], n_CO2[i], err_n_CO2[i], ref_n_CO2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2 + gamma_CO2 + n_CO2" ')