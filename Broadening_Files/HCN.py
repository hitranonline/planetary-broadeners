# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2
and n_H2 of HCN transitions in the HITRAN database as described in:
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

#--------------read HCN HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He- and H2-broadening and temperature dependence of HCN:')

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
J = np.array(total['J'])
Branch = np.array(total['branch'])

#-----------------calcuating |m| of HCN lines for H2----------------------------------
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
    if m[i]==0:
        m[i]=1
    else:
        pass
        
#-----------------calcuating |m| of HCN lines for He----------------------------------
# m stands for |m| which is related to the lower J rotational quantum number as follows:
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
    if m_He[i]==0:
        m_He[i]=2
    elif m_He[i]==1:
        m_He[i]=2
    elif m_He[i]>16:
        m_He[i]=16
    else:
        pass

#-----------------define function for gHe---------------------------------------
def gHe(x):
    a0 = -9.807238
    a1 =  9.53324
    a2 =  0.50085
    a3 =  0.31568
    b1 = 133.30485
    b2 = -13.64947
    b3 = 13.12444
    b4 = -0.22919

    ggHe = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggHe

def err_gHe(x):
    if x<16:
        err = 5
    elif x>=16 and x<30:
        err = 4
    elif x>=30:
        err = 3
    return err    
    
#-----------------define function for gH2---------------------------------------
def gH2(x):
    a0 = -2.91752
    a1 =  3.99556
    a2 = -0.42136
    a3 =  1.27061
    b1 = -4.30304
    b2 = 12.16122
    b3 =  7.01587
    b4 =  0.18831
    
    ggH2 = ((a0+a1*x+a2*x**2+a3*x**3)/(1+b1*x+b2*x**2+b3*x**3+b4*x**4))
    return ggH2

def err_gH2(x):
    if x<31:
        err = 6
    elif x>=31:
        err = 5
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
for i in range(len(m)):
    xx = m[i]
    x = m_He[i]
    gamma_He.append(float(gHe(x))) # He broadening
    err_He.append(str(err_gHe(xx))) # He uncertainty code
    ref_He.append("1496") # He Broadening Data References: Tan et al. 2022 Padé fit to the data provided by Rohart et al. 2007 https://doi.org/10.1016/j.jms.2007.09.009 and D'Eu et al. 2002 https://doi.org/10.1006/jmsp.2002.8520
    n_He.append("0.71") # He Temperature Dependence
    err_n_He.append("3") # He Temperature Dependence uncertainty code
    ref_n_He.append("1494") # He Temperature Dependence reference: Rohart et al. 2007 https://doi.org/10.1016/j.jms.2007.09.009
    gamma_H2.append(float(gH2(xx))) # H2 broadening
    err_H2.append(str(err_gH2(xx))) # H2 uncertainty code
    ref_H2.append("1498") # H2 Broadening Data References: Tan et al. 2022 Padé fit to the data provided by Charròn et al. 1980 https://doi.org/10.1063/1.440354 and Lemaire et al. 1996 https://doi.org/10.1006/jmsp.1996.0115 and Landrain et al. 1997 https://doi.org/10.1006/jmsp.1996.7223 and Mehrotra et al. 1985 https://doi.org/10.1016/0301-0104(85)85053-9 and Rohart et al. 2007 https://doi.org/10.1016/j.jms.2007.09.009
    n_H2.append("0.90") # H2 Temperature Dependence
    err_n_H2.append("3") # H2 Temperature Dependence uncertainty code
    ref_n_H2.append("1497") # H2 Temperature Dependence reference: Tan et al. 2022 averaged HCN H2-temperature dependence measurements are provided by Charròn et al. 1980 https://doi.org/10.1063/1.440354 and Rohart et al. 2007 https://doi.org/10.1016/j.jms.2007.09.009

#------------create new HITRAN format with H2 and He broadening and temperature dependence for HCN--------
with open(savepath,'w') as out:
    for i in range(len(m)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %3s, %3s, %3s \n" %(Parline[i], 
        gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i], 
        gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2" ')