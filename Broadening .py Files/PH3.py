# -*- coding: utf-8 -*-
'''
This is a Python code provided for calculating gamma_He, n_He, gamma_H2 
and n_H2 of PH3 transitions in the HITRAN database as described in:
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

#--------------read PH3 HITRAN data-------------------------------

readpath = input('input HITRAN 160 .par file to do the calculation for He- and H2-broadening and temperature dependence of PH3:')
savepath = input('output file name:')

file1 = open(readpath,'r') 
table = file1.read()
file1.close


total = ascii.read(table, Reader = ascii.FixedWidthNoHeader,
                    names = ('parline','J_upp','Ka_upp','J_low'),
                         col_starts = (0,  98, 101, 113),
                         col_ends = (159, 100, 102, 114),
                         )# This work assumes the HITRAN .par format is the input data
                           
# unused columns in next steps
Parline = np.array(total['parline'])  

# columns used in next steps
J_low = np.array(total['J_low'])    
Ka_upp = np.array(total['Ka_upp'])  
J_upp = np.array(total['J_upp'])  

#-----------------calcuating mjval for PH3 lines----------------------------------
# mjval stands for |m| which is related to the lower J" rotational quantum number as follows:
# J' = J" - 1 then m = -J" (However in this work we are using |m| so this is just J")
# J' = J" then m = J"
# J' = J" + 1 then m = J" + 1
# for these specific |m| values associated with PH3, any |m|>22 is set to a value of 22 

mjval = []
for i in range(len(J_upp)):
    if J_upp[i]==J_low[i] + 1:
        mjval.append(J_low[i]+1)
    elif J_upp[i]==J_low[i] - 1: 
        mjval.append(J_low[i])
    elif J_upp[i]==J_low[i]:
        mjval.append(J_low[i])
    if mjval[i]>22:
        mjval[i]=22
    else:
        pass

#-----------------calcuating jvalhe for PH3 lines----------------------------------
# *Note that a cutoff has been applied to J" values here for J" >= 14 then J is set to 14

jvalhe = []
for i in range(len(J_upp)):
    if J_low[i]<14:
        jvalhe.append(J_low[i])
    elif J_low[i]>=14:
        jvalhe.append(14)
    else:
        pass

#-----------------calcuating jvalh2 for PH3 lines----------------------------------
# *Note that a cutoff has been applied to J" values here for J" >= 11 then J is set to 11

jvalh2 = []
for i in range(len(J_upp)):
    if J_low[i]<11:
        jvalh2.append(J_low[i])
    elif J_low[i]>=11:
        jvalh2.append(11)
    else:
        pass

#-----------------calcuating kauppval for PH3 lines----------------------------------
# *Note that a cutoff has been applied to Ka values here for Ka > 22 then Ka is set to 22

kauppval = []
for i in range(len(J_upp)):
    if Ka_upp[i]<=22:
        kauppval.append(Ka_upp[i])
    elif Ka_upp[i]>22:
        kauppval.append(22)
    else:
        pass

#-----------------define function for gH2---------------------------------------

def gH2(mlo, Ka_upp):
    a0 =  1.134E-01
    a1 = -1.658E-03
    a2 = -1.880E-03
    a3 = -1.956E-05
    a4 = -7.558E-04
    a5 =  7.189E-04
    a6 =  1.643E-06
    a7 = -1.943E-05
    a8 = -3.443E-05
    a9 =  5.511E-05
    ggH2 = a0 + (a1*mlo) + (a2*Ka_upp) + (a3*mlo*mlo) + (a4*Ka_upp*Ka_upp) + (a5*mlo*Ka_upp) + \
    (a6*mlo*mlo*mlo) + (a7*Ka_upp*Ka_upp*Ka_upp) + (a8*mlo*mlo*Ka_upp)  + (a9*mlo*Ka_upp*Ka_upp) 
    return ggH2# mlo in this calculation stands for |m| and Ka_upp are the calculated kauppval
    
#-----------------define function for nH2---------------------------------------
def nH2(jvalh2):
    nnH2 = -0.0103*jvalh2 + 0.7247
    return nnH2
    
#-----------------define function for gHe---------------------------------------
def gHe(jvalhe):
    ggHe = -0.00104*jvalhe + 0.05915
    return ggHe
    
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
for i in range(len(Ka_upp)):
    jhe = jvalhe[i]
    jh2 = jvalh2[i]
    mj = mjval[i]
    Ka = kauppval[i]
    gamma_He.append(float(gHe(jhe)))  # He broadening
    err_He.append('3')                # He uncertainty code
    ref_He.append("1313") # He Broadening Data References: Tan et al. 2022 linear fit to data from Pickett et al. 1981 https://doi.org/10.1016/0022-4073(81)90113-8 and Sergent-Rozey et al. 1988 https://doi.org/10.1016/0022-2852(88)90107-5 and Salem et al. 2005 https://doi.org/10.1016/j.jms.2005.04.014
    n_He.append('0.3030')             # He temperature dependence
    err_n_He.append('1')              # He temperature dependence uncertainty code
    ref_n_He.append("1314") # He temperature dependence Data Reference: Levy et al. 1994 https://doi.org/10.1006/jmsp.1994.1168
    gamma_H2.append((float(gH2(mj, Ka)))) # H2 broadening
    err_H2.append('4')                # H2 uncertainty code
    ref_H2.append("1307") # H2 Broadening Data References: Tan et al. 2022 polynomial fit to data from Bouanich et al. 2004 https://doi.org/10.1016/S0022-4073(03)00143-2 and Butler et al. 2006 https://doi.org/10.1016/j.jms.2006.04.021
    n_H2.append(float(nH2(jh2)))      # H2 temperature dependence
    err_n_H2.append('3')              # H2 temperature dependence uncertainty code
    ref_n_H2.append("1309") # H2 temperature dependence Data Reference: Described in Tan et al. 2022, data from Salem et al. 2004 https://doi.org/10.1016/j.jms.2004.06.015 are linearly fit 

#------------create new HITRAN data file with He and H2 broadening and temperature dependence for PH3--------
with open(savepath,'w') as out:
    for i in range(len(Ka_upp)):
        out.write("%160s, %8.4f, %3s, %3s, %3s, %3s, %3s, %8.4f, %3s, %3s, %8.4f, %3s, %3s \n" 
        %(Parline[i], gamma_He[i], err_He[i], ref_He[i], n_He[i], err_n_He[i], ref_n_He[i],
         gamma_H2[i], err_H2[i], ref_H2[i], n_H2[i], err_n_H2[i], ref_n_H2[i]))
    else:
        print('end for calculation: output "160.par + gamma_He + n_He + gamma_H2 + n_H2" ')
