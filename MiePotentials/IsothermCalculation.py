import numpy as np
import sys
import matplotlib.pyplot as plt 
import matplotlib
#from itertools import cycle
from collections import OrderedDict

import OZ.Potentials as Potentials
import OZ.OZ_Functions as OZF
import OZ.PP_Functions  as PP_Functions
#import FigFuncs
#import HNC
#import OZ.RHNC as RHNC
from OZ.SinglePointCalculation import * 


if len(sys.argv) != 5:
    print "Usage: Identify_Tc.py m n eps_crit BRANCH"
    sys.exit()
else:
    m = int(sys.argv[1])
    n = int(sys.argv[2])
    eps_crit = float(sys.argv[3])
    branch = sys.argv[4]

    T_crit_estimate = 1.0/eps_crit


if branch == 'vap':
    rho_array = np.arange(0.01,0.80, 0.01)
elif branch == 'liq':
    rho_array = np.arange(0.80,0.01,-0.01)
else:
    print "Unknown branch. Has to be 'vap' or 'liq'"
    print branch
    sys.exit()


dr = 0.02
nr = 2048 
sig=1.0

r,k = OZF.Create_r_k(dr, nr)
sig= 1.0

Ur_ref = Potentials.WCAPotential(r,sig=sig, eps=1.0, m=m, n=n)
Ur = Potentials.LJPotential(r,sig=sig, eps=1.0, m=m, n=n)


#---------
cr_guess = None

PlotFreq=5 
i = 0

#FigSK = plt.figure()
#FigSK1 = FigSK.add_subplot(111)
#FigSK1.set_xlabel('k$\sigma$/2$\pi$', fontsize=FS)
#FigSK1.set_ylabel('S(k)', fontsize=FS)
#FigSK1.set_xscale('log', fontsize=FS)

#if eps_crit < 0.15:
#    print "Very low critical attraction (or very high Tc). Are you sure?"
#    print "If so, change the starting point of the eps_array"
#    sys.exit()

delta_eps = 0.02
#Supercritical isotherms for eps < eps_crit.  
eps_array = np.arange(eps_crit * 0.80, 2 * eps_crit, delta_eps) 
T_array = 1.0/eps_array
nT = len(T_array)

print "Epsilon array = ", eps_array
print "T_array array = ", T_array
print "n(T) = ", nT

for iT in range(nT): #starting from highest temperature
    

    T = T_array[iT]
    eps = eps_array[iT]
    print "iT = ", iT, " T = ", T, " eps = ", eps
    i = 0
    cr_guess = 0  #start with cr_guess = 0 for each isotherm. 

    #Create one file for each m-n potential. 
    OutFileH = open(branch + "-m-" + str(m) + "-n-" + str(n) + "-eps-" + str(eps) + ".dat", 'w')
    for rho in rho_array:
        i = i + 1
        try:
            #Solve setting epsilon = 1.0, but instead changing T. 
            ListHeader, ListValues , gr, cr,  Sk = SinglePointCalculation(r, k, Ur, Ur_ref, T, rho, OutFile=None, cr_guess = cr_guess)        
            cr_guess = cr

            print ListHeader
            print ListValues
            OutPuts = OrderedDict(zip(ListHeader, ListValues))    

        
            if i == 1:
                OutFileH.write(OZF.ListToTabbedStr(OutPuts.keys()))

            OutFileH.write(OZF.ListToTabbedStr(OutPuts.values()))
                         

        except TypeError:
            break #Move to next T

    OutFileH.close() #for a given epsilon (or temperature)

        

#FigSK1.legend(loc='best', ncol=2)
#plt.show()
