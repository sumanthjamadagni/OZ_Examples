import numpy as np
import sys
import matplotlib.pyplot as plt 
import matplotlib
#from itertools import cycle

import OZ.Potentials as Potentials
import OZ.OZ_Functions as OZF
import OZ.PP_Functions  as PP_Functions
#import FigFuncs
#import HNC
import OZ.RHNC as RHNC
from OZ.SinglePointCalculation import * 


dr = 0.02
nr = 2048 
sig=1.0

r,k = OZF.Create_r_k(dr, nr)

Ur_ref = Potentials.WCAPotential(r,sig=1.0, eps=1.0, m=12, n=6) 
Ur = Potentials.LJPotential(r,sig=1.0, eps=1.0, m=12, n=6)

T = float(sys.argv[1])

OutFile = sys.argv[2]
OutFileH = open(OutFile, 'w')

branch = sys.argv[3]

if branch == 'vap':
    rho_array = np.arange(0.01,0.25,0.005)
elif branch == 'liq'
    rho_array = np.arange(0.80,0.50,-0.01)
else:
    print "Unknown branch. Has to be 'vap' or 'liq'"
    sys.exit()

cr_guess = None
i = 0
for rho in rho_array:
    i = i + 1
    try:
        ListHeader, ListValues , gr, cr,  Sk = SinglePointCalculation(r, k, Ur, Ur_ref, T, rho, OutFile=None, cr_guess = cr_guess)
        cr_guess = cr
        print ListHeader
        print ListValues
        if i == 1:
            OutFileH.write(OZF.ListToTabbedStr(ListHeader))

        OutFileH.write(OZF.ListToTabbedStr(ListValues))
                         

    except TypeError:
        break

