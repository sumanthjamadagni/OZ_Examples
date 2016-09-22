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

dr = 0.01
nr = 4096
sig=1.0

r,k = OZF.Create_r_k(dr, nr)
sig= 1.0
m = 12
n = 6


Ur_ref = Potentials.WCAPotential_v2(r,sig=sig, eps=1.0, m=m, n=n)
Ur = Potentials.LJPotential_v2(r,sig=sig, eps=1.0, m=m, n=n)



fig = plt.figure()
fig1 = fig.add_subplot(111)
fig1.set_xlim([0.5,3])
fig1.set_ylim([-1,3])
fig1.plot(r,Ur_ref, 'r-d', label='Ur_ref')
fig1.plot(r,Ur, 'b-d', label='Ur')
fig1.legend(loc='upper right')



T = 0.8  #/0.90
Fig_gr  = plt.figure()
Fig1  = Fig_gr.add_subplot(121)
Fig1.set_xlim([0.5,3])



Fig2  = Fig_gr.add_subplot(122)
Fig2.set_xscale('log')


rho_array = np.arange(0.01,0.10,0.005)
cr_guess =  None
i = 0
for rho in rho_array:
        i = i + 1
        try:
            #Solve setting epsilon = 1.0, but instead changing T. 
            ListHeader, ListValues , gr, cr,  Sk = SinglePointCalculation(r, k, Ur, Ur_ref, T, rho, OutFile=None, cr_guess = cr_guess)        
            cr_guess = cr

            print ListHeader
            print ListValues
            OutPuts = OrderedDict(zip(ListHeader, ListValues))

            Fig1.plot(r,gr,label=str(rho))
            Fig2.plot(k,Sk,label=str(rho))

        
#            if i == 1:
 #               OutFileH.write(OZF.ListToTabbedStr(OutPuts.keys()))

#            OutFileH.write(OZF.ListToTabbedStr(OutPuts.values()))
                         

        except TypeError:
            break #Move to next T


Fig1.plot(r, np.exp(-Ur/T), color='black', linewidth=2.0)
Fig1.legend(loc='upper right')
plt.show()
