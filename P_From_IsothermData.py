import numpy as np
import sys
import matplotlib.pyplot as plt 
import matplotlib
#from itertools import cycle
import pandas as pd

import OZ.Potentials as Potentials
import OZ.OZ_Functions as OZF
import OZ.PP_Functions  as PP_Functions
import OZ.FigFuncs as FigFuncs
#import HNC
import OZ.RHNC as RHNC
from OZ.SinglePointCalculation import * 


DataFile = sys.argv[1]

df = pd.read_csv(DataFile, sep='\t')
print df
kappa = df['kappa'].values
rho = df['rho'].values
B2 = df['rho'].values[0]
T = df['T'].values[0]

P_kappa = PP_Functions.Calc_Pcompressibility(rho, kappa, opt='cum', offset=None, kT=T, B2=B2)
print len(kappa), len(P_kappa)

Fig = FigFuncs.CreateFig(rho, P_kappa, title=None, xlabel='$\\rho$', ylabel='Pressure', lw=2, marker='d', ls='-')
plt.show()

df.insert(6,'P_kappa', P_kappa)

df.to_csv(DataFile[:-4] + "_Pkaapa.dat", sep='\t')

