import numpy as np
import OZ.Sk_Analysis_Functions as SkF
import OZ.PP_Functions as PP
import sys
import matplotlib.pyplot as plt
import matplotlib
from OZ.OZ_Functions import ListToTabbedStr
from itertools import cycle

import pandas as pd 

colors= ['red', 'green', 'black', 'cyan', 'blue', 'magenta', 'black']
colorcycler_liq = cycle(colors)
colorcycler_vap = cycle(colors)
 

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

if len(sys.argv) == 3:
    FileName=sys.argv[1] #Data from a density scan along an isotherm 
    branch=sys.argv[2]

    print FileName
else:
    print "Usage: FitSpinodal_SingleIsotherm.py filename branch"
    sys.exit()


try:
    df = pd.read_csv(FileName, sep='\t')
except ValueError:
    print FileName
    print "Not able to read this file!!"
    #likely because this is an empty file with nothing in it. 
    sys.exit()
    
#print df

rho = df['rho'].values
kappa = df['kappa'].values
Sk0 = df['Sk0'].values

B2 = df['B2'].values[0]
print "B2 = ", B2
T = df['T'].values[0]
eps = 1.0/T


if branch == 'vap':
    P_Kappa = PP.Calc_Pcompressibility(rho, kappa, opt='cum', offset=None, kT=T, B2=B2)    
elif branch == 'liq':
    #Pressure is off by the pressure at coexistance. 
    P_Kappa = PP.Calc_Pcompressibility(rho, kappa, opt='cum', offset=None, kT=T, B2=0)
else:
    print "Error: Unknown branch"
    sys.exit()

df.insert(4, 'P_kappa', P_Kappa)
#print df
df.to_csv(FileName[:-4] + "_Pkappa.dat", sep='\t')

if len(rho)>3: #fitting 3 data points needs atleast 3 data points.
    if branch == 'vap':
        params_init_guess = (1.0, np.max(rho)*1.25, 1.5)            
        params_opt =  SkF.Fit(params_init_guess, rho, Sk0, 'left')
    elif branch == 'liq':
        params_init_guess = (1.0, np.min(rho)-0.05, 1.5)            
        params_opt =  SkF.Fit(params_init_guess, rho, Sk0, 'right')
    else:
        print "Error: Unknown branch"
        sys.exit( )
else:
    print "Not enough datapoints to fit. Need atleast 3"
    sys.exit()
        

c, rho_spinodal, exponent = params_opt
outfile = open('SpinodalPoint.dat', 'w')

Str = "#T \t eps \t branch \t rho_spinodal \n"
outfile.write(Str)
List = [T, eps , branch, rho_spinodal]
outfile.write(ListToTabbedStr(List))

T_str = str(np.round(T,2))
Fig = plt.figure(figsize=[8,6])
Fig1 = Fig.add_subplot(111)
Fig1.set_ylabel('Sk[0]', fontsize=FS)
Fig1.set_xlabel('$\\rho$', fontsize=FS)
Fig1.plot(rho, Sk0, 'r-d')
Fig1.set_title('T = ' + T_str, fontsize=FS)
Fig.savefig('Sk0-T' + T_str + ".pdf")


