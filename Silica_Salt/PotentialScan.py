import numpy as np
import matplotlib.pyplot as plt
import matplotlib

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

R = 8.314e-3 #kJ/mol/K
T = 300.0

RT = R*T

import OZ.Potentials as Pot


NaClMolWt = 58.0 #gram/mol

NaCl_WtPerc = np.arange(0.50,6.25,0.50)
NaCl_Conc = NaCl_WtPerc * 10.0 / NaClMolWt #mol/Liter 

DebyeLength_nm = 0.304/np.sqrt(NaCl_Conc)

print DebyeLength_nm
lb = 0.70 #nm, Bjerrum length, in water

ColloidDiameter = 40.0 #nm, colloid diameter
ColloidArea = np.pi * ColloidDiameter **2


#http://physics.nyu.edu/grierlab/charge6c/ - charge density of glass
#is given for a few different pH values
ChargeDensity = 1e4 #e/micron^2 !! PLACE HOLDER
ChargeDensity_nmsq = ChargeDensity/1e6 #/nm^2

ColloidCharge = np.round(ChargeDensity_nmsq * ColloidArea,2)
print "Colloid Charge = ", ColloidCharge 


#----------
sig = 1.0 
l_debye = DebyeLength_nm/ColloidDiameter
l_bjerrum = lb/ColloidDiameter

Z = ColloidCharge
#Equation 4: Bollinger, 2016 paper. 
A_array = (Z**2 * l_bjerrum)/ (1 + 0.50/l_debye)**2  #A/kT actually 

FigA = plt.figure()
FigA1 = FigA.add_subplot(111)
FigA1.set_xlabel('NaCl, wt%', fontsize=FS)
FigA1.set_ylabel('A/kT', fontsize=FS)
FigA1.plot(NaCl_WtPerc, A_array, 'r-d')

TitleStr = '$Z_{max}$ = %3.2f ; $\sigma_q = $ %2.1E e/nm$^2$' %(Z, ChargeDensity)
FigA1.set_title(TitleStr, fontsize=FS)


#-----

m = 100
n = 50

FigPot = plt.figure(figsize=[15,15])
FigB2 = plt.figure()
FigB21 = FigB2.add_subplot(111)
FigB21.set_xlabel('A', fontsize=FS)
FigB21.set_ylabel('B2', fontsize=FS)

eps_array = np.arange(0.0, 6.0, 1.0)

dr = 0.02
nr = 2048
r = np.arange(nr)*dr + dr
ifig = 1
for eps in eps_array:
    Fig1 = FigPot.add_subplot(3,2,ifig)
    Fig1.set_xlim([0.8*sig, 5.0])
    Fig1.set_ylim([-5.0, 5.0])
    Fig1.set_title('$\\epsilon$ = ' + str(eps), fontsize=FS)
    B2 = []
    for A in A_array:
        Pot_100_50 = Pot.SALRPotential(r, sig=sig, eps=eps, A=A, m=100, n=50)
        Fig1.plot(r, Pot_100_50, marker='None', linestyle='-', linewidth=2)
        B2.append(Pot.CalcB2(r, Pot_100_50, kT=RT))

    print eps, B2

    FigB21.plot(A_array, B2, linestyle='-', linewidth=2, marker='o', label='$\\epsilon=$' + str(eps))
    

    ifig = ifig + 1
        

FigB21.legend(loc='best')
plt.show()







