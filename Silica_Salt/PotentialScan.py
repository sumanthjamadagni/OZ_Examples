import numpy as np
import matplotlib.pyplot as plt
import matplotlib

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

R = 8.314e-3 #kJ/mol/K
T = 300.0

RT = R*T 
print "RT = ", RT , "kJ/mol"

import OZ.Potentials as Pot


NaClMolWt = 58.0 #gram/mol

NaCl_WtPerc = [5e-1, 1.0, 2.0, 5.0, 10.0]
NaCl_WtPerc = np.array(NaCl_WtPerc)
NaCl_Conc = NaCl_WtPerc/100.0 * 10.0 / NaClMolWt #mol/Liter
nsalt = len(NaCl_WtPerc)

DebyeLength_nm = 0.304/np.sqrt(NaCl_Conc)
lb = 0.70 #nm, Bjerrum length, in water at room temperature. 

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


print "A_array = ", A_array

FigA = plt.figure()
FigA1 = FigA.add_subplot(111)
FigA1.set_xlabel('NaCl, wt%', fontsize=FS)
FigA1.set_ylabel('A/kT', fontsize=FS)
FigA1.plot(NaCl_WtPerc, A_array, 'r-d')

TitleStr = '$Z_{max}$ = %3.2f ; $\sigma_q = $ %2.1E e/nm$^2$' %(Z, ChargeDensity_nmsq)
FigA1.set_title(TitleStr, fontsize=FS)


#-----


FigPot = plt.figure(figsize=[15,15])
FigB2 = plt.figure(figsize=[12,6])
FigB21 = FigB2.add_subplot(121)
FigB21.set_xlabel('A', fontsize=FS)
FigB21.set_ylabel('B2', fontsize=FS)

FigB22 = FigB2.add_subplot(122)
FigB22.set_xlabel('NaCl Wt %', fontsize=FS)
FigB22.set_ylabel('B2', fontsize=FS)

eps_array = np.arange(1.6, 2.6, 0.2)



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

    for isalt in range(nsalt):
        Ur = Pot.SALRPotential_v2(r, sig=sig, eps=eps, A=A_array[isalt], m=50, n=18)
        Fig1.plot(r, Ur, marker='None', linestyle='-', linewidth=2, label=str(NaCl_WtPerc[isalt]))
        B2.append(Pot.CalcB2(r, Ur, kT=1.0))

    print eps, B2
    Fig1.legend(loc='upper right')

    FigB21.plot(A_array, B2, linestyle='-', linewidth=2, marker='o', label='$\\epsilon=$' + str(eps))
    FigB22.plot(NaCl_WtPerc, B2, linestyle='-', linewidth=2, marker='o', label='$\\epsilon=$' + str(eps))
    

    ifig = ifig + 1

FigB21.legend(loc='best')
FigB22.legend(loc='best')
FigB2.tight_layout()

FigPot.tight_layout()

print "Nacl Conc = ", NaCl_Conc
print "Debye Length (nm) = ", DebyeLength_nm
plt.show()





m_list = [12, 100]
n_list = [6, 12, 24, 36, 50]
for m in m_list:
    for n in n_list:
        if m > n:            
            Pot_m_n = Pot.SALRPotential_v2(r, sig=sig,  eps=1.0, A=0, m=m, n=n)

            Pot_m_n_rep = Pot.SRLRPotential_v2(r, sig=sig, eps=1.0, A=0, m=m, n=n)
            B2 = Pot.CalcB2(r, Pot_m_n, kT=1.0) 
            B2_rep = Pot.CalcB2(r, Pot_m_n_rep, kT=1.0)
            print m, n, B2, B2_rep
            

