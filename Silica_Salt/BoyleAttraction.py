import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import newton

FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)



import OZ.Potentials as Pot

def ErrorSq(NaCl_WtPerc, r, sig=1.0, eps=1.0, m=50, n=18, Q=1e4, d=40):

    Ur, B2 = Potential(NaCl_WtPerc, r, sig=sig, eps=eps, m=m, n=n, Q=Q, d=d)

    return B2**2 

def Potential(NaCl_WtPerc, r, sig=1.0, eps=1.0, m=50, n=18, Q=1e4, d=40):

    NaClMolWt = 58.0
    NaCl_Conc = NaCl_WtPerc/100.0 * 10.0 / NaClMolWt #mol/Liter
    l_debye = (0.304/np.sqrt(NaCl_Conc))/d #non-dimensional debye length
    
    l_bjerrum = 0.70/d #non-dimensional bjerrum length
    
    ChargeDensity_nmsq = Q/1e6 #/nm^2
    Area = np.pi * d**2
    Z = ChargeDensity_nmsq * Area
    A = (Z**2 * l_bjerrum)/ (1 + 0.50/l_debye)**2
    
    Ur = Pot.SALRPotential_v2(r, sig=sig, eps=eps, A=A, m=m, n=n)
    B2 = Pot.CalcB2(r,Ur)

    return Ur, B2


FigCs = plt.figure()
FigCs1 = FigCs.add_subplot(111)
FigCs1.set_xlabel('$\epsilon$/kT', fontsize=FS)
FigCs1.set_ylabel('NaCl wt%', fontsize=FS)
FigCs1.set_yscale('log')    


d = 40.0 #nm, colloid diameter

sig=1.0
m = 50; n = 18
Q = 1e4

eps_array = np.linspace(2.0,4.0,20)


dr = 0.02
nr = 2048
r = np.arange(nr)*dr + dr

Q  = [5e3, 1e4, 2e4]
print "#Q (e/nm^2) \t  eps/kT \t cs_wt_perc \t B2"
for q in Q:
    Cs_Boyle = []
    eps_converged = []
    for eps in eps_array:
        try:
            cs_wt_perc  = newton(ErrorSq, 0.50, args=(r, sig, eps, m, n, q, d), tol=1.48e-08, maxiter=50)
            eps_converged.append(eps)
            Cs_Boyle.append(cs_wt_perc)
        except RuntimeError:
            cs_wt_perc = -1.0
        
        if cs_wt_perc > 0.0:
            Ur, B2 = Potential(cs_wt_perc, r, sig=sig, eps=eps, m=m, n=n, Q=q, d=d)
            print q/1e6, eps, cs_wt_perc, B2

    labelstr = 'Q =' + str(q/1e6) + " e/nm$^2$"
    
    FigCs1.plot(eps_converged, Cs_Boyle, marker='d', linewidth=2.0, label=labelstr)



FigCs1.legend(loc='upper right')

FigCs.savefig('Cs_Boyle.pdf')
plt.show()
    


#FigA = plt.figure()
#FigA1 = FigA.add_subplot(111)
#FigA1.set_xlabel('NaCl, wt%', fontsize=FS)
#FigA1.set_ylabel('A/kT', fontsize=FS)
#FigA1.plot(NaCl_WtPerc, A_array, 'r-d')
#
#TitleStr = '$Z_{max}$ = %3.2f ; $\sigma_q = $ %2.1E e/nm$^2$' %(Z, ChargeDensity_nmsq)
#FigA1.set_title(TitleStr, fontsize=FS)
#
#
##-----
#
#
#FigPot = plt.figure(figsize=[15,15])
#FigB2 = plt.figure(figsize=[12,6])
#FigB21 = FigB2.add_subplot(121)
#FigB21.set_xlabel('A', fontsize=FS)
#FigB21.set_ylabel('B2', fontsize=FS)
#
#FigB22 = FigB2.add_subplot(122)
#FigB22.set_xlabel('NaCl Wt %', fontsize=FS)
#FigB22.set_ylabel('B2', fontsize=FS)
#
#eps_array = np.arange(1.6, 2.6, 0.2)
#
#
#
#dr = 0.02
#nr = 2048
#r = np.arange(nr)*dr + dr
#ifig = 1
#for eps in eps_array:
#    Fig1 = FigPot.add_subplot(3,2,ifig)
#    Fig1.set_xlim([0.8*sig, 5.0])
#    Fig1.set_ylim([-5.0, 5.0])
#    Fig1.set_title('$\\epsilon$ = ' + str(eps), fontsize=FS)
#    B2 = []
#
#    for isalt in range(nsalt):
#        Ur = Pot.SALRPotential_v2(r, sig=sig, eps=eps, A=A_array[isalt], m=50, n=18)
#        Fig1.plot(r, Ur, marker='None', linestyle='-', linewidth=2, label=str(NaCl_WtPerc[isalt]))
#        B2.append(Pot.CalcB2(r, Ur, kT=1.0))
#
#    print eps, B2
#    Fig1.legend(loc='upper right')
#
#    FigB21.plot(A_array, B2, linestyle='-', linewidth=2, marker='o', label='$\\epsilon=$' + str(eps))
#    FigB22.plot(NaCl_WtPerc, B2, linestyle='-', linewidth=2, marker='o', label='$\\epsilon=$' + str(eps))
#    
#
#    ifig = ifig + 1
#
#FigB21.legend(loc='best')
#FigB22.legend(loc='best')
#
#FigB2.tight_layout()
#FigPot.tight_layout()
#
#print "Nacl Conc = ", NaCl_Conc
#print "Debye Length (nm) = ", DebyeLength_nm
#plt.show()
#
#FigPot.savefig('Potentials.pdf')
#FigB2.savefig('B2.pdf')
#FigA.savefig('A.pdf')




