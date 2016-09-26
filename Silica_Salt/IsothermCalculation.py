import numpy as np
import sys
import matplotlib.pyplot as plt 
import matplotlib
#from itertools import cycle
from collections import OrderedDict

import OZ.Potentials as Potentials
import OZ.OZ_Functions as OZF
import OZ.PP_Functions  as PP_Functions
from OZ.FigFuncs import CreateFig
from OZ.SinglePointCalculation import * 
import optparse 




FS = 18
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

parser = optparse.OptionParser()
parser.add_option('-m', help='Exponent of repulsive potential. Default=%default', default=50, dest='m', type='int')
parser.add_option('-n', help='Exponent of attractive potential. Default=%default', default=18, dest='n', type='int')
parser.add_option('-e', '--eps', help='Attractions. Default=%default', default=2.0, dest='eps', type='float')
parser.add_option('-T', '--Temp', help='Temperature. Default=%default', default=1.0, dest='T', type='float')
parser.add_option('-Q', help='Charge density, #e/micron^2. Default=%default', default=1e4, dest='ChargeDensity', type='float')
parser.add_option('-d', help='Diameter of silica particle, nm. Default=%default', default=40.0, dest='ColloidDiameter', type='float')
parser.add_option('--Cs', help='Salt concentration, wt percent. Default=%default', default=0.50, dest='NaCl_WtPerc', type='float')
parser.add_option('--branch', help='Equal to vap or liq. Default=%default', default='vap', dest='branch')

parser.add_option('--ptype', help='Standard or alternate version of LJ/WCA potentials. Default=%default', default='v1', dest='PotType')
parser.add_option('--B2', help='Print B2 and quit. Default=%default', default='False', dest='B2Flag', action='store_true')


 
(opts, args) = parser.parse_args()

print opts

if opts.NaCl_WtPerc == -1.0: #flag for switching off electrostatics
    A = 0
    l_debye = 1.0
    print "Setting A = 0"
else:
    #Calculate Electrostatic part of SALR potentisl
    Bjerrum_Length = 0.70 #nm, Bjerrum length, in water at room temperature.

    NaClMolWt = 58.0 #g/mol 
    NaCl_Conc = opts.NaCl_WtPerc/100.0 * 10.0 / NaClMolWt #mol/Liter
    DebyeLength_nm = 0.304/np.sqrt(NaCl_Conc)

    #Charge on colloid
    ColloidArea = np.pi * opts.ColloidDiameter **2
    #http://physics.nyu.edu/grierlab/charge6c/ - charge density of glass
    #is given for a few different pH values
    ChargeDensity_nmsq = opts.ChargeDensity/1e6 #/nm^2
    ColloidCharge = np.round(ChargeDensity_nmsq * ColloidArea,2)
    print "Colloid Charge = ", ColloidCharge


    #Non-dimensionalizing by colloid diameter  
    l_debye = DebyeLength_nm/opts.ColloidDiameter
    l_bjerrum = Bjerrum_Length/opts.ColloidDiameter

    Z = ColloidCharge
    #Equation 4: Bollinger, 2016 paper. 
    A = ((Z**2 * l_bjerrum)/ (1 + 0.50/l_debye)**2)
    print "A = ", A
# ------------------

if opts.branch == 'vap':
    rho_array = np.arange(0.01,0.81, 0.01)
elif opts.branch == 'liq':
    rho_array = np.arange(0.81,0.01,-0.01)
else:
    print "Unknown branch. Has to be 'vap' or 'liq'"
    print opts.branch
    sys.exit()


Params = OrderedDict()
Params['T'] = opts.T
Params['m'] = opts.m
Params['n'] = opts.n
Params['eps'] = opts.eps
Params['A'] = A
if A > 0:
    Params['ChargeDensity'] = opts.ChargeDensity
    Params['ColloidCharge'] = ColloidCharge    
    Params['l_debye'] = l_debye


FigSK = plt.figure()
FigSK1 = FigSK.add_subplot(111)
FigSK1.set_xlabel('k$\sigma$/2$\pi$', fontsize=FS)
FigSK1.set_ylabel('S(k)', fontsize=FS)
FigSK1.set_xscale('log', fontsize=FS)


FigGr = plt.figure()
FigGr1 = FigGr.add_subplot(111)
FigGr1.set_xlabel('r', fontsize=FS)
FigGr1.set_ylabel('g(r)', fontsize=FS)
FigGr1.set_xlim([0,4])


dr = 0.02
nr = 2048
sig=1.0

r,k = OZF.Create_r_k(dr, nr)
if opts.PotType == 'v1':
    Ur_ref = Potentials.SRLRPotential(r,sig=sig, eps=opts.eps, m=opts.m, n=opts.n, A=0, d=l_debye)
    Ur = Potentials.SALRPotential(r,sig=sig, eps=opts.eps, m=opts.m, n=opts.n, A=A, d=l_debye)
elif opts.PotType == 'v2':
    Ur_ref = Potentials.SRLRPotential_v2(r,sig=sig, eps=opts.eps, m=opts.m, n=opts.n, A=0, d=l_debye)
    Ur = Potentials.SALRPotential_v2(r,sig=sig, eps=opts.eps, m=opts.m, n=opts.n, A=A, d=l_debye)
else:
    print "ptype should be v1 or v2"
    sys.exit()


B2_ref = Potentials.CalcB2(r,Ur_ref,T=opts.T)
B2 = Potentials.CalcB2(r,Ur,T=opts.T)


Params['B2'] = B2
Params['B2_ref'] = B2_ref


B2_reduced = Potentials.CalcB2_reduced(r, Ur, Ur_ref, T=opts.T)
Params['B2_reduced'] = B2_reduced

sig_eff = Potentials.Sigma_Eff(r,Ur_ref,T=opts.T)
Params['Sigma_eff'] = sig_eff 

lab = "B2_ref = " + str(np.round(B2_ref,2))
FigUr_ref = CreateFig(r,Ur_ref,title=lab, xlabel='r', ylabel='Ur_ref', xleft=0.0, xright=4.0, ybottom = -1.0, ytop=5.0, lw=2, marker='d', ls='-')

lab = "B2 = " + str(np.round(B2,2))
FigUr = CreateFig(r,Ur,title=lab, xlabel='r', ylabel='Ur', xleft=0.0, xright=4.0, ybottom = np.min(Ur)*1.5, ytop=5.0, lw=2, marker='d', ls='-')

if opts.B2Flag is True:
    print Params
    print "B2 = ", B2
    sys.exit()


#---------
cr_guess = None

PlotFreq=5 
i = 0

delta_eps = 0.02
T_array = [opts.T]
nT = len(T_array)

#print "Epsilon array = ", eps_array
print "T_array array = ", T_array
print "n(T) = ", nT

for iT in range(nT): #starting from highest temperature
    

    T = T_array[iT]
    eps = 1.0/T
    print "iT = ", iT, " T = ", T, " eps = ", eps
    i = 0
    cr_guess = 0  #start with cr_guess = 0 for each isotherm. 

    #Create one file for each m-n potential.
    #OutFileH = open(opts.branch + "-m-" + str(opts.m) + "-n-" + str(opts.n) + "-eps-" + str(opts.eps) + "-Cs-" + str(opts.NaCl_WtPerc) + ".dat", 'w')
    OutFileH = open("OZCalc.dat", 'w')


    for key in Params.keys():
        StrVal = "#" + key + "\t" + str(Params[key]) + "\n"
        OutFileH.write(StrVal)

    for rho in rho_array:
        i = i + 1


        try:
            #Solve setting epsilon = 1.0, but instead changing T. 
            ListHeader, ListValues , gr, cr,  Sk = SinglePointCalculation(r, k, Ur, Ur_ref, rho, T=T, OutFile=None, cr_guess = cr_guess)        
            cr_guess = cr

            print ListHeader
            print ListValues
            OutPuts = OrderedDict(zip(ListHeader, ListValues))

            FigSK1.plot(k/(2.0*np.pi), Sk, marker='None', linewidth=0.5, label=str(rho))
            FigGr1.plot(r,gr, marker='None', label=str(rho), linewidth=0.5)


        
            if i == 1:
                OutFileH.write(OZF.ListToTabbedStr(OutPuts.keys()))

            OutFileH.write(OZF.ListToTabbedStr(OutPuts.values()))
                         

        except TypeError:
            print "Didn't converge"
            print "rho = ", rho
            rho_limit = rho 
            break #Move to next T

    OutFileH.close() #for a given epsilon (or temperature)

if rho == rho_array[-1]:
    rho_limit = rho #converged everywhere

FigGr1.set_title('$\\rho_{max}$ = ' + str(rho_limit), fontsize=FS)
FigSK1.set_title('$\\rho_{max}$ = ' + str(rho_limit), fontsize=FS)

FigSK.savefig("Sk.pdf")
FigGr.savefig('gr.pdf')

FigUr.savefig('Ur.pdf')
FigUr_ref.savefig('Ur_ref.pdf')

plt.show()
