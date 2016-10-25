import pyoz as oz
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import newton
import sys
from scipy.integrate import simps
#My own module
import Functions

FS = 15
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

dpd_att = oz.System()
r = dpd_att.r
k = dpd_att.k

#Hexane:
B2_Dortmund = -1729 #cm^3/mol
B2_nondim = Functions.CalcB2_nondim(B2_Dortmund)
print("B2_nondim = ", B2_nondim)

# --------------
a_list = [10, 15, 30, 50, 100, 200, 500, 750, 1000]

#Calculate sigma_eff as per Frenkel and Noro
sigma_eff = []
for a in a_list:
    U_rep = Functions.DPDPotential(a,r)
    sigma_eff.append( Functions.Calc_SigmaEff(r,U_rep))
    print("a = ", a, "; sigma_eff = ", sigma_eff[-1])


# --------------
form = 'II'   # POTENTIAL FORM

#Find values of epsilon such that B2(a, epsilon) = B2_nondim
#for different values of a.

eps_newton_list = []
U_min_list = []
for a in a_list:
    eps_newton, B2Error = Functions.Find_eps(a, form, B2_nondim)
    B2_newton = Functions.Calc_B2(a, eps_newton, form)
    eps_newton_list.append(eps_newton)
    U_att = Functions.CalcU_DPD_Attractive(a,r,eps_newton,form)
    Umin = np.abs(np.min(U_att))
    print("Umin = ", Umin)
    U_min_list.append(Umin)
    
#Plot epsilon(a) such that B2(a,epsilon) = B2_nondim.    
fig_eps = plt.figure()
fig_eps1 = fig_eps.add_subplot(121)
fig_eps1.set_xlabel('$a$', fontsize=FS+3)
fig_eps1.set_ylabel('$\epsilon$', fontsize=FS+3)

fig_eps2 = fig_eps.add_subplot(122)
fig_eps2.set_xlabel('$a$', fontsize=FS+3)
fig_eps2.set_ylabel('$U_{min}$', fontsize=FS+3)

fig_eps1.set_xscale('log')
fig_eps2.set_xscale('log')

fig_eps1.plot(a_list, eps_newton_list, linewidth=2, marker='o')
fig_eps2.plot(a_list, U_min_list, linewidth=2, marker='o')
fig_eps.tight_layout()

plt.show()


# ------------- NOW LET"S TRY TO CALCULTE THE SPINODAL ----------
#Plot the potential for a particular value of a in a_list.

i = 0
a = a_list[i]
eps_newton = eps_newton_list[i]

#Calculate the potential
U_att = Functions.CalcU_DPD_Attractive(a,r, eps_newton, form=form)
dpd_att.set_interaction(0,0,U_att)
#Calculate B2. Should be equal to B2_nondim
B2 = oz.properties.second_virial_coefficient(dpd_att)
assert np.abs(B2-B2_nondim) < 1e-3 #Check
B2 = np.round(B2,2)

U_min = np.round(np.min(U_att),2)
fig_ur = plt.figure()
fig_ur1 = fig_ur.add_subplot(111)
fig_ur1.plot(r, U_att, linewidth=2.0)
fig_ur1.set_xlim([0,1.1])
fig_ur1.axhline(0, color='black', linestyle='--', linewidth=2)
fig_ur1.set_title('eps = ' + str(np.round(eps_newton,2)) + '; Umin = ' + str(U_min) + "; B2 = " + str(B2))

#Mapping onto square well potential as per Noro and Frenkel:
U_rep = Functions.DPDPotential(a,r)
sigma_eff = Functions.Calc_SigmaEff(r,U_rep) #Effective hard sphere diameter. 
B2_star = Functions.Calc_B2Star(B2_nondim, sigma_eff) #Equivalent B2 for a square well
tau = Functions.Calc_tau(B2_star)
lmbd  = Functions.Calc_lambda(tau, U_att) #Effective range of potential.

#Create dictionary of stuff mapped to square-well fluid. 
SW_Map = {}
SW_Map['sigma_eff'] = sigma_eff
SW_Map['B2_star'] = B2_star
SW_Map['tau'] = tau
SW_Map['lmbd'] = lmbd
SW_Map['Tstar'] = np.abs(1.0/np.min(U_att))

print(SW_Map)

#plt.show()

#sys.exit()

#----

#Calculate g(r) and S(k) on the vapor side. 
#rho_array = np.arange(0,0.10, 0.001)
#Calculate g(r) and S(k) on the liquid side
rho_array = np.arange(10, 0.25, -0.5)

fig = plt.figure(figsize=[12,5])
fig_P = fig.add_subplot(121)
fig_U = fig.add_subplot(122)

fig_P.set_xlabel('$\\rho$', fontsize=FS)
fig_P.set_ylabel('$P$', fontsize=FS)

fig_U.set_xlabel('$\\rho$', fontsize=FS)
fig_U.set_ylabel('$U$', fontsize=FS)


rho_converged  = []
P_array = []
U_array = []

fig_grSk = plt.figure()
fig_gr = fig_grSk.add_subplot(121)
fig_gr.set_xlabel('$r$', fontsize=FS)
fig_gr.set_ylabel('$g(r)$', fontsize=FS)
fig_gr.set_xlim([0,4])

fig_Sk = fig_grSk.add_subplot(122)
fig_Sk.set_xlabel('$k/2\pi$', fontsize=FS)
fig_Sk.set_ylabel('$S(k)$', fontsize=FS)
fig_Sk.set_xlim([0,25])
