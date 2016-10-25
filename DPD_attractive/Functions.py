import pyoz as oz
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys

from scipy.optimize import newton
from scipy.integrate import simps

FS = 15
matplotlib.rc('xtick', labelsize=FS)
matplotlib.rc('ytick', labelsize=FS)

def DPDPotential(a,r):
    U_DPD = np.where(r < 1, 0.5 * a  * (1-r)**2 , 0)
    return U_DPD
        
def CalcU_Perturbation(a,r, eps, form='I'):
    a_att = a/8.0 + eps
    U_p = a_att * np.sin(np.pi * r)**2
    U_perturbation = np.where(r<1.0, U_p, 0)
    if form == 'I':
        return U_perturbation        
    elif form == 'II':
        U0 = np.where(r < 0.5, a_att, 0)
        U1 = np.where(r > 0.5, U_perturbation, 0)
        U_perturbation_WCA  = U0 + U1
        return U_perturbation_WCA

    else:
        print("Unknown form of  potential. Has to be I or II")
        return np.nan
    
def CalcU_DPD_Attractive(a,r, eps, form):    
    U_total = DPDPotential(a,r) - CalcU_Perturbation(a,r,eps, form)
    return U_total

#-------------------

def Solve_DPD_Att(rho, a, eps, form, e_r_guess = None, mix_param = 0.60, max_iter = 10000):
    dpd_att = oz.System()
    r = dpd_att.r
    U_att = CalcU_DPD_Attractive(a,r, eps, form)
    dpd_att.set_interaction(0,0,U_att)    
    g_r, c_r, e_r, H_k = dpd_att.solve(rhos=rho, closure_name='hnc', initial_e_r = e_r_guess, max_iter=max_iter, mix_param=mix_param)
    return dpd_att

def Calc_P_U(dpd_att):
    try:
        P_virial = oz.properties.pressure_virial(dpd_att)
    except TypeError:
        P_virial = -1

    try:
        U = oz.properties.internal_energy(dpd_att)
    except TypeError:
        U = -1
    return P_virial, U

#def Error(params, P_target, U_target, w_p=0.5):
#    a, eps, form = params
#    try:
#        dpd_att = Solve_DPD_Att(rho, a, eps, form)
#    except:
#        return -1 
#    
#    P_virial, U = Calc_P_U(dpd_att)
#    Error = w_p * np.abs(1.0 - P_virial/P_target) + (1.0-w_p) * np.abs(1.0 - U/U_target)
#    return Error

def Calc_B2(a,eps, form):
    dpd_att = oz.System()
    r = dpd_att.r
    U_att = CalcU_DPD_Attractive(a,r, eps,form)
    dpd_att.set_interaction(0,0,U_att)
    B2 = oz.properties.second_virial_coefficient(dpd_att)
    return B2

def Calc_B2_epsarray(a,eps_array,form):
    B2_array = []
    for eps in eps_array:
        B2 = Calc_B2(a,eps,form)
        B2_array.append(B2)

    B2_array = np.array(B2_array)
    return B2_array
    
def Calc_B2Error(eps, a, B2_target, form):    
    B2 = Calc_B2(a,eps,form)
    print("Func = B2Error; ", "eps = ", eps, ";B2 = ", B2)
    B2Error = (1.0 - B2/B2_target)**2
    return B2Error

def CalcB2_nondim(B2_Dortmund, rc=7.66):
    '''
        B2_Dortmund = cc/mol
        rc = Lengthscale of DPD, Angstroms. 
    '''
    B2_nondim = B2_Dortmund * 1e-6/(rc * 1e-10)**3 / 6.023e23 
    return B2_nondim

def Find_eps(a, form, B2_nondim):
    """
    Find the value of eps to use such that B2(a,eps) = B2_nondim. 

    #NEWTON RAPHSON DOESN"T CONVERGE UNLESS GIVEN VERY GOOD GUESS.
    #SO BIT OF WORK
    """
    print("a=",a)
    print("form=",form)
    print("B2_nondim=",B2_nondim)
    N = 1000
    eps_array = np.linspace(-a/8.0,10,N)
    B2_array = Calc_B2_epsarray(a, eps_array, form)
    idx = np.argmax(B2_array < B2_nondim)
    print("idx = ", idx)
    eps_guess = 0.5 * (eps_array[idx-1] + eps_array[idx])
    #Now use the newton solver
    eps_newton  = newton(Calc_B2Error, eps_guess, args=(a, B2_nondim, form))
    B2Error =  np.sqrt(Calc_B2Error(eps_newton, a, B2_nondim, form))

    return eps_newton, B2Error



#Noro & Frenkel Calculations. 
def Calc_SigmaEff(r, U_rep, kT=1.0):
    #Noro and Frenkel, Corresponding state for colloids, Eqn. 4.
    Integrand = (1.0 - np.exp(-U_rep/kT))
    sigma_eff = simps(Integrand, r)
    return sigma_eff

def Calc_B2Star(B2, sigma_eff):
    #Eqn. 6, Noro and Frenkel, 2000
    B2_star = B2/ (2 * np.pi * sigma_eff**3 / 3.0)
    return B2_star

def Calc_tau(B2_star):
    #Eqn. 8, Noro and Frenkel
    tau = 1.0/(4.0 * (1.0 - B2_star))
    return tau

def Calc_lambda(B2_star, U_att):
    
    Tstar = np.abs(1.0/np.min(U_att))
    #Eqn. 10, Noro and Frenkel
    lmbd = (1.0 - B2_star)/(1.0 - np.exp(1.0/Tstar)) + 1.0
    lmbd = lmbd ** (1.0/3.0)
    return lmbd


