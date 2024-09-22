import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import erf


# Compute analytical estimate from Najmabadi. It's more easily implemented from Post 1987
def Najmabadi_confinement_time(P,R, ZpFl=1):
    w_term = np.sqrt(1 + 1/(R*(ZpFl - 1/(4*P)))) #Is this 1/4P or P/4
    u_eff  = P + np.log(w_term) #same as Najmabadi's a**2

    integrandNaj = lambda t: np.exp(-t) / t
    I_term, error = quad(integrandNaj, u_eff, np.inf)
    I_term = (ZpFl + 1/4)*u_eff*np.exp(u_eff)*I_term - 1/4 

    Loss_Najmabadi = \
        np.sqrt(np.pi)/4 * u_eff*np.exp(u_eff)/I_term \
        * (np.log((w_term+1)/(w_term-1))-0.84)#0.84
    
    return Loss_Najmabadi

def Rosen_confinement_time(P,R,coeff = 0.84, ZpFl = 1):
    w_term = np.sqrt(1 + 2*P/(R*ZpFl)) 
    a_term  = np.sqrt(P + np.log(w_term)) 

    integrandRosen = lambda t: (ZpFl/2 )*np.exp(-t**2)

    I_term, error = quad(integrandRosen, a_term, np.inf)

    Loss_Rosen = \
        1/(4*ZpFl/2 / (np.log((w_term+1)/(w_term-1))- coeff) * (1-erf(a_term)))
    
    return Loss_Rosen

ephi_T = 4.468
R = 17
Z_N = 1/2
Z_dough = 1
m_e = 9.10938356e-28
Te = 940
elc_chg = 1.60217662e-19
vareps = 8.85418782e-12
density = 3e19 ##m^-3
density_cm3 = density * 1e-6

kTe = 1.6e-12 * Te

# CGS
logLambdaElc = 23.5 - np.log(np.sqrt(density_cm3) * Te**(-5/4)) - np.sqrt(1e-5 + (np.log(Te)-2)**2/16)

# SI mks
m_e = 9.10938356e-31
vth_e  = np.sqrt(2*Te*elc_chg/m_e)
print("vth_e: ", vth_e)
nu_e = 4 * np.pi / (m_e**2 * vth_e**3) * (elc_chg**2 / (4*np.pi *vareps))**2 * density * logLambdaElc

print("Collision frequency najmabadi: ", nu_e)
print("Collision period najmabadi: ", 1/nu_e)
print("logLambda Najmabadi: ", logLambdaElc)


logLambdaElc = 6.6 - 0.5 * np.log(density / 1e20) + 1.5 * np.log(Te);
vth_e  = np.sqrt(Te/m_e)
nu_e = logLambdaElc * elc_chg**4 * density / (6 * np.sqrt(2) * np.pi**(3/2) * vareps**2 * np.sqrt(m_e) * (Te*elc_chg)**(3/2))
# This is the same as egedal 2022 Nucl. Fusion


  # double nuElc = nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
  #                (6. * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(me) * pow(Te0*elc_chg, 3. / 2.));

print("Collision frequency rosen: ", nu_e)
print("Collision period rosen: ", 1/nu_e) 
print("logLambda Rosen: ", logLambdaElc)
Z_e = 1/2
Z_dough = 1

Loss_Najmabadi = Najmabadi_confinement_time(ephi_T, R, ZpFl=Z_e)
Loss_Rosen = Rosen_confinement_time(ephi_T, R, coeff=0.96, ZpFl=Z_dough)

fraction = Loss_Najmabadi / Loss_Rosen
print(fraction)
