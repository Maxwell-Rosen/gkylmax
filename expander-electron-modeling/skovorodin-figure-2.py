import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, root_scalar
import postgkyl as pg

def model_phi(R, R_wall, miome):
    term1 = np.log(R)
    term2 = np.log(np.sqrt(miome)) + 1 - (miome/R**2)**(1/3)
    ephiote = 1/(1 + np.exp( (R - np.sqrt(miome)))) * term1 
    ephiote +=1/(1 + np.exp(-(R - np.sqrt(miome)))) * term2

    # # Asymmetrical sigmoidal function using mycurvefit.com
    # ### Data
    # # R = [100, 200, 300, 400]
    # # fit factor magnitude = [0.3, 0.65, 0.85, 1]
    # # Using a point (1,0)
    # fac_in_front = 149818 + (-7.926893e-8 - 149818)/(1 + (R_wall/55.25151)**7.370451)**4.566507e-7
    # # Asymmetrical sigmoidal function using mycurvefit.com
    # ### Data
    # # R = [100, 200, 300, 400]
    # # fit factor scale = [5, 30, 45, 60]
    # # Using a point (1, 1e-2)
    # fac_scale = 17378520 + (0.06375668 - 17378520)/(1 + (R_wall/102.732)**5.375703)**4.648203e-7

    # if (R_wall == 100):
    #   a, b = 0.270096, 1.858309
    # elif (R_wall == 200):
    #   a, b = 0.646006, 15.762495
    # elif (R_wall == 300):
    #   a, b = 0.886529, 26.772415
    # elif (R_wall == 400):
    #   a, b = 1.078145, 42.958035
    # ephiote = np.log(R) - a + a/(R/b + 1)
    param_vals = np.array([[-2.90613517e-02,  1.93056686e-02, -4.32697290e-04,  4.29466477e-06,  -1.64631036e-08], 
                           [-7.33629896e-05,  1.35346122e-02, -1.36326298e-04,  6.60618528e-07,  -1.26530356e-09], 
                           [ 2.71723325e-02,  1.08481831e-02, -6.93139061e-05,  2.15473634e-07,  -2.66748298e-10], 
                           [ 6.87163743e-02,  8.49415386e-03, -3.87134748e-05,  9.13360993e-08,  -9.04080404e-11]])
    if (R_wall == 100):
        param = param_vals[0,:]
    elif (R_wall == 200):
        param = param_vals[1,:]
    elif (R_wall == 300):
        param = param_vals[2,:]
    elif (R_wall == 400):
        param = param_vals[3,:]

    
    ephiote = np.log(R) -( param[0] + param[1]*R + param[2]*R**2 + param[3]*R**3 + param[4]*R**4)
    return ephiote

def sheath_jump(R_wall, miome):
    phi_sheath = 2.6/(R_wall**(2/3))*(miome)**(1/3)
    return phi_sheath

paper_boltzmann = np.loadtxt('Boltzmann.csv', delimiter=',')
paper_line1 = np.loadtxt('Line 1.csv', delimiter=',')
paper_line2 = np.loadtxt('Line 2.csv', delimiter=',')
paper_line3 = np.loadtxt('Line 3.csv', delimiter=',')
paper_line4 = np.loadtxt('Line 4.csv', delimiter=',')

filename_phi = str('./gk_wham-field_300.gkyl')
pgData_phi = pg.GData(filename_phi)
polyOrder = 1
pgInterp_phi = pg.GInterpModal(pgData_phi, polyOrder, 'ms')
x_phi, dataOut_phi = pgInterp_phi.interpolate()

filename_bmag = str('./gk_wham-bmag.gkyl')
pgData_bmag = pg.GData(filename_bmag)
pgInterp_bmag = pg.GInterpModal(pgData_bmag, polyOrder, 'ms')
x_bmag, dataOut_bmag = pgInterp_bmag.interpolate()

filename_BiMaxwellianMoms = str('./gk_wham-ion_BiMaxwellianMoments_300.gkyl')
pgData_BiMaxwellianMoms = pg.GData(filename_BiMaxwellianMoms)
pgInterp_BiMaxwellianMoms = pg.GInterpModal(pgData_BiMaxwellianMoms, polyOrder, 'ms')
x, nOut = pgInterp_BiMaxwellianMoms.interpolate(0)
x, uOut = pgInterp_BiMaxwellianMoms.interpolate(1)
x, TparOut = pgInterp_BiMaxwellianMoms.interpolate(2)
x, TperpOut = pgInterp_BiMaxwellianMoms.interpolate(3)

bmag_shape = dataOut_bmag.shape
midpoint = int(bmag_shape[0]/2)
upperhalf = dataOut_bmag[midpoint:]
peak = np.argmax(upperhalf)
peak_idx = midpoint-peak -10

Te0 = 940 # eV
phi = dataOut_phi[:peak_idx]/ Te0
bmag = dataOut_bmag[:peak_idx]
n = nOut[:peak_idx]
u = uOut[:peak_idx]
Tpar = TparOut[:peak_idx]
Tperp = TperpOut[:peak_idx]
Teomi = Te0 * 1.602176634e-19 / (2 * 1.67262192595e-27)

K = bmag[-1] / bmag
phi = -(phi - phi[-1])
u_par_mirror = u[-1]

phi_hammett = np.zeros(len(n))
for i in range(len(n)):
    idx = i
    def model_phi_root(phi):
        return (2* phi) - np.log((K[idx]**2) * \
            (1 + 2*(Tperp[-1]*(1 - 1/K[idx]) + Teomi * phi)/(u_par_mirror**2)))
    # Find the root of the function
    phi_hammett[i] = root_scalar(model_phi_root, bracket=[0, 10]).root
print(phi_hammett)

# print("u_par_mirror = ", u_par_mirror)
# print("u = ", u)

# phi_hammett = np.log(K * u/u_par_mirror)
R = np.linspace(1, 450, 100)
boltzmann = np.log(R)
miome = 1836
R100 = np.linspace(1, 100, 100)
R200 = np.linspace(1, 200, 100)
R300 = np.linspace(1, 300, 100)
R400 = np.linspace(1, 400, 100)
phi_100 = model_phi(R100, 100, miome)
phi_200 = model_phi(R200, 200, miome)
phi_300 = model_phi(R300, 300, miome)
phi_400 = model_phi(R400, 400, miome)

phi_100 = np.append(phi_100, phi_100[-1] + sheath_jump(100, miome))
phi_200 = np.append(phi_200, phi_200[-1] + sheath_jump(200, miome))
phi_300 = np.append(phi_300, phi_300[-1] + sheath_jump(300, miome))
phi_400 = np.append(phi_400, phi_400[-1] + sheath_jump(400, miome))
R100 = np.append(R100, 100)
R200 = np.append(R200, 200)
R300 = np.append(R300, 300)
R400 = np.append(R400, 400)

plt.figure()
plt.plot(paper_boltzmann[:, 0], paper_boltzmann[:, 1], 'r', label = "S. Boltzmann")
plt.plot(paper_line1[:, 0], paper_line1[:, 1], 'b', label = "S. Line 1")
plt.plot(paper_line2[:, 0], paper_line2[:, 1], 'g', label = "S. Line 2")
plt.plot(paper_line3[:, 0], paper_line3[:, 1], 'm', label = "S. Line 3")
plt.plot(paper_line4[:, 0], paper_line4[:, 1], 'c', label = "S. Line 4")
plt.plot(R, boltzmann, 'k--', label = "Boltzmann")
plt.plot(R100, phi_100, 'b--', label = "Model $R_{wall}=100$")
plt.plot(R200, phi_200, 'g--', label = "Model $R_{wall}=200$")
plt.plot(R300, phi_300, 'm--', label = "Model $R_{wall}=300$")
plt.plot(R400, phi_400, 'c--', label = "Model $R_{wall}=400$")
plt.plot(K, phi, 'r--', label = "Gkeyll")
plt.plot(K, phi_hammett, 'y--', label = "Gkeyll Hammett")
plt.xlabel('R')
plt.ylabel('$-e \phi / T_e$')
plt.xlim([0, 450])
plt.ylim([0, 6.5])
plt.legend()
plt.savefig('skovorodin-figure-2.png')
plt.close()


R100 = np.linspace(2, 95, 100)
R200 = np.linspace(2, 190, 100)
R300 = np.linspace(2, 290, 100)
R400 = np.linspace(2, 390, 100)

boltzmann_line4 = np.log(paper_line4[:, 0])
diff_all = boltzmann_line4 - paper_line4[:, 1]
bad_cells = np.isnan(boltzmann_line4) | (diff_all < 0.0)
bad_cells[0:3] = True
bad_cells[-1] = True
diff_400 = boltzmann_line4[~bad_cells] - paper_line4[~bad_cells, 1]
R_400_diff = paper_line4[~bad_cells, 0]
# diff_400 = np.append(diff_400, np.tile(diff_400[-1],50))
# R_400_diff = np.append(R_400_diff, np.tile(R_400_diff[-1],50))

boltzmann_line3 = np.log(paper_line3[:, 0])
diff_all = boltzmann_line3 - paper_line3[:, 1]
bad_cells = np.isnan(boltzmann_line3) | (diff_all < 0.0)
bad_cells[0:3] = True
bad_cells[-1] = True
diff_300 = boltzmann_line3[~bad_cells] - paper_line3[~bad_cells, 1]
R_300_diff = paper_line3[~bad_cells, 0]
# diff_300 = np.append(diff_300, np.tile(diff_300[-1],50))
# R_300_diff = np.append(R_300_diff, np.tile(R_300_diff[-1],50))

boltzmann_line2 = np.log(paper_line2[:, 0])
diff_all = boltzmann_line2 - paper_line2[:, 1]
bad_cells = np.isnan(boltzmann_line2) | (diff_all < 0.0)
bad_cells[0:3] = True
bad_cells[-1] = True
diff_200 = boltzmann_line2[~bad_cells] - paper_line2[~bad_cells, 1]
R_200_diff = paper_line2[~bad_cells, 0]
# diff_200 = np.append(diff_200, np.tile(diff_200[-1],50))
# R_200_diff = np.append(R_200_diff, np.tile(R_200_diff[-1],50))

boltzmann_line1 = np.log(paper_line1[:, 0])
diff_all = boltzmann_line1 - paper_line1[:, 1]
bad_cells = np.isnan(boltzmann_line1) | (diff_all < 0.0)
bad_cells[0:3] = True
bad_cells[-1] = True
diff_100 = boltzmann_line1[~bad_cells] - paper_line1[~bad_cells, 1]
R_100_diff = paper_line1[~bad_cells, 0]
# diff_100 = np.append(diff_100, np.tile(diff_100[-1],50))
# R_100_diff = np.append(R_100_diff, np.tile(R_100_diff[-1],50))


diff_400 = np.append(np.tile(0,50), diff_400)
diff_300 = np.append(np.tile(0,50), diff_300)
diff_200 = np.append(np.tile(0,50), diff_200)
diff_100 = np.append(np.tile(0,50), diff_100)
R_400_diff = np.append(np.tile(1,50), R_400_diff)
R_300_diff = np.append(np.tile(1,50), R_300_diff)
R_200_diff = np.append(np.tile(1,50), R_200_diff)
R_100_diff = np.append(np.tile(1,50), R_100_diff)

plt.figure()
plt.plot(R_100_diff, diff_100, 'b--', label = "diff 100")
plt.plot(R_200_diff, diff_200, 'g--', label = "diff 200")
plt.plot(R_300_diff, diff_300, 'm--', label = "diff 300")
plt.plot(R_400_diff, diff_400, 'c--', label = "diff 400")
# expression = np.log(0.5*R400)/7
def fit_diff_func(R, a, b, c, d, e):
    # return a - a/(R/b + 1) 
    return a + b * R + c * R**2 + d * R**3 + e * R**4

fit_params = np.zeros((4, 5))
R_vals = [100, 200, 300, 400]

popt, pcov = curve_fit(fit_diff_func, R_100_diff, diff_100)
fit_params[0] = popt
plt.plot(R100, fit_diff_func(R100, *popt), 'b', label = "fits 400")

popt, pcov = curve_fit(fit_diff_func, R_200_diff, diff_200)
fit_params[1] = popt
plt.plot(R200, fit_diff_func(R200, *popt), 'g', label = "fits 200")

popt, pcov = curve_fit(fit_diff_func, R_300_diff, diff_300)
fit_params[2] = popt
plt.plot(R300, fit_diff_func(R300, *popt), 'm', label = "fits 300")

popt, pcov = curve_fit(fit_diff_func, R_400_diff, diff_400)
fit_params[3] = popt
plt.plot(R400, fit_diff_func(R400, *popt), 'c', label = "fits 400")


a_vals = fit_params[:, 0]
b_vals = fit_params[:, 1]
c_vals = fit_params[:, 2]
d_vals = fit_params[:, 3]
e_vals = fit_params[:, 4]
# print(R_vals)
# print(a_vals)
# print(b_vals)
# print(c_vals)
# print(d_vals)
# print(e_vals)

plt.xlabel('R')
plt.ylabel('$-e \phi / T_e$')
plt.legend()
plt.savefig('skovorodin-figure-2-diff.png')
plt.close()

plt.figure()
plt.subplot(2, 3, 1)

plt.plot(R_vals, a_vals, 'ko')
plt.plot([R_vals[0], R_vals[-1]], [a_vals[0], a_vals[-1]], 'k--')
plt.xlabel('R')
plt.title('Values of a')

plt.subplot(2, 3, 2)
plt.plot(R_vals, b_vals, 'ko')
plt.plot([R_vals[0], R_vals[-1]], [b_vals[0], b_vals[-1]], 'k--')
plt.xlabel('R')
plt.title('Values of b')

plt.subplot(2, 3, 3)
plt.plot(R_vals, c_vals, 'ko')
plt.plot([R_vals[0], R_vals[-1]], [c_vals[0], c_vals[-1]], 'k--')
plt.xlabel('R')
plt.title('Values of c')

plt.subplot(2, 3, 4)
plt.plot(R_vals, d_vals, 'ko')
plt.plot([R_vals[0], R_vals[-1]], [d_vals[0], d_vals[-1]], 'k--')
plt.xlabel('R')
plt.title('Values of d')

plt.subplot(2, 3, 5)
plt.plot(R_vals, e_vals, 'ko')
plt.plot([R_vals[0], R_vals[-1]], [e_vals[0], e_vals[-1]], 'k--')
plt.xlabel('R')
plt.title('Values of e')
plt.tight_layout()
plt.savefig('skovorodin-figure-2-fit-params.png')
