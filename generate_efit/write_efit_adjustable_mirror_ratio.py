import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate, interp1d
from scipy.special import ellipk, ellipe
import scipy.constants
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union

#Remove all files with the name conditioned_coil_R*.geqdsk
import os
import multiprocessing as mp

for file in os.listdir():
    if file.endswith(".geqdsk") and file.startswith("conditioned_coil_R"):
        os.remove(file)

## TODO: Control for minimum B = constant, then set Bmax = R * Bmin

def psi_f(R, Z, I_high_field_coils, I_low_field_coils, coil_dist):
    if R < 1e-17:
        R = 1e-17
    a = 0.45 #Radius of all the coils. Must be > R_max
    mu0 = scipy.constants.mu_0

    Zl = Z - coil_dist/2
    k2 = 4*a*R/((a+R)**2 + Zl**2)
    k = np.sqrt(k2)
    Aphi = -mu0*I_high_field_coils/np.pi * np.sqrt(a/R) * ((k2 - 2)/2/k*ellipk(k2) + ellipe(k2)/k)
    Zu = Z + coil_dist/2
    k2 = 4*a*R/((a+R)**2 + Zu**2)
    k = np.sqrt(k2)
    Aphi += -mu0*I_high_field_coils/np.pi * np.sqrt(a/R) * ((k2 - 2)/2/k*ellipk(k2) + ellipe(k2)/k)
    Z_middle = np.linspace(Zl,Zu,50)
    Z_middle = Z_middle[1:-1]
    k2 = 4*a*R/((a+R)**2 + Z_middle**2)
    k = np.sqrt(k2)
    Aphi += np.sum(-mu0*I_low_field_coils/np.pi * np.sqrt(a/R) * ((k2 - 2)/2/k*ellipk(k2) + ellipe(k2)/k))
    return Aphi * R


B0 = 0.0 #Domain center
R0 = .02

mirror_ratio_vector = np.array([20])#, 15, 10, 5, 2, 1.5])

I_high_field = 812998 # High field current. Optimized such that Bmax on psi_eval = aim_Bmax
I_middle_low_field = 0 # Low field current. Optimized such that Bmin on psi_eval = aim_Bmin
h = 2 

#optimizer parameters
aim_Bmax = 10.0
aim_Bmin = 0.5
I_previous_step = 0
Bmax_previous_step = 0
psi_eval = 1e-3


    #RZ box
NW = 257
NH = 257
RMIN,RMAX = .001, 0.4
ZMIN,ZMAX = -2, 2
RDIM = RMAX - RMIN
ZDIM = ZMAX - ZMIN
RLEFT = RMIN
ZMID = (ZMAX+ZMIN)/2.0
RMAXIS = 0.0
ZMAXIS = 0.0
NPSI = NW #don't write
RCENTR = 0.0
BCENTR = B0
CURRENT = 0


#Solve GS in RZ coords
Rgrid = np.linspace(RMIN,RMAX,NW)
Zgrid = np.linspace(ZMIN,ZMAX,NH)
dR = Rgrid[1] - Rgrid[0]
dZ = Zgrid[1] - Zgrid[0]
#rthetagrid = np.zeros((len(Rgrid),len(Zgrid),2))
psiRZ = np.zeros((len(Rgrid),len(Zgrid)))
BzRZ = np.zeros((len(Rgrid),len(Zgrid)))
BrRZ = np.zeros((len(Rgrid),len(Zgrid)))
BRZ = np.zeros((len(Rgrid),len(Zgrid)))

B_psi = np.zeros((len(mirror_ratio_vector),len(Zgrid)))

def minimize_function(I_high_field, I_middle, coil_distance):



    while True:
        # # From Jackson. There seems to be a typo in the k**2 in the eliptic functions. They should be k**2, not k
        print("I = %g"%I_high_field)
        def compute_psi_B(i, Ri, Zgrid, dR, dZ, I):
            psiRZ_row = np.zeros(len(Zgrid))
            BzRZ_row = np.zeros(len(Zgrid))
            BrRZ_row = np.zeros(len(Zgrid))
            BRZ_row = np.zeros(len(Zgrid))
            for j, Zj in enumerate(Zgrid):
                psiRZ_row[j] = psi_f(Ri, Zj, I)
                if Ri < 1e-17:
                    Ri = Rgrid[1]
                BzRZ_row[j] = 1/Ri**2 * (psi_f(Ri+dR, Zj, I) - psi_f(Ri-dR, Zj, I))/(2*dR)
                BrRZ_row[j] = -1/Ri * (psi_f(Ri, Zj+dZ, I) - psi_f(Ri, Zj-dZ, I))/(2*dZ)
                BRZ_row[j] = np.sqrt(BrRZ_row[j]**2 + BzRZ_row[j]**2)
            return i, psiRZ_row, BzRZ_row, BrRZ_row, BRZ_row

        with mp.Pool(mp.cpu_count()) as pool:
            results = pool.starmap(compute_psi_B, [(i, Ri, Zgrid, dR, dZ, I_high_field) for i, Ri in enumerate(Rgrid)])

        for i, psiRZ_row, BzRZ_row, BrRZ_row, BRZ_row in results:
            psiRZ[i, :] = psiRZ_row
            BzRZ[i, :] = BzRZ_row
            BrRZ[i, :] = BrRZ_row
            BRZ[i, :] = BRZ_row
        
        for j, Zj in enumerate(Zgrid):
            # For this position of Z, find which R gives psi = psi_eval
            psiR = interp1d(Rgrid, psiRZ[:,j], kind='cubic')
            R_psi = sco.brentq(lambda R: psiR(R) - psi_eval, RMIN, RMAX)
            psi_interp = psi_f(R_psi, Zj, I_high_field)
            #Interpolate the magnetic field at this point
            B_psi_interp_model = interp1d(Rgrid, BRZ[:,j], kind='cubic')
            B_psi[i_current,j] = B_psi_interp_model(R_psi)
            


        Bmax = np.max(B_psi[i_current,:])
        Bmin = B_psi[i_current,int(NH/2)]
        print("Bmax = %g at (R = %g, Z = %g)"%(Bmax, R_psi, Zgrid[np.argmax(B_psi[i_current,:])]))
        print("Bmin = %g at (R = %g, Z = %g)"%(Bmin, Rgrid[int(NW/2)], Zgrid[int(NH/2)]))
        print("Mirror ratio = %g"%(Bmax/Bmin))

        if np.isclose(Bmax, aim_Bmax, rtol=1e-2):
            a_fit = (Bmin - Bmin_previous_step)/(I_middle_low_field - I_middle_low_field_previous_step)
            b_fit = Bmin - a_fit*I_middle_low_field
            I_middle_low_field_previous_step = I_middle_low_field
            Bmin_previous_step = Bmin
            I_middle_low_field = (aim_Bmax - b_fit)/a_fit

            break

        a_fit = (Bmax - Bmax_previous_step)/(I_high_field - I_previous_step)
        b_fit = Bmax - a_fit*I_high_field
        I_previous_step = I_high_field
        Bmax_previous_step = Bmax
        I_high_field = (aim_Bmax - b_fit)/a_fit

    mirrorRatio = Bmax/Bmin
    outFileName = 'conditioned_coil_R'+str(int(Bmax/Bmin))+'.geqdsk'
    print("Making EFIT file with Bmax = %g, Bmin = %g, mirror ratio = %g"%(Bmax, Bmin, mirrorRatio))
    print("Filename: %s"%outFileName)

    Bmin = BRZ[0,int(NH/2)]

    SIMAG = psi_f(2.0,0, I_high_field) 
    SIBRY = psi_f(4.0,0, I_high_field)

    # plt.figure()
    # plt.plot(Zgrid, B_psi)
    # plt.xlabel('Z')
    # plt.ylabel('B')
    # plt.title("B calculated in RZ coords at psi = "+str(psi_eval))


    # plt.figure()
    # contour = plt.contour(Rgrid,Zgrid, psiRZ.T, levels = np.linspace(0, 1e-4, 50))
    # plt.xlabel('R')
    # plt.ylabel('Z')
    # plt.colorbar()
    # plt.title("Psi calculated in RZ coords")
    # # plt.show()

    # plt.figure()
    # plt.pcolormesh(Rgrid,Zgrid, np.log(BRZ.T))
    # plt.xlabel('R')
    # plt.ylabel('Z')
    # plt.colorbar()
    # plt.title("log(B) calculated in RZ coords")
    # # plt.show()

    # plt.figure()
    # plt.plot(Zgrid, BRZ[int(NW/2),:])
    # plt.xlabel('Z')
    # plt.ylabel('B')
    # plt.title("B calculated in RZ coords at R = "+str(Rgrid[int(NW/2)]))
    # # plt.show()

    # plt.figure()
    # plt.plot(Zgrid, BRZ[0,:])
    # plt.xlabel('Z')
    # plt.ylabel('B')
    # plt.title("B calculated in RZ coords at R = "+str(Rgrid[0]))

    # plt.figure()
    # plt.plot(Rgrid, BRZ[:,int(NH/4)], label='Z = '+str(Zgrid[int(NH/4)]))
    # plt.plot(Rgrid, BRZ[:,int(NH/4)-1], label='Z = '+str(Zgrid[int(NH/4)-1]))
    # plt.plot(Rgrid, BRZ[:,int(NH/4)+1], label='Z = '+str(Zgrid[int(NH/4)+1]))
    # plt.xlabel('R')
    # plt.ylabel('B')
    # plt.legend()
    # plt.title("B calculated in RZ coords around Z = "+str(Zgrid[int(NH/4)]))

    # plt.figure()
    # plt.plot(Rgrid, BRZ[:,int(NH/2)], label='Z = '+str(Zgrid[int(NH/2)]))
    # plt.xlabel('R')
    # plt.ylabel('B')
    # plt.legend()
    # plt.title("B calculated in RZ coords around Z = "+str(Zgrid[int(NH/2)]))

    #PSI quantities
    PSIGRID = np.linspace(SIMAG, SIBRY,NPSI)
    FPOL = (B0*R0/Rgrid)*Rgrid # F = RBphi
    FFPRIM = np.repeat(0.0, NPSI)
    PPRIME = np.repeat(-1e-6,NPSI)
    PRES = integrate.cumulative_trapezoid(PPRIME,PSIGRID,initial=0)
    PSIZR = psiRZ.T


    QPSI = np.zeros(NPSI)
    for i in range(NPSI):
        QPSI[i] = 0



    writeList = [NW, NH,                                        #3i4
                RDIM, ZDIM, RCENTR, RLEFT, ZMID,               #5E16.9
                RMAXIS, ZMAXIS, SIMAG, SIBRY, BCENTR,          #5e16.9
                CURRENT, SIMAG, 0, RMAXIS, 0,                  #5E16.9
                ZMAXIS, 0, SIBRY, 0, 0,                        #5E16.9
                FPOL, PRES, FFPRIM, PPRIME,                    #5E16.9
                PSIZR, QPSI]                                   #5E16.9

    #Header stuff
    header_fmt = "(3i4)"
    label = 'FREEGS'
    creation_date = date.today().strftime("%d/%m/%Y")
    shot = int(0)
    time = int(0)
    shot_str = f"# {shot:d}"
    time_str = f"  {time:d}ms"
    comment = f"{label:11}{creation_date:10s}   {shot_str:>8s}{time_str:16s}"
    def write_line(data: Iterable[Any], fh: TextIO, fmt: str) -> None:
        r"""
        Writes to a Fortran formatted ASCII data file. The file handle will be left on a
        newline.

        Parameters
        ---------
        data:
            The data to write.
        fh:
            File handle. Should be in a text write mode, i.e. ``open(filename, "w")``.
        fmt:
            A Fortran IO format string, such as ``'(6a8,3i3)'``.
        """
        fh.write(ff.FortranRecordWriter(fmt).write(data))
        fh.write("\n")
    #write_line((comment, 3, NW, NH), fh, header_fmt) #3 is idum


    #Now write the EFIT FILE
    with open(outFileName,'w',newline='') as f:
        ##NW,NH
        ##for i in range(2):
        ##    f.write('%d '%writeList[i])
        #f.write('FREEGS     19/06/2023        # 0  0ms              3  80  90')
        #f.write('\n')
        write_line((NW, NH), f, header_fmt) #3 is idum
        # rdim,zdim,rcentr,rleft,zmid
        for i in range(2,7):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        # rmaxis,zmaxis,simag,sibry,bcentr
        for i in range(7,12):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        # current, simag, xdum, rmaxis, xdum
        for i in range(12,17):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        #zmaxis,xdum,sibry,xdum,xdum
        for i in range(17,22):
            f.write('%16.9E'%writeList[i])
        f.write('\n')
        #FPOL,PRES,FFPRIM,PPRIME
        for i in range(22,26):
            count=0
            for j in range(0,NW):
                f.write('%16.9E'%writeList[i][j])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0
        #PSIZR
        for i in range(26,27):
            count = 0
            for j in range(0,NH):
                for k in range(0,NW):
                    f.write('%16.9E'%writeList[i][j][k])
                    count = count+1
                    if count==5:
                        f.write('\n')
                        count=0
        #QPSI
        for i in range(27,28):
            count=0
            for j in range(0,NW):
                f.write('%16.9E'%writeList[i][j])
                count = count+1
                if count==5:
                    f.write('\n')
                    count=0



for i in range(len(Imiddle_vec)):
    plt.plot(Zgrid, B_psi[i, :], label=f'I = {Imiddle_vec[i]}')
plt.xlabel('Z')
plt.ylabel('B')
plt.title("B calculated in RZ coords at psi = "+str(psi_eval))
plt.legend()
plt.tight_layout()
plt.show()
    # plt.show()
