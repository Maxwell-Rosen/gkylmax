import postgkyl as pg
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.integrate as integrate
from scipy.interpolate import pchip_interpolate
import scipy.optimize as sco
import fortranformat as ff
from datetime import date
from typing import Any, Generator, Iterable, List, TextIO, Union

# Remember to delete the top line and that Gkeyll takes the _psi file which is made from the efit file reader

# Parameters for different mirror ratios R
R_values = [3, 5, 10, 15, 22, 32, 40, 50]
mcB_inner_values = [2.005804, 2.523435, 3.500993, 4.258079, 5.131404, 6.165174, 6.879326, 7.677658]
mcB_outer_values = [2.021283, 2.541062, 3.535280, 4.311313, 5.210360, 6.278592, 7.019136, 7.848735]
gamma_inner_values = [0.454541, 0.333259, 0.226795, 0.182798, 0.149642, 0.123304, 0.109952, 0.098086]
gamma_outer_values = [0.453590, 0.333783, 0.228416, 0.184758, 0.151760, 0.125465, 0.112109, 0.100215]

B0 = 0.0
R0 = .02

L_center = 4.0
L_plugs = 2.0
L_to_wall = 1.5

def psi_f(R, Z, mcB_inner, mcB_outer, gamma_inner, gamma_outer):
  return 0.5*np.power(R,2)*mcB_inner*( 1./(np.pi*gamma_inner*(1.+((Z-L_center/2)/gamma_inner)**2)) \
                                +1./(np.pi*gamma_inner*(1.+((Z+L_center/2)/gamma_inner)**2))) \
        +0.5*np.power(R,2)*mcB_outer*( 1./(np.pi*gamma_outer*(1.+((Z-(L_center/2+L_plugs))/gamma_outer)**2)) \
                                +1./(np.pi*gamma_outer*(1.+((Z+(L_center/2+L_plugs))/gamma_outer)**2)) )

#RZ box
NW = 513
NH = 513
RMIN,RMAX = 1e-3, .2 # The _zero is for Rmin=0, but 1e-3 was used for 1x cases
ZMIN,ZMAX = -(L_center/2 + L_plugs + L_to_wall), (L_center/2 + L_plugs + L_to_wall)
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

#Header stuff
header_fmt = "(a48,3i4)"
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

# Loop over all mirror ratios
for idx, (mcB_inner, mcB_outer, gamma_inner, gamma_outer, R_mirror) in enumerate(
    zip(mcB_inner_values, mcB_outer_values, gamma_inner_values, gamma_outer_values, R_values)
):
    outFileName = f'tandem_R{R_mirror}.geqdsk'
    
    print(f"\n--- Processing R = {R_mirror} ---")
    print(
        f"mcB_inner = {mcB_inner}, mcB_outer = {mcB_outer}, "
        f"gamma_inner = {gamma_inner}, gamma_outer = {gamma_outer}"
    )
    
    SIMAG = psi_f(2.0, 0, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
    SIBRY = psi_f(4.0, 0, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
    
    print("simag = %g"%SIMAG)
    print("sibry = %g"%SIBRY)
    
    #rthetagrid = np.zeros((len(Rgrid),len(Zgrid),2))
    psiRZ = np.zeros((len(Rgrid),len(Zgrid)))
    for i,Ri in enumerate(Rgrid):
        for j,Zj in enumerate(Zgrid):
            psiRZ[i,j] = psi_f(Ri, Zj, mcB_inner, mcB_outer, gamma_inner, gamma_outer)
    
    plt.figure()
    plt.contour(Rgrid,Zgrid, psiRZ.T, levels = (np.linspace(0.001,0.004,12)))
    plt.xlabel('R')
    plt.ylabel('Z')
    plt.colorbar()
    plt.title(f"Psi calculated in RZ coords (R={R_mirror})")
    
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
    
    #Now write the EFIT FILE
    with open(outFileName,'w',newline='') as f:
        write_line((comment, 3, NW, NH), f, header_fmt) #3 is idum
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
    
    print(f"Written: {outFileName}")

plt.show()
