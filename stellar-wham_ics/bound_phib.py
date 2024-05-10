#S. Frank Realta Fusion 01/2024

import numpy as np

#physical constants
erg2keV = 1.6022e-9
clight = 2.99792458e10 #[cm/s]

def lossBound(cql,ispec,ipsi):
    # Calculates the loss boundaries at a given psi and z location along the field
    # line based on the CQL3D-m simulation data. Based on the subroutine
    # boundaries_phib in the CQL3D-m code found in lossorbm.f
    # Inputs:
    # cql   - CQL3D-m output class containing the values found in the CQL3D netCDF
    # ipsec - species index of interest
    # ipsi  - psi index of interest
    #
    # Outputs:
    # vbndry1
    # thbndry1
    # sin02thbndry1
    # vbndry2
    # thbndry2
    # vp2mn
    # zbnce2
    #

    # Read in important dimensions
    N_spec = cql.dim['gen_species_dim'] # number of species
    rDim = cql.dim['rdim'] # number of flux surfaces
    zDim = cql.dim['zdim'] # number of points along flux surface
    jx = cql.dim['xdim']   # number of points in energy grid
    iy = cql.dim['ydim']   # number of points in the pitch angle grid

    # Read in important variables
    v_norm = cql.var['vnorm']      # velocity normalization [cm/s]
    sqPsiGrid = cql.var['rya']     # normalized sqrt poloidal flux grid has dim rdim
    psiGrid = cql.var['equilpsi']  # poloidal flux at radial bin center
    mid_f = cql.var['f']/v_norm**3 # distribution f(species, psi, u, theta) [cm^-3 cm/s^-3]
    x = cql.var['x']               # normalized momentum per rest mass grid
    y0 = cql.var['y'][ipsi,:]         # pitch angle grid [radians]
    z = cql.var['z'][ipsi,:]          # value of z along the field line [cm]
    BdB0 = cql.var['bbpsi'][ipsi,:]   # value of B(z)/B(z=0) along field line
    phi = cql.var['ephiz'][-1,ipsi,:] # value of electric potential along field line [kV]
    mid_f = cql.var['f'][ispec,ipsi,:,:]/v_norm**3   # distribution [cm^-3 cm/s^-3]
    q = cql.var['bnumb'][ispec]       # charge [e]
    m = cql.var['fmass'][ispec]       # mass [g]
    
    # Derived variables
    u0 = v_norm*x
    B0dB = 1.0/BdB0
    phi0 = phi[0]
    dPhi = phi0-phi
    z0 = z[0]
    zmax = z[-1]
    qm = q/m

    #set initial values for the boundaries
    vbndry1 = u0
    thbndry2 = y0
    thbndry1 = np.zeros(jx)
    sin02thbndry1 = np.zeros(jx)

    #minimum of 2*(q/m)*(Phi(0)-Phi(z))
    vp2 = 2.*qm*erg2keV*dPhi
    vp2mn = np.amin(vp2)

    # first boundary
    # sin^2(theta) = (B/B0)(v0^2/v^2)sin^2(theta0)
    # loop over v0>0 points (note: I always assume lossmode=mirrorz here) 
    for j in range(1,jx):
        v02 = vbndry1[j]**2
        if (v02+vp2mn) > 0.0: 
            BBVV = v02*BdB0/(v02+vp2)
            BBVVmx = np.amax(BBVV)
            km = np.argmax(BBVV)
            if BBVVmx < 1.0:
                BBVVmx = 1.0
            sin02thbndry1[j] = 1.0/BBVVmx
            if (sin02thbndry1[j] <=1.0):
                thbndry1[j] = np.arcsin(np.sqrt(sin02thbndry1[j]))
        else:
            thbndry1[j] = 0.0    
            sin02thbndry1[j] = 0.0

    
    # second boundary
    # v^2 = v0^2 + vp2
    km = np.argmin(abs(vp2-vp2mn))
    zbnce2 = z[km] - z0
    if vp2mn < 0:
        vbndry2 = np.sqrt(-vp2mn)
    else:
        vbndry2 = 0.0

    itype00 = np.full((iy,jx),1)
    zbnce = np.full((iy,jx),1.0)
    for jo in range(1,jx):
        v02 = u0[jo]**2
        if (vp2mn<0.0) and (v02+vp2mn <= 0.0):
            itype00[:,jo] = -1 #trapped
            km = np.argmin(abs(vp2-vp2mn))
            zz0 = z[km] - z0
            zbnce[:,jo] = zz0
        else:
            BBVV = v02*BdB0/(v02+vp2)
            BBVVmx = np.amax(BBVV)
            km = np.argmax(BBVV)
            if(BBVVmx<1.0):
                BBVVmx=1.0
            zz0 = z[km]-z0
            for io in range(0,iy):
                sin2thet = BBVVmx*np.sin(y0[io])**2
                if (sin2thet>1.0):
                    itype00[io,jo] = -1 #trapped
                    zbnce[io,jo]=zz0 #point where particle bounced
                else:
                    itype00[io,jo] = 1 #passing
                    zbnce[io,jo] = zmax - z0 #passing particle reached edge
    return vbndry1, thbndry1, sin02thbndry1, vbndry2, thbndry2, vp2, zbnce, zbnce2, itype00


