import numpy as np

# Physical Constants
erg2kev = 1.6022e-9
clight = 2.99792458e10 #[cm/s]

def remap_f(cql,ispec,ipsi,il,vbndry1,sin02thbndry1,vbndry2,vp2,zbnce,zbnce2,itype00):
    # copy dims from the cql class
    jx = cql.dim['xdim']           # number of points in energy grid
    iy = cql.dim['ydim']           # number of points in the pitch angle grid

    # copy vars from the cql class
    v_norm = cql.var['vnorm']         # velocity normalization
    x = cql.var['x']                  # normalized momentum per rest mass grid
    y = cql.var['y'][ipsi,:]          # pitch angle grid [radians]
    z = cql.var['z'][ipsi,:]          # value of z along the field line [cm]
    BdB0 = cql.var['bbpsi'][ipsi,:]   # value of B(z)/B(z=0) along field line
    phi = cql.var['ephiz'][-1,ipsi,:] # value of electric potential along field line [kV]
    mid_f = cql.var['f']/v_norm**3    # distribution [cm^-3 cm/s^-3]
    q = cql.var['bnumb'][ispec]       # charge [e]
    
    # Derived varibles from cql3d netcdf vars
    dPhi = phi - phi[0]
    iyh = int(iy/2)
    
    # Create grids
    dtheta = np.pi/(iy-1) 
    theta  = np.linspace(0,np.pi,iy)
    u0 = x*v_norm
    
    #initialize some arrays needed later on
    itype = np.zeros((iy,jx),dtype=int)
    floc_rz = np.zeros((jx,iy))
    uloc_z = np.zeros((iy,jx))
    thet_bnd1 = np.zeros(jx)
    uloc_bnd1 = np.zeros(jx)
    wk_f = mid_f[ispec,ipsi,:,:]

    #initialize local scalar balues
    zz0 = z[il] - z[0]
    BdB0_loc = BdB0[il]
    dPhi_loc = dPhi[il]
    vp2_loc = vp2[il]
    vp2mn = np.amin(vp2)
    uloc_bnd2 = np.sqrt(max(0.0,vp2_loc-vp2mn))

    # loop over theta
    for iloc in range(iy):
        theta_loc = theta[iloc]
        sin_loc = np.sin(theta_loc)
        sin_loc2 = sin_loc*sin_loc
        
        #special treatment for v0=0 
        jloc=0
        jo = jloc
        if (vp2_loc >= 0) and (zz0 <= zbnce2):
            v_loc = np.sqrt(vp2_loc)
            sin2thet_bnd = 0.0
            thet_bnd1[jo] = 0.0 
            uloc_bnd1[jo] = v_loc 
            uloc_z[jo] = v_loc
            #the way yuri treats this now in cql3d the if here is redundant
            #but we may want it for cross code coupling
            if iloc == 0:
                itype[iloc,jloc] = 1 #passing
                io = iloc
                floc = wk_f[jo,io] #note indices reversed on our wk_f from yuri's
                floc_rz[jloc,iloc] = floc
            elif iloc == iy-1:
                itype[iloc,jloc] = 1 #passing
                io = iloc
                floc = wk_f[jo,io] 
                floc_rz[jloc,iloc] = floc
            else:
                itype[iloc,jloc] = 1 #passing
                io = iloc
                floc = wk_f[jo,io] 
                floc_rz[jloc,iloc] = floc
        else:
            uloc_z[jloc] = 0.0
            floc = 0.0
            floc_rz[jloc,iloc] = floc
            itype[iloc,jloc] = -1 #disallowed
            sin2thet_bnd=0.0
            thet_bnd1[jo]=0.0
            if vp2_loc > 0:
                uloc_bnd1[jo]=np.sqrt(vp2_loc-vp2mn)
            else:
                uloc_bnd1[jo]=0.0

        #remaining cases for v0 > 0
        for jo in range(1,jx):
            v02 = u0[jo]**2
            if ((v02+vp2_loc<0.0) or ((u0[jo]<=vbndry2) and (zz0>zbnce2))):
                jloc = jo
                uloc_z[jloc] = 0.0
                floc = 0.0
                floc_rz[jloc,iloc]=floc
                itype[iloc,jloc]=-1 #yuri's code sets this to zero at all pitch angles (not necessary as already done in loop?)
                thet_bnd1[jo] = 0.0
                if v02+vp2_loc > 0:
                    uloc_bnd1[jo]=np.sqrt(vp2_loc-vp2mn)
                else:
                    uloc_bnd1[jo]=0.0
            else:
                uloc2 = v02+vp2_loc
                v02_sin02_bndry1 = sin02thbndry1[jo]*vbndry1[jo]**2
                sin2thet_bnd = (BdB0_loc/uloc2)*v02_sin02_bndry1
                if (sin2thet_bnd>1) and (sin2thet_bnd<=1.002):
                    sin2thet_bnd = 1.0
                thet_bnd1[jo] = np.arcsin(np.sqrt(sin2thet_bnd))
                uloc_bnd1[jo] = np.sqrt(uloc2)

                v2_v02 = (v02+vp2_loc)/v02
                sinn02 = v2_v02*sin_loc2/BdB0_loc
                if sinn02>1.0:
                    v_loc = np.sqrt(v02+vp2_loc)
                    jloc = jo
                    uloc_z[jloc] = v_loc
                    itype[iloc,jloc] = 0
                    io = iyh #set to electron distribution at 90 deg
                    if q == -1: #allow locally trapped e-
                        floc = wk_f[jo,io] #figure out what iyh is
                    else: #ions
                        floc = 0.0
                    floc_rz[jloc,iloc] = floc
                else:
                    v_loc = np.sqrt(v02+vp2_loc)
                    jloc = jo
                    uloc_z[jloc] = v_loc
                    theta0 = np.arcsin(np.sqrt(sinn02))
                    io = np.absolute(y[0:iyh]-theta0).argmin()
                    if theta_loc>0.5*np.pi:
                        theta0 = np.pi - theta0
                        io = iy-io-1
                    if (itype00[io,jo]==-1) and (zz0>zbnce[io,jo]):
                        itype[iloc,jloc]=3
                    elif (theta_loc > thet_bnd1[jo]) \
                         and (theta_loc < np.pi- thet_bnd1[jo]) \
                         and (zz0 > zbnce[io,jo]):
                        itype[iloc,jloc]=3
                    else:
                        itype[iloc,jloc]=1
                    floc = wk_f[jo,io]
                    floc_rz[jloc,iloc] = floc
        #out of jx loop
    #out of iy loop
    
    return floc_rz
