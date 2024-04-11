from cql3d_utils import cql3d_ncdf
from bound_phib import lossBound
from map_f import remap_f
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import interp1d

# Move to uniform grid
regrid = True

# Physical Constants
erg2kev = 1.6022e-9
clight = 2.99792458e10 #[cm/s]

# Read in the cql3d netcdf file and copy it ot the cql class
cql = cql3d_ncdf(cql_mnemonic='WHAM_max')

# Read in important dimensions
N_spec = cql.dim['gen_species_dim'] # number of species
rDim = int(cql.dim['rdim']) # number of flux surfaces
zDim = cql.dim['zdim'] # number of points along flux surface
jx = cql.dim['xdim']   # number of points in energy grid
iy = cql.dim['ydim']   # number of points in the pitch angle grid

# Read in important variables
v_norm = cql.var['vnorm']      # velocity normalization [cm/s]
sqPsiGrid = cql.var['rya']     # normalized sqrt poloidal flux grid has dim rdim
psiGrid = cql.var['equilpsi']  # poloidal flux at radial bin center
mid_f = cql.var['f']/v_norm**3 # distribution f(species, psi, u, theta) [cm^-3 cm/s^-3]
x = cql.var['x']               # normalized momentum per rest mass grid
y = cql.var['y']               # pitch angle grid [radians]
zGrid = cql.var['z']           # value of z along the field line [cm]
BdB0 = cql.var['bbpsi']        # value of B(z)/B(z=0) along field line
B0 = cql.var['bmidplne']
phi = cql.var['ephiz'][-1,:,:]        # value of electric potential along field line [kV]
q = cql.var['bnumb']                  # charge of the particle [e]
kspeci = cql.var['kspeci'][:,0,:]     # species names
m = cql.var['fmass']                  # mass [g]

# Create offmidplane f array
f_tmp = np.zeros((N_spec,rDim,zDim,jx,iy)) # output f at each point along the field line

# Local velocity grid
u_loc = np.zeros((N_spec,rDim,zDim,jx))
v_loc = np.zeros((N_spec,rDim,zDim,jx))

# Create uniform output grid for theta
ilh = int(zDim/2)
dtheta = np.pi/(iy-1) 
theta  = np.linspace(0,np.pi,iy)

# Velocity grid 
u0 = x*v_norm
u02 = u0*u0
gamma_rel = np.sqrt(1+u02)
v0 = u0/gamma_rel

# split problem up over mpi ranks in ipsi
comm = MPI.COMM_WORLD
mpisize = int(comm.Get_size())
iam = int(comm.Get_rank())
lstep = int(rDim/mpisize)
lstep1 = lstep+1
lremain = rDim%mpisize

f_flat = np.zeros(N_spec*rDim*zDim*jx*iy)     # flattened f for use as reduction target on root
if not regrid:
    uloc_flat= np.zeros(N_spec*rDim*zDim*jx)       
    vloc_flat = np.zeros(N_spec*rDim*zDim*jx)     

# Only allocate output arrays on root proc to save memory
if iam == 0: 
    f_out = np.zeros((N_spec,rDim,zDim,jx,iy)) # output f at each point along the field line
  
# Allocate bsegment and lsegment
# based on my simple mpi partitioning routine used for FORTRAN
bsegment = np.zeros(mpisize,dtype=int)
lsegment = np.zeros(mpisize,dtype=int)
mpiassign = np.zeros(rDim,dtype=int)
for iproc in range(mpisize):
    if iproc < lremain:
        bsegment[iproc] = lstep1 * iproc 
        lsegment[iproc] = lstep1 
    else:
        bsegment[iproc] = (lremain*lstep1) + lstep * (iproc-lremain)
        lsegment[iproc] = lstep
    lower = bsegment[iproc]
    upper = bsegment[iproc] + lsegment[iproc]
    if(upper>rDim):
        upper = rDim
    mpiassign[lower:upper] = iproc

# Loop over flux surface
f_tmp[:,:,0,:,:] = mid_f
for ipsi in range(bsegment[iam],bsegment[iam]+lsegment[iam]):
    # loop over species
    for ispec in range(N_spec):
        #this is equiv to the boundaries_phib routine in CQL3D
        vbndry1, thbndry1, sin02thbndry1, vbndry2, \
            thbndry2, vp2, zbnce, zbnce2, itype00 = lossBound(cql,ispec,ipsi)
        #print('ipsi, zbnce1/2 ', ipsi, ' ', zbnce, ' ', zbnce2)
        # loop over z location
        for il in range(1,zDim):
            # This is equiv. to cfpleg routine in CQL3D
            f_tmp[ispec,ipsi,il,:,:], u_loc[ispec,ipsi,il,:] = remap_f(cql,ispec,ipsi,il,vbndry1,sin02thbndry1,vbndry2,vp2,zbnce,zbnce2,itype00)
            utmp = u_loc[ispec,ipsi,il,:]
            #numerical kludge for plotting and coupling. this can cause issues if grid does not
            #have sufficient refinement
            #set negative values to zero then find first finite index in array
            utmp0 = np.copy(utmp)
            utmp0[utmp0<0.0]=0
            indx = np.nonzero(utmp0)[0][0]-1
            mask = np.insert(np.nonzero(utmp0),0,indx)
            utmp_crop = utmp0[mask]
            
            f_crop = f_tmp[ispec,ipsi,il,mask,:]
            #find point in baseline distribution function and linearly remap based on weight
            #from local velocity grid (imperfect, but not terrible)
            mid_f = cql.var['f'][ispec,ipsi,:,:]/v_norm**3
            wght = utmp[indx+1]/(utmp[indx+1] - utmp[indx])
            #this will slightly underestimate density in most cases but it's fine for now, better approach later
            f_crop[0] = (1-wght)*mid_f[indx,:]+wght*f_crop[1,:] 
            if regrid:
                for iloc in range(iy):
                    f_interp = interp1d(utmp_crop,f_crop[:,iloc],fill_value='extrapolate')
                    f_tmp[ispec,ipsi,il,:,iloc] = f_interp(u0/v_norm)
        #out of il loop
    #out of ispec loop
    print('finished ipsi', ipsi)
#out of ipsi loop

comm.Barrier() #this might not be necessary

#really terrible reduce of solution
#do this properly eventually with gather if faster comms desired
comm.Reduce(f_tmp.flatten(),f_flat)

if not regrid:
    comm.Reduce(u_loc.flatten(),uloc_flat)
    u_loc = uloc_flat.reshape((N_spec,rDim,zDim,jx))

    #right now i am not calculating this because i wanted to save memory 
    #comm.Reduce(v_loc.flatten(),vloc_flat)
    #v_loc = vloc_flat.reshape((Nspec,rDim,zDim,jx))
    
#prepare the output file
if (iam == 0):
    f_out = f_flat.reshape((N_spec,rDim,zDim,jx,iy))
    hf = h5py.File('cql3d_f_RZ.h5', 'w')
    hf.create_dataset('v_norm',data=v_norm)    #[cm/s]
    hf.create_dataset('psiGrid',data=psiGrid)  #[Gauss cm^2]
    hf.create_dataset('zGrid',data=zGrid[0,:]) #[cm]
    if not regrid:
        hf.create_dataset('uGrid',data=u_loc)      #[normalized]
        #hf.create_dataset('vGrid',data=v_loc)      #[cm/s]
    else:
        hf.create_dataset('uGrid',data=u0)
        hf.create_dataset('vGrid',data=v0)
    hf.create_dataset('thGrid',data=theta)     #[rad]
    hf.create_dataset('charge',data=q)         #[e]
    hf.create_dataset('mass',data=m)           #[g]
    hf.create_dataset('BdB0',data=BdB0)        #[normalized]
    hf.create_dataset('B0',data=B0)            #[Gauss]
    hf.create_dataset('phi',data=phi)          #[kV]
    hf.create_dataset('f_dist',data=f_out)     #[cm^-3 s^3]
    hf.close()