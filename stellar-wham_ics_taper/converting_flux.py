from cql3d_utils import cql3d_ncdf

# Read in the cql3d netcdf file and copy it ot the cql class
cql = cql3d_ncdf(cql_mnemonic='WHAM_max')

# Read in important variables
sqPsiGrid = cql.var['rya']     # normalized sqrt poloidal flux grid has dim rdim
psiGrid = cql.var['equilpsi']  # poloidal flux at radial bin center
psiLim = cql.var['psilim']     # poloidal flux at the limiter

#grids
normPsi1 = sqPsiGrid**2
normPsi2 = psiGrid/psiLim

#unit conversions to psi(0.0) = 0 grid in MKS
psiLimMKS = abs(psiLim - psiGrid[0]) 
psiMKS = psiLimMKS*normPsi2

#print out grids
print('Normalized Psi Grid by 2 different methods')
print(normPsi1)
print(normPsi2)
print('')
print('MKS')
print('Limiter Flux')
print('psiLim ',psiLimMKS)
print('psi ', psiMKS/1e8)