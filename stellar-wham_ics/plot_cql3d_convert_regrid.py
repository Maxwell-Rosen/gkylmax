import h5py
import matplotlib.pyplot as plt
import numpy as np

iy = 400
jx = 200
lz = 128
ispec = 0
isurf = 3
ulim = 0.015

theta  = np.linspace(0,np.pi,iy)
hf = h5py.File('cql3d_f_RZ.h5','r')
v_norm = hf['v_norm']
f_dist = hf['f_dist']
u0 = hf['uGrid']

print(np.shape(f_dist))
print(np.shape(u0))

upar = np.zeros((jx,iy))
uprp = np.zeros((jx,iy))
for j in range(jx):
    upar[j,:] = u0[j]*np.cos(theta)/v_norm
    uprp[j,:] = u0[j]*np.sin(theta)/v_norm

fig, axs = plt.subplots(1, 5, figsize=(15, 3))

cont_lvls = np.linspace(np.amin(f_dist[ispec,1,1,:,:]),np.amax(f_dist[ispec,1,1,:,:]),100)

im = axs[0].contourf(upar,uprp,f_dist[ispec,1,1,:,:],levels=cont_lvls)
for i in range(len(axs)):
    axs[i].set_xlim(-ulim,ulim)
    axs[i].set_ylim(0.0,ulim)
axs[1].contourf(upar,uprp,f_dist[ispec,isurf,30,:,:],levels=cont_lvls)
axs[2].contourf(upar,uprp,f_dist[ispec,isurf,90,:,:],levels=cont_lvls)
axs[3].contourf(upar,uprp,f_dist[ispec,isurf,100,:,:],levels=cont_lvls)
axs[4].contourf(upar,uprp,f_dist[ispec,isurf,110,:,:],levels=cont_lvls)
fig.colorbar(im, ax=axs[4])
plt.show()