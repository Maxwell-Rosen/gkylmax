import numpy as np
import postgkyl as pg
import matplotlib.pyplot as plt


# simName='outputs/gk_wham_nonunif_tanmap_30vth_b1,4
simName="outputs/gk_wham_nonunif_tanmap_30vth_b1,4"

species='ion'
#pos integ moms
idata = pg.GData('%s-%s_positivity_shift_integrated_moms.gkyl'%(simName, species))
itime = idata.get_grid()[0]
vals = idata.get_values()
pm0 = vals[:,0]
pm2 = vals[:,2] + vals[:,3]
#f integ moms
idata = pg.GData('%s-%s_integrated_moms.gkyl'%(simName, species))
itime = idata.get_grid()[0]
vals = idata.get_values()
m0 = vals[:,0]
m2 = vals[:,2] + vals[:,3]
#source integ moms
# idata = pg.GData('%s-%s_source_integrated_moms.gkyl'%(simName, species))
# itime = idata.get_grid()[0]
# vals = idata.get_values()
# sm0 = vals[:,0]
# sm2 = vals[:,2] + vals[:,3]


#Now the tricky part we need to get dt from slurm at the same timestamps as integ moms
fname='slurm-443060.out'
f = open(fname,'r')
num_lines=0
for line in f:
    split_line = line.split("dt =")
    if(len(split_line))>1:
        num_lines+=1
f.close()

dt = np.zeros(num_lines)
time = np.zeros(num_lines)
i=0
f=open(fname,'r')
for line in f:
    split_line = line.split("dt =") 
    if(len(split_line))<2:
        continue
    else:
        split0 = split_line[0].split("t =")[1].split("...")
        time[i] = float(split0[0])
        dt[i] = float(split_line[1])
        i+=1
f.close()
#Now we want the corresponding indices to the integ mom data
inds = np.zeros(len(itime), dtype="int")
icoarse = 0
for ifine, t in enumerate(time):
    if icoarse == len(itime):
        break
    if t > itime[icoarse]:
        inds[icoarse] = ifine
        icoarse+=1

coarse_time = time[inds]
coarse_dt = dt[inds]





#Make the Plots
fig, ax  = plt.subplots(2,2)
# ax[0,0].plot(coarse_time, pm0/(coarse_dt)/sm0)
# ax[0,0].set_ylabel(r'$\frac{\langle \Delta f\rangle / \Delta t}{\langle S \rangle}$', fontsize=20)
# ax[1,0].plot(coarse_time, pm2/(coarse_dt)/sm2)
# ax[1,0].set_ylabel(r'$\frac{\langle v^2 \Delta f\rangle / \Delta t}{\langle v^2 S \rangle}$', fontsize=20)

dm0 = np.diff(m0)
dm2 = np.diff(m2)
dpm0 = np.diff(pm0)
dpm2 = np.diff(pm2)

ax[0,0].plot(itime[1:], dpm0)
ax[0,0].set_ylabel(r'$\frac{\langle \Delta f\rangle / \Delta t}{\langle f \rangle / \Delta t}$', fontsize=20)
ax[1,0].plot(itime[1:], dpm2)
ax[1,0].set_ylabel(r'$\frac{\langle v^2 \Delta f\rangle / \Delta t}{\langle v^2 f \rangle / \Delta t}$', fontsize=20)

ax[0,1].plot(itime[1:], dm0)
ax[0,1].set_ylabel(r'$\frac{\langle \Delta f\rangle}{\langle f \rangle}$', fontsize=20)
ax[1,1].plot(itime[1:], dm2)
ax[1,1].set_ylabel(r'$\frac{\langle v^2 \Delta f\rangle }{\langle v^2 f \rangle}$', fontsize=20)

# ax[0,1].plot(itime, pm0/m0)
# ax[0,1].set_ylabel(r'$\frac{\langle \Delta f\rangle}{\langle f \rangle}$', fontsize=20)
# ax[1,1].plot(itime, pm2/m2)
# ax[1,1].set_ylabel(r'$\frac{\langle v^2 \Delta f\rangle }{\langle v^2 f \rangle}$', fontsize=20)

for axi in ax.flatten():
    axi.set_xlabel(r'$t$', fontsize=20)

fig.suptitle("Positivity Source Rate Compared to IZ Source Rate (left) and Positivity Shift Compared to Distribution Function (right)")
fig.tight_layout()
plt.show()