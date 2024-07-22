import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import postgkyl as pg

plt.figure()
simNames = ['stepouter', 'steplower', 'stepcore', 'stepcorel', 'stepupper', 'steppflor', 'steppflol', 'stepinnerupper', 'stepinnermiddle', 'stepinnerlower', 'steppfupr', 'steppfupl']
#simNames = ['steppfupl', 'steppfupr']
for simName in simNames:
    dir = "./"
    data = pg.GData(dir+"xyz"+simName+"_nodes.gkyl")
    vals = data.get_values()
    X = vals[:,0,:,0]
    Y = vals[:,0,:,1]
    Z = vals[:,0,:,2]
    R=np.sqrt(X**2+Y**2)
    
    psid = pg.GData("step_psi.gkyl")
    #psid = pg.GData(simName+"_psi.gkyl")
    interp = pg.GInterpModal(psid,2,"mt")
    grid, psi = interp.interpolate()
    #pg.output.plot(psid, contour=True)
    
    for d in range(len(grid)):
        grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])
    
    psisep= 1.50982
    psisep = 1.5093073363672997
    psi_min = 1.3
    psi_max = psisep
    npsi = 3
    dpsi = (psi_max-psi_min)/npsi
    
    clevels = np.linspace(psi_min, psi_max, npsi)
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="k")
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[psisep], colors="r")
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[1.8], colors="g")
    
    #clevels = np.linspace(psi_min-3*dpsi, psi_max, npsi+4)
    #plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")
    
    plt.plot(R,Z,marker=".", color="k", linestyle="none")
    plt.scatter(R,Z, marker=".")
    segs1 = np.stack((R,Z), axis=2)
    segs2 = segs1.transpose(1,0,2)
    plt.gca().add_collection(LineCollection(segs1))
    plt.gca().add_collection(LineCollection(segs2))

plt.grid()
#plt.axis("tight")
plt.xlabel("R")
plt.ylabel("Z")
#plt.axis("image")
plt.show()
