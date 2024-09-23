import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import postgkyl as pg

#simName='step_outboard_shaped_plate'
#simName='step_outboard_horizontal_plate'
#simName='step_outboard_fixed_z'
#simName='stepbry'
#simName='stepcore'

plt.figure()
simNames = ['gk_wham']
for simName in simNames:
    dir = "./"
    data = pg.GData("gk_wham-nodes.gkyl")
    vals = data.get_values()
    R = vals[:,:,0]
    Z = vals[:,:,1]

    
    psid = pg.GData("data-lores/wham_psi.gkyl")
    interp = pg.GInterpModal(psid,1,"ms")
    grid, psi = interp.interpolate()
    #pg.output.plot(psid, contour=True)
    
    for d in range(len(grid)):
        grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])
    
    psisep= 0.003
    psi_min = 0.001
    psi_max =  0.003
    npsi = 3
    dpsi = (psi_max-psi_min)/npsi
    
    clevels = np.linspace(psi_min, psi_max, npsi)
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="k")
    plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[psisep], colors="r")
    
    #clevels = np.linspace(psi_min-3*dpsi, psi_max, npsi+4)
    #plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")
    
    plt.plot(R,Z,marker=".", color="k", linestyle="none")
    # plt.scatter(R,Z, marker=".", )
    segs1 = np.stack((R,Z), axis=2)
    segs1_reduced = np.stack((R[:, ::6],Z[:, ::6]), axis=2)
    segs2 = segs1_reduced.transpose(1,0,2)
    plt.gca().add_collection(LineCollection(segs1_reduced, linewidths=1))
    plt.gca().add_collection(LineCollection(segs2, linewidths=1))
    ## Set the linestyle to be thin
    plt.grid()
    #plt.axis("tight")
    #plt.axis("image")
plt.xlim(0.0, .1)
plt.ylim(-2, 2)
plt.title("WHAM nodes")
plt.xlabel("R, m")
plt.ylabel("Z, m")
plt.savefig("wham-nodes.png", dpi=1000)
plt.show()
