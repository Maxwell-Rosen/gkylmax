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
    plt.contour(grid[1], grid[0], psi[:,:,0], levels=clevels, colors="k")
    plt.contour(grid[1], grid[0], psi[:,:,0], levels=np.r_[psisep], colors="r")
    
    #clevels = np.linspace(psi_min-3*dpsi, psi_max, npsi+4)
    #plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")
    
    plt.plot(Z, R,marker=".", color="k", linestyle="none")
    # plt.scatter(R,Z, marker=".", )
    segs1 = np.stack((Z,R), axis=2)
    segs1_reduced = np.stack((Z[:, ::6],R[:, ::6]), axis=2)
    segs2 = segs1_reduced.transpose(1,0,2)
    plt.gca().add_collection(LineCollection(segs1_reduced, linewidths=1))
    plt.gca().add_collection(LineCollection(segs2, linewidths=1))
    ## Set the linestyle to be thin
    plt.grid()
    #plt.axis("tight")
    #plt.axis("image")
plt.ylim(0.0, .1)
plt.xlim(-2, 2)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.title("Simulation Node Locations", fontsize=26)
plt.ylabel("R, m", fontsize=22)
plt.xlabel("Z, m", fontsize=22)
plt.tight_layout()
plt.savefig("wham-nodes.png", dpi=1000)
plt.show()
