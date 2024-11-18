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
simNames = ['gk_mirror_uniform', 'gk_mirror_nonuniform_128']
for simName in simNames:
    dir = "./"
    data = pg.GData(str(simName + '-nodes.gkyl'))
    vals = data.get_values()
    print(np.shape(vals))
    R = vals[:,0]
    Z = vals[:,2]

    
    # psid = pg.GData("data-lores/wham_psi.gkyl")
    # interp = pg.GInterpModal(psid,1,"ms")
    # grid, psi = interp.interpolate()
    # #pg.output.plot(psid, contour=True)
    
    # for d in range(len(grid)):
    #     grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])
    
    # psisep= 0.0026
    # psi_min = 0.001
    # psi_max =  0.003
    # npsi = 3
    # dpsi = (psi_max-psi_min)/npsi
    
    # clevels = np.linspace(psi_min, psi_max, npsi)
    # plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="k")
    # plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=np.r_[psisep], colors="r")
    
    #clevels = np.linspace(psi_min-3*dpsi, psi_max, npsi+4)
    #plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")
    
    plt.plot(Z,R,marker=".", linestyle="none", label = simName)
    plt.scatter(Z,R, marker=".")
    # segs1 = np.stack((R,Z), axis=1)
    # segs2 = segs1.transpose(1,2)
    # plt.gca().add_collection(LineCollection(segs1))
    # plt.gca().add_collection(LineCollection(segs2))
    plt.grid()
    #plt.axis("tight")
    #plt.axis("image")
    plt.xlabel("Z (m)")
    plt.ylabel("R (m)")
plt.legend()
plt.savefig("node-locs.png", dpi=300)
plt.show()
