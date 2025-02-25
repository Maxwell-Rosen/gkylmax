import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import postgkyl as pg

data = pg.GData("Geometry/gk_wham-nodes.gkyl")
vals = data.get_values()
R = vals[:,0]
Z = vals[:,1]


plt.figure()
plt.grid()
plt.plot(R,Z, marker=".", linestyle='none', color="k")
# plt.scatter(R,Z, marker=".")
plt.xlabel("R")
plt.ylabel("Z")
plt.savefig("python-plots/wham-nodes-dn.png", dpi=300)
plt.close()
