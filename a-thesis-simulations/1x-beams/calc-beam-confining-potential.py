import postgkyl as pg
import numpy as np

max_field = np.max(pg.load('zzim-field_65.gkyl').interpolate())
zero_field = pg.load('zzim-field_65.gkyl').interpolate().select(z0 = 0.0)

diff = max_field - zero_field

print(diff)