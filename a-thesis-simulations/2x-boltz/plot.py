import postgkyl as pg

single_field = pg.load("../1x-beams/zzim-ion_BiMaxwellianMoments_20.gkyl").interpolate()
many_field = pg.load("zzim-ion_BiMaxwellianMoments_20.gkyl").interpolate().select(z0=-1)

pg.plot(single_field, many_field, title="1x vs 2x", legend=["1x", "2x"])