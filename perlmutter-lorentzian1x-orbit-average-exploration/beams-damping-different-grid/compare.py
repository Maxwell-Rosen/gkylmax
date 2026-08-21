import postgkyl as pg

dmap = pg.load("gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl").interpolate().map("gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl")
norm = pg.load("../beams-damping/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl").interpolate().map("../beams-damping/gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl")
pg.plot(dmap, norm)