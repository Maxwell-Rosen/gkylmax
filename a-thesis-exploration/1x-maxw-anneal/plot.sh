pgkyl zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl zzim-ion_source_integrated_moms.gkyl sel -c0 --z0 2: ev "f0 f1 / 2 *" pl --saveas "bfluxratio.png" --no_show --ymax 2 --ymin 0 -g &

pgkyl zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl zzim-ion_source_integrated_moms.gkyl sel -c0 --z0 2: ev "f0 f1 / 2 *" pl --saveas "bfluxratio_log.png" --no_show --logy -g &

pgkyl zzim-ion_integrated_moms.gkyl sel -c0 --z0 2: pl --saveas "iinteg.png" --no_show &

pgkyl zzim-ion_integrated_moms.gkyl sel -c0 --z0 2: fit -f exp_plateau -p pl -f0 --saveas "iinteg_fit.png" --no_show &

pgkyl zzim-ion_M1_80.gkyl \
  zzim-geo_int_jacobgeo.gkyl \
  zzim-geo_int_rtg33inv.gkyl \
  interp -n 1 \
  ev 'f0 f1 * f2 * 4.1637715094428562e+20 /' \
  map zzim-geo_corn_mc2nu_pos_deflated.gkyl \
  pl --xlabel 'Physical z (m)' --ylabel 'Particle rate / source' --saveas "total_fluxes.png" --no_show -g &