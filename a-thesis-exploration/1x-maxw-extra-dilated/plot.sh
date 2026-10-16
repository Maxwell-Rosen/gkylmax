pgkyl zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl zzim-ion_source_integrated_moms.gkyl sel -c0 --z0 2: ev "f0 f1 / 2 *" pl --saveas "bfluxratio.png" --no_show &

pgkyl zzim-ion_integrated_moms.gkyl sel -c0 --z0 2: pl --saveas "iinteg.png" --no_show &

pgkyl zzim-ion_integrated_moms.gkyl sel -c0 --z0 2: fit -f exp_plateau -p pl -f0 --saveas "iinteg_fit.png" --no_show &