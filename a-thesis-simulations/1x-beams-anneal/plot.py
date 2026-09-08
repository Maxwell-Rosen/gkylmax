import postgkyl as pg

# curr_int = pg.load("zzim-ion_fdot_abs_integrated_moms.gkyl")
# curr_src = pg.load("zzim-ion_source_integrated_moms.gkyl")

# div_curr = curr_int / curr_src
# div_curr.select(comp = 0, inplace = True)

# old_int = pg.load("diagnostics-before-final-anneal-20260903T084551-1722040/zzim-ion_fdot_abs_integrated_moms.gkyl")
# old_src = pg.load("diagnostics-before-final-anneal-20260903T084551-1722040/zzim-ion_source_integrated_moms.gkyl")

# div_old = old_int / old_src
# div_old.select(comp = 0, inplace = True)


# pg.plot(div_old,div_curr,
#         figure = 0,
#         logy = True,
#         xmin = 0.623333)

# curr_int_dens = pg.load("zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp = 0)
# old_int_dens = pg.load("diagnostics-before-final-anneal-20260903T084551-1722040/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp = 0)

# pg.plot(old_int_dens, curr_int_dens,
#         figure = 0,
#         logy = True,
#         grid_indices = True,
# )

# curr_int_dens = pg.load("zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp = 0)
# old_int_dens = pg.load("diagnostics-before-final-anneal-20260903T084551-1722040/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp = 0)
# master_int_dens = pg.load("../1x-beams/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp = 0)

# pg.plot(old_int_dens, curr_int_dens,master_int_dens,
#         legend = True,
#         legend_labels = ["Before Anneal", "After Anneal", "Beams"],
#         figure = 0,
#         logy = True,
#         grid_indices = True,
# )


curr_int_dens = pg.load("zzim-ion_integrated_moms.gkyl").select(comp = 0)
old_int_dens = pg.load("diagnostics-before-final-anneal-20260903T084551-1722040/zzim-ion_integrated_moms.gkyl").select(comp = 0)
master_int_dens = pg.load("../1x-beams/zzim-ion_integrated_moms.gkyl").select(comp = 0)

pg.plot(old_int_dens, curr_int_dens,master_int_dens,
        legend = True,
        legend_labels = ["Before Anneal", "After Anneal", "Beams"],
        figure = 0,
        # logy = True,
        grid_indices = True,
)