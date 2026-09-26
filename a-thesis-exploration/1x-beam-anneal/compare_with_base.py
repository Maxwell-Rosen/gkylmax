import postgkyl as pg

bflux_anneal = pg.load("zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp=0)
isrc_anneal = pg.load("zzim-ion_source_integrated_moms.gkyl").select(comp=0)
ratio_anneal = bflux_anneal / isrc_anneal * 2

bflux_base = pg.load("../1x-beam-collisionless-dilation/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp=0)
isrc_base = pg.load("../1x-beam-collisionless-dilation/zzim-ion_source_integrated_moms.gkyl").select(comp=0)
ratio_base = bflux_base / isrc_base * 2

bflux_ecycles = pg.load("../1x-beam-extra-cycles/zzim-ion_bflux_xlower_integrated_HamiltonianMoments.gkyl").select(comp=0)
isrc_ecycles = pg.load("../1x-beam-extra-cycles/zzim-ion_source_integrated_moms.gkyl").select(comp=0)
ratio_ecycles = bflux_ecycles / isrc_ecycles * 2

# pg.plot(ratio_anneal, ratio_base, grid_indices = True, ymin = 0, ymax = 4,
#         legend_labels=["Anneal", "Base"], figure = 0)

end_tail_anneal = ratio_anneal.select(z0="-1200:")
end_tail_anneal.grid[0] = end_tail_anneal.grid[0] - end_tail_anneal.grid[0][0]
fit_anneal = end_tail_anneal.fit(fit_type = "exp_plateau", print_coeffs = True)

end_tail_base = ratio_base.select(z0="-1200:")
end_tail_base.grid[0] = end_tail_base.grid[0] - end_tail_base.grid[0][0]
fit_base = end_tail_base.fit(fit_type = "exp_plateau", print_coeffs = True)

end_tail_ecycles = ratio_ecycles.select(z0="-1200:")
end_tail_ecycles.grid[0] = end_tail_ecycles.grid[0] - end_tail_ecycles.grid[0][0]
fit_ecycles = end_tail_ecycles.fit(fit_type = "linear", print_coeffs = True)

pg.plot(end_tail_anneal, fit_anneal, end_tail_base, fit_base, end_tail_ecycles, fit_ecycles,  figure = 0, legend_labels=["Anneal", "Fit, C=1.14", "Base", "Fit, C=1.18", "Extra Cycles", "Fit, C=1.14"], title = "Beam source. Boundary flux during final FDP. Fits A exp(-bx) + C", ylabel = "Bflux / Source", xlabel = "Time, s", linestyle = ["-", "--", "-", "--", "-", "--"])