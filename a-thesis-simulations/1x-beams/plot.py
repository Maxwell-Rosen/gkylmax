import postgkyl as pg

data = pg.load("zzim-ion_BiMaxwellianMoments_*.gkyl").interpolate().collect().plot()

# # new_map = pg.load("gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl").interpolate()
# # old_map = pg.load("../beams-damping/gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl").interpolate()

# # pg.plot(
# #     new_map, old_map,
# #     legend_labels=["new map", "old map"],
# #     legend=True,
# #     legend_subplot=0,
# #     legend_loc="best",
# #     title="Comparing the old and new non-uniform map",
# #     ylabel="z [m]",
# #     xlabel="z_0 computational",
# #     figsize="12,10",
# # )


# new_map_bimax = pg.load("gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl").interpolate().map("gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl")
# old_map_bimax = pg.load("../beams-damping/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl").interpolate().map("../beams-damping/gk_lorentzian_mirror-geo_corn_mc2nu_pos_deflated.gkyl")

# factor = 1.67e-27 * 2.014 / 1e-19 / 1000

# for data in (new_map_bimax, old_map_bimax):
#     data[2:4][...] *= factor

# pg.plot(
#     new_map_bimax, old_map_bimax,
#     legend_labels=["new map", "old map"],
#     legend=True,
#     legend_subplot=0,
#     legend_loc="best",
#     title="Comparing the old grid and new grid at Nz=400",
#     xlabel="z [m]",
#     figsize="12,10",
#     # split_linear_log=True,
#     # split_point=0.0,
#     # split_log_side="right",
#     # split_width_ratios=(1.0, 1.0),
#     # split_gap=0.0,
#     # split_legend_side="log",
#     # split_log_nonpositive="mask",
#     # split_linear_ylim=[
#     #     (0.0, 3e19),
#     #     (-1.4e6, 1.4e6),
#     #     (-2, 22.0),
#     #     (0.0, 45.0),
#     # ],
#     # split_log_ylim=[
#     #     (1e13, 4e19),
#     #     (1e2, 2e6),
#     #     (1e-1, 30.0),
#     #     (1e-1, 100.0),
#     # ],
#     subplot_ylabels=(
#         r"Density [$\mathrm{m}^{-3}$],"
#         r"$U_\parallel$ [m/s],"
#         r"$T_\parallel$ [keV],"
#         r"$T_\perp$ [keV]"
#     ),
# )

