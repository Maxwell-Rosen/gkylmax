import postgkyl as pg

frame = 65
single_field = pg.load("../1x-beams/zzim-ion_BiMaxwellianMoments_" + str(frame) + ".gkyl").interpolate()
many_field_inner = pg.load("zzim-ion_BiMaxwellianMoments_" + str(frame) + ".gkyl").interpolate().select(z0=0)
many_field_outer = pg.load("zzim-ion_BiMaxwellianMoments_" + str(frame) + ".gkyl").interpolate().select(z0=-1)

factor = 1.602176634e-27 * 2.014 / 1.602176634e-19 / 1000

for data in (single_field, many_field_inner, many_field_outer):
    data[..., 2:4] *= factor

pg.plot(single_field, many_field_inner, many_field_outer,
    color=["#0072B2", "#D55E00", "#059E00"],
    title="1x vs 2x",
    legend = True,
    legend_labels = ["1x", "2x inner", "2x outer"],
    legend_subplot = 0,
    legend_loc = "best",
    xlabel = "z [m]",
    figsize = "12,10",
    split_linear_log = True,
    split_point = 0.0,
    split_log_side = "right",
    split_width_ratios = (1.0, 1.0),
    split_gap = 0.0,
    split_legend_side = "log",
    split_log_nonpositive="mask",
    split_linear_ylim=[
        (0.0, 3e19),
        (-1.4e6, 1.4e6),
        (-2, 14.0),
        (0.0, 27.0),
    ],
    split_log_ylim=[
        (1e13, 4e19),
        (1e2, 2e6),
        (1e-1, 20.0),
        (1e-1, 50.0),
    ],
    subplot_ylabels=(
        r"Density [$\mathrm{m}^{-3}$],"
        r"$U_\parallel$ [m/s],"
        r"$T_\parallel$ [keV],"
        r"$T_\perp$ [keV]"
    ),
    saveas="two-vs-one-x-comparison-bimax.pdf",
)