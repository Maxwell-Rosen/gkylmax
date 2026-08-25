import string

import postgkyl as pg

bimax_filename = "zzim-ion_BiMaxwellianMoments_65.gkyl"
mc2nu_filename = "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
orig_folder = "/global/homes/m/mhrosen/scratch/gkylmax/a-thesis-simulations/1x-beams/"

with_pos = pg.load(bimax_filename).interpolate().map(mc2nu_filename)
original = pg.load(orig_folder + bimax_filename).interpolate().map(orig_folder + mc2nu_filename)

factor = 1.602176634e-27 * 2.014 / 1.602176634e-19 / 1000

for data in (with_pos, original):
    data[..., 2:4] *= factor

pg.plot(
    with_pos, original,
    # linestyle=["-", "--"],
    color=["#D55E00", "#0072B2"],
    legend_labels=[r"without $\nu_{ie}$", r"with $\nu_{ie}$"],
    legend=True,
    legend_subplot=1,
    legend_loc="best",
    title=r"Comparing POA simulations with and without ion-electron collisions",
    xlabel="z [m]",
    figsize="12,10",
    split_linear_log=True,
    split_point=0.0,
    split_log_side="right",
    split_width_ratios=(1.0, 1.0),
    split_gap=0.0,
    split_legend_side="linear",
    split_log_nonpositive="mask",
    split_linear_ylim=[
        (0.0, 5e19),
        (-1.6e6, 1.6e6),
        (-2, 30.0),
        (0.0, 60.0),
    ],
    split_log_ylim=[
        (1e13, 7e19),
        (1e2, 2e6),
        (1e-1, 35.0),
        (1e-1, 60.0),
    ],
    subplot_ylabels=(
        r"Density [$\mathrm{m}^{-3}$],"
        r"$U_\parallel$ [m/s],"
        r"$T_\parallel$ [keV],"
        r"$T_\perp$ [keV]"
    ),
    saveas="ion-electron-collisions-comparison.pdf",
)
