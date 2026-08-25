import string

import postgkyl as pg

bimax_filename = "zzim-ion_BiMaxwellianMoments_65.gkyl"
mc2nu_filename = "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
orig_folder = "/global/homes/m/mhrosen/scratch/gkylmax/a-thesis-simulations/1x-beams/"

high_alpha = pg.load(bimax_filename).interpolate().map(mc2nu_filename)
original = pg.load(orig_folder + bimax_filename).interpolate().map(orig_folder + mc2nu_filename)

factor = 1.602176634e-27 * 2.014 / 1.602176634e-19 / 1000

for data in (high_alpha, original):
    data[..., 2:4] *= factor

pg.plot(
    high_alpha, original,
    color=["#D55E00", "#0072B2"],
    legend_labels=["2e-4", "2e-5"],
    legend=True,
    legend_subplot=0,
    legend_loc="best",
    title=r"Comparing $\alpha$ = 2e-4 to $\alpha$ = 2e-5",
    xlabel="z [m]",
    figsize="12,10",
    split_linear_log=True,
    split_point=0.0,
    split_log_side="right",
    split_width_ratios=(1.0, 1.0),
    split_gap=0.0,
    split_legend_side="log",
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
    saveas="alpha-converge.pdf",
)
