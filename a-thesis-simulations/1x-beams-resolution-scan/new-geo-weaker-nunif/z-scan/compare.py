import string

import postgkyl as pg

bimax_filename = "zzim-ion_BiMaxwellianMoments_65.gkyl"
mc2nu_filename = "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
bimax_filename_gkl = "zzim-ion_BiMaxwellianMoments_65.gkyl"
mc2nu_filename_gkl = "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
orig_folder = "/global/homes/m/mhrosen/scratch/gkylmax/a-thesis-simulations/1x-beams/"

elow = pg.load("128/" + bimax_filename_gkl).interpolate().map("128/" + mc2nu_filename_gkl)
low = pg.load("192/" + bimax_filename_gkl).interpolate().map("192/" + mc2nu_filename_gkl)
# med = pg.load("256/" + bimax_filename_gkl).interpolate().map("256/" + mc2nu_filename_gkl)
# high = pg.load("320/" + bimax_filename_gkl).interpolate().map("320/" + mc2nu_filename_gkl)
original = pg.load(orig_folder + bimax_filename).interpolate().map(orig_folder + mc2nu_filename)

factor = 1.602176634e-27 * 2.014 / 1.602176634e-19 / 1000

for data in (elow, low, original):
    data[..., 2:4] *= factor

pg.plot(
    elow, low, original,
    color=["#E61F00", "#E69F00", "#009E73"],
    figure = 0,
    legend_labels=["128", "192", "Original"],
    legend_subplot=0,
    legend_loc="best",
    title=r"Resolution Scan along $z$",
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
    saveas="res-scan-z.pdf",
)
