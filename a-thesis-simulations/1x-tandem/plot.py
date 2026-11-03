import matplotlib.pyplot as plt
import postgkyl as pg

data = pg.load("zzim-ion_BiMaxwellianMoments_75.gkyl").interpolate().map("zzim-geo_corn_mc2nu_pos_deflated.gkyl")

factor = 1.67e-27 * 2.014 / 1e-19 / 1000

data[...,2:4] *= factor

single_comparison = pg.load("../1x-beams-R-scan/R-scan/R-15/zzim-ion_BiMaxwellianMoments_75.gkyl").interpolate().map("../1x-beams-R-scan/R-scan/R-15/zzim-geo_corn_mc2nu_pos_deflated.gkyl")
single_comparison[...,2:4] *= factor

simple_left = single_comparison.clone()
simple_right = single_comparison.clone()
simple_left.grid[0] -= 3.0
simple_right.grid[0] += 3.0

fig = pg.plot(
    simple_left, simple_right, data,
    legend_labels=["Simple R=15", "Simple R=15", "Tandem Mirror"],
    legend_subplot=1,
    figure = 0,
    legend_loc=(0.03, 0.68),
    title="Tandem mirror R=15 with an R=15 simple mirror",
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
        (0.0, 2.2e19),
        (-1.5e6, 1.5e6),
        (-2.0, 26.0),
        (0.0, 55.0),
    ],
    split_log_ylim=[
        (1e13, 4e19),
        (1e-1, 2e6),
        (8e-2, 30.0),
        (2e-1, 70.0),
    ],
    subplot_ylabels=(
        r"Density [$\mathrm{m}^{-3}$],"
        r"$U_\parallel$ [m/s],"
        r"$T_\parallel$ [keV],"
        r"$T_\perp$ [keV]"
    ),
    no_show=True,
)

# Each panel has a linear axis followed by a logarithmic axis.
for ax, label in zip(fig.axes[::2], ("(a)", "(b)", "(c)", "(d)")):
    ax.text(
        0.03, 0.97, label,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=14,
    )

for extension in ("pdf", "png"):
    fig.savefig(f"tandem-mirror-15.{extension}", dpi=200)

plt.show()
