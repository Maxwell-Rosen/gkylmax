import postgkyl as pg

tandem_bmag = pg.load("zzim-geo_corn_bmag.gkyl").interpolate().map("zzim-geo_corn_mc2nu_pos_deflated.gkyl")

simple_bmag = pg.load("../1x-beams-R-scan/R-scan/R-15/zzim-geo_corn_bmag.gkyl").interpolate().map("../1x-beams-R-scan/R-scan/R-15/zzim-geo_corn_mc2nu_pos_deflated.gkyl")

simple_left = simple_bmag.clone()
simple_right = simple_bmag.clone()
simple_left.grid[0] -= 3.0
simple_right.grid[0] += 3.0

pg.plot(
    simple_left, simple_right, tandem_bmag,
    legend=True,
    legend_labels=["Simple R=15", "Simple R=15", "Tandem Mirror"],
    legend_loc="best",
    title="Magnetic field magnitude for the tandem mirror R=15 compared to an R=15 simple mirror",
    xlabel="z [m]",
    figsize="12,5",
    split_linear_log=True,
    split_point=0.0,
    split_log_side="right",
    split_width_ratios=(1.0, 1.0),
    split_gap=0.0,
    split_legend_side="linear",
    split_log_nonpositive="mask",
    # split_linear_ylim=[
    #     (0.0, 1.0),
    #     (0.0, 1.0),
    # ],
    # split_log_ylim=[
    #     (1e-2, 1.0),
    #     (1e-2, 1.0),
    # ],
    ylabel=(
        r"$|B|$, T"
    ),
    saveas="tandem-mirror-15-bmag.pdf",
)