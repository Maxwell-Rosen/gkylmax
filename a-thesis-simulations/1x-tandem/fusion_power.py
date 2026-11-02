"""Plot D-T fusion and source power densities and calculate field-line gain.

The simulated ion density is interpreted as the total density of a 50:50
deuterium-tritium mixture.  The anisotropic ion temperature is reduced to
T_i = (T_parallel + 2 T_perp)/3.  CFSPOPCON supplies the Maxwellian-averaged
D-T reactivity; no fit coefficients are duplicated in this script.

The source M2 moment supplies the injected ion kinetic power density. For the
psi radial coordinate selected in sim.c, 2*pi*integral(p*jacobgeo dxi) is
dP/dpsi [W/Wb], not total device power. The 2*pi factor integrates azimuth.
Q_line = (dP_out/dpsi)/(dP_in/dpsi); Q_line > 1 means fusion output exceeds
the injected ion kinetic power on this field line.
Total device powers require radial integration of both dP/dpsi profiles;
Q_line equals device Q only with an appropriate radial-profile assumption.

Install the additional dependency into the Python environment that contains
Postgkyl with ``python -m pip install cfspopcon``.

Cite both bosch-Hale and the cfspopcon repo
"""

from pathlib import Path

import numpy as np
import postgkyl as pg

try:
    from cfspopcon.unit_handling import ureg
except ImportError:
    try:
        # Compatibility with releases that exported the registry at top level.
        from cfspopcon import ureg
    except ImportError as error:
        raise ImportError(
            "CFSPOPCON and its unit registry are required; install or update it "
            "with `python -m pip install --upgrade cfspopcon`."
        ) from error

try:
    from cfspopcon.formulas.fusion_power import DTFusionBoschHale
except ImportError:
    # Compatibility with releases that do not re-export the reaction class.
    from cfspopcon.formulas.fusion_power.fusion_data import DTFusionBoschHale


FRAME = 75
CASE_DIR = Path(__file__).resolve().parent
MOMENTS_FILE = CASE_DIR / f"zzim-ion_BiMaxwellianMoments_{FRAME}.gkyl"
MAP_FILE = CASE_DIR / "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
JACOBIAN_FILE = CASE_DIR / "zzim-geo_int_jacobgeo.gkyl"
OUTPUT_FILE = CASE_DIR / f"dt-fusion-power-{FRAME}.pdf"
SOURCE_OUTPUT_FILE = CASE_DIR / f"ion-source-power-{FRAME}.pdf"

ELEMENTARY_CHARGE = 1.602176634e-19  # J/eV
PROTON_MASS = 1.67262192369e-27  # kg
ION_MASS = 2.014 * PROTON_MASS  # deuterium mass used in sim.c
DT_ENERGY = 17.6e6 * ELEMENTARY_CHARGE  # energy released per reaction [J]
DT_REACTION = DTFusionBoschHale()


def dt_reactivity_cfspopcon(temperature_kev):
    """Return CFSPOPCON's Maxwellian-averaged D-T reactivity in m^3/s."""
    temperature_kev = np.asarray(temperature_kev, dtype=float)
    reactivity = np.zeros_like(temperature_kev)
    positive = temperature_kev > 0.0
    if np.any(positive):
        rate = DT_REACTION.calc_rate_coefficient(
            temperature_kev[positive] * ureg.keV
        )
        if hasattr(rate, "pint"):
            # CFSPOPCON normally returns a Pint-aware xarray.DataArray.
            rate_values = rate.pint.to("m^3 / s").pint.magnitude
        else:
            # Retain compatibility with versions returning a Pint Quantity.
            rate_values = rate.to("m^3 / s").magnitude
        reactivity[positive] = np.asarray(rate_values)
    return reactivity


def fusion_power_density(moments):
    """Return D-T power density [W/m^3] from evaluated ion moments."""
    ion_density_z = moments[..., 0]
    temperature_factor = ION_MASS / (ELEMENTARY_CHARGE * 1.0e3)  # keV
    temperature_parallel_z = moments[..., 2] * temperature_factor
    temperature_perpendicular_z = moments[..., 3] * temperature_factor
    temperature_ion_z = (
        temperature_parallel_z + 2.0 * temperature_perpendicular_z
    ) / 3.0

    # Interpret the local simulated density as an equimolar D-T mixture.
    # Nonphysical negative interpolated values contribute no fusion power.
    physical_ion_density_z = np.where(ion_density_z > 0.0, ion_density_z, 0.0)
    deuterium_density_z = 0.5 * physical_ion_density_z
    tritium_density_z = 0.5 * physical_ion_density_z
    reactivity_z = dt_reactivity_cfspopcon(temperature_ion_z)

    return (
        deuterium_density_z
        * tritium_density_z
        * reactivity_z
        * DT_ENERGY
    )


def field_line_powers(moments, source, jacobian, num_quad=4):
    """Return (dP_in/dpsi, dP_out/dpsi) in W/Wb, integrated over azimuth.

    Gauss integration stays on the computational grid: jacobgeo already
    includes the position-map Jacobian, so mapped Z widths must not be used.
    sim.c selects PSI_CART_Z; a sqrt(psi) radial coordinate would require
    an additional division by 2*sqrt(psi) to obtain powers per unit psi.
    The source is time independent in sim.c. Its M2 moment excludes the
    electrostatic potential energy present in the Hamiltonian diagnostic.
    """
    for data in (moments, source, jacobian):
        if data.num_dims != 1:
            raise ValueError("Field-line integration requires 1D configuration data.")
        if not np.array_equal(data.grid[0], moments.grid[0]):
            raise ValueError("Moments, source, and Jacobian must share a grid.")

    num_cells = len(moments.grid[0]) - 1
    # Postgkyl stores each component's quadrature points consecutively.
    moments_q = moments.represent(to="quad", num_quad=num_quad).values.reshape(
        num_cells, 4, num_quad
    ).transpose(0, 2, 1)
    source_q = source.represent(to="quad", num_quad=num_quad).values.reshape(
        num_cells, num_quad
    )
    jacobian_q = jacobian.represent(to="quad", num_quad=num_quad).values.reshape(
        num_cells, num_quad
    )
    _, weights = np.polynomial.legendre.leggauss(num_quad)
    # dV/dpsi = 2*pi*jacobgeo*dxi for this axisymmetric psi-coordinate case.
    volume_per_flux_weights = (
        2.0 * np.pi
        * 0.5 * np.diff(moments.grid[0])[:, np.newaxis] * weights * jacobian_q
    )
    source_power_q = 0.5 * source.ctx["mass"] * source_q
    dp_in_dpsi = float(np.sum(source_power_q * volume_per_flux_weights))
    dp_out_dpsi = float(np.sum(fusion_power_density(moments_q) * volume_per_flux_weights))
    if not (np.isfinite(dp_in_dpsi) and np.isfinite(dp_out_dpsi)
            and dp_in_dpsi > 0.0 and dp_out_dpsi > 0.0):
        raise ValueError("Gain requires finite, positive input and fusion powers.")
    return dp_in_dpsi, dp_out_dpsi


def main():
    raw_moments = pg.load(str(MOMENTS_FILE))
    # A static source is written at initialization and again on restart.
    source_files = [
        path for path in CASE_DIR.glob("zzim-ion_source_M2_*.gkyl")
        if int(path.stem.rsplit("_", 1)[1]) <= FRAME
    ]
    if not source_files:
        raise FileNotFoundError(f"No source M2 diagnostic at or before frame {FRAME}.")
    source_file = max(source_files, key=lambda path: int(path.stem.rsplit("_", 1)[1]))
    source = pg.load(str(source_file))
    jacobian = pg.load(str(JACOBIAN_FILE))
    dp_in_dpsi, dp_out_dpsi = field_line_powers(raw_moments, source, jacobian)
    gain = dp_out_dpsi / dp_in_dpsi
    print(f"Frame {FRAME}, t = {raw_moments.ctx['time']:.8g} s")
    print(f"Source diagnostic: {source_file.name}")
    print("Axisymmetric powers per unit psi, including the 2*pi azimuthal integral:")
    print(f"dP_in/dpsi  (ion source kinetic power) = {dp_in_dpsi:.8g} W/Wb")
    print(f"dP_out/dpsi (D-T fusion power)         = {dp_out_dpsi:.8g} W/Wb")
    print(f"Q_line = (dP_out/dpsi) / (dP_in/dpsi) = {gain:.8g}")
    print("Total device powers and Q require radial profiles or an explicit radial assumption.")

    moments = raw_moments.interpolate().map(str(MAP_FILE))
    fusion_power_values_z = fusion_power_density(moments)

    # Reuse the mapped Postgkyl grid and replace the four moment components
    # with the single derived fusion-power-density component.
    fusion_power = moments.clone(metadata_only=True)
    fusion_power.push(list(moments.grid), fusion_power_values_z[..., np.newaxis])
    fusion_power.label = "D-T fusion power density"

    source_power = source.interpolate().map(str(MAP_FILE)) * (0.5 * source.ctx["mass"])
    source_power.label = "Ion source power density"

    pg.plot(
        fusion_power / 1000,
        legend_labels=["D-T fusion"],
        xlabel=r"$Z$ [m]",
        ylabel=r"D-T fusion power density [kW m$^{-3}$]",
        title="",
        figsize="10,6",
        saveas=str(OUTPUT_FILE),
    )

    pg.plot(
        source_power / 1000,
        legend_labels=["Ion beam source"],
        color="tab:orange",
        xlabel=r"$Z$ [m]",
        ylabel=r"Ion source power density [kW m$^{-3}$]",
        title="",
        figsize="10,6",
        saveas=str(SOURCE_OUTPUT_FILE),
    )


if __name__ == "__main__":
    main()
