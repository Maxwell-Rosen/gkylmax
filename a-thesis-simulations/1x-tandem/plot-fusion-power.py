"""Plot the local D-T fusion power density for the tandem-mirror case.

The simulated ion density is interpreted as the total density of a 50:50
deuterium-tritium mixture.  The anisotropic ion temperature is reduced to
T_i = (T_parallel + 2 T_perp)/3.  CFSPOPCON supplies the Maxwellian-averaged
D-T reactivity; no fit coefficients are duplicated in this script.

Install the additional dependency into the Python environment that contains
Postgkyl with ``python -m pip install cfspopcon``.
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


FRAME = 130
CASE_DIR = Path(__file__).resolve().parent
MOMENTS_FILE = CASE_DIR / f"zzim-ion_BiMaxwellianMoments_{FRAME}.gkyl"
MAP_FILE = CASE_DIR / "zzim-geo_corn_mc2nu_pos_deflated.gkyl"
OUTPUT_FILE = CASE_DIR / f"dt-fusion-power-{FRAME}.pdf"

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


def main():
    moments = pg.load(str(MOMENTS_FILE)).interpolate().map(str(MAP_FILE))

    # Every array below is evaluated pointwise on the mapped Z grid.
    ion_density_z = moments[..., 0]
    temperature_factor = ION_MASS / (ELEMENTARY_CHARGE * 1.0e3)
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

    fusion_power_values_z = (
        deuterium_density_z
        * tritium_density_z
        * reactivity_z
        * DT_ENERGY
    )

    # Reuse the mapped Postgkyl grid and replace the four moment components
    # with the single derived fusion-power-density component.
    fusion_power = moments.clone(data=False)
    fusion_power.push(list(moments.grid), fusion_power_values_z[..., np.newaxis])
    fusion_power.label = "D-T fusion power density"

    pg.plot(
        fusion_power,
        xlabel=r"$Z$ [m]",
        ylabel=r"D-T fusion power density [W m$^{-3}$]",
        title="Tandem mirror D-T fusion power",
        figsize="10,6",
        saveas=str(OUTPUT_FILE),
    )


if __name__ == "__main__":
    main()
