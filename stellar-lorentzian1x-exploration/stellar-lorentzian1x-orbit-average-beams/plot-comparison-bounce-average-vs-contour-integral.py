"""
Beam source comparison: spatial vs orbit-averaged vs contour-integrated sources.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from multiprocessing import Pool, cpu_count

def _contour_worker(args):
  """Worker function for parallel contour integration."""
  z, vpar, mu, obj = args
  return obj.contour_integral_beam_source(z, vpar, mu)

class BeamSourceComparison:
  """
  Class to compare different beam source models in a mirror machine.
  
  Models:
  1. Spatial beam source (Dorf et al. 2025): Gaussian in z, beam distribution in v
  2. Orbit-averaged source: Maps velocities to midplane using energy conservation
  3. Contour-integrated source: Averages spatial source over energy contours
  """
  
  # Physical constants
  eV = 1.602176634e-19
  mp = 1.67262192369e-27
  mi = 2.014 * mp  # Deuteron mass
  
  def __init__(self, 
         # Geometry parameters
         mcB=6.51292,
         gamma_geo=0.124904,
         Z_m=0.98,
         B_p=0.53,
         # Beam parameters
         gamma_0=4174.226,
         gamma_oavg=1109.5135852642,
         T_beam_eV=200,
         E_beam_eV=25000,
         Lb=0.2):
    """
    Initialize beam source comparison.
    
    Parameters:
      mcB : float - Magnetic field parameter
      gamma_geo : float - Lorentzian width parameter
      Z_m : float - Mirror throat location (m)
      B_p : float - Reference magnetic field (T)
      gamma_0 : float - Spatial source normalization
      gamma_oavg : float - Orbit-averaged source normalization
      T_beam_eV : float - Beam temperature (eV)
      E_beam_eV : float - Beam energy (eV)
      Lb : float - Spatial source width (m)
    """
    # Geometry
    self.mcB = mcB
    self.gamma_geo = gamma_geo
    self.Z_m = Z_m
    self.B_p = B_p
    
    # Beam parameters
    self.gamma_0 = gamma_0
    self.gamma_oavg = gamma_oavg
    self.T_beam = T_beam_eV * self.eV
    self.E_beam = E_beam_eV * self.eV
    self.v_beam = np.sqrt(self.E_beam / self.mi)
    self.sigma_beam = 2 * self.T_beam / self.mi
    self.Lb = Lb
    
    # Precompute bmag spline for fast evaluation in contour integrals
    self._z_grid = np.linspace(-1.0, 1.0, 2001)
    self._bmag_grid = self._bmag_raw(self._z_grid)
    self._bmag_spline = CubicSpline(self._z_grid, self._bmag_grid)
  
  def _bmag_raw(self, Z):
    """Raw magnetic field calculation (used for spline construction)."""
    Z = np.atleast_1d(Z)
    L1 = 1.0 / (1 + ((Z - self.Z_m) / self.gamma_geo)**2)
    L2 = 1.0 / (1 + ((Z + self.Z_m) / self.gamma_geo)**2)
    Bmag = self.mcB / (np.pi * self.gamma_geo) * (L1 + L2)
    return Bmag.item() if Bmag.size == 1 else Bmag
  
  def bmag(self, Z, use_spline=False):
    """Compute magnetic field magnitude at position Z.
    
    Args:
      Z: Position(s) to evaluate
      use_spline: If True, use precomputed spline (faster for many evaluations)
    """
    if use_spline:
      return self._bmag_spline(Z)
    return self._bmag_raw(Z)
  
  def dbmag_dz(self, Z):
    """Compute derivative of magnetic field dB/dZ."""
    Z = np.atleast_1d(Z)
    # d/dZ of 1/(1 + ((Z - Z_m)/gamma)^2) = -2(Z - Z_m) / (gamma^2 * (1 + ((Z-Z_m)/gamma)^2)^2)
    dL1 = -2 * (Z - self.Z_m) / (self.gamma_geo**2 * (1 + ((Z - self.Z_m) / self.gamma_geo)**2)**2)
    dL2 = -2 * (Z + self.Z_m) / (self.gamma_geo**2 * (1 + ((Z + self.Z_m) / self.gamma_geo)**2)**2)
    dBdZ = self.mcB / (np.pi * self.gamma_geo) * (dL1 + dL2)
    return dBdZ.item() if dBdZ.size == 1 else dBdZ
  
  def midplane_beam_source(self, vpar, mu):
    """Beam source evaluated at midplane velocities."""
    vperp = np.sqrt(2.0 * mu * self.B_p / self.mi)
    return self.gamma_oavg * np.exp(-((np.abs(vpar) - self.v_beam)**2 + 
                      (vperp - self.v_beam)**2) / self.sigma_beam)
  
  def spatial_beam_source(self, z, vpar, mu):
    """Spatial beam source with Gaussian z-dependence (Dorf et al. 2025)."""
    zdep = np.exp(-(z / self.Lb)**2)
    vdep = self.midplane_beam_source(vpar, mu) / self.gamma_oavg
    return self.gamma_0 * zdep * vdep
  
  def orbit_avg_beam_source(self, z, vpar, mu):
    """Orbit-averaged source: map vpar to midplane using energy conservation."""
    # Set source to zero outside mirror throat region
    z = np.atleast_1d(z)
    result = self.midplane_beam_source(np.sqrt(vpar**2 + 2 * mu * (self.bmag(z) - self.bmag(0)) / self.mi), mu)
    result = np.where(np.abs(z) > self.Z_m, 0.0, result)
    return result.item() if result.size == 1 else result
  
  def find_bounce_points(self, E, mu, z_min=-0.98, z_max=0.98):
    if mu <= 0:
      return z_min, z_max
    
    B_bounce = E / mu
    
    # Function whose roots give bounce points
    def f(zp):
      return self.bmag(zp) - B_bounce
    
    # Check if particle can access any region
    B_min = self.bmag(0)  # B is minimum at midplane
    if B_bounce < B_min:
      return None  # Particle cannot exist here
    
    # Find left bounce point (between z_min and 0)
    f_left = f(z_min)
    f_mid = f(0)
    if f_left * f_mid < 0:
      z_left = brentq(f, z_min, 0, xtol=1e-12)
    elif f_left <= 0:
      z_left = z_min  # Accessible all the way to z_min
    else:
      return None  # No accessible region on left
    
    # Find right bounce point (between 0 and z_max)
    f_right = f(z_max)
    if f_mid * f_right < 0:
      z_right = brentq(f, 0, z_max, xtol=1e-12)
    elif f_right <= 0:
      z_right = z_max  # Accessible all the way to z_max
    else:
      return None  # No accessible region on right
    
    return z_left, z_right
  
  def contour_integral_beam_source(self, z, vpar, mu, z_min=-0.98, z_max=0.98, use_fast=True):
    """
    Compute orbit-averaged source by integrating along energy contours.
    Uses time weighting: <S> = (∫ S dt) / (∫ dt) = (∫ S dz/|v_par|) / (∫ dz/|v_par|)
    
    Args:
      use_fast: If True, use fixed-order Gaussian quadrature with spline bmag (much faster)
    """
    E = 0.5 * self.mi * vpar**2 + mu * self.bmag(z)
    
    # Find bounce points
    bounce_pts = self.find_bounce_points(E, mu, z_min, z_max)
    if bounce_pts is None:
      return 0.0
    z_left, z_right = bounce_pts
    
    if use_fast:
      result = self._contour_integral_fast(E, mu, z_left, z_right)
    else:
      result = self._contour_integral_adaptive(E, mu, z_left, z_right)

    result = np.where(np.abs(z) > self.Z_m, 0.0, result)
    return result
  
  def _contour_integral_fast(self, E, mu, z_left, z_right, n_quad=32):
    """
    Fast contour integral using fixed-order Gaussian quadrature.
    
    Uses change of variables to handle sqrt singularity at bounce points:
    z = z_mid + (z_half) * sin(theta), theta in [-pi/2, pi/2]
    dz = z_half * cos(theta) d(theta)
    
    This removes the 1/sqrt singularity since vpar ~ sqrt(cos(theta)) near boundaries.
    """
    z_mid = 0.5 * (z_left + z_right)
    z_half = 0.5 * (z_right - z_left)
    
    two_over_mi = 2.0 / self.mi
    
    def integrand_transformed(theta):
      """Integrand in transformed coordinates (vectorized)."""
      theta = np.atleast_1d(theta)
      cos_th = np.cos(theta)
      z = z_mid + z_half * np.sin(theta)
      
      # Use spline for fast bmag evaluation
      B = self._bmag_spline(z)
      kinetic = E - mu * B
      
      # vpar * dz/dtheta = sqrt(2/m * kinetic) * z_half * cos(theta)
      # The 1/vpar factor becomes: 1/sqrt(2/m * kinetic)
      # Combined with dz = z_half * cos(theta) dtheta:
      # dz/vpar = z_half * cos(theta) / sqrt(2/m * kinetic)
      
      # Near bounce points, kinetic -> 0, but cos(theta) -> 0 too
      # The product cos(theta)/sqrt(kinetic) stays finite
      
      # For numerical stability, compute directly
      vpar_sq = two_over_mi * np.maximum(kinetic, 1e-30)
      vpar = np.sqrt(vpar_sq)
      jacobian = z_half * cos_th  # dz/dtheta
      
      # Source term (vectorized)
      zdep = np.exp(-(z / self.Lb)**2)
      vperp_sq = 2.0 * mu * self.B_p / self.mi
      vperp = np.sqrt(vperp_sq) if vperp_sq > 0 else 0.0
      vdep = np.exp(-((np.abs(vpar) - self.v_beam)**2 + 
               (vperp - self.v_beam)**2) / self.sigma_beam)
      source = self.gamma_0 * zdep * vdep
      
      # Return numerator and denominator integrands
      weight = jacobian / vpar
      return source * weight, weight
    
    # Use fixed_quad which is much faster than adaptive quad
    # Integrate from -pi/2 to pi/2
    def num_integrand(theta):
      num, _ = integrand_transformed(theta)
      return num
    
    def den_integrand(theta):
      _, den = integrand_transformed(theta)
      return den
    
    numerator, _ = fixed_quad(num_integrand, -np.pi/2, np.pi/2, n=n_quad)
    denominator, _ = fixed_quad(den_integrand, -np.pi/2, np.pi/2, n=n_quad)
    
    result = numerator / denominator if denominator > 1e-30 else 0.0
    # Set source to zero outside mirror throat region

    return result
  
  def _contour_integral_adaptive(self, E, mu, z_left, z_right):
    """Original adaptive quadrature method (slower but more accurate)."""
    def vpar_at_z(zp):
      return np.sqrt(max((2.0 / self.mi) * (E - mu * self.bmag(zp)), 1e-30))
    
    def integrand_num(zp):
      vp = vpar_at_z(zp)
      return self.spatial_beam_source(zp, vp, mu) / vp
    
    def integrand_den(zp):
      return 1.0 / vpar_at_z(zp)
    
    # Suppress integration warnings - near bounce points 1/vpar diverges but integral is still finite
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      numerator, _ = quad(integrand_num, z_left, z_right, limit=200, epsabs=1e-8, epsrel=1e-6)
      denominator, _ = quad(integrand_den, z_left, z_right, limit=200, epsabs=1e-8, epsrel=1e-6)
    
    return numerator / denominator if denominator > 1e-30 else 0.0
  
  def _compute_contour_parallel(self, Z, VPAR, MU, n_workers=None):
    """Compute contour integral in parallel using multiprocessing."""
    shape = Z.shape
    args_list = [(z, vpar, mu, self) for z, vpar, mu in 
           zip(Z.flatten(), VPAR.flatten(), MU.flatten())]
    
    if n_workers is None:
      n_workers = cpu_count()
    print(f"Using multiprocessing with {n_workers} workers for {len(args_list)} integrals...")
    with Pool(n_workers) as pool:
      results = pool.map(_contour_worker, args_list)
    print("Done.")
    
    return np.array(results).reshape(shape)
  
  def compute_all_sources_on_grid(self, z_coords, vpar_coords, mu_coords, n_workers=None):
    """
    Compute all three source models on a grid.
    
    Returns:
      dict with keys 'spatial', 'orbit_avg', 'contour' containing 3D arrays
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    
    spatial = self.spatial_beam_source(Z, VPAR, MU)
    orbit_avg = self.orbit_avg_beam_source(Z, VPAR, MU)
    
    # Parallel computation for contour integral
    contour = self._compute_contour_parallel(Z, VPAR, MU, n_workers)
    
    return {
      'Z': Z, 'VPAR': VPAR, 'MU': MU,
      'spatial': spatial,
      'orbit_avg': orbit_avg,
      'contour': contour
    }

  def compute_M0_fluxes(self, z_coords, vpar_coords, mu_coords, results=None, n_workers=None):
    """Compute integrated M0 fluxes for each source model.
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    Z = results['Z']
    jacobian = 2 * np.pi / self.mi
    
    fluxes = {}
    for key in ['spatial', 'orbit_avg', 'contour']:
      vals = results[key] * jacobian
      # Integrate over vpar, then z, then mu
      vals = np.trapz(vals, vpar_coords, axis=1)
      vals = np.trapz(vals, z_coords, axis=0)
      fluxes[key] = np.trapz(vals, mu_coords, axis=0)
    
    print("\nTotal M0 fluxes (particles per second per m^3):")
    print(f"  Spatial source:    {fluxes['spatial']:.3e} s^-1 m^-3")
    print(f"  Contour integral:  {fluxes['contour']:.3e} s^-1 m^-3")
    print(f"  Orbit-averaged:    {fluxes['orbit_avg']:.3e} s^-1 m^-3")

    # pgkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp ev "f[0] f[1] *" integ 0 info
    # pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel -c0 info

    return fluxes
  
  def compute_power(self, z_coords, vpar_coords, mu_coords, results=None, n_workers=None):
    """Compute total power (energy per second) deposited by each source model.
    
    Power = ∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    
    Returns:
      dict with keys 'spatial', 'orbit_avg', 'contour' containing power in Watts/m^3
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    Z = results['Z']
    VPAR = results['VPAR']
    MU = results['MU']
    
    # Particle energy: E = 0.5 * m * vpar^2 + mu * B(z)
    energy = 0.5 * self.mi * VPAR**2 + MU * self.bmag(Z)
    
    jacobian = 2 * np.pi / self.mi
    
    power = {}
    for key in ['spatial', 'orbit_avg', 'contour']:
      # Source * energy * jacobian
      vals = results[key] * energy * jacobian
      # Integrate over vpar, then z, then mu
      vals = np.trapz(vals, vpar_coords, axis=1)
      vals = np.trapz(vals, z_coords, axis=0)
      power[key] = np.trapz(vals, mu_coords, axis=0)
    
    print("\nTotal power deposited (energy per second per m^3):")
    print(f"  Spatial source:    {power['spatial']:.3e} W/m^3 = {power['spatial']/self.eV:.3e} eV/s/m^3")
    print(f"  Contour integral:  {power['contour']:.3e} W/m^3 = {power['contour']/self.eV:.3e} eV/s/m^3")
    print(f"  Orbit-averaged:    {power['orbit_avg']:.3e} W/m^3 = {power['orbit_avg']/self.eV:.3e} eV/s/m^3")
    
    # Also compute average energy per particle
    fluxes = self.compute_M0_fluxes(z_coords, vpar_coords, mu_coords, results, n_workers)
    print("\nAverage energy per injected particle:")
    for key in ['spatial', 'orbit_avg', 'contour']:
      if fluxes[key] > 0:
        avg_E = power[key] / fluxes[key]
        print(f"  {key:16s}: {avg_E/self.eV:.1f} eV")
    
    return power
  
  def plot_comparison_2d(self, z_coords, vpar_coords, mu_coords, results=None, mu_idx=0, n_workers=None):
    """Plot 2D comparison of all three source models.
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    
    Z = results['Z'][:,:,mu_idx]
    VPAR = results['VPAR'][:,:,mu_idx] / 1e6
    
    im1 = ax1.pcolormesh(Z, VPAR, results['spatial'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im1, ax=ax1)
    ax1.set_xlabel('z (m)')
    ax1.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax1.set_title('Spatial Beam Source')
    
    im2 = ax2.pcolormesh(Z, VPAR, results['contour'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im2, ax=ax2)
    ax2.set_xlabel('z (m)')
    ax2.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax2.set_title('Contour-Integrated Source')
    
    im3 = ax3.pcolormesh(Z, VPAR, results['orbit_avg'][:,:,mu_idx], 
               shading='auto', cmap='inferno')
    plt.colorbar(im3, ax=ax3)
    ax3.set_xlabel('z (m)')
    ax3.set_ylabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax3.set_title('Midplane-Mapped Source')
    
    plt.tight_layout()
    return fig
  
  def plot_vpar_cut(self, mu, z_cut=0.0, n_workers=None):
    """Plot 1D comparison at fixed z, with x-axis in energy (keV)."""
    vpar_coords = np.linspace(0, 2e6, 500)
    orbit_avg = self.orbit_avg_beam_source(z_cut, vpar_coords, mu)
    
    # Use parallel computation for contour integral
    Z = np.full_like(vpar_coords, z_cut)
    MU = np.full_like(vpar_coords, mu)
    contour = self._compute_contour_parallel(
      Z.reshape(-1, 1, 1), vpar_coords.reshape(-1, 1, 1), 
      MU.reshape(-1, 1, 1), n_workers).flatten()

    # Compute the integral for both sources
    contour_integral = np.trapz(contour, vpar_coords)
    orbit_avg_integral = np.trapz(orbit_avg, vpar_coords)
    print(f"At z={z_cut} m, mu={mu:.3e}:")
    print(f"  Contour integral M0 flux: {contour_integral:.3e}")
    print(f"  Orbit-averaged M0 flux:   {orbit_avg_integral:.3e}")
    
    # Convert vpar to energy in keV: E = 0.5 * m * vpar^2
    # Use signed energy to preserve direction information
    energy_keV = np.sign(vpar_coords) * 0.5 * self.mi * vpar_coords**2 / (1000 * self.eV) + mu * self.bmag(z_cut) / (1000 * self.eV)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(energy_keV, contour, label='Contour Integral', linestyle='-')
    ax.plot(energy_keV, orbit_avg, label='Midplane Mapped', linestyle='--')
    ax.set_xlabel(r'$E =\frac{1}{2} m v_\parallel^2 + \mu B $ (keV)')
    ax.set_ylabel('Source')
    ax.set_title(f'Source comparison at z = {z_cut} m')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return fig

# ============================================================================
# Script execution
# ============================================================================

if __name__ == "__main__":
  bsc = BeamSourceComparison()
  
  # z_coords = np.linspace(-0.6, 0.6, 100)
  # vpar_coords = np.linspace(-1.5e6, 1.5e6, 100)
  z_coords = np.linspace(-0.8, 0.8, 100)
  vpar_coords = np.linspace(-2e6, 2e6, 50)
  mu_val = 0.5 * bsc.mi * bsc.v_beam**2 / bsc.B_p
  mu_coords = np.array([mu_val])
  mu_coords_flux = np.linspace(0.0, 2.0 * mu_val, 50)
  
  # results = bsc.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords_flux)
  # fluxes = bsc.compute_M0_fluxes(z_coords, vpar_coords, mu_coords_flux, results)
  # power = bsc.compute_power(z_coords, vpar_coords, mu_coords_flux, results)
  
  # # Plot 2D comparison
  # fig1 = bsc.plot_comparison_2d(z_coords, vpar_coords, mu_coords)
  # plt.show()
  
  # Plot 1D vpar cut at z=0
  fig2 = bsc.plot_vpar_cut(mu_val, z_cut=0.0)
  plt.show()