"""
Beam source comparison: spatial vs orbit-averaged vs contour-integrated sources.
"""
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize
from multiprocessing import Pool, cpu_count
import matplotlib

matplotlib.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.size': 12,
    'axes.titlesize': 18,
    'axes.labelsize': 22,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
})

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
         gamma_spatial=5476,
         gamma_oavg=1448,
         T_beam_eV_spatial=200,
         T_beam_eV_oavg=200,
         E_beam_eV_spatial=25000,
         E_beam_eV_oavg=25000,
         Lb=0.2):
    """
    Initialize beam source comparison.
    
    Parameters:
      mcB : float - Magnetic field parameter
      gamma_geo : float - Lorentzian width parameter
      Z_m : float - Mirror throat location (m)
      B_p : float - Reference magnetic field (T)
      gamma_spatial : float - Spatial source normalization
      gamma_oavg : float - Orbit-averaged source normalization
      T_beam_eV_spatial : float - Beam temperature for spatial source (eV)
      T_beam_eV_oavg : float - Beam temperature for orbit-averaged source (eV)
      E_beam_eV_spatial : float - Beam energy for spatial source (eV)
      E_beam_eV_oavg : float - Beam energy for orbit-averaged source (eV)
      Lb : float - Spatial source width (m)
    """
    # Geometry
    self.mcB = mcB
    self.gamma_geo = gamma_geo
    self.Z_m = Z_m
    self.B_p = B_p
    
    # Beam parameters
    self.gamma_spatial = gamma_spatial  # Same normalization for spatial source
    self.gamma_oavg = gamma_oavg
    self.T_beam_spatial = T_beam_eV_spatial * self.eV
    self.T_beam_oavg = T_beam_eV_oavg * self.eV
    self.E_beam_spatial = E_beam_eV_spatial * self.eV
    self.E_beam_oavg = E_beam_eV_oavg * self.eV
    self.v_beam_spatial = np.sqrt(self.E_beam_spatial / self.mi)
    self.v_beam_oavg = np.sqrt(self.E_beam_oavg / self.mi)
    self.sigma_beam_spatial = 2 * self.T_beam_spatial / self.mi
    self.sigma_beam_oavg = 2 * self.T_beam_oavg / self.mi
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
  
  def midplane_beam_source(self, vpar, mu, gamma_beam, v_beam, sigma_beam):
    """Beam source evaluated at midplane velocities."""
    vperp = np.sqrt(2.0 * mu * self.B_p / self.mi)
    return gamma_beam * np.exp(-((np.abs(vpar) - v_beam)**2 + 
                      (vperp - v_beam)**2) / sigma_beam)
  
  def spatial_beam_source(self, z, vpar, mu):
    """Spatial beam source with Gaussian z-dependence (Dorf et al. 2025)."""
    zdep = np.exp(-(z / self.Lb)**2)
    vdep = self.midplane_beam_source(vpar, mu, self.gamma_spatial, self.v_beam_spatial, self.sigma_beam_spatial)
    return zdep * vdep
  
  def orbit_avg_beam_source(self, z, vpar, mu):
    """Orbit-averaged source: map vpar to midplane using energy conservation."""
    # Set source to zero outside mirror throat region
    z = np.atleast_1d(z)
    result = self.midplane_beam_source(np.sqrt(vpar**2 + 2 * mu * (self.bmag(z) - self.bmag(0)) / self.mi), mu, 
                                      self.gamma_oavg, self.v_beam_oavg, self.sigma_beam_oavg)
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
      vpar_sq = two_over_mi * np.maximum(kinetic, 1e-30)
      vpar = np.sqrt(vpar_sq)
      jacobian = z_half * cos_th  # dz/dtheta
      source = self.spatial_beam_source(z, vpar, mu)
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

  def _integrate_M0_from_3d(self, source_3d, z_coords, vpar_coords, mu_coords):
    """Helper: compute M0 flux from precomputed 3D source array.
    
    M0 = ∫∫∫ S(z, vpar, mu) * jacobian dz dvpar dmu
    """
    jacobian = 2 * np.pi / self.mi
    vals = source_3d * jacobian
    # Integrate over vpar, then z, then mu
    vals = np.trapz(vals, vpar_coords, axis=1)
    vals = np.trapz(vals, z_coords, axis=0)
    return np.trapz(vals, mu_coords, axis=0)
  
  def _integrate_power_from_3d(self, source_3d, Z_grid, VPAR_grid, MU_grid, 
                               z_coords, vpar_coords, mu_coords):
    """Helper: compute power from precomputed 3D source array.
    
    Power = ∫∫∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    """
    energy = 0.5 * self.mi * VPAR_grid**2 + MU_grid * self.bmag(Z_grid)
    jacobian = 2 * np.pi / self.mi
    vals = source_3d * energy * jacobian
    # Integrate over vpar, then z, then mu
    vals = np.trapz(vals, vpar_coords, axis=1)
    vals = np.trapz(vals, z_coords, axis=0)
    return np.trapz(vals, mu_coords, axis=0)

  def compute_M0_orbit_avg(self, z_coords, vpar_coords, mu_coords):
    """Compute integrated M0 flux for orbit-averaged source only (3D integral).
    
    M0 = ∫∫∫ S(z, vpar, mu) * jacobian dz dvpar dmu
    
    This is a lightweight version that only computes the orbit-averaged source,
    useful for optimization loops where we don't need all three sources.
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    orbit_avg = self.orbit_avg_beam_source(Z, VPAR, MU)
    return self._integrate_M0_from_3d(orbit_avg, z_coords, vpar_coords, mu_coords)
  
  def compute_power_orbit_avg(self, z_coords, vpar_coords, mu_coords):
    """Compute total power for orbit-averaged source only (3D integral).
    
    Power = ∫∫∫ S(z, vpar, mu) * E(vpar, mu, z) * jacobian dz dvpar dmu
    where E = 0.5 * m * vpar^2 + mu * B(z) is the particle kinetic energy.
    
    This is a lightweight version that only computes the orbit-averaged source,
    useful for optimization loops where we don't need all three sources.
    """
    Z, VPAR, MU = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    orbit_avg = self.orbit_avg_beam_source(Z, VPAR, MU)
    return self._integrate_power_from_3d(orbit_avg, Z, VPAR, MU, z_coords, vpar_coords, mu_coords)

  def compute_M0_fluxes(self, z_coords, vpar_coords, mu_coords, results=None, n_workers=None):
    """Compute integrated M0 fluxes for each source model (3D integral over z, vpar, mu).
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    fluxes = {}
    for key in ['spatial', 'orbit_avg', 'contour']:
      fluxes[key] = self._integrate_M0_from_3d(results[key], z_coords, vpar_coords, mu_coords)
    
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
    
    power = {}
    for key in ['spatial', 'orbit_avg', 'contour']:
      power[key] = self._integrate_power_from_3d(results[key], Z, VPAR, MU, 
                                                  z_coords, vpar_coords, mu_coords)
    
    print("\nTotal power deposited (energy per second per m^3):")
    print(f"  Spatial source:    {power['spatial']:.3e} W/m^3 = {power['spatial']/self.eV:.3e} eV/s/m^3")
    print(f"  Contour integral:  {power['contour']:.3e} W/m^3 = {power['contour']/self.eV:.3e} eV/s/m^3")
    print(f"  Orbit-averaged:    {power['orbit_avg']:.3e} W/m^3 = {power['orbit_avg']/self.eV:.3e} eV/s/m^3")
    
    # Also compute average energy per particle
    fluxes = self.compute_M0_fluxes(z_coords, vpar_coords, mu_coords, results, n_workers)
    print("\nAverage energy per injected particle:")
    print(f"  Spatial source:    {power['spatial']/fluxes['spatial']:.3e} J = {(power['spatial']/fluxes['spatial'])/self.eV:.3e} eV")
    print(f"  Contour integral:  {power['contour']/fluxes['contour']:.3e} J = {(power['contour']/fluxes['contour'])/self.eV:.3e} eV")
    print(f"  Orbit-averaged:    {power['orbit_avg']/fluxes['orbit_avg']:.3e} J = {(power['orbit_avg']/fluxes['orbit_avg'])/self.eV:.3e} eV")
    
    return power
  
  def plot_comparison_2d(self, z_coords, vpar_coords, mu_coords, results=None, mu_idx=0, n_workers=None):
    """Plot 2D comparison of all three source models.
    
    If results is None, compute_on_grid will be called internally.
    """
    if results is None:
      results = self.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords, n_workers)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 6))
    
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
  
  def optimize_orbit_avg_to_match_contour(self, z_coords, vpar_coords, mu_coords,
                                          n_workers=None, verbose=True):
    """
    Optimize orbit_avg_beam_source parameters to match contour_integral_beam_source.
    
    Creates a NEW BeamSourceComparison object with optimized parameters (E_beam, 
    gamma_oavg, T_beam) such that orbit_avg_beam_source best matches the 
    contour integral of the ORIGINAL (self) spatial source.
    
    This optimizes for the FULL 3D integrals over z, vpar, and mu, matching
    the total M0 flux and power computed by compute_M0_fluxes() and compute_power().
    
    The original object (self) is NOT modified.
    
    Args:
      z_coords: Array of z coordinates for the 3D grid
      vpar_coords: Array of vpar coordinates for the 3D grid
      mu_coords: Array of mu coordinates for the 3D grid
      n_workers: Number of parallel workers for contour integral
      verbose: Print optimization progress
    
    Returns:
      BeamSourceComparison: New object with optimized parameters
    """
    n_z = len(z_coords)
    n_vpar = len(vpar_coords)
    n_mu = len(mu_coords)
    
    if verbose:
      print(f"Computing contour integral on 3D grid ({n_z} x {n_vpar} x {n_mu} = {n_z*n_vpar*n_mu} points)...")
    
    # Create 3D grids
    Z_grid, VPAR_grid, MU_grid = np.meshgrid(z_coords, vpar_coords, mu_coords, indexing='ij')
    
    # Compute contour integral for all points (this is expensive but only done once)
    contour_target_3d = self._compute_contour_parallel(
      Z_grid.reshape(-1, 1, 1), 
      VPAR_grid.reshape(-1, 1, 1), 
      MU_grid.reshape(-1, 1, 1), 
      n_workers
    ).reshape(n_z, n_vpar, n_mu)
    
    # Precompute target M0 and power using class helper methods (consistent with compute_M0_fluxes/compute_power)
    M0_contour_total = self._integrate_M0_from_3d(contour_target_3d, z_coords, vpar_coords, mu_coords)
    power_contour_total = self._integrate_power_from_3d(contour_target_3d, Z_grid, VPAR_grid, MU_grid,
                                                        z_coords, vpar_coords, mu_coords)
    
    if verbose:
      print(f"  Target M0 flux (3D integral): {M0_contour_total:.3e}")
      print(f"  Target power (3D integral):   {power_contour_total:.3e} W = {power_contour_total/self.eV:.3e} eV/s")
    
    # Initial parameter guesses in physical units (eV)
    E_beam_init_eV = self.E_beam_oavg / self.eV
    T_beam_init_eV = self.T_beam_oavg / self.eV
    x0 = np.array([E_beam_init_eV, self.gamma_oavg, T_beam_init_eV])
    
    # Parameter bounds
    bounds = [
      (0.5 * E_beam_init_eV, 2.0 * E_beam_init_eV),
      (0.1 * self.gamma_oavg, 10.0 * self.gamma_oavg),
      (0.1 * T_beam_init_eV, 10.0 * T_beam_init_eV),
    ]
    
    def objective(params):
      """Objective function: combined residual of shape, M0 flux, and power (3D integrals)."""
      E_beam_test_eV, gamma_oavg_test, T_beam_test_eV = params
      
      # Create temporary object with test parameters
      temp_obj = BeamSourceComparison(
        mcB=self.mcB, gamma_geo=self.gamma_geo, Z_m=self.Z_m, B_p=self.B_p,
        gamma_oavg=gamma_oavg_test,
        T_beam_eV_oavg=T_beam_test_eV,
        E_beam_eV_oavg=E_beam_test_eV
      )
      
      # Compute orbit_avg on 3D grid (this is fast - no contour integration needed)
      orbit_avg_3d = temp_obj.orbit_avg_beam_source(Z_grid, VPAR_grid, MU_grid)
      
      # 1. Shape residual (normalized by max value)
      max_val = max(np.max(np.abs(contour_target_3d)), 1e-30)
      shape_residual = np.sum(((orbit_avg_3d - contour_target_3d) / max_val)**2) / (n_z * n_vpar * n_mu)
      
      # 2. M0 flux residual using compute_M0_orbit_avg (consistent with compute_M0_fluxes)
      M0_orbit_avg_total = temp_obj.compute_M0_orbit_avg(z_coords, vpar_coords, mu_coords)
      M0_ref = max(abs(M0_contour_total), 1e-30)
      M0_residual = ((M0_orbit_avg_total - M0_contour_total) / M0_ref)**2
      
      # 3. Power residual using compute_power_orbit_avg (consistent with compute_power)
      power_orbit_avg_total = temp_obj.compute_power_orbit_avg(z_coords, vpar_coords, mu_coords)
      power_ref = max(abs(power_contour_total), 1e-30)
      power_residual = ((power_orbit_avg_total - power_contour_total) / power_ref)**2
      
      total_residual = shape_residual + M0_residual + power_residual
      
      return total_residual, shape_residual, M0_residual, power_residual
    
    def objective_scalar(params):
      return objective(params)[0]
    
    if verbose:
      print("Optimizing orbit-averaged parameters...")
      print(f"  Initial: E_beam_oavg={x0[0]:.1f} eV, gamma_oavg={x0[1]:.3e}, T_beam_oavg={x0[2]:.1f} eV")
      init_total, init_shape, init_M0, init_power = objective(x0)
      print(f"  Initial residuals: shape={init_shape:.3e}, M0={init_M0:.3e}, power={init_power:.3e}, total={init_total:.3e}")
    
    result = minimize(objective_scalar, x0, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': 1000, 'ftol': 1e-16, 'gtol': 1e-16})
    
    E_beam_opt_eV, gamma_oavg_opt, T_beam_opt_eV = result.x
    
    if verbose:
      print(f"  Optimized: E_beam_oavg={E_beam_opt_eV:.1f} eV, gamma_oavg={gamma_oavg_opt:.3e}, T_beam_oavg={T_beam_opt_eV:.1f} eV")
      print(f"  Optimization converged: {result.success}, iterations: {result.nit}, message: {result.message}")
      final_total, final_shape, final_M0, final_power = objective(result.x)
      print(f"  Final residuals: shape={final_shape:.3e}, M0={final_M0:.3e}, power={final_power:.3e}, total={final_total:.3e}")
    
    if verbose:
      print(f"\nOptimized orbit-averaged beam parameters:")
      print(f"  E_beam_oavg: {self.E_beam_oavg/self.eV:.1f} eV -> {E_beam_opt_eV:.1f} eV")
      print(f"  T_beam_oavg: {self.T_beam_oavg/self.eV:.1f} eV -> {T_beam_opt_eV:.1f} eV")
      print(f"  gamma_oavg: {self.gamma_oavg:.3e} -> {gamma_oavg_opt:.3e}")
      print(f"\nSpatial source parameters (unchanged):")
      print(f"  E_beam_spatial: {self.E_beam_spatial/self.eV:.1f} eV")
      print(f"  T_beam_spatial: {self.T_beam_spatial/self.eV:.1f} eV")
      print(f"  gamma_spatial: {self.gamma_spatial:.3e}")
    
    # Create new object with optimized parameters
    optimized = BeamSourceComparison(
      mcB=self.mcB,
      gamma_geo=self.gamma_geo,
      Z_m=self.Z_m,
      B_p=self.B_p,
      gamma_spatial=self.gamma_spatial,
      gamma_oavg=gamma_oavg_opt,
      T_beam_eV_spatial=self.T_beam_spatial / self.eV,
      T_beam_eV_oavg=T_beam_opt_eV,
      E_beam_eV_spatial=self.E_beam_spatial / self.eV,
      E_beam_eV_oavg=E_beam_opt_eV,
      Lb=self.Lb
    )
    
    return optimized
  
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
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(energy_keV, contour, label='Contour Integral', linestyle='-')
    ax.plot(energy_keV, orbit_avg, label='Midplane Mapped', linestyle='--')
    ax.set_xlabel(r'$E =\frac{1}{2} m v_\parallel^2 + \mu B $ (keV)')
    ax.set_ylabel('Source')
    ax.set_title(f'Source comparison at z = {z_cut} m')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return fig

  def plot_z_slice(self, z_cut=0.0, n_vpar=100, n_mu=100, n_workers=None):
    """Plot 2D comparison at fixed z: contour integral, orbit-averaged, and difference.
    
    Args:
      z_cut: z position for the slice
      n_vpar: Number of vpar grid points
      n_mu: Number of mu grid points
      n_workers: Number of parallel workers for contour integral
    
    Returns:
      fig: Matplotlib figure with 3 panels
    """
    # Create vpar and mu grids
    vpar_coords = np.linspace(0, 2e6, n_vpar)
    mu_max = 2.0 * 0.5 * self.mi * self.v_beam_spatial**2 / self.B_p
    mu_coords = np.linspace(0, mu_max, n_mu)
    
    VPAR, MU = np.meshgrid(vpar_coords, mu_coords, indexing='ij')
    Z = np.full_like(VPAR, z_cut)
    
    # Compute orbit-averaged source (fast, vectorized)
    orbit_avg = self.orbit_avg_beam_source(Z, VPAR, MU)
    
    # Compute contour integral (parallel)
    print(f"Computing contour integral at z={z_cut} m...")
    contour = self._compute_contour_parallel(
      Z.reshape(-1, 1, 1), 
      VPAR.reshape(-1, 1, 1), 
      MU.reshape(-1, 1, 1), 
      n_workers
    ).reshape(n_vpar, n_mu)
    
    # Convert vpar to 10^6 m/s for plotting
    vpar_plot = VPAR / 1e6
    # Convert mu to mu*B in keV for y-axis
    B_at_z = self.bmag(z_cut)
    mu_plot = MU * B_at_z / (1000 * self.eV)
    
    # Compute difference
    diff = orbit_avg - contour
    
    # Create figure with 3 panels
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14, 4))
    
    # Common color scale for contour and orbit_avg
    vmax = max(np.max(contour), np.max(orbit_avg))
    vmin = 0
    
    # Left panel: Contour integral
    im1 = ax1.pcolormesh(vpar_plot, mu_plot, contour, shading='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(im1, ax=ax1, label='Source')
    ax1.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax1.set_ylabel(r'$\mu B$ (keV)')
    ax1.set_title(f'Contour Integral at z = {z_cut} m')
    
    # Middle panel: Orbit-averaged
    im2 = ax2.pcolormesh(vpar_plot, mu_plot, orbit_avg, shading='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(im2, ax=ax2, label='Source')
    ax2.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax2.set_ylabel(r'$\mu B$ (keV)')
    ax2.set_title(f'Orbit-Averaged at z = {z_cut} m')
    
    # Right panel: Difference (orbit_avg - contour)
    diff_max = max(abs(np.min(diff)), abs(np.max(diff)))
    im3 = ax3.pcolormesh(vpar_plot, mu_plot, diff, shading='auto', cmap='RdBu_r', 
                         vmin=-diff_max, vmax=diff_max)
    plt.colorbar(im3, ax=ax3, label='Difference')
    ax3.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
    ax3.set_ylabel(r'$\mu B$ (keV)')
    ax3.set_title(f'Orbit-Averaged − Contour at z = {z_cut} m')
    
    plt.tight_layout()
    return fig

# ============================================================================
# Script execution
# ============================================================================

if __name__ == "__main__":
  bsc = BeamSourceComparison()
  
  # z_coords = np.linspace(-0.6, 0.6, 100)
  # vpar_coords = np.linspace(-1.5e6, 1.5e6, 100)
  z_coords = np.linspace(-0.8, 0.8, 50)
  vpar_coords = np.linspace(-2e6, 2e6, 50)
  mu_val = 0.5 * bsc.mi * bsc.v_beam_spatial**2 / bsc.B_p
  mu_coords = np.linspace(0.0, 2.0 * mu_val, 50)
  
  # results = bsc.compute_all_sources_on_grid(z_coords, vpar_coords, mu_coords)
  # power = bsc.compute_power(z_coords, vpar_coords, mu_coords, results)
  # fluxes = bsc.compute_M0_fluxes(z_coords, vpar_coords, mu_coords, results)
  
  # # Plot 2D comparison
  # fig1 = bsc.plot_comparison_2d(z_coords, vpar_coords, mu_val)
  # plt.show()
  
  # Plot 1D vpar cut at z=0 (original)
  # fig2 = bsc.plot_vpar_cut(mu_val, z_cut=0.0)
  # fig2.suptitle('Before optimization')
  
  # Plot 2D slice at z=0 (before optimization)
  fig_slice = bsc.plot_z_slice(z_cut=0.0, n_vpar=100, n_mu=100)
  fig_slice.suptitle('Before optimization')
  
  # Optimize using full 3D integrals
  bsc_optimized = bsc.optimize_orbit_avg_to_match_contour(z_coords, vpar_coords, mu_coords)
  
  # fig3 = bsc_optimized.plot_vpar_cut(mu_val, z_cut=0.0)

  bsc_optimized.compute_power(z_coords, vpar_coords, mu_coords)
  fig3 = bsc_optimized.plot_z_slice(z_cut=0.0, n_vpar=100, n_mu=100)
  fig3.suptitle('After optimization')
  plt.show()