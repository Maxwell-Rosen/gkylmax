"""
Beam source comparison: spatial vs orbit-averaged vs contour-integrated sources.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from multiprocessing import Pool, cpu_count


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
    
    def bmag(self, Z):
        """Compute magnetic field magnitude at position Z."""
        Z = np.atleast_1d(Z)
        L1 = 1.0 / (1 + ((Z - self.Z_m) / self.gamma_geo)**2)
        L2 = 1.0 / (1 + ((Z + self.Z_m) / self.gamma_geo)**2)
        Bmag = self.mcB / (np.pi * self.gamma_geo) * (L1 + L2)
        return Bmag.item() if Bmag.size == 1 else Bmag
    
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
        vpar_midplane = np.sqrt(vpar**2 + 2 * mu * (self.bmag(z) - self.bmag(0)) / self.mi)
        return self.midplane_beam_source(vpar_midplane, mu)
    
    def contour_integral_beam_source(self, z, vpar, mu, z_min=-0.98, z_max=0.98):
        """
        Compute orbit-averaged source by integrating along energy contours.
        <S> = (∫ S dz/|vpar|) / (∫ dz/|vpar|)
        """
        E = 0.5 * self.mi * vpar**2 + mu * self.bmag(z)
        
        # Find bounce points
        if mu > 0:
            B_bounce = E / mu
            z_search = np.linspace(z_min, z_max, 1000)
            B_search = self.bmag(z_search)
            accessible = np.where(B_search <= B_bounce)[0]
            if len(accessible) == 0:
                return 0.0
            z_left = z_search[accessible[0]]
            z_right = z_search[accessible[-1]]
        else:
            z_left, z_right = z_min, z_max
        
        def vpar_at_z(zp):
            return np.sqrt(max((2.0 / self.mi) * (E - mu * self.bmag(zp)), 1e-30))
        
        def integrand_num(zp):
            vp = vpar_at_z(zp)
            return self.spatial_beam_source(zp, vp, mu) / vp
        
        def integrand_den(zp):
            return 1.0 / vpar_at_z(zp)
        
        numerator, _ = quad(integrand_num, z_left, z_right, limit=100, epsabs=1e-10, epsrel=1e-10)
        denominator, _ = quad(integrand_den, z_left, z_right, limit=100, epsabs=1e-10, epsrel=1e-10)
        
        return numerator / denominator if denominator > 1e-30 else 0.0
    
    def compute_on_grid(self, z_coords, vpar_coords, mu_coords, n_workers=None):
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
    
    def _compute_contour_parallel(self, Z, VPAR, MU, n_workers=None):
        """Compute contour integral in parallel."""
        if n_workers is None:
            n_workers = cpu_count()
        
        shape = Z.shape
        args_list = [(z, vpar, mu, self) for z, vpar, mu in 
                     zip(Z.flatten(), VPAR.flatten(), MU.flatten())]
        
        print(f"Computing contour integrals with {n_workers} workers...")
        with Pool(n_workers) as pool:
            results = pool.map(_contour_worker, args_list)
        print("Done.")
        
        return np.array(results).reshape(shape)
    
    def compute_fluxes(self, z_coords, vpar_coords, results):
        """Compute integrated M0 fluxes for each source model."""
        Z = results['Z'][:,:,0]
        jacobian = self.bmag(Z)
        
        fluxes = {}
        for key in ['spatial', 'orbit_avg', 'contour']:
            vals = results[key][:,:,0] * jacobian
            fluxes[key] = np.trapz(np.trapz(vals, vpar_coords, axis=1), z_coords, axis=0)
        
        return fluxes
    
    def plot_comparison_2d(self, z_coords, vpar_coords, results, mu_idx=0):
        """Plot 2D comparison of all three source models."""
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
    
    def plot_vpar_cut(self, vpar_coords, mu, z_cut=0.0):
        """Plot 1D comparison at fixed z."""
        spatial = self.spatial_beam_source(z_cut, vpar_coords, mu)
        orbit_avg = self.orbit_avg_beam_source(z_cut, vpar_coords, mu)
        contour = np.array([self.contour_integral_beam_source(z_cut, vp, mu) 
                           for vp in vpar_coords])
        
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(vpar_coords/1e6, contour, label='Contour Integral', linestyle='-')
        ax.plot(vpar_coords/1e6, orbit_avg, label='Midplane Mapped', linestyle='--')
        ax.plot(vpar_coords/1e6, spatial, label='Spatial Source', linestyle=':')
        ax.set_xlabel(r'$v_\parallel$ ($10^6$ m/s)')
        ax.set_ylabel('Source')
        ax.set_title(f'Source comparison at z = {z_cut} m')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        return fig


def _contour_worker(args):
    """Worker function for parallel contour integration."""
    z, vpar, mu, obj = args
    return obj.contour_integral_beam_source(z, vpar, mu)


# ============================================================================
# Script execution
# ============================================================================

if __name__ == "__main__":
    # Create comparison object with default parameters
    bsc = BeamSourceComparison()
    
    # Define grid
    z_coords = np.linspace(-0.6, 0.6, 100)
    vpar_coords = np.linspace(-1.5e6, 1.5e6, 100)
    mu_val = 0.5 * bsc.mi * bsc.v_beam**2 / bsc.B_p
    mu_coords = np.array([mu_val])
    
    # Compute on grid
    results = bsc.compute_on_grid(z_coords, vpar_coords, mu_coords)
    
    # Compute and print fluxes
    fluxes = bsc.compute_fluxes(z_coords, vpar_coords, results)
    print("\nTotal M0 fluxes:")
    print(f"  Spatial source:    {fluxes['spatial']:.3e}")
    print(f"  Contour integral:  {fluxes['contour']:.3e}")
    print(f"  Orbit-averaged:    {fluxes['orbit_avg']:.3e}")
    
    # Plot 2D comparison
    fig1 = bsc.plot_comparison_2d(z_coords, vpar_coords, results)
    plt.show()
    
    # Plot 1D vpar cut at z=0
    vpar_fine = np.linspace(-1.5e6, 1.5e6, 500)
    fig2 = bsc.plot_vpar_cut(vpar_fine, mu_val, z_cut=0.0)
    plt.show()
