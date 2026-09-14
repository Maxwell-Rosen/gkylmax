The orbit-average factor of approximately 12 in the Maxwellian boundary/source plot is dominated by internal MPI interfaces included in the integrated boundary diagnostic. It is not a factor of 12 in the actual wall particle flux.

This audit uses the outputs present on September 10, 2026. No input, solver, or simulation output was changed. The current Maxwellian run is distinct from the older outputs now in `1x-maxw/no-crazy-time-dilation`.

At frame 44, within the fifth orbit-average phase, the reconstruction gives:

| Run | Saved diagnostic, 2 lower / S | Actual lower wall, doubled / S | Contribution from internal lower faces / S |
|---|---:|---:|---:|
| Beam | 0.321472 | 1.085579 | -0.764108 |
| Maxwellian | 12.447665 | 1.140327 | +11.307338 |
| Maxwellian, no extra FDP dilation | 13.085587 | 1.683993 | +11.401594 |

The convention of doubling the lower flux matches `1x-maxw/plot.sh`. It assumes approximate reflection symmetry; these are not independent measurements of both physical walls. The uppermost physical surface is not included in the saved arrays, which store each cell's lower surface. Summing the existing lower and upper integrated diagnostics does not remove their respective internal-interface contamination.

![Saved diagnostic and actual lower wall](wall_vs_diagnostic.png)

The red curves are the existing integrated diagnostics; the blue points reconstruct the actual lower wall from saved surface flux snapshots. Gray intervals are OAPs. The horizontal axis is the diagnostic sample index, matching the original plots, rather than physical elapsed time. FDP curves are sampled more coarsely by the saved distribution frames; transient features between blue points are not reconstructed.

The implementation chain is explicit:

1. [`gk_species_bflux.c:662`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gk_species_bflux.c:662) sets boundary ranges from each rank's **local** skin and ghost cells, so it calculates internal subdomain faces as well as external walls.
2. [`gyrokinetic.c:254`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gyrokinetic.c:254) makes `comm_plane[0]` the full spatial communicator in 1D.
3. [`gk_species_bflux.c:325`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gk_species_bflux.c:325) sums the local lower-face moments across that communicator. Rank zero then writes this sum to the lower boundary diagnostic. The upper diagnostic has the analogous behavior.

For the eight ranks, each containing 50 of the 400 spatial cells, the lower diagnostic therefore measures

\[
R_{\rm saved}=\frac{2}{S}\sum_{r=0}^{7}\Gamma_{r,\rm lower}
=\frac{2\Gamma_{\rm wall,lower}}{S}
+\frac{2}{S}\sum_{r=1}^{7}\Gamma_{r,\rm lower}.
\]

The second term is an unintended contribution if this file is interpreted as wall loss. It sums one oriented face per subdomain, not both sides of each internal interface. There is no general identity making that sum vanish. Cancellation in these symmetric simulations is only approximate.

The reconstruction uses components 0 through 5 of each saved `ion_collisionless_surf_flux` array, which are the six nodal values on the lower spatial surface. The DG boundary kernel has

\[
G_0=\sum_{a=0}^{5}w_a F_a,\qquad
w=(5/18,4/9,5/18,5/18,4/9,5/18).
\]

Combining that kernel, the M0 moment kernel, and configuration integration gives the outward lower particle flux

\[
\Gamma_{j,\rm lower}=-\frac{\pi\Delta v_c\Delta\mu_c}{m_i}
\sum_{k_v,k_\mu}G_0(j,k_v,k_\mu).
\]

Here the velocity widths are computational widths; the saved flux already includes the mapped-velocity Jacobians. See [`dg_gyrokinetic_boundary_surfx_1x2v_ser_p1.c`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/ker/dg_gyrokinetic/dg_gyrokinetic_boundary_surfx_1x2v_ser_p1.c:24) and [`gyrokinetic_mom_1x2v_ser_p1.c`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/ker/gyrokinetic_mom/gyrokinetic_mom_1x2v_ser_p1.c:123).

Taking only spatial index zero gives the physical lower wall; summing indices `0, 50, 100, ..., 350` reproduces the saved diagnostic. Across the checked interior OAP snapshots, the largest absolute discrepancy in the source-normalized ratio is `2.33e-8`. At frame 44 the differences are below `2e-10`. Small differences near switches or during FDPs are expected: the saved surface array is an instantaneous RK-stage quantity, while the integrated diagnostic is RK-weighted.

The source dependence is visible directly in the internal-face terms. At frame 44, the normalized contributions `2 Gamma_lower / S` at zero-based spatial indices 150 and 250 are:

| Run | Face 150 | Face 250 | Sum of these two terms |
|---|---:|---:|---:|
| Beam | +17.3116 | -18.1272 | -0.8156 |
| Maxwellian | -238.7826 | +250.0339 | +11.2513 |
| Maxwellian, no extra FDP dilation | -240.2254 | +251.5442 | +11.3187 |

These are internal streaming fluxes. Their large, oppositely signed values leave a residual that nearly explains the entire plotted discrepancy. The Maxwellian changes both their magnitude and direction. A roughly five-percent imperfect cancellation of terms of size 250 leaves a term of size 11. The beam's corresponding residual is much smaller and negative, which also explains its low and sometimes negative OAP diagnostic values. The log plot cannot display negative values, so its downward features should not automatically be interpreted as tiny positive wall losses.

Both sources have total particle injection approximately `3.51344e20` in the code's integrated normalization, and the Maxwellian was tuned to match the beam's energy injection as well. Matching these two integrals does not match pitch-angle distribution, spatial source profile, pressure anisotropy, or internal streaming response. The beam is narrow around nonzero parallel and perpendicular speed; the broad Maxwellian populates low magnetic moment and the passing region as well.

This distinction matters especially for the implemented OAP equation:

\[
\partial_t F=M_T\{\alpha A[F]+C[F]+S_F\},\qquad\alpha=2\times10^{-5},
\]

where `M_T` is one in cells classified as trapped and zero otherwise. OAP damping is disabled. The whole RHS, including collisions and source, is multiplied by this mask after boundary moments are calculated. The mask is cellwise: any passing node makes the whole cell inactive. The nominal source diagnostic integrates `S_F`, not `M_T S_F`.

At frame 44, applying the saved mask to the zeroth modal coefficient of the saved source retains approximately 100% of beam injection and 93.62% of Maxwellian injection (93.69% for the older Maxwellian). This is evidence of a real source/mask difference, **not** a factor-of-ten source normalization error. The six-percent suppressed fraction cannot by itself explain the plotted twelvefold ratio.

Internal *unscaled* streaming fluxes during this artificial evolution do not obey the same balance as equilibrium wall losses. Streaming is strongly suppressed in the evolved equation, and cells in the passing region are frozen; collisions and source reshape the active trapped population. The source's velocity-space structure changes that population and hence its internal streaming fluxes. The table above measures the resulting sensitivity directly; it does not assume an analytical prediction for their exact amplitudes.

Even after restricting diagnostics to the wall, an OAP wall/source ratio is not an instantaneous particle-conservation check. In the saved Maxwellian frames 10, 11, 12, and 14, the entire lower skin-cell distribution is identical, the boundary mask is zero everywhere, and its saved spatial surface flux is unchanged. Its true doubled wall ratio stays at `0.6720576`, whereas the integrated diagnostic goes from about `0.67849` to `12.16661`. The wall distribution is being carried forward from the preceding FDP. A ratio near one at that frozen wall does not establish equilibrium of the trapped population.

The relevant code ordering is [`gk_species.c:119`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gk_species.c:119), [`gyrokinetic.c:2376`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gyrokinetic.c:2376), and [`loss_cone_mask_gyrokinetic.c:108`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/zero/loss_cone_mask_gyrokinetic.c:108). The Boltzmann field uses local boundary moments and broadcasts the physical endpoint values, rather than consuming this all-reduced diagnostic file; see [`gk_field_boltzmann.c:25`](/global/homes/m/mhrosen/gkeyll/gyrokinetic/apps/gk_field_boltzmann.c:25).

The extra FDP dilation is a separate operation. It can change the state inherited by the following OAP, so histories need not coincide, but the internal-face contamination occurs with either FDP configuration. The old Maxwellian run explicitly confirms this. The OAP always resets to the loss-cone multiplier and the same collisionless scaling.

For a wall diagnostic, the targeted correction is to include only global physical boundary faces in the reduction. In 1D, each endpoint has one owning spatial rank; its local integrated flux is sufficient. An alternative is to zero contributions on non-endpoint ranks before a full-communicator sum. Any corresponding time-integrated diagnostic should receive the same treatment. A general communicator change needs care because `comm_plane` has other consumers. This audit does not change Gkeyll or claim to resolve other possible discretization/conservation errors.

The build metadata are `53d3093b2ca3` for the beam, `93ec1a63ac59` for the current Maxwellian, and `8a3d2ccb776b` for the older Maxwellian. The problematic 1D communicator construction exists in all three revisions.

Reproduce the CSV and figures with:

```bash
/global/homes/m/mhrosen/.conda/envs/pgkyl/bin/python a-thesis-simulations/bflux-diagnostic-audit/audit.py
```

Detailed per-rank values and reconstruction differences are in [flux_comparison.csv](flux_comparison.csv). The plot is also available as [PDF](wall_vs_diagnostic.pdf).
