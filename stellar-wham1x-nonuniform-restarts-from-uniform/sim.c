#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_mirror_geo.h>
#include <gkyl_math.h>

#include <rt_arg_parse.h>

// Neccisary for the boltzmann electron density reading
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_position_map.h>
#include <gkyl_velocity_map.h>

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
  int cdim, vdim; // Dimensionality.
  // Plasma parameters
  double mi;
  double qi;
  double me;
  double qe;
  double Te0;
  double n0;
  double B_p;
  double beta;
  double tau;
  double Ti0;
  double kperpRhos;
  // Parameters controlling initial conditions.
  double alim;
  double nuFrac;
  // Electron-electron collision freq.
  double logLambdaElc;
  double nuElc;
  double elc_nuFrac;
  // Ion-ion collision freq.
  double logLambdaIon;
  double nuIon;
  // Thermal speeds.
  double vti;
  double vte;
  double c_s;
  // Gyrofrequencies and gyroradii.
  double omega_ci;
  double rho_s;
  double kperp; // Perpendicular wavenumber in SI units.
  double RatZeq0; // Radius of the field line at Z=0.
  // Axial coordinate Z extents. Endure that Z=0 is not on
  double z_min;
  double z_max;
  double psi_min;
  double psi_eval;
  double psi_max;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double t_end;
  int num_frames;
  double write_phase_freq; // Frequency of writing phase-space data.
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude;
  double ion_source_sigma;
  double ion_source_temp;

  // Boltzmann electron reading
  struct gkyl_array *f_ion;
  struct gkyl_array *jacobtot;
  struct gkyl_array *ion_jacobvel;
  struct gkyl_rect_grid f_ion_grid;
  struct gkyl_basis phase_basis;
  struct gkyl_range phase_local;
  struct gkyl_range phase_local_ext;
  struct gkyl_position_map *position_map;
  double target_z_fa;
  double tartget_vpar_fa;
  double target_mu_fa;
};

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./6.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_ion*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_ion*tan(cvpar*b)/tan(b);
  }
  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

// Evaluate collision frequencies
void
evalNuIon(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuIon;
}

void
load_ion_donor(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  struct gkyl_rect_grid mc2nu_pos_grid, f_ion_grid, jacobtot_grid, ion_jacobvel_grid;
  struct gkyl_array *mc2nu_pos, *f_ion, *jacobtot, *ion_jacobvel;

  mc2nu_pos = gkyl_grid_array_new_from_file(&mc2nu_pos_grid, 
    "../initial-conditions/boltz-elc-288z-nu2000/gk_wham-mc2nu_pos.gkyl");
  f_ion     = gkyl_grid_array_new_from_file(&f_ion_grid, 
    "../initial-conditions/boltz-elc-288z-nu2000/gk_wham-ion_400.gkyl");
  jacobtot = gkyl_grid_array_new_from_file(&jacobtot_grid,
    "../initial-conditions/boltz-elc-288z-nu2000/gk_wham-jacobtot.gkyl");
  ion_jacobvel = gkyl_grid_array_new_from_file(&ion_jacobvel_grid,
    "../initial-conditions/boltz-elc-288z-nu2000/gk_wham-ion_jacobvel.gkyl");
  

  app->f_ion = f_ion;
  app->f_ion_grid = f_ion_grid;
  app->jacobtot = jacobtot;
  app->ion_jacobvel = ion_jacobvel;

  int lower_cell[] = {1};
  int upper_cell[] = {mc2nu_pos_grid.cells[0]};

  int poly_order = 1;
  int cdim = 1;
  int vdim = 2;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  struct gkyl_basis phase_basis;
  gkyl_cart_modal_gkhybrid(&phase_basis, cdim, vdim);
  app->phase_basis = phase_basis;

  int nghost[] = { 0, 0, 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&mc2nu_pos_grid, nghost, &local, &local_ext);

  struct gkyl_range phase_local, phase_local_ext;
  gkyl_create_grid_ranges(&f_ion_grid, nghost, &phase_local, &phase_local_ext);
  app->phase_local = phase_local;
  app->phase_local_ext = phase_local_ext;

  // Create a position map object
  struct gkyl_position_map_inp pmap_inp = { };
  // Potential future issue by using the M0 ranges and basis for the position map
  struct gkyl_position_map *gpm = gkyl_position_map_new(pmap_inp, 
    mc2nu_pos_grid, local, local_ext, local, local_ext, basis);
  gkyl_position_map_set_mc2nu(gpm, mc2nu_pos);
  app->position_map = gpm;
}

static double
invert_position_map_func(double x, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double xc[] = {x};
  double xfa[3];
  gkyl_position_map_eval_mc2nu(app->position_map, xc, xfa);
  return xfa[0] - app->target_z_fa;
}

static double
invert_velocity_map_vpar_func(double vpar, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double xc[] = {vpar, 0.0};
  double xfa[3];
  mapc2p_vel_ion(0, xc, xfa, ctx);
  return xfa[0] - app->tartget_vpar_fa;
}

static double
invert_velocity_map_mu_func(double mu, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double xc[] = {0.0, mu};
  double xfa[3];
  mapc2p_vel_ion(0, xc, xfa, ctx);
  return xfa[1] - app->target_mu_fa;
}

void
read_ion_distf(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;

  // Must convert the point xn from cartesian to polar
  double z_field_aligned = xn[0];
  double vpar_field_aligned = xn[1];
  double mu_field_aligned = xn[2];

  // I have concerns about the non-uniform velocity grid
  // Must use the same one in the input file
  // Safe to just assume that it's comming from a uniform space grid
  // Assume same non-uniform velocity grid, so I can use the functions in this file rather than the velocity map.
  // The grid of the distribution function will have the -1 to 1 limits, so we must map from these coordinates to non-uniform coordinates
  
  // First we must determine the computational coordinate in non-uniform space
  app->target_z_fa = z_field_aligned;
  app->tartget_vpar_fa = vpar_field_aligned;
  app->target_mu_fa = mu_field_aligned;

  double interval_lower = -M_PI;
  double interval_upper = M_PI;
  double interval_lower_eval = invert_position_map_func(interval_lower, ctx);
  double interval_upper_eval = invert_position_map_func(interval_upper, ctx);
  struct gkyl_qr_res res = gkyl_ridders(invert_position_map_func, ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double z_computational = res.res;

  interval_lower = -1.0;
  interval_upper = 1.0;
  interval_lower_eval = invert_velocity_map_vpar_func(interval_lower, ctx);
  interval_upper_eval = invert_velocity_map_vpar_func(interval_upper, ctx);
  res = gkyl_ridders(invert_velocity_map_vpar_func, ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double vpar_computational = res.res;

  interval_lower = 0.0;
  interval_upper = 1.0;
  interval_lower_eval = invert_velocity_map_mu_func(interval_lower, ctx);
  interval_upper_eval = invert_velocity_map_mu_func(interval_upper, ctx);
  res = gkyl_ridders(invert_velocity_map_mu_func, ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double mu_computational = res.res;


  // Now we calculate the value of f at this computational coordinate
  // I'm limiting myself to 1x
  int pdim = app->cdim + app->vdim;
  double xc[3] = {z_computational, vpar_computational, mu_computational};
  int pidx[3];
  struct gkyl_basis basis = app->phase_basis;
  struct gkyl_rect_grid grid = app->f_ion_grid;
  struct gkyl_range local = app->phase_local;

  for (int i = 0; i < pdim; i++) {
    int idx_temp = local.lower[i] + (int) floor((xc[i] - grid.lower[i]) / grid.dx[i]);
    idx_temp = GKYL_MAX2(local.lower[i], GKYL_MIN2(local.upper[i], idx_temp));
    pidx[i] = idx_temp; 
  }
  long lidx = gkyl_range_idx(&local, pidx);

  const double *f_c = gkyl_array_cfetch(app->f_ion, lidx);
  const double *jacobvel_c = gkyl_array_cfetch(app->ion_jacobvel, lidx);

  int cidx[1] = {pidx[0]};
  long config_lidx = gkyl_range_idx(&app->position_map->local, cidx);
  const double *jacobtot_c = gkyl_array_cfetch(app->jacobtot, config_lidx);

  double cxc[3];
  gkyl_rect_grid_cell_center(&grid, pidx, cxc);
  for (int i = 0; i < pdim; i++) {
    xc[i] = (xc[i] - cxc[i]) / (grid.dx[i]*0.5);
  }
  double f_val = basis.eval_expand(xc, f_c);
  double config_xc[1] = {xc[0]};
  double jacobtot_val = app->position_map->basis.eval_expand(xc, jacobtot_c);
  fout[0] = f_val / jacobvel_c[0] / jacobtot_val;
}

void
free_ion_donor(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  gkyl_array_release(app->f_ion);
  gkyl_position_map_release(app->position_map);
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double src_amp = app->ion_source_amplitude;
  double z_src = 0.0;
  double src_sigma = app->ion_source_sigma;
  double src_amp_floor = src_amp*1e-2;
  if (fabs(z) <= 1.0)
  {
    fout[0] = fmax(src_amp_floor, (src_amp / sqrt(2.0 * M_PI * pow(src_sigma, 2))) *
      exp(-1 * pow((z - z_src), 2) / (2.0 * pow(src_sigma, 2))));
  }
  else
  {
    fout[0] = 1e-16;
  }
}

void
eval_upar_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double TSrc0 = app->ion_source_temp;
  double Tfloor = TSrc0*1e-2;
  if (fabs(z) <= 1.0)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

void
output_diagnostics(struct gk_mirror_ctx ctx, void *app_inp, struct gkyl_app_args app_args)
{
  struct gkyl_gk *app = app_inp;
  int my_rank = 0;
  int comm_sz = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi){
    gkyl_comm_get_rank(app->parallelism.comm, &my_rank);
    int comm_sz;
    gkyl_comm_get_size(app->parallelism.comm, &comm_sz);
  }
#endif
  if (my_rank == 0) {
    if (app->cdim == 1){printf("Grid size = %d in Z\n", app->cells[0]);}
    else if (app->cdim == 2){printf("Grid size = %d in psi, %d in Z\n", app->cells[0], app->cells[1]);}
    else if (app->cdim == 3){printf("Grid size = %d in psi, %d in Z, %d in theta\n", app->cells[0], app->cells[1], app->cells[2]);}

    printf("Velocity grid is %d in vpar and %d in mu \n", ctx.Nvpar, ctx.Nmu);
    if (app_args.use_mpi)
      printf("Number of MPI ranks: %d\n", app_args.cuts[0]);
    if (app_args.use_gpu)
      printf("Number of GPUs: %d\n", app_args.cuts[0]);
    printf("psi_eval = %g, psi_min = %g, psi_max = %g\n", ctx.psi_eval, ctx.psi_min, ctx.psi_max);
    printf("z_min = %g, z_max = %g\n", ctx.z_min, ctx.z_max);
    printf("vpar_max_ion/vti = %g, mu_max_ion/mu_ti = %g\n", ctx.vpar_max_ion/ctx.vti, sqrt(ctx.mu_max_ion/ctx.mi*2.0*ctx.B_p)/ctx.vti);
    printf("vpar_max_elc/vte = %g, mu_max_elc/mu_te = %g\n", ctx.vpar_max_elc/ctx.vte, sqrt(ctx.mu_max_elc/ctx.me*2.0*ctx.B_p)/ctx.vte);
    printf("vti = %.4e, vte = %.4e, c_s = %.4e, mu_ti = %.4e, mu_te = %.4e\n", ctx.vti, ctx.vte, ctx.c_s, ctx.mi * pow(ctx.vti, 2.) / (2. * ctx.B_p),
     ctx.me * pow(ctx.vte, 2.) / (2. * ctx.B_p));
    printf("omega_ci = %.4e, rho_s = %.4e, kperp = %.4e\n", ctx.omega_ci, ctx.rho_s, ctx.kperp);
    printf("nu_elc = %.4e, nu_ion = %.4e\n", ctx.nuElc, ctx.nuIon);
    printf("1/nuElc = %.4e, 1/nuIon = %.4e\n", 1./ctx.nuElc, 1./ctx.nuIon);
    printf("App name = %s\n", app->name);
    printf("Using positivity = %d\n", app->enforce_positivity);
    printf("Nonuniform mapping fraction = %g\n", app->geometry.position_map_info.map_strength);
    // Print the clock time
    time_t now;
    time(&now);
    printf("Date and time: %s", ctime(&now));
  }
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 1, vdim = 2; // Dimensionality.

  // Universal constant parameters.
  double eps0 = GKYL_EPSILON0;
  double mu0 = GKYL_MU0; // Not sure if this is right
  double eV = GKYL_ELEMENTARY_CHARGE;
  double mp = GKYL_PROTON_MASS; // ion mass
  double me = GKYL_ELECTRON_MASS;
  double qi = eV;  // ion charge
  double qe = -eV; // electron charge

  // Plasma parameters.
  double mi = 2.014 * mp;
  double Te0 = 940 * eV;
  double n0 = 3e19;
  double B_p = 0.53;
  double beta = 0.4;
  double tau = pow(B_p, 2.) * beta / (2.0 * mu0 * n0 * Te0) - 1.;
  double Ti0 = tau * Te0;
  double kperpRhos = 0.1;

  // Parameters controlling initial conditions.
  double alim = 0.125;
  double alphaIC0 = 2;
  double alphaIC1 = 10;

  double nuFrac = 1.0;
  double elc_nuFrac = 1/5.489216862238348;
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Te0 / eV);
  double nuElc = elc_nuFrac * nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
                 (6. * sqrt(2.) * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(me) * pow(Te0, 3. / 2.));
  // Ion-ion collision freq.
  double logLambdaIon = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Ti0 / eV);
  double nuIon = nuFrac * logLambdaIon * pow(eV, 4.) * n0 /
                 (12 * pow(M_PI, 3. / 2.) * pow(eps0, 2.) * sqrt(mi) * pow(Ti0, 3. / 2.));

  // Thermal speeds.
  double vti = sqrt(Ti0 / mi);
  double vte = sqrt(Te0 / me);
  double c_s = sqrt(Te0 / mi);

  // Gyrofrequencies and gyroradii.
  double omega_ci = eV * B_p / mi;
  double rho_s = c_s / omega_ci;

  // Perpendicular wavenumber in SI units:
  double kperp = kperpRhos / rho_s;

  // Geometry parameters.
  double z_min = -M_PI + 1e-2;
  double z_max = M_PI - 1e-2;
  double psi_min = 1e-6; // Go smaller. 1e-4 might be too small
  double psi_eval= 1e-3;
  double psi_max = 3e-3; // aim for 2e-2

  // Grid parameters
  double vpar_max_elc = 30 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 30 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nx = 16;
  int Nz = 216;
  int Nvpar = 32; // 96 uniform
  int Nmu = 32;  // 192 uniform
  int poly_order = 1;
  double t_end = 4e-3;//100e-6;
  int num_frames = 4e2;
  double write_phase_freq = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude = 7.39462473347548e22;
  double ion_source_sigma = 0.1;
  double ion_source_temp = 5000. * eV;

  struct gk_mirror_ctx ctx = {
    .cdim = cdim,
    .vdim = vdim,
    .mi = mi,
    .qi = qi,
    .me = me,
    .qe = qe,
    .Te0 = Te0,
    .n0 = n0,
    .B_p = B_p,
    .beta = beta,
    .tau = tau,
    .Ti0 = Ti0,
    .kperpRhos = kperpRhos,
    .alim = alim,
    .nuFrac = nuFrac,
    .logLambdaElc = logLambdaElc,
    .nuElc = nuElc,
    .elc_nuFrac = elc_nuFrac,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .kperp = kperp, 
    .z_min = z_min,
    .z_max = z_max,
    .psi_min = psi_min,
    .psi_eval = psi_eval,
    .psi_max = psi_max,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .ion_source_amplitude = ion_source_amplitude,
    .ion_source_sigma = ion_source_sigma,
    .ion_source_temp = ion_source_temp,
  };
  load_ion_donor(&ctx);
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_calc)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
  }

  bool trig_now_phase = gkyl_tm_trigger_check_and_bump(iot_phase, t_curr);
  if (trig_now_phase || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;

    gkyl_gyrokinetic_app_write_phase(app, t_curr, frame);
  }
}

int main(int argc, char **argv)
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Init(&argc, &argv);
  }
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct gk_mirror_ctx ctx = create_ctx(); // context for init functions
  int cells_x[ctx.cdim], cells_v[ctx.vdim];
  for (int d=0; d<ctx.cdim; d++)
    cells_x[d] = APP_ARGS_CHOOSE(app_args.xcells[d], ctx.cells[d]);
  for (int d=0; d<ctx.vdim; d++)
    cells_v[d] = APP_ARGS_CHOOSE(app_args.vcells[d], ctx.cells[ctx.cdim+d]);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, stderr);
  
  int my_rank = 0;
  int comm_sz = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi){
    gkyl_comm_get_rank(comm, &my_rank);
    int comm_sz;
    gkyl_comm_get_size(comm, &comm_sz);
  }
#endif

  struct gkyl_gyrokinetic_projection ion_ic = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = read_ion_distf,
      .ctx_func = &ctx,  
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .no_by = true,
    .projection = ion_ic,
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .evolve = true,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
      .diagnostics = { 
        .num_diag_moments = 6,
        .diag_moments = { "BiMaxwellianMoments", "M0", "M1", "M2", "M2par", "M2perp" },
      },
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Ti0,
      .ctx = &ctx,
      .self_nu = evalNuIon, 
      .write_diagnostics = true, 
      .nuFrac = 2000.0,
    },
    .num_diag_moments = 8,
    .diag_moments = { "BiMaxwellianMoments", "M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp" },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { "FourMoments" },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { "FourMoments" },
    },
  };
  
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
  };

struct gkyl_efit_inp efit_inp = {
    .filepath = "/home/mr1884/scratch/gkylmax/eqdsk/wham.geqdsk",
    .rz_poly_order = 2,                     // polynomial order for psi(R,Z) used for field line tracing
    .flux_poly_order = 1,                   // polynomial order for fpol(psi)
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary 
  };
  
  struct gkyl_gk app_inp = {  // GK app
    .name = "gk_wham",
    .cdim = ctx.cdim ,  .vdim = ctx.vdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .enforce_positivity = true,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {ctx.psi_eval, 0.0},
      .efit_info = efit_inp,
      .mirror_grid_info = grid_inp,
      .position_map_info = {
        .id = GKYL_PMAP_CONSTANT_DB_NUMERIC,
        .map_strength = 0.25,
        .maximum_slope_at_max_B = 1.0,
        .maximum_slope_at_min_B = 4.0,
      },
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 1,
    .species = {ion},
    .field = field,
    .parallelism = {
      .use_gpu = app_args.use_gpu,
      .cuts = { app_args.cuts[0] },
      .comm = comm,
    },
  };
  
  // Create app object.

  output_diagnostics(ctx, &app_inp, app_args);

  clock_t start_time = clock();
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);
  clock_t end_time = clock();
  if (my_rank == 0)
    printf("Time to create app object: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

  // Initial and final simulation times.
  start_time = clock();
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.t_end;
  // Initialize simulation.
  if (app_args.is_restart) {
    struct gkyl_app_restart_status status = gkyl_gyrokinetic_app_read_from_frame(app, app_args.restart_frame);

    if (status.io_status != GKYL_ARRAY_RIO_SUCCESS) {
      gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read restart file! (%s)\n",
        gkyl_array_rio_status_msg(status.io_status));
      goto freeresources;
    }

    frame_curr = status.frame;
    t_curr = status.stime;

    gkyl_gyrokinetic_app_cout(app, stdout, "Restarting from frame %d", frame_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " at time = %g\n", t_curr);
  }
  else {
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
  }  
  end_time = clock();
  if (my_rank == 0)
    printf("Time to load initial conditions: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
  start_time = clock();

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = t_end/(ctx.write_phase_freq*num_frames), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1;
  end_time = clock();
  if (my_rank == 0)
    printf("Time to write diagnostics: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
  start_time = clock();
  double init_time = t_curr;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);    
    if (step % 1000 == 0 || step == 1) {
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g ... ", status.dt_actual);
      end_time = clock();
      double time_per_hour = (t_curr + status.dt_actual - init_time) / ((double)(end_time - start_time) / CLOCKS_PER_SEC) * 3600.0;
      gkyl_gyrokinetic_app_cout(app, stdout, "will cover %g s in 1 hour\n", time_per_hour);
    }

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
          break;
      }

    t_curr += status.dt_actual;
      dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, t_curr > t_end);

    if (dt_init < 0.0) {
      dt_init = status.dt_actual;
    }
    else if (status.dt_actual < dt_failure_tol * dt_init) {
      num_failures += 1;

      gkyl_gyrokinetic_app_cout(app, stdout, "WARNING: Time-step dt = %g", status.dt_actual);
      gkyl_gyrokinetic_app_cout(app, stdout, " is below %g*dt_init ...", dt_failure_tol);
      gkyl_gyrokinetic_app_cout(app, stdout, " num_failures = %d\n", num_failures);
      if (num_failures >= num_failures_max) {
        gkyl_gyrokinetic_app_cout(app, stdout, "ERROR: Time-step was below %g*dt_init ", dt_failure_tol);
        gkyl_gyrokinetic_app_cout(app, stdout, "%d consecutive times. Aborting simulation ....\n", num_failures_max);
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, true);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

      step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);
  
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);
  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0)
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  free_ion_donor(&ctx);
  gkyl_gyrokinetic_comms_release(comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
