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
  double elc_source_amplitude;
  double elc_source_sigma;
  double elc_source_temp;

  // Boltzmann electron reading
  struct gkyl_array *field;
  struct gkyl_array *ion_M0;
  double ni_sheath;
  struct gkyl_position_map *position_map;
  double target_z_fa;
};

void
scaleIC(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0/3.0; // Scale factor for initial conditions
}

// Evaluate collision frequencies
void
evalNuElc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->nuElc;
}

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
  struct gkyl_rect_grid field_grid, mc2nu_pos_grid, M0_grid;
  struct gkyl_array *field, *mc2nu_pos, *M0;

  field     = gkyl_grid_array_new_from_file(&field_grid, 
    "/home/mr1884/scratch/gkylmax/initial-conditions/boltz-elc-288z-nu2000/gk_wham-field_0.gkyl");
  mc2nu_pos = gkyl_grid_array_new_from_file(&mc2nu_pos_grid, 
    "/home/mr1884/scratch/gkylmax/initial-conditions/boltz-elc-288z-nu2000/gk_wham-mc2nu_pos.gkyl");
  M0 = gkyl_grid_array_new_from_file(&M0_grid, 
    "/home/mr1884/scratch/gkylmax/initial-conditions/boltz-elc-288z-nu2000/gk_wham-ion_M0_0.gkyl");

  app->field = field;
  app->ion_M0 = M0;

  int lower_cell[] = {1};
  int upper_cell[] = {M0_grid.cells[0]};

  int poly_order = 1;
  int cdim = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  int nghost[] = { 0, 0, 0 };
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&M0_grid, nghost, &local, &local_ext);

  // Create a position map object
  struct gkyl_position_map_inp pmap_inp = { };
  // Potential future issue by using the M0 ranges and basis for the position map
  struct gkyl_position_map *gpm = gkyl_position_map_new(pmap_inp, 
    mc2nu_pos_grid, local, local_ext, local, local_ext, basis);
  gkyl_position_map_set_mc2nu(gpm, mc2nu_pos);
  app->position_map = gpm;
}

void
free_ion_donor(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  gkyl_array_release(app->field);
  gkyl_array_release(app->ion_M0);
  gkyl_position_map_release(app->position_map);
}

static double
invert_position_map_func(double x, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double xc[] = {x};
  double xfa[3];
  gkyl_position_map_eval_mc2nu(app->position_map, xc, xfa);
  // printf("xc = %g maps to xfa = %g \n", x, xfa[0]);
  return xfa[0] - app->target_z_fa;
}

void
botlzmann_elc_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Since boltzmann electrons are particle density and gyrokinetic sims take guiding center density
  // Boltzmann electron simulations have n_i = n_e, and we must determine the ion density
  // using the polarization density.

  struct gk_mirror_ctx *app = ctx;
  double z_field_aligned = xn[0];

  // First we must determine the computational coordinate in non-uniform space
  app->target_z_fa = z_field_aligned;
  double interval_lower = -M_PI;
  double interval_upper = M_PI;
  double interval_lower_eval = invert_position_map_func(interval_lower, ctx);
  double interval_upper_eval = invert_position_map_func(interval_upper, ctx);
  struct gkyl_qr_res res = gkyl_ridders(invert_position_map_func, ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double z_computational = res.res;

  // Now we calculate the value of the field at this computational coordinate
  int cdim = 1;
  struct gkyl_basis basis = app->position_map->basis;
  struct gkyl_rect_grid grid = app->position_map->grid;
  struct gkyl_range local = app->position_map->local;

  // I'm limiting myself to 1x
  int idx_temp = local.lower[0] + (int) floor((z_computational - grid.lower[0]) / grid.dx[0]);
  idx_temp = GKYL_MAX2(local.lower[0], GKYL_MIN2(local.upper[0], idx_temp));
  long lidx = gkyl_range_idx(&local, &idx_temp);
  const double *field_coeffs = gkyl_array_cfetch(app->ion_M0, lidx);
  double cxc[3];
  gkyl_rect_grid_cell_center(&grid, &idx_temp, cxc);
  double x_log = (z_computational - cxc[0]) / (grid.dx[0]*0.5);
  double M0_val = basis.eval_expand(&x_log, field_coeffs);

  fout[0] = M0_val;
}

void
boltzmann_elc_upar(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
boltzmann_elc_T(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  fout[0] = app->Te0;
}


void
botlzmann_elc_field(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  // There is no position map in the polarization solve presently. It needs that
  double z_field_aligned = xn[0];

  // First we must determine the computational coordinate in non-uniform space
  app->target_z_fa = z_field_aligned;
  double interval_lower = -M_PI;
  double interval_upper = M_PI;
  double interval_lower_eval = invert_position_map_func(interval_lower, ctx);
  double interval_upper_eval = invert_position_map_func(interval_upper, ctx);
  struct gkyl_qr_res res = gkyl_ridders(invert_position_map_func, ctx,
    interval_lower, interval_upper, interval_lower_eval, interval_upper_eval, 10, 1e-6);
  double z_computational = res.res;

  // Now we calculate the value of the field at this computational coordinate
  int cdim = 1;
  struct gkyl_basis basis = app->position_map->basis;
  struct gkyl_rect_grid grid = app->position_map->grid;
  struct gkyl_range local = app->position_map->local;

  // I'm limiting myself to 1x
  int idx_temp = local.lower[0] + (int) floor((z_computational - grid.lower[0]) / grid.dx[0]);
  idx_temp = GKYL_MAX2(local.lower[0], GKYL_MIN2(local.upper[0], idx_temp));
  long lidx = gkyl_range_idx(&local, &idx_temp);
  const double *field_coeffs = gkyl_array_cfetch(app->field, lidx);
  double cxc[3];
  gkyl_rect_grid_cell_center(&grid, &idx_temp, cxc);
  double x_log = (z_computational - cxc[0]) / (grid.dx[0]*0.5);
  double field_val = basis.eval_expand(&x_log, field_coeffs);

  if (field_val < 10.0)
   field_val = 10.0;
  fout[0] = field_val;
}

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

void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  double b = 1.45;
  double linear_velocity_threshold = 1./6.;
  double frac_linear = 1/b*atan(linear_velocity_threshold*tan(b));
  if (fabs(cvpar) < frac_linear) {
    double func_frac = tan(frac_linear*b) / tan(b);
    vp[0] = vpar_max_elc*func_frac*cvpar/frac_linear;
  }
  else {
    vp[0] = vpar_max_elc*tan(cvpar*b)/tan(b);
  }
  // vp[0] = vc[0] * vpar_max_elc;
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
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
  double z_min = -2.0;
  double z_max =  2.0;
  double psi_min = 1e-6; // Go smaller. 1e-4 might be too small
  double psi_eval= 1e-3;
  double psi_max = 3e-3; // aim for 2e-2

  // Grid parameters
  double vpar_max_elc = 30 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 30 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nx = 16;
  int Nz = 288;
  int Nvpar = 32; // 96 uniform
  int Nmu = 32;  // 192 uniform
  int poly_order = 1;
  double t_end = 1e-9;//100e-6;
  int num_frames = 10;
  double write_phase_freq = 1;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude = 7.39462473347548e22;
  double ion_source_sigma = 0.1;
  double ion_source_temp = 5000. * eV;
  double elc_source_amplitude = 7.39462473347548e22;
  double elc_source_sigma = 0.1;
  double elc_source_temp = Te0;

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
    .elc_source_amplitude = elc_source_amplitude,
    .elc_source_sigma = elc_source_sigma,
    .elc_source_temp = elc_source_temp,
  };
  load_ion_donor(&ctx);
  return ctx;
}

void
calc_integrated_diagnostics(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app,
  double t_curr, bool is_restart_IC, bool force_calc, double dt)
{
  if (!is_restart_IC && (gkyl_tm_trigger_check_and_bump(iot, t_curr) || force_calc)) {
    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);

    if ( !(dt < 0.0) )
      gkyl_gyrokinetic_app_save_dt(app, t_curr, dt);
  }
}

void
write_data(struct gkyl_tm_trigger* iot_conf, struct gkyl_tm_trigger* iot_phase,
  gkyl_gyrokinetic_app* app, double t_curr, bool is_restart_IC, bool force_write)
{
  bool trig_now_conf = gkyl_tm_trigger_check_and_bump(iot_conf, t_curr);
  if (trig_now_conf || force_write) {
    int frame = (!trig_now_conf) && force_write? iot_conf->curr : iot_conf->curr-1;
    gkyl_gyrokinetic_app_write_conf(app, t_curr, frame);

    if (!is_restart_IC) {
      gkyl_gyrokinetic_app_write_field_energy(app);
      gkyl_gyrokinetic_app_write_integrated_mom(app);
      gkyl_gyrokinetic_app_write_dt(app);
    }
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

  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .density = botlzmann_elc_density,
    .ctx_density = &ctx,
    .temp = boltzmann_elc_T,
    .ctx_temp = &ctx,
    .upar = boltzmann_elc_upar,
    .ctx_upar = &ctx,
  };

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .no_by = true,
    .projection = elc_ic,
    // .init_from_file = {
    //   .type = GKYL_IC_IMPORT_F,
    //   .file_name = "gk_wham-elc_0.gkyl",
    // },
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Te0,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = {"ion"},
      .write_diagnostics = true,
      .nuFrac = ctx.elc_nuFrac * 2000.0,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .num_adapt_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
        .center_gauss = {0.0},
        .sigma_gauss = {0.2},
        .particle = 0.0,
        .energy = 0.0,
        .temp_max = 10000*GKYL_ELEMENTARY_CHARGE,
      },
      .adapt[0] = {
        .adapt_species_name = "ion",
        .adapt_particle = true,
        .adapt_energy = true,
        .num_boundaries = 2,
        .dir = {0, 0},
        .edge = {GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
      },
      .diagnostics = {
        .num_diag_moments = 6,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      },
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
  };

  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .scale_with_polarization = true,
    .no_by = true,
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .init_from_file = {
      .type = GKYL_IC_IMPORT_F,
      // .conf_scale = &scaleIC,
      // .conf_scale_ctx = &ctx,
      .file_name = "../initial-conditions/kinet-elc-288z-nu2000/gk_wham-ion_0.gkyl",
    },
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Ti0,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = {"elc"},
      .write_diagnostics = true,
      .nuFrac = 2000.0,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .num_adapt_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_GAUSSIAN,
        .center_gauss = {0.0},
        .sigma_gauss = {0.2},
        .particle = 0.0,
        .energy = 0.0,
        .temp_max = 10000*GKYL_ELEMENTARY_CHARGE,
      },
      .adapt[0] = {
        .adapt_species_name = "ion",
        .adapt_particle = true,
        .adapt_energy = true,
        .num_boundaries = 2,
        .dir = {0, 0},
        .edge = {GKYL_LOWER_EDGE, GKYL_UPPER_EDGE},
      },
      .diagnostics = {
        .num_diag_moments = 5,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP },
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
      }
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_HAMILTONIAN },
    },
  };

  struct gkyl_gyrokinetic_field field = {
    .polarization_bmag = ctx.B_p,
    .kperpSq = pow(ctx.kperp, 2.),
    .polarization_potential = botlzmann_elc_field,
    .polarization_potential_ctx = &ctx,
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "/home/mr1884/gkylzero/data/unit/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_MIRROR_GRID_GEN_SQRT_PSI_CART_Z, // coordinate system for psi grid
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
      .mirror_grid_info = grid_inp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 2,
    .species = {ion, elc},
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
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, app_args.is_restart, false, -1.0);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, app_args.is_restart, false);

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

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, t_curr > t_end, status.dt_actual);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, t_curr > t_end);

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
        calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false, true, status.dt_actual);
        write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false, true);
        break;
      }
    }
    else {
      num_failures = 0;
    }

      step += 1;
  }

  gkyl_gyrokinetic_app_stat_write(app);

  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app); // fetch simulation statistics
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
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld\n", stat.n_io);
  gkyl_gyrokinetic_app_print_timings(app, stdout);

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
