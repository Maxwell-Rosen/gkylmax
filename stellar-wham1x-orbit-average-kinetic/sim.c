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
  double ion_source_temp;
  double elc_source_amplitude;
  double elc_source_temp;

  // Boltzmann electron reading
  struct gkyl_array *field;
  struct gkyl_array *ion_M0;
  double ni_sheath;
  struct gkyl_rect_grid grid;
  struct gkyl_range local; // Local range of the grid.
  struct gkyl_basis basis; // Basis for the grid.
  double target_z_fa;

  double alpha_ion; // Multirate factor.
  double alpha_elc; // Multirate factor.
  bool static_field; // Whether the potential is constant in time.
  double I_loss_ion; // Value of the mask in the loss region.
  double I_loss_elc; // Value of the mask in the loss region.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_ion; // Type of df/dt multiplier.
  enum gkyl_gyrokinetic_fdot_multiplier_type fdot_mult_type_elc; // Type of df/dt multiplier.
  bool use_positivity_hack; // Whether to use the positivity hack.
  bool static_ion; // Whether the ion distribution is static.
  bool static_elc; // Whether the electron distribution is static.
};

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
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double src_amp = app->ion_source_amplitude;
  if (fabs(z) <= 1.0)
    fout[0] = src_amp * (1 - pow(fabs(z), 6));
  else
    fout[0] = 1e-16;
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
    fout[0] = TSrc0;
  else
    fout[0] = Tfloor;
}


void
eval_density_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double src_amp = app->elc_source_amplitude;
  if (fabs(z) <= 1.0)
    fout[0] = src_amp * (1 - pow(fabs(z), 6));
  else
    fout[0] = 1e-16;
}

void
eval_upar_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  fout[0] = 0.0;
}

void
eval_temp_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double z = xn[0];
  double TSrc0 = app->ion_source_temp;
  double Tfloor = TSrc0*1e-2;
  if (fabs(z) <= 1.0)
    fout[0] = TSrc0;
  else
    fout[0] = Tfloor;
}

void
load_ion_donor(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  struct gkyl_rect_grid field_grid, mc2nu_pos_grid, M0_grid;
  struct gkyl_array *field, *mc2nu_pos, *M0;

  field = gkyl_grid_array_new_from_file(&field_grid,
    "initial-condition/gk_wham-field_500.gkyl");
  M0 = gkyl_grid_array_new_from_file(&M0_grid,
    "initial-condition/gk_wham-ion_M0_500.gkyl");

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

  app->grid = M0_grid;
  app->local = local;
  app->basis = basis;
}

void
free_ion_donor(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  gkyl_array_release(app->field);
  gkyl_array_release(app->ion_M0);
}

void
botlzmann_elc_density(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  // Since boltzmann electrons are particle density and gyrokinetic sims take guiding center density
  // Boltzmann electron simulations have n_i = n_e, and we must determine the ion density
  // using the polarization density.

  struct gk_mirror_ctx *app = ctx;
  double z_computational = xn[0];

  // Now we calculate the value of the field at this computational coordinate
  int cdim = 1;
  struct gkyl_basis basis = app->basis;
  struct gkyl_rect_grid grid = app->grid;
  struct gkyl_range local = app->local;

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
  double z_computational = xn[0];

  // Now we calculate the value of the field at this computational coordinate
  int cdim = 1;
  struct gkyl_basis basis = app->basis;
  struct gkyl_rect_grid grid = app->grid;
  struct gkyl_range local = app->local;

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
  double elc_nuFrac = 1/5.489216862238348; // Factor from Pastukhov calculation with Dougherty collisions
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
  // double t_end = 1.5e-3;//100e-6;
  // int num_frames = 1500;
  double write_phase_freq = 1;
  // int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

  // Source parameters
  double ion_source_amplitude = 1.e20;
  double ion_source_temp = 5000. * eV;
  double elc_source_amplitude = 1.e20;
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
    // .t_end = t_end,
    // .num_frames = num_frames,
    .write_phase_freq = write_phase_freq,
    // .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
    .ion_source_amplitude = ion_source_amplitude,
    .ion_source_temp = ion_source_temp,
    .elc_source_amplitude = elc_source_amplitude,
    .elc_source_temp = elc_source_temp,
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
    gkyl_gyrokinetic_app_write_field_energy(app);
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
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
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

  // Extract variables from command line arguments.
// "alpha_ion=$,alpha_elc=$,t_end=$,num_frames=$,static_field=$,I_loss_ion=$,I_loss_elc=$,fdot_mult_type_ion=$,fdot_mult_type_elc=$,use_positivity_hack=$,static_ion=$,static_elc=$"

  // sscanf(app_args.opt_args, "alpha_ion=%lf,alpha_elc=%lf,t_end=%lf,num_frames=%d,static_field=%d,I_loss_ion=%lf,I_loss_elc=%lf,fdot_mult_type_ion=%d,fdot_mult_type_elc=%d,use_positivity_hack=%d,static_ion=%d,static_elc=%d",
    // &ctx.alpha_ion, &ctx.alpha_elc, &ctx.t_end, &ctx.num_frames, &ctx.static_field, &ctx.I_loss_ion, &ctx.I_loss_elc, &ctx.fdot_mult_type_ion, &ctx.fdot_mult_type_elc, &ctx.use_positivity_hack, &ctx.static_ion, &ctx.static_elc);

// alpha_ion=1.0,alpha_elc=1.0,t_end=.001,num_frames=5,static_field=0,I_loss_ion=1.0,I_loss_elc=1.0,fdot_mult_type_ion=0,fdot_mult_type_elc=0,use_positivity_hack=1,static_ion=1,static_elc=0
  ctx.alpha_ion = 1.0;
  ctx.alpha_elc = 1.0;
  ctx.t_end = 0.001; // in seconds
  ctx.num_frames = 5;
  ctx.static_field = false;
  ctx.I_loss_ion = 1.0; // Value of the mask in the loss region.
  ctx.I_loss_elc = 1.0; // Value of the mask in the loss region.
  ctx.fdot_mult_type_ion = 0; // Type of df/dt multiplier.
  ctx.fdot_mult_type_elc = 0; // Type of df/dt multiplier.
  ctx.use_positivity_hack = false; // Whether to use the positivity hack.
  ctx.static_ion = true; // Whether the ion distribution is static.
  ctx.static_elc = false; // Whether the electron distribution is static. 

  // Convert t_end to seconds.
  ctx.t_end = ctx.t_end*1e-6;
  if (my_rank == 0) {
    printf("Using command line arguments:\n");
    printf("  alpha = %.9e\n", ctx.alpha_ion);
    printf("  alpha_elc = %.9e\n", ctx.alpha_elc);
    printf("  t_end = %.9e s\n", ctx.t_end);
    printf("  num_frames = %d\n", ctx.num_frames);
    printf("  static_field = %d\n", ctx.static_field);
    printf("  I_loss_ion = %.9e\n", ctx.I_loss_ion);
    printf("  I_loss_elc = %.9e\n", ctx.I_loss_elc);
    printf("  fdot_mult_type_ion = %d\n", ctx.fdot_mult_type_ion);
    printf("  fdot_mult_type_elc = %d\n", ctx.fdot_mult_type_elc);
    printf("  use_positivity_hack = %d\n", ctx.use_positivity_hack);
    printf("  static_ion = %d\n", ctx.static_ion);
    printf("  static_elc = %d\n", ctx.static_elc);
  }
  ctx.int_diag_calc_num = ctx.num_frames*100;
  struct gkyl_gyrokinetic_projection elc_ic = {
    .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
    .density = botlzmann_elc_density,
    .ctx_density = &ctx,
    .upar = boltzmann_elc_upar,
    .ctx_upar = &ctx,
    .temp = boltzmann_elc_T,
    .ctx_temp = &ctx,
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
    .is_static = ctx.static_elc,
    .projection = elc_ic,
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },

    // .collisionless_scale_factor = ctx.alpha_elc,

    // .time_rate_multiplier = {
    //   .type = ctx.fdot_mult_type_elc,
    // },

    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .normNu = true,
      .n_ref = ctx.n0,
      .T_ref = ctx.Te0,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = {"ion"}, // There is an issue with the cross-species collisions. If we comment this out, the simulation runs fine.
      .write_diagnostics = true,
      .nuFrac = ctx.elc_nuFrac,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM,
        .density = eval_density_elc_source,
        .ctx_density = &ctx,
        .upar = eval_upar_elc_source,
        .ctx_upar = &ctx,
        .temp = eval_temp_elc_source,
        .ctx_temp = &ctx,
      },
      .diagnostics = {
        .num_diag_moments = 6,
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
      },
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
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
    .is_static = ctx.static_ion,
    .init_from_file = {
      .type = GKYL_IC_IMPORT_F,
      .file_name = "initial-condition/gk_wham-ion_500.gkyl",
    },
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },

    .collisionless_scale_factor = ctx.alpha_ion,

    .time_rate_multiplier = {
      .type = ctx.fdot_mult_type_ion,
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
      .num_cross_collisions = 1,
      .collide_with = {"elc"},
      .write_diagnostics = true,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .num_sources = 1,
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
        .diag_moments = { GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_HAMILTONIAN},
        .num_integrated_diag_moments = 1,
        .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
      },
    },
    .write_omega_cfl = true,
    .num_diag_moments = 8,
    .diag_moments = {GKYL_F_MOMENT_BIMAXWELLIAN, GKYL_F_MOMENT_M0, GKYL_F_MOMENT_M1, GKYL_F_MOMENT_M2, GKYL_F_MOMENT_M2PAR, GKYL_F_MOMENT_M2PERP, GKYL_F_MOMENT_M3PAR, GKYL_F_MOMENT_M3PERP },
    .num_integrated_diag_moments = 1,
    .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP },
    .time_rate_diagnostics = true,

    .boundary_flux_diagnostics = {
      .num_integrated_diag_moments = 1,
      .integrated_diag_moments = { GKYL_F_MOMENT_M0M1M2PARM2PERP},
    },
  };

  struct gkyl_gyrokinetic_field field = {
    .polarization_bmag = ctx.B_p,
    .kperpSq = pow(ctx.kperp, 2.),
    .polarization_potential = botlzmann_elc_field,
    .polarization_potential_ctx = &ctx,
    .is_static = ctx.static_field,
  };

  struct gkyl_mirror_geo_grid_inp grid_inp = {
    .filename_psi = "../eqdsk/wham_hires.geqdsk_psi.gkyl", // psi file to use
    .rclose = 0.2, // closest R to region of interest
    .zmin = -2.0,  // Z of lower boundary
    .zmax =  2.0,  // Z of upper boundary
    .include_axis = false, // Include R=0 axis in grid
    .fl_coord = GKYL_GEOMETRY_MIRROR_GRID_GEN_PSI_CART_Z, // coordinate system for psi grid
  };

  struct gkyl_gk app_inp = {  // GK app
    .name = "gk_wham",
    .cdim = ctx.cdim ,  .vdim = ctx.vdim,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = { cells_x[0] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .enforce_positivity = ctx.use_positivity_hack,
    .geometry = {
      .geometry_id = GKYL_GEOMETRY_MIRROR,
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
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  // Initial and final simulation times.
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

  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write_conf = { .dt = (t_end-t_curr)/(num_frames-frame_curr), .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_write_phase = { .dt = (t_end-t_curr)/(ctx.write_phase_freq*(num_frames-frame_curr)), .tcurr = t_curr, .curr = frame_curr};
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = (t_end-t_curr)/GKYL_MAX2(num_frames-frame_curr, (num_int_diag_calc/num_frames)*(num_frames-frame_curr)),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write_conf, &trig_write_phase, app, t_curr, false);

  // start, end and initial time-step
  double dt = t_end-t_curr;
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;

  long step = 1, num_steps = app_args.num_steps;
  while ((t_curr < t_end) && (step <= num_steps)) {
    // if (step == 1 || step % 100 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %i at t = %g ...", step, t_curr);

    dt = t_end - t_curr;
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);

    // if (step == 1 || step % 100 == 0)
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
      break;
    }
    t_curr += status.dt_actual;
    dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr >= t_end);
    write_data(&trig_write_conf, &trig_write_phase, app, t_curr, t_curr >= t_end);

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

  // fetch simulation statistics
  struct gkyl_gyrokinetic_stat stat = gkyl_gyrokinetic_app_stat(app);

  gkyl_gyrokinetic_app_cout(app, stdout, "\n");
  gkyl_gyrokinetic_app_cout(app, stdout, "Simulation completed at t = %g s\n", t_curr);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of update calls %ld\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of forward-Euler calls %ld\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-2 failures %ld\n", stat.nstage_2_fail);
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.n_io);
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
