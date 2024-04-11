#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fem_parproj.h>
#include <gkyl_fem_poisson_bctype.h>
#include <gkyl_gyrokinetic.h>
#include <gkyl_math.h>
#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

#include <rt_arg_parse.h>
#include <gkyl_mirror_geo.h>

// Define the context of the simulation. This is basically all the globals
struct gk_mirror_ctx
{
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
  double alphaIC0;
  double alphaIC1;
  double nuFrac;
  // Electron-electron collision freq.
  double logLambdaElc;
  double nuElc;
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
  double psi_eval;
  // Magnetic equilibrium model.
  double mcB;
  double gamma;
  double Z_m;
  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double B_bt;
  double R_bt;
  double Z_bt;
  double z_bt;
  double R_m;
  double B_m;
  double z_m;
double Z_m_computational;
  // Physics parameters at mirror throat
  double n_m;
  double Te_m;
  double Ti_m;
double Ti_perp0;
  double Ti_par0;
  double Ti_perp_m;
  double Ti_par_m;
  double Te_perp0;
  double Te_par0;
  double Te_perp_m;
  double Te_par_m;
  double cs_m;
  // Source parameters
  double NSrcIon;
  double lineLengthSrcIon;
  double sigSrcIon;
  double NSrcFloorIon;
  double TSrc0Ion;
  double TSrcFloorIon;
  double NSrcElc;
  double lineLengthSrcElc;
  double sigSrcElc;
  double NSrcFloorElc;
  double TSrc0Elc;
  double TSrcFloorElc;
  // Grid parameters
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int num_cell_vpar;
  int num_cell_mu;
  int num_cell_z;
  int poly_order;
  double final_time;
  int num_frames;
  double psi_in;
  double z_in;
  // For non-uniform mapping
  double diff_dz;
  double psi_in_diff;
  int mapping_order_center;
  int mapping_order_expander;
  double mapping_frac;
  void *mapc2p;
};


struct gkyl_mirror_geo_efit_inp inp = {
  // psiRZ and related inputs
  .filepath = "../eqdsk/wham.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.48,
  .zmax =  2.48,
  .write_node_coord_array = true,
  .node_file_nm = "wham_nodes.gkyl",
  .nonuniform_mapping_fraction = 0.0,
};

// -- Source functions.
void
eval_density_elc_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double NSrc = app->NSrcElc;
  double zSrc = app->lineLengthSrcElc;
  double sigSrc = app->sigSrcElc;
  double NSrcFloor = app->NSrcFloorElc;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2.))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2.))));
  }
  else
  {
    fout[0] = 1e-16;
  }
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
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double sigSrc = app->sigSrcElc;
  double TSrc0 = app->TSrc0Elc;
  double Tfloor = app->TSrcFloorElc;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

void
eval_density_ion_source(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double NSrc = app->NSrcIon;
  double zSrc = app->lineLengthSrcIon;
  double sigSrc = app->sigSrcIon;
  double NSrcFloor = app->NSrcFloorIon;
  if (fabs(z) <= app->z_m)
  {
    fout[0] = fmax(NSrcFloor, (NSrc / sqrt(2.0 * M_PI * pow(sigSrc, 2))) *
                                  exp(-1 * pow((z - zSrc), 2) / (2.0 * pow(sigSrc, 2))));
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
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double sigSrc = app->sigSrcIon;
  double TSrc0 = app->TSrc0Ion;
  double Tfloor = app->TSrcFloorIon;
  if (fabs(z) <= 2.0 * sigSrc)
  {
    fout[0] = TSrc0;
  }
  else
  {
    fout[0] = Tfloor;
  }
}

// Electrons initial conditions
void
eval_density_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}

void
eval_temp_par_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double Te_par0 = app->Te_par0;
  double Te_par_m = app->Te_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Te_par_m+(Te_par0-Te_par_m)*tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Te_par_m;
  }
}

void
eval_temp_perp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double Te_perp0 = app->Te_perp0;
  double Te_perp_m = app->Te_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Te_perp_m - Te_perp0*tanh(3.*fabs(z_m-fabs(z)));
  }
  else
  {
    fout[0] = Te_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

void
eval_temp_elc(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  eval_temp_par_elc(t, xn, fout, ctx);
  double Tpar = fout[0];
  eval_temp_perp_elc(t, xn, fout, ctx);
  double Tperp = fout[0];
  fout[0] = (Tpar + 2 * Tperp) / 3;
}


// Ion initial conditions
void
eval_density_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double sigma = 0.9*z_m;
  if (fabs(z) <= sigma)
  {
    fout[0] = 0.5*app->n0*(1. + tanh(10. * sigma * fabs(sigma - fabs(z))));
  }
  else
  {
    fout[0] = 0.5*app->n0*exp(-5 * (fabs(sigma - fabs(z))));
  }
}

void
eval_upar_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double cs_m = app->cs_m;
  double z_m = app->z_m;
  double z_max = app->z_max;
  if (fabs(z) <= z_m)
  {
    fout[0] = 0.0;
  }
  else
  {
    fout[0] = fabs(z) / z * cs_m * tanh(3 * (z_max - z_m) * fabs(fabs(z) - z_m)); // Maybe put a 5 here
  }
}

void
eval_temp_par_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double Ti_par0 = app->Ti_par0;
  double Ti_par_m = app->Ti_par_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_par_m + (Ti_par0 - Ti_par_m) * tanh(4 * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_par_m * GKYL_MAX2(1.e-2, 4 * log(fabs(fabs(z) - z_m) + 1));
  }
}

void
eval_temp_perp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double psi = app->psi_eval; // Magnetic flux function psi of field line.
  double z = xn[0];
  double z_m = app->z_m;
  double Ti_perp0 = app->Ti_perp0;
  double Ti_perp_m = app->Ti_perp_m;
  if (fabs(z) <= z_m)
  {
    fout[0] = Ti_perp_m + (Ti_perp0 - Ti_perp_m) * tanh(3. * fabs(z_m - fabs(z)));
  }
  else
  {
    fout[0] = Ti_perp_m * GKYL_MAX2(1.e-3, exp(-5. * (fabs(z_m - fabs(z)))));
  }
}

void
eval_temp_ion(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  eval_temp_par_ion(t, xn, fout, ctx);
  double Tpar = fout[0];
  eval_temp_perp_ion(t, xn, fout, ctx);
  double Tperp = fout[0];
  fout[0] = (Tpar + 2 * Tperp) / 3;
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

struct gk_mirror_ctx
create_ctx(void)
{
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
  // Electron-electron collision freq.
  double logLambdaElc = 6.6 - 0.5 * log(n0 / 1e20) + 1.5 * log(Te0 / eV);
  double nuElc = nuFrac * logLambdaElc * pow(eV, 4.) * n0 /
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
  double RatZeq0 = 0.10; // Radius of the field line at Z=0.
  double z_min = -M_PI + 1e-1;
  double z_max = M_PI - 1e-1;
  double psi_eval = 0.0026530898059565;

  // Parameters controlling the magnetic equilibrium model.
  double mcB = 6.51292;
  double gamma = 0.124904;
  double Z_m = 0.98;

  // Source parameters
  double NSrcIon = 3.1715e23 / 8.0;
  double lineLengthSrcIon = 0.0;
  double sigSrcIon = Z_m / 4.0;
  double NSrcFloorIon = 0.05 * NSrcIon;
  double TSrc0Ion = Ti0 * 1.25;
  double TSrcFloorIon = TSrc0Ion / 8.0;
  double NSrcElc = NSrcIon;
  double lineLengthSrcElc = lineLengthSrcIon;
  double sigSrcElc = sigSrcIon;
  double NSrcFloorElc = NSrcFloorIon;
  double TSrc0Elc = TSrc0Ion / tau;
  double TSrcFloorElc = TSrcFloorIon / tau;

  // Grid parameters
  double vpar_max_elc = 20 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 20 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int num_cell_vpar = 128; // Number of cells in the paralell velocity direction 96
  int num_cell_mu = 192;  // Number of cells in the mu direction 192
  int num_cell_z = 240;
  int poly_order = 1;
  double final_time = 100e-6;
  int num_frames = 20;

  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double B_bt = 1.058278;
  double R_bt = 0.071022;
  double Z_bt = 0.467101;
  double z_bt = 0.468243;

  // Physics parameters at mirror throat
  double n_m = 1.105617e19;
  double Te_m = 346.426583 * eV;
  double Ti_m = 3081.437703 * eV;
  double cs_m = 4.037740e5;

// Initial conditions parameters
  double Ti_perp0 = 10000 * eV;
  double Ti_perp_m = 15000 * eV;
  double Ti_par0 = 7500 * eV;
  double Ti_par_m = 1000 * eV;

  double Te_par0 = 1800 * eV;  
  double Te_par_m = 300 * eV;
  double Te_perp0 = 2000 * eV;
  double Te_perp_m = 3000 * eV;

  struct gk_mirror_ctx ctx = {
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
    .alphaIC0 = alphaIC0,
    .alphaIC1 = alphaIC1,
    .nuFrac = nuFrac,
    .logLambdaElc = logLambdaElc,
    .nuElc = nuElc,
    .logLambdaIon = logLambdaIon,
    .nuIon = nuIon,
    .vti = vti,
    .vte = vte,
    .c_s = c_s,
    .omega_ci = omega_ci,
    .rho_s = rho_s,
    .kperp = kperp, 
    .RatZeq0 = RatZeq0,
    .z_min = z_min,
    .z_max = z_max,
    .psi_eval = psi_eval,
    .mcB = mcB,
    .gamma = gamma,
    .Z_m = Z_m,
    .B_bt = B_bt,
    .R_bt = R_bt,
    .Z_bt = Z_bt,
    .z_bt = z_bt,
    .n_m = n_m,
    .Te_m = Te_m,
    .Ti_m = Ti_m,
    .Ti_perp0 = Ti_perp0,
    .Ti_par0 = Ti_par0,
    .Ti_perp_m = Ti_perp_m,
    .Ti_par_m = Ti_par_m,
    .Te_par0 = Te_par0,
    .Te_par_m = Te_par_m,
    .Te_perp0 = Te_perp0,
    .Te_perp_m = Te_perp_m,
    .cs_m = cs_m,
    .NSrcIon = NSrcIon,
    .lineLengthSrcIon = lineLengthSrcIon,
    .sigSrcIon = sigSrcIon,
    .NSrcFloorIon = NSrcFloorIon,
    .TSrc0Ion = TSrc0Ion,
    .TSrcFloorIon = TSrcFloorIon,
    .NSrcElc = NSrcElc,
    .lineLengthSrcElc = lineLengthSrcElc,
    .sigSrcElc = sigSrcElc,
    .NSrcFloorElc = NSrcFloorElc,
    .TSrc0Elc = TSrc0Elc,
    .TSrcFloorElc = TSrcFloorElc,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .num_cell_z = num_cell_z,
    .num_cell_vpar = num_cell_vpar,
    .num_cell_mu = num_cell_mu,
    .poly_order = poly_order,
    .final_time = final_time,
    .num_frames = num_frames,
  };
  ctx.z_m = 1;
  return ctx;
}

void
write_data(struct gkyl_tm_trigger *iot, gkyl_gyrokinetic_app *app, double tcurr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, tcurr))
  {
    gkyl_gyrokinetic_app_write(app, tcurr, iot->curr - 1);
    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, tcurr, iot->curr - 1);
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
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.num_cell_z);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.num_cell_vpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.num_cell_mu);

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int ccells[] = { NZ };
  int cdim = sizeof(ccells)/sizeof(ccells[0]);
  struct gkyl_range cglobal_r;
  gkyl_create_global_range(cdim, ccells, &cglobal_r);

  // create decomposition
  int cuts[cdim];
#ifdef GKYL_HAVE_MPI  
  for (int d=0; d<cdim; d++)
    cuts[d] = app_args.use_mpi? app_args.cuts[d] : 1;
#else
  for (int d=0; d<cdim; d++) cuts[d] = 1;
#endif  
    
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(cdim, cuts, &cglobal_r);

  // construct communcator for use in app
  struct gkyl_comm *comm;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_gpu && app_args.use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
#else
    printf("Using -g and -M together requires NCCL.\n");
    assert( 0 == 1);
#endif
  } else if (app_args.use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
        .decomp = decomp
      }
    );
  } else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = decomp,
        .use_gpu = app_args.use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp,
      .use_gpu = app_args.use_gpu
    }
  );
#endif

  int my_rank, comm_sz;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &comm_sz);

  int ncuts = 1;
  for (int d=0; d<cdim; d++) ncuts *= cuts[d];
  if (ncuts != comm_sz) {
    if (my_rank == 0)
      fprintf(stderr, "*** Number of ranks, %d, do not match total cuts, %d!\n", comm_sz, ncuts);
    goto mpifinalize;
  }  

  for (int d=0; d<cdim-1; d++) {
    if (cuts[d] > 1) {
      if (my_rank == 0)
        fprintf(stderr, "*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n", cuts[d], d);
      goto mpifinalize;
    }
  }

  if (my_rank == 0) printf("Grid size = %d in Z, %d in Vpar, %d in mu\n", NZ, NV, NMU);
struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-ctx.vpar_max_elc, 0.0},
    .upper = {ctx.vpar_max_elc, ctx.mu_max_elc},
    .cells = {NV, NMU},
    .polarization_density = ctx.n0,
    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density_elc,
      .ctx_upar = &ctx,
      .upar= eval_upar_elc,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_elc,
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_elc,     
    },
    .collisions =  {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = { "ion" },
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .num_sources = 1,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_elc_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_elc_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_elc_source,      
      }, 
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .num_diag_moments = 7, // Copied from GKsoloviev, but
    .diag_moments = {"M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"},
  };
    struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-ctx.vpar_max_ion, 0.0},
    .upper = {ctx.vpar_max_ion, ctx.mu_max_ion},
    .cells = {NV, NMU},
    .polarization_density = ctx.n0,
    .projection = {
      .proj_id = GKYL_PROJ_BIMAXWELLIAN, 
      .ctx_density = &ctx,
      .density = eval_density_ion,
      .ctx_upar = &ctx,
      .upar= eval_upar_ion,
      .ctx_temppar = &ctx,
      .temppar = eval_temp_par_ion,      
      .ctx_tempperp = &ctx,
      .tempperp = eval_temp_perp_ion,    
    },
    .bcx = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    }, 
   .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
    },
    .source = {
      .source_id = GKYL_PROJ_SOURCE,
      .write_source = true,
      .projection[0] = {
        .proj_id = GKYL_PROJ_MAXWELLIAN_PRIM, 
        .ctx_density = &ctx,
        .density = eval_density_ion_source,
        .ctx_upar = &ctx,
        .upar= eval_upar_ion_source,
        .ctx_temp = &ctx,
        .temp = eval_temp_ion_source,      
      }, 
    },
    .num_diag_moments = 7,
    .diag_moments = {"M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"},
  };
  struct gkyl_gyrokinetic_field field = {
    .gkfield_id = GKYL_GK_FIELD_BOLTZMANN,
    .electron_mass = ctx.me,
    .electron_charge = ctx.qe,
    .electron_temp = ctx.Te0,
    .bmag_fac = ctx.B_p, // Issue here. B0 from soloviev, so not sure what to do. Ours is not constant
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
      };
  struct gkyl_gk gk = {  // GK app
    .name = "gk_wham_unif",
    .cdim = 1,
    .vdim = 2,
    .lower = {ctx.z_min},
    .upper = {ctx.z_max},
    .cells = {NZ},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .world = {ctx.psi_eval, 0.0},
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
    .num_species = 1,
    .species = {ion},
    .field = field,
    .use_gpu = app_args.use_gpu,
    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };
  if (my_rank == 0) printf("Creating app object ...\n");
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);  // create app object
  double tcurr = 0.0, tend = ctx.final_time; // start, end and initial time-step
  double dt = tend - tcurr;
  int nframe = ctx.num_frames;
  struct gkyl_tm_trigger io_trig = {.dt = tend / nframe}; // create trigger for IO
  if (my_rank == 0) printf("Applying initial conditions ...\n");
  gkyl_gyrokinetic_app_apply_ic(app, tcurr);  // initialize simulation
  if (my_rank == 0) printf("Computing initial diagnostics ...\n");
  write_data(&io_trig, app, tcurr);
  if (my_rank == 0) printf("Computing initial field energy ...\n");
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  if (my_rank == 0) printf("Starting main loop ...\n");
  long step = 1, num_steps = app_args.num_steps;
  while ((tcurr < tend) && (step <= num_steps))
  {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step at t = %g ...", tcurr);
      struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
      if (step % 100 == 0)
      {
      gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
      }
    if (!status.success)
      {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
          break;
      }
    tcurr += status.dt_actual;
      dt = status.dt_suggested;
      write_data(&io_trig, app, tcurr);
      step += 1;
  }
  if (my_rank == 0) printf(" ... finished\n");
  gkyl_gyrokinetic_app_calc_field_energy(app, tcurr);
  gkyl_gyrokinetic_app_write_field_energy(app);
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
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Updates took %g secs\n", stat.total_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  // simulation complete, free app
  gkyl_gyrokinetic_app_release(app);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  
  mpifinalize:
  ;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Finalize();
#endif
  return 0;
}
