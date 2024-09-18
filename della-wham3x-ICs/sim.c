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
  double psi_max;
  // Magnetic equilibrium model.
  double mcB;
  double gamma;
  double Z_m;
  // Bananna tip info. Hardcoad to avoid dependency on ctx
  double z_m;
  double Z_m_computational;
  // Physics parameters at mirror throat
  double vpar_max_ion;
  double vpar_max_elc;
  double mu_max_ion;
  double mu_max_elc;
  int Nx;
  int Ny;
  int Nz;
  int Nvpar;
  int Nmu;
  int cells[GKYL_MAX_DIM]; // Number of cells in all directions.
  int poly_order;
  double t_end;
  int num_frames;
  int int_diag_calc_num; // Number of integrated diagnostics computations (=INT_MAX for every step).
  double dt_failure_tol; // Minimum allowable fraction of initial time-step.
  int num_failures_max; // Maximum allowable number of consecutive small time-steps.

  
  // Initial conditions reading
  double *f_dist_ion;
  double *f_dist_elc;
  double *phi_vals;
  double *psi_grid;
  double *z_grid;
  double *v_grid;
  double *theta_grid;
  double *B_grid;
  int *dims;
  int rank;

  //Nonuniform grid
  void *mirror_geo_c2fa_ctx;

};


struct gkyl_mirror_geo_efit_inp inp = {
  // psiRZ and related inputs
  // .filepath = "../eqdsk/wham_dia_hires.geqdsk",
  .filepath = "../eqdsk/wham_vac_hires.geqdsk",
  .rzpoly_order = 2,
  .fluxpoly_order = 1,
  .plate_spec = false,
  .quad_param = {  .eps = 1e-10 }
};


struct gkyl_mirror_geo_grid_inp ginp = {
  .rclose = 0.2,
  .zmin = -2.0,
  .zmax =  2.0,
  .write_node_coord_array = true,
  .node_file_nm = "wham_nodes.gkyl",
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
  // vp[0] = vc[0];
  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
}

struct gk_mirror_ctx
create_ctx(void)
{
  int cdim = 3, vdim = 2; // Dimensionality.

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

  // Geometry parameters.
  double z_min = -M_PI + 1e-1;
  double z_max = M_PI - 1e-1;
  double psi_min = 5e-4; // Go smaller. 1e-4 might be too small
  double psi_max = 3e-3; // aim for 2e-2

  // Grid parameters
  double vpar_max_elc = 30 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 30 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int Nvpar = 16; // 96 uniform
  int Nmu = 16;  // 192 uniform
  int Nz = 16;
  int Ny = 4;
  int Nx = 4;
  int poly_order = 1;
  double t_end = 500e-6;
  int num_frames = 5000;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

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
    .z_min = z_min,
    .z_max = z_max,
    .psi_min = psi_min,
    .psi_max = psi_max,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .Nx = Nx,
    .Ny = Ny,
    .Nz = Nz,
    .Nvpar = Nvpar,
    .Nmu = Nmu,
    .cells = {Nx, Ny, Nz, Nvpar, Nmu},
    .poly_order = poly_order,
    .t_end = t_end,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
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
write_data(struct gkyl_tm_trigger* iot, gkyl_gyrokinetic_app* app, double t_curr, bool force_write)
{
  bool trig_now = gkyl_tm_trigger_check_and_bump(iot, t_curr);
  if (trig_now || force_write) {
    int frame = (!trig_now) && force_write? iot->curr : iot->curr-1;

    gkyl_gyrokinetic_app_write(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_mom(app);
    gkyl_gyrokinetic_app_write_mom(app, t_curr, frame);
    gkyl_gyrokinetic_app_write_source_mom(app, t_curr, frame);

    gkyl_gyrokinetic_app_calc_field_energy(app, t_curr);
    gkyl_gyrokinetic_app_write_field_energy(app);

    gkyl_gyrokinetic_app_calc_integrated_mom(app, t_curr);
    gkyl_gyrokinetic_app_write_integrated_mom(app);
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

  // Create decomposition.
  struct gkyl_rect_decomp *decomp = gkyl_gyrokinetic_comms_decomp_new(ctx.cdim, cells_x, app_args.cuts, app_args.use_mpi, stderr);

  // Construct communicator for use in app.
  struct gkyl_comm *comm = gkyl_gyrokinetic_comms_new(app_args.use_mpi, app_args.use_gpu, decomp, stderr);

  int my_rank = 0;
  int comm_sz = 1;
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi){
    gkyl_comm_get_rank(comm, &my_rank);
    int comm_sz;
    gkyl_comm_get_size(comm, &comm_sz);
  }
#endif

  if (my_rank == 0) {
    printf("Grid size = %d in psi, %d in theta, %d in Z, %d in Vpar, %d in mu\n", cells_x[0], cells_x[1], cells_x[2], cells_v[0], cells_v[1]);
    if (app_args.use_mpi)
      printf("Number of MPI ranks: %d\n", decomp->ndecomp);
    if (app_args.use_gpu)
      printf("Number of GPUs: %d\n", decomp->ndecomp);
    printf("psi_min = %g, psi_max = %g\n", ctx.psi_min, ctx.psi_max);
    printf("z_min = %g, z_max = %g\n", ctx.z_min, ctx.z_max);
    printf("vpar_max_ion/vti = %g, mu_max_ion/mu_ti = %g\n", ctx.vpar_max_ion/ctx.vti, sqrt(ctx.mu_max_ion/ctx.mi*2.0*ctx.B_p)/ctx.vti);
    printf("vpar_max_elc/vte = %g, mu_max_elc/mu_te = %g\n", ctx.vpar_max_elc/ctx.vte, sqrt(ctx.mu_max_elc/ctx.me*2.0*ctx.B_p)/ctx.vte);
    printf("vti = %.4e, vte = %.4e, c_s = %.4e, mu_ti = %.4e, mu_te = %.4e\n", ctx.vti, ctx.vte, ctx.c_s, ctx.mi * pow(ctx.vti, 2.) / (2. * ctx.B_p),
     ctx.me * pow(ctx.vte, 2.) / (2. * ctx.B_p));
    printf("omega_ci = %.4e, rho_s = %.4e, kperp = %.4e\n", ctx.omega_ci, ctx.rho_s, ctx.kperp);
    printf("1/nuElc = %.4e, 1/nuIon = %.4e\n", 1./ctx.nuElc, 1./ctx.nuIon);
  }

  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-1.0, 0.0},
    .upper = { 1.0, 1.0},
    .cells = { cells_v[0], cells_v[1]},
    .polarization_density = ctx.n0,
    .no_by = true,
    .enforce_positivity = true,
    .init_from_file = {
      .type = GKYL_IC_IMPORT_F,
      .file_name = "gk_wham2x-elc_0_lores.gkyl",
      // .conf_scale = ic_conf_fac,
      // .conf_scale_ctx = &ctx,
    },
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },
    .bcx = {
      .lower = {.type = GKYL_SPECIES_REFLECT,},
      .upper = {.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuElc,
      .num_cross_collisions = 1,
      .collide_with = {"ion"},
    },
    .num_diag_moments = 1,
    .diag_moments = {"BiMaxwellianMoments"},
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
    .enforce_positivity = true,
    .init_from_file = {
      .type = GKYL_IC_IMPORT_F,
      .file_name = "gk_wham2x-ion_0_lores.gkyl",
      // .conf_scale = ic_conf_fac,
      // .conf_scale_ctx = &ctx,
    },
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    .bcx = {
      .lower = {.type = GKYL_SPECIES_REFLECT},
      .upper = {.type = GKYL_SPECIES_ZERO_FLUX},
    },
    .bcz = {
      .lower={.type = GKYL_SPECIES_GK_SHEATH,},
      .upper={.type = GKYL_SPECIES_GK_SHEATH,},
    },    
    .collisions = {
      .collision_id = GKYL_LBO_COLLISIONS,
      .ctx = &ctx,
      .self_nu = evalNuIon,
      .num_cross_collisions = 1,
      .collide_with = {"elc"},
    },
    .num_diag_moments = 1,
    .diag_moments = {"BiMaxwellianMoments"},
  };
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .poisson_bcs = {
      .lo_type = {GKYL_POISSON_NEUMANN,   GKYL_POISSON_PERIODIC},
      .up_type = {GKYL_POISSON_DIRICHLET, GKYL_POISSON_PERIODIC},
      .lo_value = {0.0, 0.0},
      .up_value = {0.0, 0.0},
    },
  };
  struct gkyl_gk app_inp = {  // GK app
    .name = "gk_wham",
    .cdim = 3,  .vdim = 2,
    .lower = {ctx.psi_min, -M_PI, ctx.z_min},
    .upper = {ctx.psi_max,  M_PI, ctx.z_max},
    .cells = { cells_x[0], cells_x[1], cells_x[2] },
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = {
      .geometry_id = GKYL_MIRROR,
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
    },
    .num_periodic_dir = 1,
    .periodic_dirs = {1},
    .num_species = 2,
    .species = {elc, ion},
    .field = field,
    .use_gpu = app_args.use_gpu,
    .has_low_inp = true,
    .low_inp = {
      .local_range = decomp->ranges[my_rank],
      .comm = comm
    }
  };

  // start timer
  clock_t start_time = clock();
  
  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&app_inp);

  clock_t end_time = clock();
  printf("Time to create app object: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

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
    start_time = clock();
    gkyl_gyrokinetic_app_apply_ic(app, t_curr);
    end_time = clock();
    printf("Time to apply IC: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);
  }  

  start_time = clock();
  // Create triggers for IO.
  int num_frames = ctx.num_frames, num_int_diag_calc = ctx.int_diag_calc_num;
  struct gkyl_tm_trigger trig_write = { .dt = t_end/num_frames, .tcurr = t_curr, .curr = frame_curr };
  struct gkyl_tm_trigger trig_calc_intdiag = { .dt = t_end/GKYL_MAX2(num_frames, num_int_diag_calc),
    .tcurr = t_curr, .curr = frame_curr };

  // Write out ICs (if restart, it overwrites the restart frame).
  calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, false);
  write_data(&trig_write, app, t_curr, false);

  double dt = t_end-t_curr; // Initial time step.
  // Initialize small time-step check.
  double dt_init = -1.0, dt_failure_tol = ctx.dt_failure_tol;
  int num_failures = 0, num_failures_max = ctx.num_failures_max;
  end_time = clock();
  printf("Time to write IC: %g\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);    
    // if (step % 1000 == 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
    gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);
    // }

    if (!status.success) {
      gkyl_gyrokinetic_app_cout(app, stdout, "** Update method failed! Aborting simulation ....\n");
          break;
      }

    t_curr += status.dt_actual;
      dt = status.dt_suggested;

    calc_integrated_diagnostics(&trig_calc_intdiag, app, t_curr, t_curr > t_end);
    write_data(&trig_write, app, t_curr, t_curr > t_end);

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
        write_data(&trig_write, app, t_curr, true);
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
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  gkyl_gyrokinetic_comms_release(decomp, comm);

#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  return 0;
}
