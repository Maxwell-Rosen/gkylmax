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
  int num_cell_vpar;
  int num_cell_mu;
  int num_cell_z;
  int unif_z_cells;
  int num_cell_psi;
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

int get_lower_index(const int nx, const double *x, const double xP)
  //  nx: number of cells in the dimension
  //  x: grid points in the dimension
  //  xP: point to interpolate at
{
  if(xP < x[0])
      return 0;
  else if(xP > x[nx-1])
      return nx-2;
  else
  {
    // Return the rounded down index of the point we are closest to in x
    // use midpoint rule to find the right index
    int i = 0;
    int j = nx-1;
    while(j-i > 1)
    {
      int k = (i+j)/2;
      if(xP < x[k])
        j = k;
      else
        i = k;
    }
    return i;
  }
}

double LI_2D(const int *ncells,  // array of number of cells
             const double *x,            // x grid values
             const double *y,            // y grid values
             const double *f,            // field of function values over the x, y grid (1-D row major)
             const double *pt           // point to interpolate at
             )
{
    int i = get_lower_index(ncells[0],x,pt[0]);
    int j = get_lower_index(ncells[1],y,pt[1]);
    int ip = i+1;
    int jp = j+1;
    double xWta = fmax(0, fmin(1, (pt[0] - x[i])/(x[i+1]-x[i])));
    double yWta = fmax(0, fmin(1, (pt[1] - y[j])/(y[j+1]-y[j])));
    double xWtb = 1-xWta;
    double yWtb = 1-yWta;
    return (f[i *ncells[1] +j ]*xWtb +                    // return f[i ][j ]*xWtb*yWtb +
            f[ip*ncells[1] +j ]*xWta) * yWtb +            //        f[ip][j ]*xWta*yWtb +
           (f[i *ncells[1] +jp]*xWtb +                    //        f[i ][jp]*xWtb*yWta +
            f[ip*ncells[1] +jp]*xWta) * yWta;             //        f[ip][jp]*xWta*yWta;
}

double LI_4D(const int *ncells, // array of number of cells
             const double *x, // psi grid for data
             const double *y, // z grid for data
             const double *z, // vpar grid for data
             const double *w, // mu grid for data
             const double *f, // distribution function at grid locations
             const double *pt // point to interpolate at
             )
{
  // Heavily inspired by  yet adapted from the algorythm at 
  // https://github.com/BYUignite/multilinear_interpolation
  int i = get_lower_index(ncells[0],x,pt[0]); //16
  int j = get_lower_index(ncells[1],y,pt[1]);
  int k = get_lower_index(ncells[2],z,pt[2]);
  int l = get_lower_index(ncells[3],w,pt[3]);

  // Linear interpolate
  int ip = i+1;
  int jp = j+1;
  int kp = k+1;
  int lp = l+1;
  double xWta = fmax(0, fmin(1, (pt[0] - x[i])/(x[i+1]-x[i])));
  double yWta = fmax(0, fmin(1, (pt[1] - y[j])/(y[j+1]-y[j])));
  double zWta = fmax(0, fmin(1, (pt[2] - z[k])/(z[k+1]-z[k])));
  double wWta = fmax(0, fmin(1, (pt[3] - w[l])/(w[l+1]-w[l])));
  double xWtb = 1-xWta;
  double yWtb = 1-yWta;
  double zWtb = 1-zWta;
  double wWtb = 1-wWta;
  int nwzy = ncells[3]*ncells[2]*ncells[1];
  int nwz = ncells[3]*ncells[2];
  int nw = ncells[3];
  return (((f[i *nwzy + j *nwz + k *nw + l ]*xWtb +                                // return f[i ][j ][k ][l ]*xWtb*yWtb*zWtb*wWtb +
            f[ip*nwzy + j *nwz + k *nw + l ]*xWta) * yWtb +                        //        f[ip][j ][k ][l ]*xWta*yWtb*zWtb*wWtb +
            (f[i *nwzy + jp*nwz + k *nw + l ]*xWtb +                                //        f[i ][jp][k ][l ]*xWtb*yWta*zWtb*wWtb +
            f[ip*nwzy + jp*nwz + k *nw + l ]*xWta) * yWta) * zWtb +                //        f[ip][jp][k ][l ]*xWta*yWta*zWtb*wWtb +
          ((f[i *nwzy + j *nwz + kp*nw + l ]*xWtb +                                //        f[i ][j ][kp][l ]*xWtb*yWtb*zWta*wWtb +
            f[ip*nwzy + j *nwz + kp*nw + l ]*xWta) * yWtb +                        //        f[ip][j ][kp][l ]*xWta*yWtb*zWta*wWtb +    
            (f[i *nwzy + jp*nwz + kp*nw + l ]*xWtb +                                //        f[i ][jp][kp][l ]*xWtb*yWta*zWta*wWtb +
            f[ip*nwzy + jp*nwz + kp*nw + l ]*xWta) * yWta) * zWta) * wWtb +        //        f[ip][jp][kp][l ]*xWta*yWta*zWta*wWtb +
          (((f[i *nwzy + j *nwz + k *nw + lp]*xWtb +                                //        f[i ][j ][k ][lp]*xWtb*yWtb*zWtb*wWta +
            f[ip*nwzy + j *nwz + k *nw + lp]*xWta) * yWtb +                        //        f[ip][j ][k ][lp]*xWta*yWtb*zWtb*wWta +
            (f[i *nwzy + jp*nwz + k *nw + lp]*xWtb +                                //        f[i ][jp][k ][lp]*xWtb*yWta*zWtb*wWta +
            f[ip*nwzy + jp*nwz + k *nw + lp]*xWta) * yWta) * zWtb +                //        f[ip][jp][k ][lp]*xWta*yWta*zWtb*wWta +
          ((f[i *nwzy + j *nwz + kp*nw + lp]*xWtb +                                //        f[i ][j ][kp][lp]*xWtb*yWtb*zWta*wWta +
            f[ip*nwzy + j *nwz + kp*nw + lp]*xWta) * yWtb +                        //        f[ip][j ][kp][lp]*xWta*yWtb*zWta*wWta +    
            (f[i *nwzy + jp*nwz + kp*nw + lp]*xWtb +                                //        f[i ][jp][kp][lp]*xWtb*yWta*zWta*wWta +
            f[ip*nwzy + jp*nwz + kp*nw + lp]*xWta) * yWta) * zWta) * wWta;         //        f[ip][jp][kp][lp]*xWta*yWta*zWta*wWta;
}


// Function to calculate the linear index from multi-dimensional indices
size_t calculate_index(const size_t* indices, const int* dims, int num_dims) {
    size_t index = 0;
    size_t multiplier = 1; // Multiplier to calculate the linear index
    for (int i = num_dims - 1; i >= 0; --i) {
        index += indices[i] * multiplier;
        multiplier *= dims[i];
    }
    return index;
}

double* load_binary_file(const char* filename, size_t* num_elements) {
    FILE* file = fopen(filename, "rb");
    // Determine the size of the file
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);
    *num_elements = file_size / sizeof(double);
    double* data = (double*)malloc(file_size);
    size_t elements_read = fread(data, sizeof(double), *num_elements, file);
    fclose(file);
    return data;
}

void
read_ion_distf(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx app = *(struct gk_mirror_ctx*)ctx;
  double* f_dist = app.f_dist_ion;
  double* psi_grid = app.psi_grid;
  double* z_grid = app.z_grid;
  double* v_grid = app.v_grid;
  double* theta_grid = app.theta_grid;
  int* dims = app.dims;

  // FROM GEOMETRY INPUT HARDCOPY
  double z_min_geo = -2.48;
  double z_max_geo =  2.48;
  double tmin = app.z_min;
  double tmax = app.z_max;
  double x_fa[2];
  gkyl_mirror_geo_comp2fieldalligned_advance(t, xn, x_fa, app.mirror_geo_c2fa_ctx);
  double t_norm_cord = (x_fa[1] + tmin) / (tmax - tmin);
  double z_cord = fabs(-z_min_geo + t_norm_cord * (z_max_geo - z_min_geo));

  // Calculate magnetic field at this point
  double interp_pt[4];
  interp_pt[0] = xn[0];
  interp_pt[1] = z_cord;
  double B_val = LI_2D(dims, psi_grid, z_grid, app.B_grid, interp_pt);

  // Must convert the point xn from cartesian to polar
  double vpar = xn[2];
  double mu = xn[3];
  double vperp = sqrt((2.0 * B_val * mu) / app.mi);
  double v = sqrt(pow(vpar, 2) + pow(vperp, 2));
  double theta = atan2(vperp, vpar);

  interp_pt[2] = v;
  interp_pt[3] = theta;

  double interp_val = LI_4D(dims, psi_grid, z_grid, v_grid, theta_grid, f_dist, interp_pt);
  fout[0] = interp_val;
}

void
read_elc_distf(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx app = *(struct gk_mirror_ctx*)ctx;
  double* f_dist = app.f_dist_elc;
  double* psi_grid = app.psi_grid;
  double* z_grid = app.z_grid;
  double* v_grid = app.v_grid;
  double* theta_grid = app.theta_grid;
  int* dims = app.dims;

  // FROM GEOMETRY INPUT HARDCOPY
  double z_min_geo = -2.48;
  double z_max_geo =  2.48;
  double tmin = app.z_min;
  double tmax = app.z_max;
  double x_fa[2];
  gkyl_mirror_geo_comp2fieldalligned_advance(t, xn, x_fa, app.mirror_geo_c2fa_ctx);
  double t_norm_cord = (x_fa[1] + tmin) / (tmax - tmin);
  double z_cord = fabs(-z_min_geo + t_norm_cord * (z_max_geo - z_min_geo));

  double interp_pt[4];
  interp_pt[0] = xn[0];
  interp_pt[1] = z_cord;
  double B_val = LI_2D(dims, psi_grid, z_grid, app.B_grid, interp_pt);

  // Must convert the point xn from cartesian to polar
  double vpar = xn[2];
  double mu = xn[3];
  double vperp = sqrt((2.0 * B_val * mu) / app.me);
  double v = sqrt(pow(vpar, 2) + pow(vperp, 2));
  double theta = atan2(vperp, vpar);

  interp_pt[2] = v;
  interp_pt[3] = theta;

  double interp_val = LI_4D(dims, psi_grid, z_grid, v_grid, theta_grid, f_dist, interp_pt);
  fout[0] = interp_val;
}

void
read_phi(double t, const double *GKYL_RESTRICT xn, double *GKYL_RESTRICT fout, void *ctx)
{
  struct gk_mirror_ctx app = *(struct gk_mirror_ctx*)ctx;
  double* psi_grid = app.psi_grid;
  double* z_grid = app.z_grid;
  double* phi_vals = app.phi_vals;
  int* dims = app.dims;

  double interp_pt[2];
  interp_pt[0] = xn[0];

  // FROM GEOMETRY INPUT HARDCOPY
  double z_min_geo = -2.48;
  double z_max_geo =  2.48;
  double tmin = app.z_min;
  double tmax = app.z_max;
  double x_fa[2];
  gkyl_mirror_geo_comp2fieldalligned_advance(t, xn, x_fa, app.mirror_geo_c2fa_ctx);
  double t_norm_cord = (x_fa[1] + tmin) / (tmax - tmin);
  double z_cord = fabs(-z_min_geo + t_norm_cord * (z_max_geo - z_min_geo));
  interp_pt[1] = z_cord;

  double interp_val = LI_2D(dims, psi_grid, z_grid, phi_vals, interp_pt);
  fout[0] = interp_val * (1 - pow((xn[0] - app.psi_min)/(app.psi_max - app.psi_min),2));
}

void
load_wham_distf(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;

  const char* filename_f_dist_ion = "../binary_files/f_dist_ion.bin";
  size_t num_elements_f_dist_ion;
  double* f_dist_ion = load_binary_file(filename_f_dist_ion, &num_elements_f_dist_ion);
  app->f_dist_ion = f_dist_ion;

  const char* filename_f_dist_elc = "../binary_files/f_dist_elc.bin";
  size_t num_elements_f_dist_elc;
  double* f_dist_elc = load_binary_file(filename_f_dist_elc, &num_elements_f_dist_elc);
  app->f_dist_elc = f_dist_elc;

  const char *filename_phiVals = "../binary_files/phi.bin";
  size_t num_elements_phiVals;
  double *phi_vals = load_binary_file(filename_phiVals, &num_elements_phiVals);
  app->phi_vals = phi_vals; // Need to check units. I think this is kV

  const char *filename_psiGrid = "../binary_files/psiGrid.bin";
  size_t num_elements_psiGrid;
  double *psi_grid = load_binary_file(filename_psiGrid, &num_elements_psiGrid);
  app->psi_grid = psi_grid;

  const char *filename_zGrid = "../binary_files/zGrid.bin";
  size_t num_elements_zGrid;
  double *z_grid = load_binary_file(filename_zGrid, &num_elements_zGrid);
  app->z_grid = z_grid;

  const char *filename_uGrid = "../binary_files/uGrid.bin";
  size_t num_elements_uGrid;
  double *u_grid = load_binary_file(filename_uGrid, &num_elements_uGrid);

  const char *filename_v_norm = "../binary_files/v_norm.bin";
  size_t num_elements_v_norm;
  double *v_norm = load_binary_file(filename_v_norm, &num_elements_v_norm);

  // multiply every element of u_grid with v_norm and call it v_grid
  double *v_grid = (double*)malloc(num_elements_uGrid * sizeof(double));
  for (int i = 0; i < num_elements_uGrid; i++) {
    v_grid[i] = u_grid[i] * v_norm[0];
  }
  free(v_norm); // free v_norm (not needed anymore
  free(u_grid); // free u_grid (not needed anymore)
  app->v_grid = v_grid;
  size_t num_elements_vGrid = num_elements_uGrid;

  const char *filename_theta = "../binary_files/theta.bin";
  size_t num_elements_theta;
  double *theta_grid = load_binary_file(filename_theta, &num_elements_theta);
  app->theta_grid = theta_grid;

  const char *filename_BGrid = "../binary_files/BGrid.bin";
  size_t num_elements_BGrid;
  double *B_grid = load_binary_file(filename_BGrid, &num_elements_BGrid);
  app->B_grid = B_grid;

  int rank = 4;
  int *dims = (int*)malloc(rank * sizeof(int));
  dims[0] = num_elements_psiGrid;
  dims[1] = num_elements_zGrid;
  dims[2] = num_elements_vGrid;
  dims[3] = num_elements_theta;
  app->dims = dims;
  app->rank = rank;
}

void
free_wham_distf(void* ctx)
{
  struct gk_mirror_ctx *app = ctx;
  free(app->f_dist_ion);
  free(app->f_dist_elc);
  free(app->phi_vals);
  free(app->psi_grid);
  free(app->z_grid);
  free(app->v_grid);
  free(app->theta_grid);
  free(app->B_grid);
  free(app->dims);
}

void mapc2p_vel_ion(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_ion = app->vpar_max_ion;
  double mu_max_ion = app->mu_max_ion;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_ion*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_ion*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_ion*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_ion*pow(cmu,2);
}

void mapc2p_vel_elc(double t, const double *vc, double* GKYL_RESTRICT vp, void *ctx)
{
  struct gk_mirror_ctx *app = ctx;
  double vpar_max_elc = app->vpar_max_elc;
  double mu_max_elc = app->mu_max_elc;

  double cvpar = vc[0], cmu = vc[1];
  // Linear map up to vpar_max/2, then quadratic.
  if (fabs(cvpar) <= 0.5)
    vp[0] = vpar_max_elc*cvpar;
  else if (cvpar < -0.5)
    vp[0] = -vpar_max_elc*2.0*pow(cvpar,2);
  else
    vp[0] =  vpar_max_elc*2.0*pow(cvpar,2);

  // Quadratic map in mu.
  vp[1] = mu_max_elc*pow(cmu,2);
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
  double z_min = -M_PI + 1e-1;
  double z_max = M_PI - 1e-1;
  double psi_min = 1e-3; // Go smaller. 1e-4 might be too small
  double psi_max = 1e-2; // aim for 2e-2

  // Grid parameters
  double vpar_max_elc = 20 * vte;
  double mu_max_elc = me * pow(3. * vte, 2.) / (2. * B_p);
  double vpar_max_ion = 20 * vti;
  double mu_max_ion = mi * pow(3. * vti, 2.) / (2. * B_p);
  int num_cell_vpar = 64; // 96 uniform
  int num_cell_mu = 96;  // 192 uniform
  int num_cell_z = 144;
  int unif_z_cells = 288;
  int num_cell_psi = 16;
  int poly_order = 1;
  double final_time = 3e-6;
  int num_frames = 60;
  int int_diag_calc_num = num_frames*100;
  double dt_failure_tol = 1.0e-4; // Minimum allowable fraction of initial time-step.
  int num_failures_max = 20; // Maximum allowable number of consecutive small time-steps.

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
    .z_min = z_min,
    .z_max = z_max,
    .psi_min = psi_min,
    .psi_max = psi_max,
    .vpar_max_ion = vpar_max_ion,
    .vpar_max_elc = vpar_max_elc,
    .mu_max_ion = mu_max_ion,
    .mu_max_elc = mu_max_elc,
    .num_cell_psi = num_cell_psi,
    .num_cell_z = num_cell_z,
    .unif_z_cells = unif_z_cells,
    .num_cell_vpar = num_cell_vpar,
    .num_cell_mu = num_cell_mu,
    .poly_order = poly_order,
    .final_time = final_time,
    .num_frames = num_frames,
    .int_diag_calc_num = int_diag_calc_num,
    .dt_failure_tol = dt_failure_tol,
    .num_failures_max = num_failures_max,
  };
  load_wham_distf(&ctx);
  ctx.mirror_geo_c2fa_ctx = gkyl_malloc(sizeof(struct gkyl_mirror_geo_c2fa_ctx));
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
  if (app_args.use_mpi)
    MPI_Init(&argc, &argv);
#endif

  if (app_args.trace_mem)
  {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  struct gk_mirror_ctx ctx = create_ctx(); // context for init functions
  int NPSI = APP_ARGS_CHOOSE(app_args.xcells[0], ctx.num_cell_psi);
  int NZ = APP_ARGS_CHOOSE(app_args.xcells[1], ctx.num_cell_z);
  int NV = APP_ARGS_CHOOSE(app_args.vcells[0], ctx.num_cell_vpar);
  int NMU = APP_ARGS_CHOOSE(app_args.vcells[1], ctx.num_cell_mu);

  int nrank = 1; // number of processors in simulation
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi)
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
#endif  

  // create global range
  int ccells[] = { NPSI, NZ };
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

  if (my_rank == 0) printf("Grid size = %d in psi, %d in Z, %d in Vpar, %d in mu\n", NPSI, NZ, NV, NMU);

  struct gkyl_gyrokinetic_projection elc_ic = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = read_elc_distf,
      .ctx_func = &ctx,  
  };
  struct gkyl_gyrokinetic_species elc = {
    .name = "elc",
    .charge = ctx.qe,
    .mass = ctx.me,
    .lower = {-1.0/sqrt(2.), 0.0},
    .upper = { 1.0/sqrt(2.), 1.0},
    .cells = {NV, NMU},
    .polarization_density = ctx.n0,
    .projection = elc_ic,
    .mapc2p = {
      .mapping = mapc2p_vel_elc,
      .ctx = &ctx,
    },
    .bcx = {
      .lower = {.type = GKYL_SPECIES_REFLECT,},
      .upper = {.type = GKYL_SPECIES_ZERO_FLUX,},
    },
    .bcy = {
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
    .num_diag_moments = 7,
    .diag_moments = {"M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"},
  };
  struct gkyl_gyrokinetic_projection ion_ic = {
      .proj_id = GKYL_PROJ_FUNC,
      .func = read_ion_distf,
      .ctx_func = &ctx,  
  };
  struct gkyl_gyrokinetic_species ion = {
    .name = "ion",
    .charge = ctx.qi,
    .mass = ctx.mi,
    .lower = {-1.0/sqrt(2.), 0.0},
    .upper = { 1.0/sqrt(2.), 1.0},
    .cells = {NV, NMU},
    .polarization_density = ctx.n0,
    .no_by = true,
    .nG_from_npol = true,
    .projection = ion_ic,
    .mapc2p = {
      .mapping = mapc2p_vel_ion,
      .ctx = &ctx,
    },
    .bcx = {
      .lower = {.type = GKYL_SPECIES_REFLECT},
      .upper = {.type = GKYL_SPECIES_ZERO_FLUX},
    },
    .bcy = {
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
    .num_diag_moments = 7,
    .diag_moments = {"M0", "M1", "M2", "M2par", "M2perp", "M3par", "M3perp"},
  };
  struct gkyl_gyrokinetic_field field = {
    .fem_parbc = GKYL_FEM_PARPROJ_NONE,
    .poisson_bcs = {
      .lo_type = {GKYL_POISSON_NEUMANN},
      .up_type = {GKYL_POISSON_DIRICHLET},
      .lo_value = {0.0},
      .up_value = {0.0},
    },
    .polarization_phi = read_phi,
    .polarization_phi_ctx = &ctx,
  };
  struct gkyl_gk gk = {  // GK app
    .name = "gk_wham",
    .cdim = 2,
    .vdim = 2,
    .lower = {ctx.psi_min, ctx.z_min},
    .upper = {ctx.psi_max, ctx.z_max},
    .cells = {NPSI, NZ},
    .poly_order = ctx.poly_order,
    .basis_type = app_args.basis_type,
    .geometry = {
      // .geometry_id = GKYL_GEOMETRY_FROMFILE
      .geometry_id = GKYL_MIRROR,
      .world = {0.0},
      .mirror_efit_info = &inp,
      .mirror_grid_info = &ginp,
      .mirror_geo_c2fa_ctx = ctx.mirror_geo_c2fa_ctx,
      .nonuniform_mapping_fraction = (1 - (double)NZ / ctx.unif_z_cells),
    },
    .num_periodic_dir = 0,
    .periodic_dirs = {},
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
  
  // Create app object.
  gkyl_gyrokinetic_app *app = gkyl_gyrokinetic_app_new(&gk);

  // Initial and final simulation times.
  int frame_curr = 0;
  double t_curr = 0.0, t_end = ctx.final_time;
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

  long step = 1;
  while ((t_curr < t_end) && (step <= app_args.num_steps)) {
    gkyl_gyrokinetic_app_cout(app, stdout, "Taking time-step %ld at t = %g ...", step, t_curr);
      struct gkyl_update_status status = gkyl_gyrokinetic_update(app, dt);
      gkyl_gyrokinetic_app_cout(app, stdout, " dt = %g\n", status.dt_actual);

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
  if (stat.nstage_2_fail > 0) {
    gkyl_gyrokinetic_app_cout(app, stdout, "  Max rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[1]);
    gkyl_gyrokinetic_app_cout(app, stdout, "  Min rel dt diff for RK stage-2 failures %g\n", stat.stage_2_dt_diff[0]);
  }
  gkyl_gyrokinetic_app_cout(app, stdout, "Number of RK stage-3 failures %ld\n", stat.nstage_3_fail);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species RHS calc took %g secs\n", stat.species_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisions RHS calc took %g secs\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Field RHS calc took %g secs\n", stat.field_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Species collisional moments took %g secs\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, stdout, "Total updates took %g secs\n", stat.total_tm);

  gkyl_gyrokinetic_app_cout(app, stdout, "Number of write calls %ld,\n", stat.nio);
  gkyl_gyrokinetic_app_cout(app, stdout, "IO time took %g secs \n", stat.io_tm);

// Crude and should be moved to a lower layer. Not sure how to access the context there though. It's not in app
// I think I would at least need to add a property to app->gk_geom with the geometry_id and the arcl_map_ctx
  freeresources:
  // Free resources after simulation completion.
  gkyl_gyrokinetic_app_release(app);
  free_wham_distf(&ctx);
  struct gkyl_mirror_geo_c2fa_ctx *c2fa_release = ctx.mirror_geo_c2fa_ctx;
  gkyl_array_release(c2fa_release->c2fa);
  gkyl_array_release(c2fa_release->c2fa_deflate);
  gkyl_free(c2fa_release);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  mpifinalize:
#ifdef GKYL_HAVE_MPI
  if (app_args.use_mpi) {
    MPI_Finalize();
  }
#endif
  
  return 0;
}
