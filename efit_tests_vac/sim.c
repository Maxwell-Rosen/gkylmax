#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>
#include <gkyl_basis.h>

#include <gkyl_efit.h>


void test_wham(){
  char* filepath = "../eqdsk/wham_vac_hires.geqdsk";
  int rzpoly_order = 2;
  int fluxpoly_order = 1;
  struct gkyl_efit* efit = gkyl_efit_new(filepath,rzpoly_order,  GKYL_BASIS_MODAL_SERENDIPITY, fluxpoly_order, false);

  printf( "rdim=%g zdim=%g rcentr=%g rleft=%g zmid=%g  rmaxis=%g zmaxis=%g simag=%1.16e sibry=%1.16e bcentr=%g  current=%g simag=%g rmaxis=%g   zmaxis=%g sibry=%g \n", efit->rdim, efit->zdim, efit->rcentr, efit->rleft, efit->zmid, efit->rmaxis, efit->zmaxis, efit->simag, efit->sibry, efit->bcentr, efit-> current, efit->simag, efit->rmaxis, efit-> zmaxis, efit->sibry);
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psizr, "wham_psi.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyrzr, "wham_psibyr.gkyl");
  gkyl_grid_sub_array_write(efit->rzgrid, efit->rzlocal, 0, efit->psibyr2zr, "wham_psibyr2.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->fpolflux, "wham_fpol.gkyl");
  gkyl_grid_sub_array_write(efit->fluxgrid, efit->fluxlocal, 0, efit->qflux, "wham_q.gkyl");

  gkyl_efit_release(efit);

}

void main(int argc, char **argv) {
  test_wham();
}