#include "qrca_kernels.h"
#include "qrca_utils.h"


void split_dlarfb_hier_task(int cm, int cn, double *C, int ldim_C, double *U, int ldim_U, int dimt, double *T, int ldim_T, int fm, double *F, double *work)
{
	double d_mone=-1.00;
	double d_one=1.00;

  	NoFLA_Copy(fm, cn, F, fm, work, fm);

	dtrmm_("Left", "Upper", "Transpose", "Non-Unit",
		&fm, &cn, &d_one,
		T, &ldim_T,
		work, &fm);	

	//printf("cm %i cn %i fm %i ldim_U %i ldim_C %i\n", cm, cn, fm, ldim_U, ldim_C);
	dgemm_("No Transpose", "No Transpose", 
		&cm, &cn, &fm, 
		&d_mone, U, &ldim_U, 
		work, &fm, 
		&d_one, C, &ldim_C);
}
