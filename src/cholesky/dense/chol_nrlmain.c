#include "chol_nrlmain.h"
#include "tasks_nested_gemm.h"
#include "tasks_nested_syrk.h"
#include "tasks_nested_trsm.h"
#include "tasks_nested_potrf.h"


#include "fptype.h"


int CHOL_NRL(int mt, int sb, int b, int t, fp_t **Ah) 
{
	int k;
	for ( k = 0; k < mt ;k++ ) {

		NTASK_POTRF(sb, b, t, Ah[k* mt + k]);

		int i;
		for (i = k+1; i < mt; i++) {
			NTASK_TRSM(sb, b, t, Ah[k * mt + k], Ah[k * mt + i]);
		}

		for (i = k+1; i < mt; i++) {
			int j;
			for (j=k+1; j<i; j++) {
				NTASK_GEMM(sb, b, t, Ah[k * mt + i], Ah[k * mt + j], Ah[j * mt + i]);
			}
			NTASK_SYRK(sb, b, t, Ah[k * mt + i], Ah[i * mt + i]);
		}
	}

	return 0;
}
