#include "qrca_kernels.h"


#include<math.h>


extern double **Ah;
extern double **Th;
extern double **Sh;


void qrcab(int diagl, int tr, int tc, int bs, int mt, int nt) {
	int rk=0;
	int k;
	for (k=0; k<diagl; k++ ) {
		
		int d=1;
		int mtl=mt-k;
		int mtlodd=mtl & 1;
		int treedepth=log2(mtl+mtlodd);
		printf("\nk %i depth %i mtl %i\n",k,treedepth,mtl);
		int l;
		for (l=0; l<treedepth; l++) {
			int jd=d<<1;
			int jodd=mtlodd?!(l&1):0;
			int j;
			for(j=k; j<mt-jodd; j+=jd) {
				printf("%i & %i\n",j,j+d);
			}
			if(jodd) {

				dgeqt2_task( bs, mt, nt, 0, 
        			double * A, double * T, double * S); 
				printf("odd %i\n",j);
			}

			d=jd;
		}
		
#if 0
		for (j=k; k<nt; k+=2) {
		int row=rk / tr;
		int rskip=rk % tr;
		dgeqt2_task( bs, tr, tc, rskip,  
			Ah[ mt*k+row ], 
			Th[ mt*k+row ],
			Sh[ mt*k+row ]);


		int j;
		for( j = k+1; j < nt; j++ ) {
			dlarfb_task( bs, tr, tc, tc, rskip, 
				Ah[ mt*k+row ], 
				Sh[ mt*k+row ], 
				Ah[ mt*j+row ]);
		}

		int i;
		for( i = row+1; i < mt; i++ ) {
			NoFLA_Compute_td_QR_var31a( bs, tr, tc, rskip, 
				Ah[ mt*k+row ], /*ts,*/
			    	Ah[ mt*k+i ], /*ts,*/
			    	Th[ mt*k+i ],
                                    Sh[ mt*k+i ] /*,bs*/ );

			int j;
			for( j = k+1; j < nt; j++ ) {
				NoFLA_Apply_td_QT_var31a( bs, tr, tc, rskip, 
						Ah[ mt*k+i ], 
						Sh[ mt*k+i ], 
						Ah[ mt*j+row ], 
						Ah[ mt*j+i ]);
			}
		}
		rk+=tc;
#endif
	}
}
