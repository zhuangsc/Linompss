#include "qrca_kernels.h"
#include "qrca_utils.h"


extern double **gl_Ah;
extern double **gl_Th;
extern double **gl_Sh;

#if 0
static inline __attribute__((always_inline)) unsigned int log2( unsigned int x ) {
  unsigned int ans = 0 ;
  while( x>>=1 ) ans++;
  return ans ;
}
#endif

void qrcabin(int diagl, int tr, int tc, int bs, int mt, int nt) {
	printf("mt %i nt %i\n",mt,nt);
	int rk=0;
	int k;
	for( k = 0; k < diagl; k++ ) {

		unsigned int cmt = mt==1?0:mt;
		unsigned int mtr = cmt - k;
		int jump = 1;

		//printf("k=%i jmp=%i mtr=%i \n",k,jump,mtr);

		int b;
		for ( b=k; b<mt; b++) {
			//printf("dgeqr (%i,%i)\n",b,k);
			dgeqt2_task( bs, tr, tc, 0,  
				gl_Ah[ mt*k+b ], 
				gl_Th[ mt*k+b ],
				gl_Sh[ mt*k+b ]);

			int j;
			for( j = k+1; j<nt; j++ ) {
	//			printf("dlarfb (%i,%i)\n",b,j);
				dlarfb_task( bs, tr, tc, tc, 0, 
					gl_Ah[ mt*k+b ], 
					gl_Sh[ mt*k+b ], 
					gl_Ah[ mt*j+b ]);
			}
		}

		int pc = mt>>1;
		int l = cmt - k;
		while ( l > 1 )	{
			int djump = jump<<1;
			int pc = l >> 1; 
			int pl = l & 1;
			int cpc = pc==0?1:pc;
			//printf("cpc=%i l=%i\n",cpc,l);

			int p = k; 
			int pi;
			for ( pi=0; pi<cpc; pi++, p+=djump ) {
				//printf("\t%i %i\n",p,p+jump);

#if 1
	//			printf("dgeqr_split %i %i,%i\n",k,p,p+jump);
				NoFLA_Compute_td_QR_var31b( bs, tr, tc, 0, 
					gl_Ah[ mt*k+p ], /*ts,*/
			    		gl_Ah[ mt*k+p+jump ], /*ts,*/
			    		gl_Th[ mt*k+p+jump ],
                                    	gl_Sh[ mt*k+p+jump ] /*,bs*/ );
#endif
#if 1
				int j;
				for( j = k+1; j < nt; j++ ) {
	//				printf("dlarfb_split %i %i,%i\n",j,p,p+jump);
					NoFLA_Apply_td_QT_var31a( bs, tr, tc, 0, 
						gl_Ah[ mt*k+p+jump ], 
						gl_Sh[ mt*k+p+jump ], 
						gl_Ah[ mt*j+p ], 
						gl_Ah[ mt*j+p+jump ]);
				}
#endif
			}

			

			l = l - cpc;
			jump = jump<<1;
			pc = pc>>1;
		} 
		//printf("\n");
	}
}
