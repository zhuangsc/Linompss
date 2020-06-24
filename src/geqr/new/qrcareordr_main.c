#include "qrca_kernels.h"

extern double **gl_Ah;
extern double **gl_Sh;
extern double **gl_Th;


void qrcareordr(int diagl, int tr, int tc, int bs, int mt, int nt) {
	int d;
	for(d=0;d<diagl;d++) {

		// Process [d,0:d-1]
		int j;
		for(j=0;j<d;j++) {
			// Update element j in row d
			int i;
			for(i=0;i<j;i++) {
				NoFLA_Apply_td_QT_var31a( bs, tr, tc, 0, 
					gl_Ah[ mt*i+d ], 
                                    		gl_Sh[ mt*i+d ], 
                                    		gl_Ah[ mt*j+i ], 
                                    		gl_Ah[ mt*j+d ]);
					//printf("%i split dlarfb %i|%i,%i\n",d,i,d,j);
			}

			NoFLA_Compute_td_QR_var31a( bs, tr, tc, 0, 
					gl_Ah[ mt*j+j ], /*ts,*/
			    		gl_Ah[ mt*j+d ], /*ts,*/
			    		gl_Th[ mt*j+d ],
			    		gl_Sh[ mt*j+d ] /*,bs*/ );

                        	//printf("%i split dgeqt2 %i|%i,%i\n",d,d,j,j);
		}


		// Process [0:d-1,d]
		int i;
		for(i=0;i<d;i++) {
			// Update element [i,d]

			int j;
			for(j=0;j<i;j++) {
				NoFLA_Apply_td_QT_var31a( bs, tr, tc, 0, 
                                    	gl_Ah[ mt*j+i ], 
                                    	gl_Sh[ mt*j+i ], 
                                    	gl_Ah[ mt*d+j ], 
                                    	gl_Ah[ mt*d+i ]);
			}

			dlarfb_task( bs, tr, tc, tc, 0, 
					gl_Ah[ mt*i+i ], 
					gl_Sh[ mt*i+i ], 
					gl_Ah[ mt*d+i ]);
		}


		int dd;
		for(dd=0;dd<d;dd++) {
			NoFLA_Apply_td_QT_var31a( bs, tr, tc, 0, 
		    			gl_Ah[ mt*dd+d ], 
			    		gl_Sh[ mt*dd+d ], 
			    		gl_Ah[ mt*d+dd ], 
			    		gl_Ah[ mt*d+d ]);
		}

		dgeqt2_task( bs, tr, tc, 0,  
				gl_Ah[ mt*d+d ], 
				gl_Th[ mt*d+d ],
				gl_Sh[ mt*d+d ]);

     			//           printf("%i split dlarfb  (%i,%i) %i\n",d,d,d,task++);
      			//          printf("%i dgeqt2 (%i,%i) %i\n\n",d,d,d,task++);
	}

	for(d=diagl;d<mt;d++) {
		int j;
		for(j=0;j<nt;j++) {
			int i;
			for(i=0;i<j;i++) {
				NoFLA_Apply_td_QT_var31a( bs, tr, tc, 0, 
                                    		gl_Ah[ mt*i+d ], 
                                    		gl_Sh[ mt*i+d ], 
                                    		gl_Ah[ mt*j+i ], 
                                    		gl_Ah[ mt*j+d ]);
			}

			NoFLA_Compute_td_QR_var31a( bs, tr, tc, 0, 
					gl_Ah[ mt*j+j ], /*ts,*/
					gl_Ah[ mt*j+d ], /*ts,*/
					gl_Th[ mt*j+d ],
					gl_Sh[ mt*j+d ] /*,bs*/ );
		}
	}
}
