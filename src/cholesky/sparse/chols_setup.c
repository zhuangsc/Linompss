#include "chol_setup.h"

#include <time.h>
#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"
#include "chols_warm.h"
#include "symfac.h"
#include "matfprint.h"


extern void *Acsr;
extern void *Acsc;
extern void *Ahb;
extern int format;
extern int *work;
extern long nzA;


int chol_setup(int check, int m, int mr, int ts, int bs, int tm, int mleft) {
	struct timeval start, stop;

	gettimeofday(&start, NULL);

	one2zero(Ahb);
	hbmat_t *lAcsr, *Ahhb, *Atrans;
	int *etree_ptr;

	if ( format == MAT_CSC ) {
		lAcsr = malloc(sizeof(hbmat_t));
		hb_init(lAcsr, Ahb);
		hb_csrcsc(lAcsr, Ahb, CSC2CSR);
		etree_ptr = etree(lAcsr);
		Ahhb = hb2hbh_sym_etree(Ahb, ts, etree_ptr);
		Atrans = malloc(sizeof(hbmat_t));
		hbh_init(Atrans, Ahhb);
		hbh_csrcsc(Atrans, Ahhb, CSC2CSR);
	} else{
		etree_ptr = etree(Ahb);
		//Ahhb = hb2hbh_sym_etree_csr_p(Ahb, ts, etree_ptr);
		Ahhb = hb2hbh_symcsr(Ahb, ts, etree_ptr, 1);
		Atrans = malloc(sizeof(hbmat_t));
		hbh_init(Atrans, Ahhb);
		hbh_csrcsc(Ahhb, Atrans, CSR2CSC);
	}

	hbmat_t* Z = Ahb;
	nzA = Z->elemc;
	gettimeofday(&stop, NULL);

	unsigned long elapsed = stop.tv_usec - start.tv_usec;
	elapsed += (stop.tv_sec - start.tv_sec) * 1000000;
	printf("info: symbolic fac %lu us\n", elapsed);

	if ( !format ) {
		hb_free(lAcsr);
	}

	if ( format==MAT_CSC ) {
		Acsc = Ahhb;
		Acsr = Atrans;
	} else {
		Acsc = Atrans;
		Acsr = Ahhb;
	}

	//fprint_csc("Aorig.mat", Ahb, 0);
	work = malloc(mr* sizeof(int));
	chol_warmup();

	return 0;
}

void chol_shutdown() {
	if ( !format ) {
		hb_free(Acsc);
	}
	else {
		hbh_free2(Acsr);
	}

	free(work);
	free(Ahb);
//	free(etree_g);
//	free(vptr_pool);
//	free(vpos_pool);
//	free(vval_pool);
}
