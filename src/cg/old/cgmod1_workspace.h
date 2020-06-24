#ifndef __CGMOD1_WORKSPACE_H__
#define __CGMOD1_WORKSPACE_H__


#include <stdlib.h>

#include "fptype.h"


static inline __attribute__((always_inline)) int cgmod1_malloc(size_t n, int s, int dupl, fp_t *work, unsigned long works, fp_t ***p, fp_t ***v, fp_t ***r, fp_t ***z,\
	fp_t ***gamma, fp_t ***sigma, fp_t ***delta) {
	unsigned long ns = n * s;
	unsigned long nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	unsigned long salgn = ((s + 0x7f ) >> 7 ) << 7;

	fp_t **ptmp = malloc(8 * dupl * sizeof(fp_t*));
	if ( ptmp == NULL ) {
		return 1;
	}

	//printf("workspace %i\n", works);

	int offs = 0;
	/* ns */
	fp_t **lp = *p = &ptmp[offs];
	offs += dupl;
	fp_t **lv = *v = &ptmp[offs];
	offs += dupl;
	fp_t **lr = *r = &ptmp[offs];
	offs += dupl;
	fp_t **lz = *z = &ptmp[offs];
	offs += dupl;
	/* s */
	fp_t **lgamma = *gamma = &ptmp[offs];
	offs += dupl;
	fp_t **lsigma = *sigma = &ptmp[offs];
	offs += dupl;
	fp_t **ldelta = *delta = &ptmp[offs];

	unsigned long allocsize = 0;
	int d;
	for ( d=0; d<dupl; ++d ) {
		lp[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}
	for ( d=0; d<dupl; ++d ) {
		lv[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}
	for ( d=0; d<dupl; ++d ) {
		lr[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}
	for ( d=0; d<dupl; ++d ) {
		lz[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}

	for ( d=0; d<dupl; ++d ) {
		lgamma[d] = work;
		work += salgn;
		lsigma[d] = work;
		work += salgn;
		ldelta[d] = work;
		work += salgn;
		allocsize += salgn * 3; 
	}

	if ( allocsize != works ) {
		fprintf(stderr, "err: insufficient workspace (avail %i req %i)\n", works, allocsize);
		return 1;
	}

	return 0;
}


#endif // __CGMOD1_WORKSPACE_H__
