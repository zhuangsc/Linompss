#ifndef __CGAS_WORKSPACE_H__
#define __CGAS_WORKSPACE_H__


#include <stdlib.h>

#include "fptype.h"


static inline __attribute__((always_inline)) int cgas_malloc(size_t n, int s, int dupl, fp_t *work, unsigned long works, fp_t ***p, fp_t ***tmp, fp_t ***r, fp_t ***alpha1, fp_t ***alpha2) 
{
	unsigned long ns = n * s;
	unsigned long nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	unsigned long salgn = ((s + 0x7f ) >> 7 ) << 7;

	fp_t **ptmp = malloc(5 * dupl * sizeof(fp_t*));
	if ( ptmp == NULL ) {
		return 1;
	}

	int offs = 0;
	/* ns */
	fp_t **lp = *p = &ptmp[offs];
	offs += dupl;
	fp_t **lr = *r = &ptmp[offs];
	offs += dupl;
	fp_t **ltmp = *tmp = &ptmp[offs];
	offs += dupl;
	/* s */
	fp_t **lalpha1 = *alpha1 = &ptmp[offs];
	offs += dupl;
	fp_t **lalpha2 = *alpha2 = &ptmp[offs];

	unsigned long allocsize = 0;
	int d;
	for ( d=0; d<dupl; ++d ) {
		lp[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}
	for ( d=0; d<dupl; ++d ) {
		lr[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}
	for ( d=0; d<dupl; ++d ) {
		ltmp[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}

	for ( d=0; d<dupl; ++d ) {
		lalpha1[d] = work;
		work += salgn;
		lalpha2[d] = work;
		work += salgn;
		allocsize += salgn<<1; 
	}

	if ( allocsize > works ) {
		fprintf(stderr, "err: insufficient workspace (avail %i req %i)\n", works, allocsize);
		return 1;
	}

	return 0;
}


#endif // __CGAS_WORKSPACE_H__
