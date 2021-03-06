#ifndef __CGPIPE_WORKSPACE_H__
#define __CGPIPE_WORKSPACE_H__


#include <stdlib.h>

#include "fptype.h"


static inline __attribute__((always_inline)) int cgpipe_malloc(size_t n, int s, int dupl, fp_t *work, unsigned long works, fp_t **p, fp_t **w, fp_t **q, fp_t **z, fp_t **ss, fp_t ***r,\
	fp_t ***alpha, fp_t ***gamma, fp_t ***delta) {
	unsigned long ns = n * s;
	unsigned long nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	unsigned long salgn = ((s + 0x7f ) >> 7 ) << 7;

	fp_t **ptmp = malloc(4 * dupl * sizeof(fp_t*));
	if ( ptmp == NULL ) {
		return 1;
	}

	int offs = 0;
	/* ns */
	fp_t **lr = *r = &ptmp[offs];
	offs += dupl;

	/* s */
	fp_t **lalpha = *alpha = &ptmp[offs];
	offs += dupl;
	fp_t **lgamma = *gamma = &ptmp[offs];
	offs += dupl;
	fp_t **ldelta = *delta = &ptmp[offs];

	unsigned long allocsize = 0;
	int d;

	*p = work;
	work += nsalgn;
	allocsize += nsalgn; 

	*w = work;
	work += nsalgn;
	allocsize += nsalgn; 

	*q = work;
	work += nsalgn;
	allocsize += nsalgn; 

	*z = work;
	work += nsalgn;
	allocsize += nsalgn; 

	*ss = work;
	work += nsalgn;
	allocsize += nsalgn; 

	for ( d=0; d<dupl; ++d ) {
		lr[d] = work;
		work += nsalgn;
		allocsize += nsalgn; 
	}

	for ( d=0; d<dupl; ++d ) {
		lalpha[d] = work;
		work += salgn;
		lgamma[d] = work;
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


#endif // __CGPIPE_WORKSPACE_H__
