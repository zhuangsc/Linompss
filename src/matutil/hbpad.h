#ifndef __HBPAD_H__


#include "hb.h"


#ifdef SINGLE_PRECISION

#define csc_pad			csc_pad_single

#else

#define csc_pad			csc_pad_double

#endif


hbmat_t* csc_pad_single(hbmat_t *A, int mr, int offs);
hbmat_t* csc_pad_double(hbmat_t *A, int mr, int offs);


#endif // __HBPAD_H__
