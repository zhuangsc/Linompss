#ifndef __CHOLS_LLMAIN_H__
#define __CHOLS_LLMAIN_H__


#include "hb.h"
#include "hbconvrt.h"
#include "matfprint.h"

#ifdef SINGLE_PRECISION

#define chols_ll			schols_ll
#define chols_ll_upper		schols_ll_upper

#endif

#ifdef DOUBLE_PRECISION

#define chols_ll			dchols_ll
#define chols_ll_upper		dchols_ll_upper

#endif


int schols_ll(hbmat_t* Acsc, hbmat_t *Acsr, int* work);
int schols_ll_upper(hbmat_t* Acsr, hbmat_t *Acsc, int* work);
int dchols_ll(hbmat_t* Acsc, hbmat_t *Acsr, int* work);
int dchols_ll_upper(hbmat_t* Acsr, hbmat_t *Acsc, int* work);



#endif // __CHOLS_LLMAIN_H__
