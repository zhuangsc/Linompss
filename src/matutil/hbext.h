#ifndef __HBEXT_H__
#define __HBEXT_H__


#include "hb.h"


void hb_markudiag(hbmat_t *A);
void one2zero(hbmat_t* in_matrix);
int* etree(hbmat_t* in_matrix);
void hb_setdiag(hbmat_t *A);


#endif // __HBEXT_H__
