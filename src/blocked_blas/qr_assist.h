#ifndef __ASSIST_QR_H__
#define __ASSIST_QR_H__


#include "fptype.h"


int sapply_td_BlockH_var1( int m, int n, int k,
        float * buff_U2, int ldim_U2,
        float * buff_S,  int ldim_S,
        float * buff_W,  int ldim_W,
        float * buff_C1, int ldim_C1,
        float * buff_C2, int ldim_C2 ); 


int dapply_td_BlockH_var1( int m, int n, int k,
        double * buff_U2, int ldim_U2,
        double * buff_S,  int ldim_S,
        double * buff_W,  int ldim_W,
        double * buff_C1, int ldim_C1,
        double * buff_C2, int ldim_C2 ); 


#endif // __ASSIST_QR_H__
