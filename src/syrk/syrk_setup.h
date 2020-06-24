#ifndef __SYRK_SETUP_H__
#define __SYRK_SETUP_H__


int syrk_setup(int n, int k, int b, void **A, void **C, void **Cc); 
int syrk_shutdown(int n, int k, int b, void *A, void *C, void *Cc);


#endif // __SYRK_SETUP_H__
