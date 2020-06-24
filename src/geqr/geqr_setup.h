#ifndef __QRCA_SETUP__
#define __QRCA_SETUP__


int geqr_setup(int check, int m, int mr, int n, int nr, int tr, int tc, int bs, int mt, int nt, int mleft, int nleft, double ***Ah, double ***Th, double ***Sh, double **Aorig);
void geqr_shutdown();


#endif // __QRCA_SETUP__
