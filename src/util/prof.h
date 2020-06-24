#ifndef __PROF_H__
#define __PROF_H__


#include <math.h>


typedef struct {
	float en;
	int i;
} profel_t;


typedef struct proft {
	int s;
	int dupl;
	int profc;
	int writec;
	int read;
	int readc;
	int write;
	float resolution;
	float *p;
	int *accept;
	//float min;
	//int max;
	//int *prior;
} prof_t;


#include <float.h>
#include <stdlib.h>


int prof_compar(const void *a, const void *b);
void prof_prioritize(prof_t *prof);
void prof_print(prof_t *prof);


static inline __attribute__((always_inline)) int prof_malloc(prof_t *p, float resolution, int n, int dupl) 
{
	int nd = n * dupl;

	p->dupl = dupl;

	p->s = n;
	p->read = 0;
	p->write = 0;
	p->readc = 0;
	p->writec = 0;
	p->profc = 0;

	p->resolution = resolution;

	p->p = calloc((2 + nd), sizeof(float));
	p->p[0] = -1.0;
	p->p[1] = -1.0;
	p->p += 2;
	p->accept = calloc(n, sizeof(int));
}


static inline __attribute__((always_inline)) int prof_free(prof_t *p) 
{
	if ( p->s != 0 ) {
		
		float *lp = p->p - 2;
		free(lp);
		int *la = p->accept;
		free(la);
		//free(p->prior);
	}
}


static inline __attribute__((always_inline)) int prof_add(prof_t *p, int i, float en) 
{
	int profs = p->s;
	if ( profs ) {
		int w = p->write;
		int wc = p->writec + 1;

		float *prof = p->p + w * profs;
		prof[i] = en;

		int init = wc == 1;
		float *pM = &p->p[-w-1]; /* save new maximum */ 
		float M = *pM;
		float lM = en > M || init ? en : M;
		*pM = lM;

		//printf("%i %i i %i en %e M %e\n", w, wc, i, en, lM);

		int done = wc == p->s;
		p->writec = done ? 0 : wc;
		p->profc += done ? 1 : 0;
		w += done ? 1 : 0 ;
		p->write = w & 0x1;

#if 0
		if ( done && p->readc ) {
			printf("prof changed\n");
		}
#endif
	}
}


static inline __attribute__((always_inline)) int prof_priority(prof_t *p, int i) 
{
	int profs = p->s;
	if ( profs ) {
		int r = (p->profc + 1) & 1;
		//int rc = p->readc + 1;

		float *prof = p->p + r * profs;
		float M = p->p[-r-1];

		float scal = (1 / M * p->resolution);
		float len = fabs(prof[i]) * scal;
//		printf("prof[%d]: %.4f scal: %.4f len: %.4f\n", i, prof[i], scal, len);

		len = ( M < 0 ) ? 1 : len;

		//printf("%i %i %i using M %e scale %e: %e %f\n", r, rc, i, M, scal, len, prof[i]);

		//int done = rc == profs;
		//p->readc = done ? 0 : rc ;

		return (int) len;
	} else { 
		return 0;
	}
}


static inline __attribute__((always_inline)) void prof_dump(prof_t *p)
{
	int profs = p->s;
	float sum = 0.0;
	float sum_acc = 0.0;

	int r = (p->profc+1) & 1;
	float *prof = p->p + r * profs;
	float M = p->p[-r-1];
	int *acc = p->accept;
	int i;
	for( i = 0; i < profs; ++i ) {
		sum += fabs(prof[i]);
		if ( acc[i] )
			sum_acc+=fabs(prof[i]);
	}
	for ( i = 0; i < profs; ++i ) {
		if (acc[i])
			printf("%.4f  ", fabs(prof[i]/sum));
		else
			printf("%.4f(*)  ", fabs(prof[i]/sum));
	}
	printf("\n");
	printf("acc/tol: %.4f\n", fabs(sum_acc)/fabs(sum));

}

static inline __attribute__((always_inline)) void prof_dump_fake(prof_t *p, int fake)
{
	int profs = p->s;
	float sum = 0.0;
	float sum_acc = 0.0;

	int r = (p->profc+1) & 1;
	float *prof = p->p + r * profs;
	float M = p->p[-r-1];
	int *acc = p->accept;
	int i;
	for( i = 0; i < profs; ++i ) {
		sum += prof[i];
		if ( acc[i] )
			sum_acc+=prof[i];
	}
	for ( i = 0; i < profs; ++i ) {
		if (acc[i])
			printf("%.4f  ", prof[i]/sum);
		else
			printf("%.4f(*)  ", prof[i]/sum);
	}
	printf("\n");
	if ( !fake )
		printf("acc/tol: %.4f\n", sum_acc/sum);
	else
		printf("acc/tol: 1.0000\n");
}


static inline __attribute__((alwasys_inline)) void prof_reg(prof_t *p, int idx, int accept)
{
	int profs = p->s;
	int r = (p->profc+1) & 1;
	int *acc = p->accept;
	acc[idx] = accept;
}

#endif // __PROF_H__
