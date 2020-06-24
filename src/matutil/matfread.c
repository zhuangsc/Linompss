#include "matfread.h"

#include "fptype.h"


#if DOUBLE_PRECISION
#define __read_mm2dense		read_mm2ddense
#else 
#define __read_mm2dense		read_mm2sdense
#endif


int __read_mm2dense(FILE *f, int m, int n, fp_t *A) 
{
	char buf[1024];
	fp_t el;

	int i = 0;
	while ( fgets(buf, sizeof(buf), f) != NULL && i < m) {
    	if ( buf[0] != '#' ) {
			sscanf(buf, FP_SCANSPEC, &el);
			*A++ = el;
			++i;
		}
	}

	return 0;
}
