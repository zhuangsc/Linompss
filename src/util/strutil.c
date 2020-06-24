#include "strutil.h"


#include <string.h>


int str2int(const char *str, const char *smap[]) {
	int i;
	for ( i=0; i<8; ++i ) {
		if ( strcmp(str, smap[i]) == 0 ) {
			return i;
		}
	}

	return -1;
}
