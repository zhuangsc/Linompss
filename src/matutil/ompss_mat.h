#ifndef __OMPSS_MAT_H__
#define __OMPSS_MAT_H__


typedef enum {
	GENMATCONF_SPD, 
	GENMATCONF_COVAR, 
	GENMATCONF_RAND, 
	GENMATCONF_EYE, 
	GENMATCONF_DIAGDOM, 
	GENMATCONF_ONES, 
	GENMATCONF_ZERO,  
	GENMATCONF_COSTAS 
} ompssmat_t;


#endif // __OMPSS_MAT_H__

