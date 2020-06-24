#ifndef __OMPSS_APPS_H__
#define __OMPSS_APPS_H__


typedef enum {
	OMPSSAPP_CG, 
	OMPSSAPP_CGMOD1, 
	OMPSSAPP_CGPROF, 
	OMPSSAPP_POTRF, 
	OMPSSAPP_LU,
	OMPSSAPP_QR
} ompssapp_t;


#define LIBCG_EXPORT __attribute__((__visibility__("default")))


extern const char *ompssapp_map[];

extern const char *ompssapp_CG;
extern const char *ompssapp_CGMOD1;
extern const char *ompssapp_CGPROF;
extern const char *ompssapp_POTRF;
extern const char *ompssapp_LU;
extern const char *ompssapp_QR;
extern const char *ompssapp_ITREF;


extern LIBCG_EXPORT int ompssapp_readid(char *app, char *id, ompssapp_t deflt); 


#endif // __OMPSS_APPS_H__
