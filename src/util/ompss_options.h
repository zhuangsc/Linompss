#ifndef __OMPSS_OPTIONS_H__
#define __OMPSS_OPTIONS_H__


typedef enum {
	OMPSSOPT_NO, 
	OMPSSOPT_YES 
} ompssopt_t;


#define LIBCG_EXPORT __attribute__((__visibility__("default")))


extern const char *ompssopt_map[];

extern const char *ompssopt_NO;
extern const char *ompssopt_YES;


extern LIBCG_EXPORT int ompssopt_read(char *app, char *opt, ompssopt_t deflt); 


#endif // __OMPSS_OPTIONS_H__
