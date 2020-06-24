#ifndef __STRUTIL_H__
#define __STRUTIL_H__


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


extern LIBBBLAS_EXPORT int str2int(const char *str, const char *smap[]);


#endif // __STRUTIL_H__
