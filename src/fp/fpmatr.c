#include "fpmatr.h"


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif



#if USE_MKL

#include "mkl.h"

ompssblas_t OMPSSBLAS_ROWMAJOR 			= {.e = CblasRowMajor};
ompssblas_t OMPSSBLAS_COLMAJOR			= {.e = CblasColMajor};
ompssblas_t OMPSSBLAS_TRANSP			= {.e = CblasTrans};
ompssblas_t OMPSSBLAS_NTRANSP			= {.e = CblasNoTrans};
ompssblas_t OMPSSBLAS_UPPERTRIANG		= {.e = CblasUpper};
ompssblas_t OMPSSBLAS_LOWERTRIANG		= {.e = CblasLower};
ompssblas_t OMPSSBLAS_NTRIANG			= {.e = 0};
ompssblas_t OMPSSBLAS_DIAGUNIT			= {.e = CblasUnit};
ompssblas_t OMPSSBLAS_NDIAGUNIT			= {.e = CblasNonUnit};
ompssblas_t OMPSSBLAS_LEFT				= {.e = CblasLeft};
ompssblas_t OMPSSBLAS_RIGHT				= {.e = CblasRight};

ompsslapack_t OMPSSLAPACK_UPPERTRIANG		= {.e = 'U'};
ompsslapack_t OMPSSLAPACK_LOWERTRIANG		= {.e = 'L'};


#else

#include "blas.h"

ompssblas_t OMPSSBLAS_ROWMAJOR;			
ompssblas_t OMPSSBLAS_COLMAJOR;
ompssblas_t OMPSSBLAS_TRANSP 			= {.s = "T"};
ompssblas_t OMPSSBLAS_NTRANSP			= {.s = "N"};
ompssblas_t OMPSSBLAS_UPPERTRIANG		= {.s = "U"};
ompssblas_t OMPSSBLAS_LOWERTRIANG		= {.s = "L"};
ompssblas_t OMPSSBLAS_NTRIANG			= {.s = "x"};
ompssblas_t OMPSSBLAS_DIAGUNIT			= {.s = "U"};
ompssblas_t OMPSSBLAS_NDIAGUNIT			= {.s = "N"};
ompssblas_t OMPSSBLAS_LEFT				= {.s = "L"};
ompssblas_t OMPSSBLAS_RIGHT				= {.s = "R"};

ompsslapack_t OMPSSLAPACK_UPPERTRIANG		= {.s = "U"};
ompsslapack_t OMPSSLAPACK_LOWERTRIANG		= {.s = "L"};

#endif
