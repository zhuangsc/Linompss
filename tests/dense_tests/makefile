CC = mcc --ompss
CFLAGS=-Wall -fdiagnostics-color
LDFLAGS=-lm

#include paths
LO_INC = -I../../build/include
MKL_INC = -I${MKLROOT}/include

#link paths
LO_LIB = ../../build/lib

#actual link line
MKL_LINK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
GEMM_LINK = ${LDFLAGS} ${MKL_LINK} -L${LO_LIB} -lompssmm
SYRK_LINK = ${LDFLAGS} ${MKL_LINK} -L${LO_LIB} -lompsssyrk
TRSM_LINK = ${LDFLAGS} ${MKL_LINK} -L${LO_LIB} -lompsstrsm
LU_LINK	  = ${LDFLAGS} ${MKL_LINK} -L${LO_LIB} -lompsslu
POTRF_LINK = ${LDFLAGS} ${MKL_LINK} -L${LO_LIB} -lompssdchol -lompssutil -lbblas

############################################################################################

SOURCE = gemm_test.c syrk_test.c trsm_test.c lu_test.c potrf_test.c

BIN = dgemm_test sgemm_test \
	  dsyrk_test ssyrk_test \
	  dtrsm_test strsm_test \
	  dlu_test slu_test \
	  dpotrf_test spotrf_test

OBJ = dgemm_test.o sgemm_test.o \
	  dsyrk_test.o ssyrk_test.o \
	  dtrsm_test.o strsm_test.o \
	  dlu_test.o slu_test.o \
	  dpotrf_test.o spotrf_test.o

############################################################################################

all : $(BIN)

spotrf_test.o : potrf_test.c
	$(CC) $(CFLAGS) -DSINGLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dpotrf_test.o : potrf_test.c
	$(CC) $(CFLAGS) -DDOUBLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dpotrf_test : dpotrf_test.o
	$(CC) $(CFLAGS) -o $@ $< $(POTRF_LINK)

spotrf_test : spotrf_test.o
	$(CC) $(CFLAGS) -o $@ $< $(POTRF_LINK)


slu_test.o : lu_test.c
	$(CC) $(CFLAGS) -DSINGLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dlu_test.o : lu_test.c
	$(CC) $(CFLAGS) -DDOUBLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dlu_test : dlu_test.o
	$(CC) $(CFLAGS) -o $@ $< $(LU_LINK)

slu_test : slu_test.o
	$(CC) $(CFLAGS) -o $@ $< $(LU_LINK)


strsm_test.o : trsm_test.c
	$(CC) $(CFLAGS) -DSINGLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dtrsm_test.o : trsm_test.c
	$(CC) $(CFLAGS) -DDOUBLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dtrsm_test : dtrsm_test.o
	$(CC) $(CFLAGS) -o $@ $< $(TRSM_LINK)

strsm_test : strsm_test.o
	$(CC) $(CFLAGS) -o $@ $< $(TRSM_LINK)


ssyrk_test.o : syrk_test.c
	$(CC) $(CFLAGS) -DSINGLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dsyrk_test.o : syrk_test.c
	$(CC) $(CFLAGS) -DDOUBLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dsyrk_test : dsyrk_test.o
	$(CC) $(CFLAGS) -o $@ $< $(SYRK_LINK)

ssyrk_test : ssyrk_test.o
	$(CC) $(CFLAGS) -o $@ $< $(SYRK_LINK)


sgemm_test.o : gemm_test.c
	$(CC) $(CFLAGS) -DSINGLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dgemm_test.o : gemm_test.c
	$(CC) $(CFLAGS) -DDOUBLE_PRECISION -DUSE_MKL -c -o $@ $< $(LO_INC) $(MKL_INC)

dgemm_test : dgemm_test.o
	$(CC) $(CFLAGS) -o $@ $< $(GEMM_LINK)

sgemm_test : sgemm_test.o
	$(CC) $(CFLAGS) -o $@ $< $(GEMM_LINK)

.PHONY : clean

clean : 
	rm -f $(OBJ) $(BIN)
