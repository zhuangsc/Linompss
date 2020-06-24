CC = gcc -Wall -g
LIB_PATH = ../../build/lib
INC_PATH = ../../build/include

MKL_LINK = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
MKL_INC = -I${MKLROOT}/include
LIN_LINK = -L${LIB_PATH} -lmatutil -lompssjacobi -lompssgs -lsbblas -lm
#LIN_INC = -I${INC_PATH}
LIN_INC = -I${INC_PATH} -I../../src/matutil -I../../src/util -I../../src/blocked_sblas -I../../src/blocked_blas \
		  -I../../src/fp -I../../src/blas -I../../src/jacobi -I../../src/gauss_seidel


JAC_SOURCE = jacobi_test.c 
JAC_DOBJ = djacobi_test.o
JAC_SOBJ = sjacobi_test.o
JAC_DBIN = djacobi_test
JAC_SBIN = sjacobi_test

GS_SOURCE = gs_test.c 
GS_DOBJ = dgs_test.o
GS_SOBJ = sgs_test.o
GS_DBIN = dgs_test
GS_SBIN = sgs_test

OBJ=$(JAC_DOBJ) $(JAC_SOBJ) $(GS_DOBJ) $(GS_SOBJ)
BIN=$(JAC_DBIN) $(JAC_SBIN) $(GS_DBIN) $(GS_SBIN)

#%.o : %.c 
#	$(CC) -DDOUBLE_PRECISION -c -o $@ $< $(LIN_INC)
#	$(CC) -DSINGLE_PRECISION -c -o $@ $< $(LIN_INC)

#$(BIN) : $(OBJ)
#	$(CC) -o $@ $< $(MKL_LINK) $(LIN_LINK)

.PHONY : all
all : $(BIN)

$(JAC_SOBJ) : $(JAC_SOURCE)
	$(CC) -DSINGLE_PRECISION -c -o $(JAC_SOBJ) $(JAC_SOURCE) $(LIN_INC)

$(JAC_DOBJ) : $(JAC_SOURCE)
	$(CC) -DDOUBLE_PRECISION -c -o $(JAC_DOBJ) $(JAC_SOURCE) $(LIN_INC)

$(JAC_SBIN) : $(JAC_SOBJ)
	$(CC) -o $(JAC_SBIN) $(JAC_SOBJ) $(MKL_LINK) $(LIN_LINK)

$(JAC_DBIN) : $(JAC_DOBJ)
	$(CC) -o $(JAC_DBIN) $(JAC_DOBJ) $(MKL_LINK) $(LIN_LINK)

### GS
$(GS_SOBJ) : $(GS_SOURCE)
	$(CC) -DSINGLE_PRECISION -c -o $(GS_SOBJ) $(GS_SOURCE) $(LIN_INC)

$(GS_DOBJ) : $(GS_SOURCE)
	$(CC) -DDOUBLE_PRECISION -c -o $(GS_DOBJ) $(GS_SOURCE) $(LIN_INC)

$(GS_SBIN) : $(GS_SOBJ)
	$(CC) -o $(GS_SBIN) $(GS_SOBJ) $(MKL_LINK) $(LIN_LINK)

$(GS_DBIN) : $(GS_DOBJ)
	$(CC) -o $(GS_DBIN) $(GS_DOBJ) $(MKL_LINK) $(LIN_LINK)


.PHONY : clean

clean : 
	rm -f $(OBJ) $(BIN) *.log
