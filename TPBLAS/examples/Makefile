PRGM = test_dot test_copy test_axpy

all: $(PRGM)

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -L$(LIB_DIR) $(LIBST)  
OPTIONS_LINK_DYN  = -L$(LIB_DIR) $(LIBDYN)  

#test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

#test_dot_dyn: test_dot.o flop.o
#	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN)

test_dot: test_dot.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_copy: test_copy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_axpy: test_axpy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_dot.o: test_dot.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc $(OPTIONS_COMPIL) -c $<

test_copy.o: test_copy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc $(OPTIONS_COMPIL) -c $<

test_axpy.o: test_axpy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc $(OPTIONS_COMPIL) -c $<

clean:
	rm -f *.o $(PRGM) *~
