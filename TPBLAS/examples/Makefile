TEST_PRGM = test_dot test_copy test_axpy test_swap

PERF_PRGM = perf_copy perf_dot perf_axpy perf_gemv perf_gemm

all: $(TEST_PRGM) $(PERF_PRGM)

INC_DIR=../include

LIB_DIR=../src

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -fPIC -I$(INC_DIR)

OPTIONS_LINK_STATIC  =   -L$(LIB_DIR) $(LIBST)  
OPTIONS_LINK_DYN  = -L$(LIB_DIR) $(LIBDYN)  

#test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
#	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

#test_dot_dyn: test_dot.o flop.o
#	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN)

#Test Programs
test_dot: test_dot.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_copy: test_copy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_axpy: test_axpy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

test_swap: test_swap.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

#Perf programs

perf_copy: perf_copy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

perf_dot: perf_dot.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

perf_axpy: perf_axpy.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

perf_gemv: perf_gemv.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

perf_gemm: perf_gemm.o $(LIB_DIR)/libmnblas.a
	gcc $(OPTIONS_COMPIL) -o $@ $< $(OPTIONS_LINK_STATIC)

%.o: %.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc $(OPTIONS_COMPIL) -c $<

clean:
	rm -f *.o $(TEST_PRGM) $(PERF_PRGM) *~
