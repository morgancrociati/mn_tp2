# une bibliotheque statique (libmnblas.a)
# une bibliotheque dynamique (libmnblasdyn.so)

LIB_DIR=../lib
INC_DIR=../include

FUNC = complexe.o copy.o swap.o dot.o axpy.o asum.o iamax.o iamin.o nrm2.o gemv.o gemm.o blas_test.o

OPTIONS=-Wall -fPIC $(OPTIM) -I$(INC_DIR)

all: libmnblas.a libmnblasdyn.so

# Bibliothèque statique

libmnblas.a: $(FUNC)
	rm -f libmnblas.a $(LIB_DIR)/libmnblas.a
	ar -r libmnblas.a $(FUNC)
	cp libmnblas.a $(LIB_DIR)

# Bibliothèque dynamique

libmnblasdyn.so: $(FUNC)
	rm -f libmnblasdyn.so $(LIB_DIR)/libmnblasdyn.so
	gcc -shared -o libmnblasdyn.so $(FUNC)
	cp libmnblasdyn.so $(LIB_DIR)

# Complexe

complexe.o: complexe.c $(INC_DIR)/complexe.h
	gcc $(OPTIONS) -c $<

# BLAS 1

copy.o : copy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc $(OPTIONS) -c $<

swap.o : swap.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

dot.o : dot.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

axpy.o : axpy.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc $(OPTIONS) -c $<

asum.o : asum.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

iamax.o : iamax.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

iamin.o : iamin.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

nrm2.o : nrm2.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $< -lm

# BLAS 2

gemv.o : gemv.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

# BLAS 3

gemm.o : gemm.c $(INC_DIR)/mnblas.h $(INC_DIR)/complexe.h
	gcc  $(OPTIONS) -c $<

# Tests

blas_test.o : blas_test.c $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc  $(OPTIONS) -c $<

clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a 


