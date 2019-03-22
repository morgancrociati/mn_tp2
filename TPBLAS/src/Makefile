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


%.o: %.c $(INC_DIR)/complexe.h $(INC_DIR)/blas_test.h
	gcc  $(OPTIONS) -c $<

clean:
	rm -f *.o *.a *.so *~  $(INC_DIR)/*~ $(LIB_DIR)/libmnblasdyn.so  $(LIB_DIR)/libmnblas.a 


