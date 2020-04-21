CC  = gcc
CXX = g++

ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib

ifeq ($(BLISLAB_USE_BLAS),true)
CFLAGS = -$(COMPILER_OPT_LEVEL) -I/usr/include/x86_64-linux-gnu -march=armv8-a+sve -fPIC -DUSE_BLAS
LDLIBS = $(LIBBLISLAB) -lpthread -lm -lblas                                 #sudo apt-get install libblas-dev (Netlib BLAS)
#LDLIBS = $(LIBBLISLAB) ${BLAS_DIR}/lib/libblis.a -lpthread -lm              #blis: https://github.com/flame/blis/
#LDLIBS = $(LIBBLISLAB) ${BLAS_DIR}/lib/libopenblas.a -lpthread -lm          #openblas: https://github.com/xianyi/OpenBLAS
#LDLIBS = $(LIBBLISLAB) ${BLAS_DIR}/lib/*.a -lpthread -lm                     #*/General for any static link file: .a
else
CFLAGS = -$(COMPILER_OPT_LEVEL) -march=armv8-a+sve -fPIC
LDLIBS = $(LIBBLISLAB) -lpthread -lm
endif

LIBBLISLAB = $(BLISLAB_DIR)/lib/libblislab.a
SHAREDLIBBLISLAB = $(BLISLAB_DIR)/lib/libblislab.so

LDFLAGS = -I$(BLISLAB_DIR)/include -I$(BLISLAB_DIR)/kernels -I${BLAS_DIR}/include

OS_NAME = $(shell uname -s)
ifeq ($(OS_NAME),Linux)
LDLIBS += -lrt
endif

