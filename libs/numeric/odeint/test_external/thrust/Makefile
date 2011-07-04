CUDA_ROOT = /usr/local/cuda

CC = g++
CXX = g++
NVCC = $(CUDA_ROOT)/bin/nvcc

INCLUDES += -I$(BOOST_ROOT) -I$(THRUST_ROOT) -I$(CUDA_ROOT)/include -I../../../../..

NVCCFLAGS = -O3 $(INCLUDES) --compiler-bindir=/usr/bin/g++

LDLIBS = -lcudart
LDFLAGS = -L$(CUDA_ROOT)/lib64

%.co : %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<


all : check_thrust


check_thrust : check_thrust.co
	$(CC) -o check_thrust $(LDFLAGS) $(LDLIBS) check_thrust.co
check_thrust.co : check_thrust.cu

clean :
	-rm *~ *.o *.co check_thrust

