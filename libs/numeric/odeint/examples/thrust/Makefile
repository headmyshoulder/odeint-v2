CUDA_ROOT = /usr/local/cuda
ARCH = sm_13

CC = g++
CXX = g++
NVCC = $(CUDA_ROOT)/bin/nvcc

INCLUDES += -I$(BOOST_ROOT) -I$(THRUST_ROOT) -I$(CUDA_ROOT)/include -I../../../../..

NVCCFLAGS = -O3 $(INCLUDES) -arch $(ARCH) --compiler-bindir=/usr/bin/g++

LDLIBS = -lcudart
LDFLAGS = -L$(CUDA_ROOT)/lib64

%.co : %.cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<


all : phase_osc_chain


phase_osc_chain : phase_osc_chain.co
	$(CC) -o phase_osc_chain $(LDFLAGS) $(LDLIBS) phase_osc_chain.co
phase_osc_chain.co : phase_osc_chain.cu

clean :
	-rm *~ *.o *.co phase_osc_chain