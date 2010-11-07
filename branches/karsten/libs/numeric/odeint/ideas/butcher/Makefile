CC = g++

RATIO_ROOT = $(HOME)/boost/chrono

CXXFLAGS = -O3 -I$(BOOST_ROOT) -I$(ODEINT_ROOT) -I$(RATIO_ROOT)

HEADERS = algebra.hpp convert_value.hpp diagnostic.hpp explicit_runge_kutta.hpp predefined_steppers.hpp

all : butcher_performance butcher_lorenz butcher_vdp

butcher_performance : butcher_performance.o
butcher_performance.o : butcher_performance.cpp $(HEADERS)

butcher_vdp : butcher_vdp.o
butcher_vdp.o : butcher_vdp.cpp $(HEADERS)

butcher_lorenz : butcher_lorenz.o
butcher_lorenz.o : butcher_lorenz.cpp $(HEADERS)


clean :
	-rm *.o *~ butcher_performance butcher_vdp butcher_lorenz