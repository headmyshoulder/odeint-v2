# Copyright 2014 Mario Mulansky
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or
# copy at http://www.boost.org/LICENSE_1_0.txt)

# make sure BOOST_ROOT is pointing to your boost directory
# otherwise, set it here:
# BOOST_ROOT = /path/to/boost
# you also need NT2s SIMD libary available set the include path here:
# SIMD_INCLUDE = /path/to/simd/include

INCLUDES += -I../../include/ -I$(BOOST_ROOT) -I${SIMD_INCLUDE}
#CXXFLAGS = -O3 -ffast-math -march=core-avx-i -std=c++0x -DNDEBUG -I${ODEINT_ROOT} -I${BOOST_ROOT} -I${SIMD_INCLUDE}
# -mtune=native -march=native

CXXFLAGS = -Ofast -xHost -fno-alias -ip -inline-forceinline -std=c++0x -DNDEBUG ${INCLUDES}
CXX = icpc
CC = icpc
