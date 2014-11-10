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

# change this if you want to cross-compile
ARCH = Host
# ARCH = AVX
# ARCH = SSE4.2

CXX = icpc
CC = icpc

INCLUDES += -I../../include/ -I$(BOOST_ROOT) -I${SIMD_INCLUDE}
CXXFLAGS = -Ofast -x${ARCH} -fno-alias -ip -inline-forceinline -std=c++0x -DNDEBUG ${INCLUDES}
