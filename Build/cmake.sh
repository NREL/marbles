#!/bin/bash
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=clang++ \
      -DCMAKE_C_COMPILER:STRING=clang \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DMARBLES_ENABLE_MPI:BOOL=OFF \
      -DMARBLES_ENABLE_CPPCHECK:BOOL=OFF \
      -DMARBLES_ENABLE_CLANG_TIDY:BOOL=OFF \
      .. 
nice cmake --build . --parallel $(sysctl -n hw.ncpu)
#ctest
