#!/bin/bash

function clean(){
    rm -f CMakeCache.txt
    rm -f CTestCustom.cmake
    rm -f CTestTestfile.cmake
    rm -f DartConfiguration.tcl
    rm -f Makefile
    rm -f MarblesConfig.cmake
    rm -f marbles
    rm -f cmake_install.cmake
    rm -f compile_commands.json
    rm -rf CMakeFiles
    rm -rf Source
    rm -rf Submodules
    rm -rf Testing
    rm -rf Tests
    rm -rf marbles_obj
}

CLEAN='false'

while getopts :c flag
do
    case "${flag}" in
        c)
            CLEAN='true'
            ;;
        '?')
            echo "INVALID OPTION -- ${OPTARG}" >&2
            exit 1
            ;;
        ':')
            echo "MISSING ARGUMENT for option -- ${OPTARG}" >&2
            exit 1
            ;;
        *)
            echo "UNIMPLEMENTED OPTION -- ${flag}" >&2
            exit 1
            ;;
    esac
done

if ${CLEAN}; then
    echo "Cleaning build"
    clean
fi

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=clang++ \
      -DCMAKE_C_COMPILER:STRING=clang \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DMARBLES_DIM:STRING=3 \
      -DMARBLES_ENABLE_MPI:BOOL=OFF \
      -DMARBLES_ENABLE_CPPCHECK:BOOL=OFF \
      -DMARBLES_ENABLE_CLANG_TIDY:BOOL=OFF \
      ..
nice cmake --build . --parallel $(sysctl -n hw.ncpu)
#ctest
