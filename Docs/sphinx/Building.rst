.. _Building:

Building
--------

`ltBoltz` can be built in two ways: GNU Make build system and CMake.

GNU Make
~~~~~~~~

Assuming one uses the submodules, e.g., ``git clone --recursive`` or ``git submodule init; git submodule update``, then building `ltBoltz` is as simple as ::

  cd Build
  make -j

One can edit the ``GNUMakefile`` to customize compile options.

CMake
~~~~~

To use the CMake option, one executes the CMake configure command in the build directory::

  cd Build
  cmake -DCMAKE_BUILD_TYPE:STRING=Release \
        -DLTBOLTZ_ENABLE_MPI:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DCMAKE_C_COMPILER:STRING=mpicc \
        .. && make
