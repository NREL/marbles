MARBLES: Multi-scale Adaptively Refined Boltzmann LatticE Solver
----------------------------------------------------------------

|CI Badge| |Documentation Badge| |License Badge| |AMReX Badge| |C++ Badge|

.. |CI Badge| image:: https://github.com/NREL/marbles/workflows/MARBLES-CI/badge.svg
   :target: https://github.com/NREL/marbles/actions

.. |Documentation Badge| image:: https://github.com/NREL/marbles/workflows/MARBLES-Docs/badge.svg
   :target: https://https://nrel.github.io/marbles

.. |License Badge| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
   :target: https://amrex-codes.github.io/amrex/

.. |C++ Badge| image:: https://img.shields.io/badge/language-C%2B%2B17-blue
   :target: https://isocpp.org/



Getting Started
~~~~~~~~~~~~~~~

To compile and run `MARBLES`, one needs a C++ compiler that supports the C++17 standard, and then execute ::

    $ git clone --recursive git@github.com:NREL/marbles.git
    $ cd Build
    $ make realclean && make -j
    $ ./marbles3d.xxx.yyy.ex example.inp

.. note::
   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least version 7, and MPI should be at least of standard version 3).
   B. The example file can be any file from the `Tests/test_files` directories.
   C. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with AMReX's `Amrvis <https://amrex-codes.github.io/amrex/docs_html/Visualization.html>`_ or `VisIt <https://visit-dav.github.io/visit-website/>`_:

      1. In VisIt, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder.
      2. With Amrvis, `$ amrvis3d plt00030`, for example.

Alternatively, one can use cmake to build the code ::

    $ cd Build
    $ ./cmake.sh
    $ ./marbles example.inp

Dependencies
~~~~~~~~~~~~

`MARBLES` is built on the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library.


Documentation
~~~~~~~~~~~~~

The full documentation for `MARBLES` exists in the Docs directory; at present this is maintained inline using `Sphinx <http://www.sphinx-doc.org>`_. To build the documentation ::

    $ cd Docs && mkdir build && cd build && sphinx-build -M html ../sphinx .

Or, using cmake ::

    $ cd Build && cmake -B build-docs ../Docs && cmake --build build-docs
