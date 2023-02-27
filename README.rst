MARBLES: Multi-scale Adaptively Refined Boltzmann LatticE Solvers
-----------------------------------------------------------------

|Documentation Badge| |License Badge| |AMReX Badge| |C++ Badge|

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
   :target: https://amrex-codes.github.io/amrex/

.. |C++ Badge| image:: https://img.shields.io/badge/language-C%2B%2B17-blue
   :target: https://isocpp.org/

.. |License Badge| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. |Documentation Badge| image:: https://img.shields.io/badge/Docs-github.io-blue.svg
   :target: https://marbles.github.io


Getting Started
~~~~~~~~~~~~~~~

* To compile and run `MARBLES`, one needs a C++ compiler that supports the C++17 standard:

1. Have `MARBLES` use the default submodules for AMReX ::

    git clone --recursive git@github.com:NREL/marbles.git
    cd Build
    make realclean && make -j
    ./marbles.xxx.yyy.ex example.inp

* Notes

   A. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least version 7, and MPI should be at least of standard version 3).
   B. The example file can be any file from the `Tests/test_files` directories.
   C. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with CCSE's Amrvis (<https://ccse.lbl.gov/Downloads/downloadAmrvis.html>) or Vis-It (<http://vis.lbl.gov/NERSC/Software/visit/>):

      1. In VisIt, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder..
      2. With Amrvis, "amrvis3d plt00030", for example.


Dependencies
~~~~~~~~~~~~

`MARBLES` is built on the `AMReX` library.


Documentation
~~~~~~~~~~~~~

The full documentation for `MARBLES` exists in the Docs directory; at present this is maintained inline using Sphinx  `Sphinx <http://www.sphinx-doc.org>`_.

    cd Docs && mkdir build && cd build && sphinx-build -M html ../sphinx .


Acknowledgment
~~~~~~~~~~~~~~

This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.