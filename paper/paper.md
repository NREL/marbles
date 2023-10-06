---
title: 'MARBLES: Multi-scale Adaptively Refined Boltzmann LatticE Solver'
tags:
  - C++
  - adaptive mesh refinement
  - CFD
  - lattice boltzmann method
authors:
  - firstname: Marc
    surname: Henry de Frahan
    orcid: 0000-0001-7742-1565
    affiliation: 1
  - firstname: Ethan
    surname: Young
    orcid: 0000-0003-1106-7406
    affiliation: 1
  - firstname: Hariswaran
    surname: Sitaraman
    orcid: 0000-0001-5304-1664
    affiliation: 1
  - firstname: Ross
    surname: Larsen
    orcid: 0000-0002-2928-9835
    affiliation: 1
  - firstname: Jon
    surname: Rood
    orcid: 0000-0002-7513-3225
    affiliation: 1
affiliations:
  - name: High Performance Algorithms and Complex Fluids, National Renewable Energy Laboratory, USA
    index: 1
date: 5 October 2023
bibliography: paper.bib
---

# Summary

MARBLES is a computational fluid dynamics solver that leverages the Lattice Boltzmann method and block-structured adaptive mesh refinement (AMR) to simulate flows in complex media. The solver leverages the AMReX [@AMReX] library, a library that provides underlying data structures and programming models to run on massively parallel computing architectures.

MARBLES (insert technical detail)

MARBLES uses the Embedded Boundary (EB) formulation in AMReX to represent complex geometries. In this approach, an arbitrary surface is defined by the user, using either compositions of simple shapes or an STL file. This is used to intersect the Cartesian mesh and defined cells that are covered by the geometry (i.e., inside the body). Bounce back conditions are imposed on these surfaces to capture the effect of the geometry on the fluid flow. This leads to a very robust handling of very complex geometry, including the flow through porous media.

MARBLES is written in C++ and is built upon the AMReX library. This library implements the parallel paradigms, grid structures, domain decomposition, and portability framework. MARBLES uses a MPI+X approach: the AMR grid patches are distributed on different CPU ranks using MPI. Each grid can be either (i) logically tiled and solved on different threads on multi-core CPU machine using OpenMP, or (ii) solved on GPU threads on GPU nodes using CUDA, HIP, and SYCL.

# Statement of Need

Several software tools that implement the Lattice Boltzmann method can found online, including a PyTorch based solver, Lettuce [@bedrunka2021lettuce], the commonly used OpenLB [@krause2021openlb], Palabos [@latt2021palabos], waLBerla [@bauer2021walberla; @godenschwager2013framework], and others [@schmieschek2017lb3d; @mora2020concise; @pastewka2019hpc].

In contrast with these solvers, MARBLES is

Its unique features consist in combining

To the authors' knowledge, MARBLES is the only AMReX-based Lattice Boltzmann code that leverages the modern programming models implemented in that library. Because of this, MARBLES is easily extensible to other capabilities and will automatically perform well on emerging architectures.

MARBLES, through its use of vendor-agnostic programming models implemented in AMReX, achieves high performances from a small desktop station to the world's largest supercomputer.

MARBLES is intended for students, researchers and engineers interested in simulating mesoscale and macroscale flows with the Lattice Boltzmann method on modern high performance computing hardware. With performance portability and being vendor agnostic foremost considerations in the design of the code, MARBLES can leverage the computational resources available on the latest heterogeneous exascale platforms, include graphics processing units. To achieve (some statement of the current world's energy needs). In this context, MARBLES is a valuable tool to study (application examples)

# Acknowledgments

This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.

# References
