MARBLES_HOME = ../
AMREX_HOME ?= $(MARBLES_HOME)/Submodules/AMReX

EBASE = marbles

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

TINY_PROFILE = TRUE
USE_PARTICLES = FALSE

BL_NO_FORT = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(MARBLES_HOME)/Source/Make.package

INCLUDE_LOCATIONS += $(MARBLES_HOME)/Source
VPATH_LOCATIONS   += $(MARBLES_HOME)/Source

include $(AMREX_HOME)/Tools/GNUMake/Make.rules


