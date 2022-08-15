AMREX_HOME ?= ../../submods/amrex
LBM_HOME = ../

EBASE = ltBoltz

DEBUG	= FALSE

DIM	= 3

COMP    = gnu

USE_MPI   = TRUE 
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_EB    = TRUE

TINY_PROFILE = TRUE
USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/EB/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(LBM_HOME)/Source/Make.package

INCLUDE_LOCATIONS += $(LBM_HOME)/Source
VPATH_LOCATIONS   += $(LBM_HOME)/Source

BL_NO_FORT = TRUE
include $(AMREX_HOME)/Tools/GNUMake/Make.rules


