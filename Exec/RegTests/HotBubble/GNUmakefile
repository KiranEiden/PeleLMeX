TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro


# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX

# PelePhysics
Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Constant

USE_SUNDIALS_PP = TRUE

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
