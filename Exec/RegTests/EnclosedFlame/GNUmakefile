TOP = ../../../Submodules
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/..
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro

# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX

FUEGO_GAS       = FALSE
# PelePhysics
ifeq ($(FUEGO_GAS), TRUE)
	Chemistry_Model = posf10325
	Eos_Model	= Fuego
	Transport_Model = Simple
else
	USE_LIBTORCH    = FALSE
	LIBTORCH_HOME = /home/bperry/libtorch
	DEFINES += -D_GLIBCXX_USE_CXX11_ABI=0
	Eos_Model = Manifold
	Chemistry_Model = Manifold2d
	Transport_Model = Manifold
endif

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
