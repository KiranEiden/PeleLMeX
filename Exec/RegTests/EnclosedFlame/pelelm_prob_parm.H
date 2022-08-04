#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real standoff = 0.0_rt;
   amrex::Real pertmag  = 0.0004_rt;
   amrex::Real flame_radius = 0.036_rt;
   int  meanFlowDir = 1;

   // stuff for HIT init
   int hitIC = 0;
   int input_resolution = 0;
   amrex::Real uin_norm = 1.0;
   amrex::Real lscale = 1.0;
   amrex::Real offset = 0.0;
   amrex::Real urms0 = 0.0;
   amrex::Real Linput = 0.0;
  amrex::Real buffer_boundary = 1.0/6.0;
  amrex::Real buffer_thickness = 1.0/24.0;
  /*
#if (AMREX_SPACEDIM == 2)
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_lo = {0.0, 0.0};
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_hi = {0.0, 0.0};
#elif (AMREX_SPACEDIM == 3)
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_lo = {0.0, 0.0, 0.0};
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_hi = {0.0, 0.0, 0.0};
#endif
   amrex::Real win_slope = 1.0;
  */
   amrex::Real* d_uinput = nullptr;
   amrex::Real* d_vinput = nullptr;
   amrex::Real* d_winput = nullptr;
   amrex::Real* d_xarray = nullptr;
   amrex::Real* d_xdiff  = nullptr;
};
#endif
