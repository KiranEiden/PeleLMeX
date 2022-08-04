#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data,
                     const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>* eosparm)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    /*AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];); */
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + 0.5 * Lx;,
                 const amrex::Real yc = prob_lo[1] + 0.5 * Ly;,
                 const amrex::Real zc = prob_lo[2] + 0.5 * Lz;);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos(eosparm);
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert = 0.0;

    amrex::Real radius = std::sqrt(AMREX_D_TERM(  (x-xc)*(x-xc),
                                                + (y-yc)*(y-yc),
                                                + (z-zc)*(z-zc)));

    if (prob_parm.pertmag > 0.0)
    {

#if ( AMREX_SPACEDIM == 2 )
       amrex::Real angle = std::atan2((y-yc),(x-xc));
       pert = prob_parm.pertmag *
             (0.7 * std::sin(2*angle) +
              //1.023 * std::sin(2 * (angle - 0.004598) / radius2) +
              //0.945 * std::sin(3 * (angle - 0.00712435) / radius2) +
              0.9 * std::sin(8 * angle ) +
              1.4 * std::sin(4 * angle  ));
    }

    amrex::Real y1 = (prob_parm.flame_radius - (radius + prob_parm.standoff - 0.5*dx[1] + pert))*100;
    amrex::Real y2 = (prob_parm.flame_radius - (radius + prob_parm.standoff + 0.5*dx[1] + pert))*100;

#elif ( AMREX_SPACEDIM == 3 )

       amrex::Real angle = std::atan2((y-yc),(x-xc));
       amrex::Real angle2 = std::atan2((z-zc),(y-yc));
       pert = prob_parm.pertmag *
             (0.5 * std::sin(2*angle) +
              0.7 * std::sin(8 * angle ) +
              1.1 * std::sin(4 * angle  )) *
             (1.04 * std::sin(2*angle2) +
              0.3 * std::sin(8 * angle2 ) +
              0.85 * std::sin(16 * angle2  ));
    }

amrex::Real y1 = (prob_parm.flame_radius - (radius + prob_parm.standoff - 1.42*dx[1] + pert))*100;
amrex::Real y2 = (prob_parm.flame_radius - (radius + prob_parm.standoff + 1.42*dx[1] + pert))*100;
#endif

    amrex::Real ymean = 0.5 * (y1+y2);
    pele::physics::PMF::pmf(pmf_data,ymean, ymean, pmf_vals);

    state(i,j,k,TEMP) = pmf_vals[0];;

    state(i,j,k,VELX) = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    state(i,j,k,VELY) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 )
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = pmf_vals[1]*1e-2;
#endif

#ifdef USE_MANIFOLD_EOS
    for (int n = 0; n < NUM_SPECIES-1; n++){
      massfrac[n] = pmf_vals[3 + n];
    }
    amrex::Real rho_temp, tmp1, tmp2;
    eos.PYT2R(tmp1, massfrac, tmp2, rho_temp);
    rho_temp *= 1.0e3;            // CGS -> MKS conversion
    state(i,j,k,DENSITY) = rho_temp;
    state(i,j,k,RHOH) = 0.0; // No RhoH transport in manifold model

    for (int n = 0; n < NUM_SPECIES-1; n++) {
      state(i,j,k,FIRSTSPEC+n) = massfrac[n] * rho_temp;
    }
    state(i,j,k,FIRSTSPEC+NUM_SPECIES-1) = rho_temp;
#else
    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);
    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
#endif

    // Read in turbulent fluctuations from MarcHDF
    if (prob_parm.hitIC) {
      amrex::Real xtf[3] = {x, y, z};
      amrex::Real u[3] = {0.0};
      amrex::Real mod[3] = {0.0};
      int idx[3] = {0};
      int idxp1[3] = {0};
      amrex::Real slp[3] = {0.0};
      for (int cnt = 0; cnt < 3; cnt++) {
        mod[cnt] = std::fmod(xtf[cnt], prob_parm.Linput);
        locate(prob_parm.d_xarray, prob_parm.input_resolution, mod[cnt], idx[cnt]);
        idxp1[cnt] = (idx[cnt] + 1) % prob_parm.input_resolution;
        slp[cnt] = (mod[cnt] - prob_parm.d_xarray[idx[cnt]]) / prob_parm.d_xdiff[idx[cnt]];
      }
      int inSize = prob_parm.input_resolution;

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
      const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
      const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
      const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
      const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
      const amrex::Real f7 = slp[0] * slp[1] * slp[2];
      u[0] =
        prob_parm.d_uinput
            [idx[0] + inSize * (idx[1] + inSize * idx[2])] *
          f0 +
        prob_parm.d_uinput
            [idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_uinput
            [idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_uinput
            [idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_uinput
            [idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_uinput
            [idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_uinput
            [idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_uinput
            [idxp1[0] +
             inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;
      u[1] =
        prob_parm.d_vinput
            [idx[0] + inSize * (idx[1] + inSize * idx[2])] *
          f0 +
        prob_parm.d_vinput
            [idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_vinput
            [idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_vinput
            [idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_vinput
            [idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_vinput
            [idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_vinput
            [idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_vinput
            [idxp1[0] +
             inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;
      u[2] =
        prob_parm.d_winput
            [idx[0] + inSize * (idx[1] + inSize * idx[2])] *
          f0 +
        prob_parm.d_winput
            [idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_winput
            [idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_winput
            [idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_winput
            [idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_winput
            [idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_winput
            [idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_winput
            [idxp1[0] +
             inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;
      AMREX_D_TERM(
      const amrex::Real decayx =
        (0.5 *
         (std::tanh((xtf[0] - (prob_lo[0] + prob_parm.buffer_boundary*Lx)) / (prob_parm.buffer_thickness*Lx)) + 1)) *
        (0.5 *
         (std::tanh(-(xtf[0] - (prob_hi[0] - prob_parm.buffer_boundary*Lx)) / (prob_parm.buffer_thickness*Lx)) + 1));,
      const amrex::Real decayy =
        (0.5 *
         (std::tanh((xtf[1] - (prob_lo[1] + prob_parm.buffer_boundary*Ly)) / (prob_parm.buffer_thickness*Ly)) + 1)) *
        (0.5 *
         (std::tanh(-(xtf[1] - (prob_hi[1] - prob_parm.buffer_boundary*Ly))  / (prob_parm.buffer_thickness*Ly)) + 1));,
      const amrex::Real decayz =
        (0.5 *
         (std::tanh((xtf[2] - (prob_lo[2] + prob_parm.buffer_boundary*Lz)) / (prob_parm.buffer_thickness*Lz)) + 1)) *
        (0.5 *
         (std::tanh(-(xtf[2] - (prob_hi[2] - prob_parm.buffer_boundary*Lz)) / (prob_parm.buffer_thickness*Lz)) + 1));)
      u[0] *= AMREX_D_TERM(decayx, * decayy, * decayz);
      u[1] *= AMREX_D_TERM(decayx, * decayy, * decayz);
      u[2] *= AMREX_D_TERM(decayx, * decayy, * decayz);
      AMREX_D_TERM(state(i,j,k,VELX) = u[0];,
                   state(i,j,k,VELY) = u[1];,
                   state(i,j,k,VELZ) = u[2];)
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>* eosparm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
