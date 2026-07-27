#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <simd.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>

#include "aster_utils.hxx"
#include "eigenvalues.hxx"
#include "fluxes.hxx"
#include "reconstruct.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace EOSX;
using namespace ReconX;
using namespace AsterUtils;

enum class flux_t { LxF, HLLE };
enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class rec_var_t { v_vec, z_vec, s_vec };

/* Everything the per-face flux computation needs from its caller.
 *
 * CalcFlux builds one instance per call and the loop lambda captures it by
 * value -- the same values it captured individually before the per-face code
 * was extracted into a helper. Members carry the same names as the
 * corresponding locals of CalcFlux (and the DECLARE_CCTK_ARGUMENTSX /
 * DECLARE_CCTK_PARAMETERS variables), so the helper can unpack them into
 * identically-named references and the extracted body stays verbatim. */
template <typename EOSType> struct FluxContext {
  /* grid functions for fluxes */
  vec<GF3D2<CCTK_REAL>, dim> fluxdenss;
  vec<GF3D2<CCTK_REAL>, dim> fluxDEnts;
  vec<GF3D2<CCTK_REAL>, dim> fluxmomxs;
  vec<GF3D2<CCTK_REAL>, dim> fluxmomys;
  vec<GF3D2<CCTK_REAL>, dim> fluxmomzs;
  vec<GF3D2<CCTK_REAL>, dim> fluxtaus;
  vec<GF3D2<CCTK_REAL>, dim> fluxDYes;
  /* off-diagonal B-fluxes, indexed by the flux direction i: flux of B_{dir_j}
   * resp. B_{dir_k} in direction i (the diagonal flux of B_{dir_i} in
   * direction dir_i vanishes identically and has no grid function) */
  vec<GF3D2<CCTK_REAL>, dim> fluxB_j;
  vec<GF3D2<CCTK_REAL>, dim> fluxB_k;

  /* grid functions */
  vec<GF3D2<const CCTK_REAL>, dim> gf_vels;
  vec<GF3D2<const CCTK_REAL>, dim> gf_zvec;
  vec<GF3D2<const CCTK_REAL>, dim> gf_svec;
  vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs;
  vec<GF3D2<const CCTK_REAL>, dim> gf_dBstags;
  vec<GF3D2<const CCTK_REAL>, dim> gf_beta;
  smat<GF3D2<const CCTK_REAL>, dim> gf_g;

  /* grid functions for Upwind CT */
  vec<GF3D2<CCTK_REAL>, dim> vbar_j;
  vec<GF3D2<CCTK_REAL>, dim> vbar_k;
  vec<GF3D2<CCTK_REAL>, dim> ap_face;
  vec<GF3D2<CCTK_REAL>, dim> am_face;

  /* grid functions for PP flux limiter */
  vec<GF3D2<CCTK_REAL>, dim> gf_theta;

  /* cell-centred scalar grid functions (DECLARE_CCTK_ARGUMENTSX names) */
  GF3D2<const CCTK_REAL> alp;
  GF3D2<const CCTK_REAL> rho;
  GF3D2<const CCTK_REAL> press;
  GF3D2<const CCTK_REAL> eps;
  GF3D2<const CCTK_REAL> entropy;
  GF3D2<const CCTK_REAL> Ye;
  GF3D2<const CCTK_REAL> temperature;
  GF3D2<const CCTK_REAL> LOflag;
  GF3D2<const CCTK_REAL> dens;
  GF3D2<const CCTK_REAL> DYe;

  /* EOS object and solver settings (CalcFlux arguments) */
  EOSType *eos_3p;
  rec_var_t rec_var;
  reconstruction_t reconstruction;
  reconstruction_t reconstruction_LO;
  reconstruct_params_t reconstruct_params;
  flux_t fluxtype;

  /* flag for tabulated EOS */
  bool istab;

  /* velocity limit from Con2PrimFactory parameters */
  CCTK_REAL v_lim;

  /* parameters (DECLARE_CCTK_PARAMETERS names) */
  CCTK_REAL r_atmo;
  CCTK_REAL rho_abs_min;
  CCTK_REAL n_rho_atmo;
  CCTK_REAL recon_thresh;
  CCTK_REAL p_atmo;
  CCTK_REAL n_press_atmo;
  CCTK_REAL t_atmo;
  CCTK_REAL n_temp_atmo;
  CCTK_REAL Ye_atmo;
  CCTK_REAL atmo_tol;
  CCTK_INT use_press_atmo;
  CCTK_INT reconstruct_with_temperature;
  CCTK_INT use_pplim;
  CCTK_INT loworder_flux;

  /* flesh scalars (DECLARE_CCTK_ARGUMENTSX names; CCTK_DELTA_TIME expands
   * to cctk_delta_time / cctk_timefac) */
  CCTK_REAL cctk_delta_time;
  int cctk_timefac;
  int cctk_iteration; // used by the CCTK_DEBUG diagnostics only
};

/* Compute the fluxes in direction `dir_i` through the cell face at p.I
 * (the lower face in `dir_i` of cell p.I): interpolate the metric to the
 * face, reconstruct the primitives from the neighbouring cells (with
 * atmosphere handling), assemble the conserved variables and their fluxes,
 * solve the Riemann problem, apply the positivity-preserving flux limiter,
 * and store the upwindCT auxiliaries. Extracted verbatim from the main
 * loop of CalcFlux.
 *
 * `face_X` are the coordinates of the face; they are passed in explicitly
 * (instead of reading p.X) so that this function does not depend on the
 * centering of the calling loop. p is otherwise used for indexing only:
 * p.I, p.DI and p.DX are centering-independent. */
template <int dir_i, bool uct, bool pplim, typename EOSType>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
CalcFluxAtFace(const FluxContext<EOSType> &fx, const PointDesc &p,
               const vect<CCTK_REAL, dim> &face_X) {
  static_assert(dir_i >= 0 && dir_i < 3, "");

  /* Unpack the context into references named exactly as the corresponding
   * locals of CalcFlux, so that the code below stays a verbatim copy of
   * the loop body it was extracted from (see FluxContext). */
  const auto &fluxdenss = fx.fluxdenss;
  const auto &fluxDEnts = fx.fluxDEnts;
  const auto &fluxmomxs = fx.fluxmomxs;
  const auto &fluxmomys = fx.fluxmomys;
  const auto &fluxmomzs = fx.fluxmomzs;
  const auto &fluxtaus = fx.fluxtaus;
  const auto &fluxDYes = fx.fluxDYes;
  const auto &fluxB_j = fx.fluxB_j;
  const auto &fluxB_k = fx.fluxB_k;
  const auto &gf_vels = fx.gf_vels;
  const auto &gf_zvec = fx.gf_zvec;
  const auto &gf_svec = fx.gf_svec;
  const auto &gf_Bvecs = fx.gf_Bvecs;
  const auto &gf_dBstags = fx.gf_dBstags;
  const auto &gf_beta = fx.gf_beta;
  const auto &gf_g = fx.gf_g;
  const auto &vbar_j = fx.vbar_j;
  const auto &vbar_k = fx.vbar_k;
  const auto &ap_face = fx.ap_face;
  const auto &am_face = fx.am_face;
  const auto &gf_theta = fx.gf_theta;
  const auto &alp = fx.alp;
  const auto &rho = fx.rho;
  const auto &press = fx.press;
  const auto &eps = fx.eps;
  const auto &entropy = fx.entropy;
  const auto &Ye = fx.Ye;
  const auto &temperature = fx.temperature;
  const auto &LOflag = fx.LOflag;
  const auto &dens = fx.dens;
  const auto &DYe = fx.DYe;
  const auto &eos_3p = fx.eos_3p;
  const auto &rec_var = fx.rec_var;
  const auto &reconstruction = fx.reconstruction;
  const auto &reconstruction_LO = fx.reconstruction_LO;
  const auto &reconstruct_params = fx.reconstruct_params;
  const auto &fluxtype = fx.fluxtype;
  const auto &istab = fx.istab;
  const auto &v_lim = fx.v_lim;
  const auto &r_atmo = fx.r_atmo;
  const auto &rho_abs_min = fx.rho_abs_min;
  const auto &n_rho_atmo = fx.n_rho_atmo;
  const auto &recon_thresh = fx.recon_thresh;
  const auto &p_atmo = fx.p_atmo;
  const auto &n_press_atmo = fx.n_press_atmo;
  const auto &t_atmo = fx.t_atmo;
  const auto &n_temp_atmo = fx.n_temp_atmo;
  const auto &Ye_atmo = fx.Ye_atmo;
  const auto &atmo_tol = fx.atmo_tol;
  const auto &use_press_atmo = fx.use_press_atmo;
  const auto &reconstruct_with_temperature = fx.reconstruct_with_temperature;
  const auto &use_pplim = fx.use_pplim;
  const auto &loworder_flux = fx.loworder_flux;
  const auto &cctk_delta_time = fx.cctk_delta_time;
  const auto &cctk_timefac = fx.cctk_timefac;
  const auto &cctk_iteration = fx.cctk_iteration; // CCTK_DEBUG only

// Prebind the velocity slice for this direction once
const auto gf_vel_dir_i = gf_vels(dir_i);

const auto reconstruct_pt = [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var,
                                            const PointDesc &p,
                                            bool gf_is_rho,
                                            bool gf_is_press) {
  return reconstruct<vec<CCTK_REAL, 2> >(var, p, reconstruction, dir_i,
                                         gf_is_rho, gf_is_press, press,
                                         gf_vel_dir_i, reconstruct_params);
};
const auto reconstruct_loworder =
    [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p,
                    bool gf_is_rho, bool gf_is_press) {
      return reconstruct<vec<CCTK_REAL, 2> >(
          var, p, reconstruction_LO, dir_i, gf_is_rho, gf_is_press, press,
          gf_vel_dir_i, reconstruct_params);
    };

// PROBE (Step 2a -- eigenvalue collapse): the numerical-flux solvers need only
// the global wavespeed bounds charmax = max(0, all eigenvalues) and
// charmin = min(0, all eigenvalues), NOT the full (degenerate) lambda. Passing
// those two scalars is bit-identical to the old laxf/hlle (max/min/fabs are
// exact reductions) and lets lambda die right after it is computed instead of
// staying live to the UCT block.
const auto laxf_cc =
    [=] CCTK_DEVICE(CCTK_REAL charmax, CCTK_REAL charmin, vec<CCTK_REAL, 2> var,
                    vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      // c = max(0, |all eigenvalues|) = max(charmax, -charmin)
      const CCTK_REAL c = fmax(charmax, -charmin);
      return 0.5 * ((flux(0) + flux(1)) - c * (var(1) - var(0)));
    };
const auto calcflux =
    [=] CCTK_DEVICE(CCTK_REAL charmax, CCTK_REAL charmin, vec<CCTK_REAL, 2> var,
                    vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      CCTK_REAL flx;
      switch (fluxtype) {
      case flux_t::LxF: {
        flx = laxf_cc(charmax, charmin, var, flux);
        break;
      }
      case flux_t::HLLE: {
        // charpm = charmax - charmin
        flx = (charmax * flux(0) - charmin * flux(1) +
               charmax * charmin * (var(1) - var(0))) /
              (charmax - charmin);
        break;
      }
      default:
        assert(0);
      }
      return flx;
    };

constexpr int dir_j = (dir_i == 0) ? 1 : ((dir_i == 1) ? 2 : 0);
constexpr int dir_k = (dir_i == 0) ? 2 : ((dir_i == 1) ? 0 : 1);

  /* Reconstruct primitives from the cells on left (indice 0) and right
   * (indice 1) side of this face rc = reconstructed variables or
   * computed from reconstructed variables */

  /* Interpolate metric components from vertices to faces */
  const CCTK_REAL alp_avg = calc_avg_v2f<dir_i>(alp, p);
  const vec<CCTK_REAL, 3> betas_avg(
      [&](int i) ARITH_INLINE { return calc_avg_v2f<dir_i>(gf_beta(i), p); });
  const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
    return calc_avg_v2f<dir_i>(gf_g(i, j), p);
  });

  /* determinant of spatial metric */
  const CCTK_REAL detg_avg = calc_det(g_avg);
  const CCTK_REAL sqrtg = sqrt(detg_avg);

  // Booleans controlling reconstruction fallbacks
  bool useLO = false;
  bool resetL = false;
  bool resetR = false;

  // Reconstruct density
  auto rho_rc = reconstruct_pt(rho, p, true, true);

  // Reconstruct entropy
  auto entropy_rc = reconstruct_pt(entropy, p, false, false);

  // Reconstruct Ye
  auto Ye_rc = reconstruct_pt(Ye, p, false, false);

  // Initialize variables for eps, pressure, and temperature
  vec<CCTK_REAL, 2> eps_rc;
  vec<CCTK_REAL, 2> press_rc;
  vec<CCTK_REAL, 2> temp_rc;

  // Setting up atmosphere for two neighboring cell centers
  vec<CCTK_REAL, 2> r_atm;
  vec<CCTK_REAL, 2> r2_atm = {0.0, 0.0};
  vec<CCTK_REAL, 2> rho_atm;
  vec<CCTK_REAL, 2> rho_cut;
  vec<CCTK_REAL, 2> press_atm;
  vec<CCTK_REAL, 2> eps_atm;
  vec<CCTK_REAL, 2> temp_atm;

  // Get coordinates at neighboring cell centers
  for (int ii = 0; ii < 3; ii++) {
    r2_atm(0) += (face_X[ii] - (ii == dir_i) * 0.5 * (p.DX[dir_i])) *
                 (face_X[ii] - (ii == dir_i) * 0.5 * (p.DX[dir_i]));
    r2_atm(1) += (face_X[ii] + (ii == dir_i) * 0.5 * (p.DX[dir_i])) *
                 (face_X[ii] + (ii == dir_i) * 0.5 * (p.DX[dir_i]));
  }
  r_atm(0) = sqrt(r2_atm(0));
  r_atm(1) = sqrt(r2_atm(1));

  // Grading rho
  rho_atm(0) = (r_atm(0) > r_atmo)
                   ? (rho_abs_min * pow((r_atmo / r_atm(0)), n_rho_atmo))
                   : rho_abs_min;
  rho_atm(0) = std::max(eos_3p->rgrho.min, rho_atm(0));
  rho_cut(0) = rho_atm(0) * recon_thresh;

  rho_atm(1) = (r_atm(1) > r_atmo)
                   ? (rho_abs_min * pow((r_atmo / r_atm(1)), n_rho_atmo))
                   : rho_abs_min;
  rho_atm(1) = std::max(eos_3p->rgrho.min, rho_atm(1));
  rho_cut(1) = rho_atm(1) * recon_thresh;

  // Grading temperature or pressure
  if (use_press_atmo) {
    press_atm(0) = (r_atm(0) > r_atmo)
                       ? (p_atmo * pow(r_atmo / r_atm(0), n_press_atmo))
                       : p_atmo;
    press_atm(0) = std::max(eos_3p->press_from_rho_temp_ye(
                                rho_atm(0), eos_3p->rgtemp.min, Ye_atmo),
                            press_atm(0));
    press_atm(1) = (r_atm(1) > r_atmo)
                       ? (p_atmo * pow(r_atmo / r_atm(1), n_press_atmo))
                       : p_atmo;
    press_atm(1) = std::max(eos_3p->press_from_rho_temp_ye(
                                rho_atm(1), eos_3p->rgtemp.min, Ye_atmo),
                            press_atm(1));
    eps_atm(0) =
        eos_3p->eps_from_rho_press_ye(rho_atm(0), press_atm(0), Ye_atmo);
    eps_atm(1) =
        eos_3p->eps_from_rho_press_ye(rho_atm(1), press_atm(1), Ye_atmo);
    temp_atm(0) =
        eos_3p->temp_from_rho_eps_ye(rho_atm(0), eps_atm(0), Ye_atmo);
    temp_atm(1) =
        eos_3p->temp_from_rho_eps_ye(rho_atm(1), eps_atm(1), Ye_atmo);
  } else {
    temp_atm(0) = (r_atm(0) > r_atmo)
                      ? (t_atmo * pow(r_atmo / r_atm(0), n_temp_atmo))
                      : t_atmo;
    temp_atm(0) = std::max(eos_3p->rgtemp.min, temp_atm(0));

    temp_atm(1) = (r_atm(1) > r_atmo)
                      ? (t_atmo * pow(r_atmo / r_atm(1), n_temp_atmo))
                      : t_atmo;
    temp_atm(1) = std::max(eos_3p->rgtemp.min, temp_atm(1));
    press_atm(0) =
        eos_3p->press_from_rho_temp_ye(rho_atm(0), temp_atm(0), Ye_atmo);
    press_atm(1) =
        eos_3p->press_from_rho_temp_ye(rho_atm(1), temp_atm(1), Ye_atmo);
    eps_atm(0) =
        eos_3p->eps_from_rho_temp_ye(rho_atm(0), temp_atm(0), Ye_atmo);
    eps_atm(1) =
        eos_3p->eps_from_rho_temp_ye(rho_atm(1), temp_atm(1), Ye_atmo);
  }
  // End atmosphere

  // Check shock detection flag
  if (LOflag(p.I) || LOflag(p.I - p.DI[dir_i]))
    useLO = true;

  if (reconstruct_with_temperature) {

    // Reconstruct temperature
    temp_rc = reconstruct_pt(temperature, p, false, false);

    // Use lower-order if reconstructed rho, entropy, Ye or T is <= 0
    if ((rho_rc(0) <= rho_cut(0)) || (entropy_rc(0) <= 0.0) ||
        (Ye_rc(0) <= 0.0) || (temp_rc(0) <= 0.0) ||
        (rho_rc(1) <= rho_cut(1)) || (entropy_rc(1) <= 0.0) ||
        (Ye_rc(1) <= 0.0) || (temp_rc(1) <= 0.0) || useLO) {

      useLO = true;

      rho_rc = reconstruct_loworder(rho, p, true, true);
      entropy_rc = reconstruct_loworder(entropy, p, false, false);
      Ye_rc = reconstruct_loworder(Ye, p, false, false);
      temp_rc = reconstruct_loworder(temperature, p, false, false);
    }

    // If reconstructed rho is still <= atmo, flag for reset
    if (rho_rc(0) <= rho_cut(0)) {
      resetL = true;
      rho_rc(0) = rho_atm(0);
      entropy_rc(0) =
          eos_3p->kappa_from_rho_eps_ye(rho_atm(0), eps_atm(0), Ye_atmo);
      temp_rc(0) = temp_atm(0);
      Ye_rc(0) = Ye_atmo;
    }
    if (rho_rc(1) <= rho_cut(1)) {
      resetR = true;
      rho_rc(1) = rho_atm(1);
      entropy_rc(1) =
          eos_3p->kappa_from_rho_eps_ye(rho_atm(1), eps_atm(1), Ye_atmo);
      temp_rc(1) = temp_atm(1);
      Ye_rc(1) = Ye_atmo;
    }
    // End lower-order

    // Compute eps_rc and press_rc using lambdas
    for (int f = 0; f < 2; ++f) {
      eps_rc(f) =
          eos_3p->eps_from_rho_temp_ye(rho_rc(f), temp_rc(f), Ye_rc(f));
      press_rc(f) =
          eos_3p->press_from_rho_temp_ye(rho_rc(f), temp_rc(f), Ye_rc(f));
    }

  } else {

    // Reconstruct pressure
    press_rc = reconstruct_pt(press, p, false, true);

    // Use lower-order if reconstructed rho, entropy, Ye or press is <= 0
    if ((rho_rc(0) <= rho_cut(0)) || (entropy_rc(0) <= 0.0) ||
        (Ye_rc(0) <= 0.0) || (press_rc(0) <= 0.0) ||
        (rho_rc(1) <= rho_cut(1)) || (entropy_rc(1) <= 0.0) ||
        (Ye_rc(1) <= 0.0) || (press_rc(1) <= 0.0) || useLO) {

      useLO = true;

      rho_rc = reconstruct_loworder(rho, p, true, true);
      entropy_rc = reconstruct_loworder(entropy, p, false, false);
      Ye_rc = reconstruct_loworder(Ye, p, false, false);
      press_rc = reconstruct_loworder(press, p, false, true);
    }

    // If reconstructed rho is still <= atmo, flag for reset
    if (rho_rc(0) <= rho_cut(0)) {
      resetL = true;
      rho_rc(0) = rho_atm(0);
      entropy_rc(0) =
          eos_3p->kappa_from_rho_eps_ye(rho_atm(0), eps_atm(0), Ye_atmo);
      press_rc(0) = press_atm(0);
      Ye_rc(0) = Ye_atmo;
    }
    if (rho_rc(1) <= rho_cut(1)) {
      resetR = true;
      rho_rc(1) = rho_atm(1);
      entropy_rc(1) =
          eos_3p->kappa_from_rho_eps_ye(rho_atm(1), eps_atm(1), Ye_atmo);
      press_rc(1) = press_atm(1);
      Ye_rc(1) = Ye_atmo;
    }
    // End lower-order

    // Compute eps_rc and temp_rc using lambdas
    for (int f = 0; f < 2; ++f) {
      eps_rc(f) =
          eos_3p->eps_from_rho_press_ye(rho_rc(f), press_rc(f), Ye_rc(f));
      temp_rc(f) =
          eos_3p->temp_from_rho_eps_ye(rho_rc(f), eps_rc(f), Ye_rc(f));
    }
  }

  const vec<CCTK_REAL, 2> rhoh_rc([&](int f) ARITH_INLINE {
    return rho_rc(f) + rho_rc(f) * eps_rc(f) + press_rc(f);
  });

  // Introduce reconstructed Bs
  // Use staggered dB for i == dir_i
  vec<vec<CCTK_REAL, 2>, 3> Bs_rc;

  // Assign the value for the primary direction
  const CCTK_REAL val = gf_dBstags(dir_i)(p.I) / sqrtg;
  Bs_rc(dir_i)(0) = val;
  Bs_rc(dir_i)(1) = val;

  // Lambda to assign the reconstructed values
  auto assign_reconstructed = [&](int d) {
    auto tmp = useLO ? reconstruct_loworder(gf_Bvecs(d), p, false, false)
                     : reconstruct_pt(gf_Bvecs(d), p, false, false);
    Bs_rc(d)(0) = tmp(0);
    Bs_rc(d)(1) = tmp(1);
  };

  // Assign reconstructed values for the two perpendicular directions
  assign_reconstructed(dir_j);
  assign_reconstructed(dir_k);
  // End of setting Bs

  vec<vec<CCTK_REAL, 2>, 3> vels_rc;
  vec<vec<CCTK_REAL, 2>, 3> vlows_rc;
  vec<CCTK_REAL, 2> w_lorentz_rc;
  switch (rec_var) {
  case rec_var_t::v_vec: {

    if (useLO) {

      for (int i = 0; i <= 2; ++i) { // loop over components
        vels_rc(i) = reconstruct_loworder(gf_vels(i), p, false, false);
      }
    } else {

      for (int i = 0; i <= 2; ++i) { // loop over components
        vels_rc(i) = reconstruct_pt(gf_vels(i), p, false, false);
      }
    }

    /* co-velocity measured by Eulerian observer: v_j */
    vlows_rc = calc_contraction(g_avg, vels_rc);
    auto v2_rc = calc_contraction(vlows_rc, vels_rc);

    /* Lower-order if above limit */
    const CCTK_REAL v2_lim = v_lim * v_lim;
    if (!(useLO) && (v2_rc(0) >= v2_lim || v2_rc(1) >= v2_lim)) {

      for (int i = 0; i <= 2; ++i) { // loop over components
        vels_rc(i) = reconstruct_loworder(gf_vels(i), p, false, false);
      }

      vlows_rc = calc_contraction(g_avg, vels_rc);
      v2_rc = calc_contraction(vlows_rc, vels_rc);
    }

    /* Last resort if lower-order is also above limit:
     * Rescale */
    if (v2_rc(0) >= v2_lim) {
      CCTK_REAL f = v_lim / sqrt(v2_rc(0));
      for (int i = 0; i <= 2; ++i) {
        vels_rc(i)(0) *= f;
        vlows_rc(i)(0) *= f;
      }
      v2_rc(0) = v2_lim;
    }

    if (v2_rc(1) >= v2_lim) {
      CCTK_REAL f = v_lim / sqrt(v2_rc(1));
      for (int i = 0; i <= 2; ++i) {
        vels_rc(i)(1) *= f;
        vlows_rc(i)(1) *= f;
      }
      v2_rc(1) = v2_lim;
    }

    /* Lorentz factor: W = 1 / sqrt(1 - v^2) */
    w_lorentz_rc(0) = 1 / sqrt(1 - v2_rc(0));
    w_lorentz_rc(1) = 1 / sqrt(1 - v2_rc(1));
    break;
  };
  case rec_var_t::z_vec: {

    vec<vec<CCTK_REAL, 2>, 3> zvec_rc([&](int i) ARITH_INLINE {
      return reconstruct_pt(gf_zvec(i), p, false, false);
    });

    // Lower-order
    if (useLO) {

      for (int i = 0; i <= 2; ++i) { // loop over components
        zvec_rc(i) = reconstruct_loworder(gf_zvec(i), p, false, false);
      }
    }
    // End lower-order

    const vec<vec<CCTK_REAL, 2>, 3> zveclow_rc =
        calc_contraction(g_avg, zvec_rc);
    const auto z2_rc = calc_contraction(zveclow_rc, zvec_rc);

    w_lorentz_rc(0) = sqrt(1 + z2_rc(0));
    w_lorentz_rc(1) = sqrt(1 + z2_rc(1));

    for (int i = 0; i <= 2; ++i) {   // loop over components
      for (int j = 0; j <= 1; ++j) { // loop over left and right state
        vels_rc(i)(j) = zvec_rc(i)(j) / w_lorentz_rc(j);
        vlows_rc(i)(j) = zveclow_rc(i)(j) / w_lorentz_rc(j);
      }
    }
    break;
  };
  case rec_var_t::s_vec: {

    vec<vec<CCTK_REAL, 2>, 3> svec_rc([&](int i) ARITH_INLINE {
      return reconstruct_pt(gf_svec(i), p, false, false);
    });

    // Lower-order
    if (useLO) {

      for (int i = 0; i <= 2; ++i) { // loop over components
        svec_rc(i) = reconstruct_loworder(gf_svec(i), p, false, false);
      }
    }
    // End lower-order

    const vec<vec<CCTK_REAL, 2>, 3> sveclow_rc =
        calc_contraction(g_avg, svec_rc);
    const auto s2_rc = calc_contraction(sveclow_rc, svec_rc);

    w_lorentz_rc(0) =
        sqrt(0.5 + sqrt(0.25 + s2_rc(0) / rhoh_rc(0) / rhoh_rc(0)));
    w_lorentz_rc(1) =
        sqrt(0.5 + sqrt(0.25 + s2_rc(1) / rhoh_rc(1) / rhoh_rc(1)));

    for (int i = 0; i <= 2; ++i) {   // loop over components
      for (int j = 0; j <= 1; ++j) { // loop over left and right state
        vels_rc(i)(j) =
            svec_rc(i)(j) / w_lorentz_rc(j) / w_lorentz_rc(j) / rhoh_rc(j);
        vlows_rc(i)(j) =
            sveclow_rc(i)(j) / w_lorentz_rc(j) / w_lorentz_rc(j) / rhoh_rc(j);
      }
    }
    break;
  };
  }

  // Set vels to 0 if flagged for atm reset
  if (resetL) {
    w_lorentz_rc(0) = 1.0;
    for (int i = 0; i <= 2; ++i) {
      vels_rc(i)(0) = 0.0;
      vlows_rc(i)(0) = 0.0;
    }
  }
  if (resetR) {
    w_lorentz_rc(1) = 1.0;
    for (int i = 0; i <= 2; ++i) {
      vels_rc(i)(1) = 0.0;
      vlows_rc(i)(1) = 0.0;
    }
  }
  /* END RECONSTRUCTION */

  /* vtilde^i = alpha * v^i - beta^i */
  const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
    return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
      return alp_avg * vels_rc(i)(f) - betas_avg(i);
    });
  });

  /* alpha * b0 = W * B^i * v_i */
  const vec<CCTK_REAL, 2> alp_b0_rc([&](int f) ARITH_INLINE {
    return w_lorentz_rc(f) * calc_contraction(Bs_rc, vlows_rc)(f);
  });
  /* covariant magnetic field measured by the Eulerian observer */
  const vec<vec<CCTK_REAL, 2>, 3> Blows_rc = calc_contraction(g_avg, Bs_rc);
  /* B^2 = B^i * B_i */
  const vec<CCTK_REAL, 2> B2_rc = calc_contraction(Bs_rc, Blows_rc);
  /* covariant magnetic field measured by the comoving observer:
   *  b_i = B_i/W + alpha*b^0*v_i */
  const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
    return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
      return Blows_rc(i)(f) / w_lorentz_rc(f) + alp_b0_rc(f) * vlows_rc(i)(f);
    });
  });
  /* b^2 = b^{\mu} * b_{\mu} */
  const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
    return (B2_rc(f) + pow2(alp_b0_rc(f))) / pow2(w_lorentz_rc(f));
  });

  /* componets correspond to the dir_i we are considering */
  const CCTK_REAL beta_avg = betas_avg(dir_i);
  const vec<CCTK_REAL, 2> vel_rc{vels_rc(dir_i)};
  const vec<CCTK_REAL, 2> B_rc{Bs_rc(dir_i)};
  const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir_i)};

  // TODO: Compute pressure based on user-specified EOS.
  // Currently, computing press for classical ideal gas from reconstructed
  // vars

  const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
    return eos_3p->csnd_from_rho_temp_ye(rho_rc(f), temp_rc(f), Ye_rc(f)) *
           eos_3p->csnd_from_rho_temp_ye(rho_rc(f), temp_rc(f), Ye_rc(f));
  });

  const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
    return 1 + eps_rc(f) + press_rc(f) / rho_rc(f);
  });

  /* Computing conservatives from primitives: */

  /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
  const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
    return sqrtg * rho_rc(f) * w_lorentz_rc(f);
  });

  /* DEnt = sqrt(g) * D * s  = sqrt(g) * (rho * W) * s */
  /*    s = entropy */
  const vec<CCTK_REAL, 2> DEnt_rc([&](int f) ARITH_INLINE {
    return sqrtg * rho_rc(f) * w_lorentz_rc(f) * entropy_rc(f);
  });

  /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
  const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
    return dens_rc(f) * h_rc(f) * w_lorentz_rc(f);
  });
  /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
  const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc =
      dens_h_W_rc + sqrtg * (pow2(alp_b0_rc) + B2_rc);
  /* auxiliary: (pgas + pmag) */
  const vec<CCTK_REAL, 2> press_plus_pmag_rc = press_rc + 0.5 * bsq_rc;

  /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
  const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
    return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
      return dens_h_W_plus_sqrtg_W2b2_rc(f) * vlows_rc(i)(f) -
             sqrtg * alp_b0_rc(f) * blows_rc(i)(f);
    });
  });

  /* tau = sqrt(g)*t =
   *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
  const vec<CCTK_REAL, 2> tau_rc =
      dens_h_W_rc - dens_rc + sqrtg * (B2_rc - press_plus_pmag_rc);

  /* Computing fluxes of conserved variables: */

  /* auxiliary: unit in 'dir_i' */
  const vec<CCTK_REAL, 3> unit_dir_i{vec<int, 3>::unit(dir_i)};
  /* auxiliary: alpha * sqrt(g) */
  const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;
  /* auxiliary: B^i / W */
  const vec<CCTK_REAL, 2> B_over_w_lorentz_rc(
      [&](int f) ARITH_INLINE { return B_rc(f) / w_lorentz_rc(f); });

  /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
  const vec<CCTK_REAL, 2> flux_dens(
      [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

  /* flux(DEnt) = sqrt(g) * D * s * vtilde^i = sqrt(g) * rho * W * s *
   * vtilde^i */
  const vec<CCTK_REAL, 2> flux_DEnt(
      [&](int f) ARITH_INLINE { return DEnt_rc(f) * vtilde_rc(f); });

  /* flux(mom_j)^i = sqrt(g)*(
   *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
  const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
    return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
      return moms_rc(j)(f) * vtilde_rc(f) +
             alp_sqrtg * (press_plus_pmag_rc(f) * unit_dir_i(j) -
                          blows_rc(j)(f) * B_over_w_lorentz_rc(f));
    });
  });

  /* flux(tau) = sqrt(g)*(
   *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
  const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
    return tau_rc(f) * vtilde_rc(f) +
           alp_sqrtg * (press_plus_pmag_rc(f) * vel_rc(f) -
                        alp_b0_rc(f) * B_over_w_lorentz_rc(f));
  });

  /* flux(DYe) = sqrt(g) * (D * Ye * vtilde^i) */
  const vec<CCTK_REAL, 2> DYe_rc(
      [&](int f) ARITH_INLINE { return dens_rc(f) * Ye_rc(f); });
  const vec<CCTK_REAL, 2> flux_DYe(
      [&](int f) ARITH_INLINE { return DYe_rc(f) * vtilde_rc(f); });

  /* Calculate eigenvalues: */

  /* variable for either g^xx, g^yy or g^zz depending on the direction */
  const CCTK_REAL u_avg = calc_inv(g_avg, detg_avg)(dir_i, dir_i);
  /* eigenvalues -- collapsed immediately to the only wavespeed bounds any
   * consumer needs (charmax/charmin); the full lambda_tmp dies here instead of
   * living all the way to the UCT block. fmax/fmin reductions match the old
   * hlle/laxf/maxspeeds_from_lambdas exactly -> bit-identical. */
  const vec<vec<CCTK_REAL, 4>, 2> lambda_tmp =
      eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
                  w_lorentz_rc, h_rc, bsq_rc);
  CCTK_REAL charmax = 0, charmin = 0;
  for (int s = 0; s < 2; ++s)
    for (int m = 0; m < 4; ++m) {
      charmax = fmax(charmax, lambda_tmp(s)(m));
      charmin = fmin(charmin, lambda_tmp(s)(m));
    }

  /* Calculate numerical fluxes */
  if (!useLO || !loworder_flux) {
    fluxdenss(dir_i)(p.I) = calcflux(charmax, charmin, dens_rc, flux_dens);
    fluxDEnts(dir_i)(p.I) = calcflux(charmax, charmin, DEnt_rc, flux_DEnt);
    fluxmomxs(dir_i)(p.I) = calcflux(charmax, charmin, moms_rc(0), flux_moms(0));
    fluxmomys(dir_i)(p.I) = calcflux(charmax, charmin, moms_rc(1), flux_moms(1));
    fluxmomzs(dir_i)(p.I) = calcflux(charmax, charmin, moms_rc(2), flux_moms(2));
    fluxtaus(dir_i)(p.I) = calcflux(charmax, charmin, tau_rc, flux_tau);
    fluxDYes(dir_i)(p.I) = calcflux(charmax, charmin, DYe_rc, flux_DYe);
  } else {
    fluxdenss(dir_i)(p.I) = laxf_cc(charmax, charmin, dens_rc, flux_dens);
    fluxDEnts(dir_i)(p.I) = laxf_cc(charmax, charmin, DEnt_rc, flux_DEnt);
    fluxmomxs(dir_i)(p.I) = laxf_cc(charmax, charmin, moms_rc(0), flux_moms(0));
    fluxmomys(dir_i)(p.I) = laxf_cc(charmax, charmin, moms_rc(1), flux_moms(1));
    fluxmomzs(dir_i)(p.I) = laxf_cc(charmax, charmin, moms_rc(2), flux_moms(2));
    fluxtaus(dir_i)(p.I) = laxf_cc(charmax, charmin, tau_rc, flux_tau);
    fluxDYes(dir_i)(p.I) = laxf_cc(charmax, charmin, DYe_rc, flux_DYe);
  }

  if constexpr (!uct) {
    // flux-CT only: off-diagonal induction fluxes. Btildes_rc / Es_rc /
    // flux_Btildes are computed here, physically inside the compile-time
    // CT-scheme guard, so the upwind-CT instantiation carries none of them in
    // its live-set (not relying on DCE to prune them). Values are identical to
    // the pre-hoist function-scope computation -> bit-identical per config.
    /* Btildes^i = sqrt(g) * B^i */
    const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
        [&](int i) ARITH_INLINE { return sqrtg * Bs_rc(i); });
    /* electric field E_i = \tilde\epsilon_{ijk} Btilde_j * vtilde_k */
    const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
        calc_cross_product(Btildes_rc, vtildes_rc);
    /* flux(Btildes) = {{0, -Ez, Ey}, {Ez, 0, -Ex}, {-Ey, Ex, 0}} */
    const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
        calc_cross_product(unit_dir_i, Es_rc);
    if (!useLO || !loworder_flux) {
      fluxB_j(dir_i)(p.I) =
          calcflux(charmax, charmin, Btildes_rc(dir_j), flux_Btildes(dir_j));
      fluxB_k(dir_i)(p.I) =
          calcflux(charmax, charmin, Btildes_rc(dir_k), flux_Btildes(dir_k));
    } else {
      fluxB_j(dir_i)(p.I) =
          laxf_cc(charmax, charmin, Btildes_rc(dir_j), flux_Btildes(dir_j));
      fluxB_k(dir_i)(p.I) =
          laxf_cc(charmax, charmin, Btildes_rc(dir_k), flux_Btildes(dir_k));
    }
  }

  /* Positivity Preserving Limiter */
  // First, check if both cells are in the atmosphere. If so, and the
  // atmosphere is graded, the PP limiter can be spuriously activated before
  // the cell is reset to atmosphere after the evolution step, making
  // debugging difficult. At face Ip, Ip refers to the right cell and Im
  // refers to the left.
  // The entire PP-limiter block is gated at compile time on the CT-independent
  // `pplim` template parameter (production runs use_pplim=no). For the
  // use_pplim=no instantiation the whole `_ppl` recomputation live-set compiles
  // out, shedding its register/AGPR pressure. Bit-identical per config: the
  // dispatch keys the template on the runtime use_pplim, so `pplim` true <=>
  // use_pplim true, matching the old runtime `if (use_pplim && !ppl_atmo)`.
  if constexpr (pplim) {
  const auto Ip = p.I;
  const auto Im = p.I - p.DI[dir_i];

  vec<CCTK_REAL, 2> rho_ppl = {rho(Im), rho(Ip)};
  const bool ppl_atmo = ((rho_ppl(0) <= rho_atm(0) * (1 + atmo_tol)) &&
                         (rho_ppl(1) <= rho_atm(1) * (1 + atmo_tol)));

  // `use_pplim` here is a guaranteed-true redundant term inside the pplim
  // guard (dispatch keys the template on it), kept so this condition is
  // byte-identical to the pre-template code and `use_pplim` stays referenced.
  if (use_pplim && !ppl_atmo) {

    /* BEGIN REPEATED RECONSTRUCTION CODE */
    // Copy fluid quantities from cell center to face
    vec<CCTK_REAL, 2> entropy_ppl = {entropy(Im), entropy(Ip)};
    vec<CCTK_REAL, 2> Ye_ppl = {Ye(Im), Ye(Ip)};
    vec<CCTK_REAL, 2> eps_ppl = {eps(Im), eps(Ip)};
    vec<CCTK_REAL, 2> press_ppl = {press(Im), press(Ip)};
    vec<CCTK_REAL, 2> temp_ppl = {temperature(Im), temperature(Ip)};

    const vec<CCTK_REAL, 2> rhoh_ppl([&](int f) ARITH_INLINE {
      return rho_ppl(f) + rho_ppl(f) * eps_ppl(f) + press_ppl(f);
    });

    // Introduce reconstructed Bs
    // Use staggered dB for i == dir_i
    vec<vec<CCTK_REAL, 2>, 3> Bs_ppl;

    // Assign the value for the primary direction
    const CCTK_REAL val = gf_dBstags(dir_i)(p.I) / sqrtg;
    Bs_ppl(dir_i)(0) = val;
    Bs_ppl(dir_i)(1) = val;

    // Lambda to assign the reconstructed values
    auto assign_reconstructed = [&](int d) {
      vec<CCTK_REAL, 2> tmp = {gf_Bvecs(d)(Im), gf_Bvecs(d)(Ip)};
      Bs_ppl(d)(0) = tmp(0);
      Bs_ppl(d)(1) = tmp(1);
    };

    // Assign reconstructed values for the two perpendicular directions
    assign_reconstructed(dir_j);
    assign_reconstructed(dir_k);
    // End of setting Bs

    // Reconstruct zvec
    vec<vec<CCTK_REAL, 2>, 3> vels_ppl;
    vec<vec<CCTK_REAL, 2>, 3> vlows_ppl;
    vec<CCTK_REAL, 2> w_lorentz_ppl;

    vec<vec<CCTK_REAL, 2>, 3> zvec_ppl([&](int i) ARITH_INLINE {
      vec<CCTK_REAL, 2> tmp = {gf_zvec(i)(Im), gf_zvec(i)(Ip)};
      return tmp;
    });

    const vec<vec<CCTK_REAL, 2>, 3> zveclow_ppl =
        calc_contraction(g_avg, zvec_ppl);
    const auto z2_ppl = calc_contraction(zveclow_ppl, zvec_ppl);

    w_lorentz_ppl(0) = sqrt(1 + z2_ppl(0));
    w_lorentz_ppl(1) = sqrt(1 + z2_ppl(1));

    for (int i = 0; i <= 2; ++i) {   // loop over components
      for (int j = 0; j <= 1; ++j) { // loop over left and right state
        vels_ppl(i)(j) = zvec_ppl(i)(j) / w_lorentz_ppl(j);
        vlows_ppl(i)(j) = zveclow_ppl(i)(j) / w_lorentz_ppl(j);
      }
    }

    /* vtilde^i = alpha * v^i - beta^i */
    const vec<vec<CCTK_REAL, 2>, 3> vtildes_ppl([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return alp_avg * vels_ppl(i)(f) - betas_avg(i);
      });
    });

    /* alpha * b0 = W * B^i * v_i */
    const vec<CCTK_REAL, 2> alp_b0_ppl([&](int f) ARITH_INLINE {
      return w_lorentz_ppl(f) * calc_contraction(Bs_ppl, vlows_ppl)(f);
    });
    /* covariant magnetic field measured by the Eulerian observer */
    const vec<vec<CCTK_REAL, 2>, 3> Blows_ppl =
        calc_contraction(g_avg, Bs_ppl);
    /* B^2 = B^i * B_i */
    const vec<CCTK_REAL, 2> B2_ppl = calc_contraction(Bs_ppl, Blows_ppl);
    /* covariant magnetic field measured by the comoving observer:
     *  b_i = B_i/W + alpha*b^0*v_i */
    const vec<vec<CCTK_REAL, 2>, 3> blows_ppl([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return Blows_ppl(i)(f) / w_lorentz_ppl(f) +
               alp_b0_ppl(f) * vlows_ppl(i)(f);
      });
    });
    /* b^2 = b^{\mu} * b_{\mu} */
    const vec<CCTK_REAL, 2> bsq_ppl([&](int f) ARITH_INLINE {
      return (B2_ppl(f) + pow2(alp_b0_ppl(f))) / pow2(w_lorentz_ppl(f));
    });

    /* componets correspond to the dir_i we are considering */
    const vec<CCTK_REAL, 2> vel_ppl{vels_ppl(dir_i)};
    const vec<CCTK_REAL, 2> B_ppl{Bs_ppl(dir_i)};
    const vec<CCTK_REAL, 2> vtilde_ppl{vtildes_ppl(dir_i)};

    const vec<CCTK_REAL, 2> h_ppl([&](int f) ARITH_INLINE {
      return 1 + eps_ppl(f) + press_ppl(f) / rho_ppl(f);
    });

    /* Computing conservatives from primitives: */

    /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
    const vec<CCTK_REAL, 2> dens_ppl([&](int f) ARITH_INLINE {
      return sqrtg * rho_ppl(f) * w_lorentz_ppl(f);
    });

    /* DEnt = sqrt(g) * D * s  = sqrt(g) * (rho * W) * s */
    /*    s = entropy */
    const vec<CCTK_REAL, 2> DEnt_ppl([&](int f) ARITH_INLINE {
      return sqrtg * rho_ppl(f) * w_lorentz_ppl(f) * entropy_ppl(f);
    });

    /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
    const vec<CCTK_REAL, 2> dens_h_W_ppl([&](int f) ARITH_INLINE {
      return dens_ppl(f) * h_ppl(f) * w_lorentz_ppl(f);
    });
    /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
    const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_ppl =
        dens_h_W_ppl + sqrtg * (pow2(alp_b0_ppl) + B2_ppl);
    /* auxiliary: (pgas + pmag) */
    const vec<CCTK_REAL, 2> press_plus_pmag_ppl = press_ppl + 0.5 * bsq_ppl;

    /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
    const vec<vec<CCTK_REAL, 2>, 3> moms_ppl([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return dens_h_W_plus_sqrtg_W2b2_ppl(f) * vlows_ppl(i)(f) -
               sqrtg * alp_b0_ppl(f) * blows_ppl(i)(f);
      });
    });

    /* tau = sqrt(g)*t =
     *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
    const vec<CCTK_REAL, 2> tau_ppl =
        dens_h_W_ppl - dens_ppl + sqrtg * (B2_ppl - press_plus_pmag_ppl);

    /* Btildes^i = sqrt(g) * B^i */
    const vec<vec<CCTK_REAL, 2>, 3> Btildes_ppl(
        [&](int i) ARITH_INLINE { return sqrtg * Bs_ppl(i); });

    /* Computing fluxes of conserved variables: */

    /* auxiliary: B^i / W */
    const vec<CCTK_REAL, 2> B_over_w_lorentz_ppl(
        [&](int f) ARITH_INLINE { return B_ppl(f) / w_lorentz_ppl(f); });

    /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
    const vec<CCTK_REAL, 2> flux_dens_ppl(
        [&](int f) ARITH_INLINE { return dens_ppl(f) * vtilde_ppl(f); });

    /* flux(DEnt) = sqrt(g) * D * s * vtilde^i = sqrt(g) * rho * W * s *
     * vtilde^i */
    const vec<CCTK_REAL, 2> flux_DEnt_ppl(
        [&](int f) ARITH_INLINE { return DEnt_ppl(f) * vtilde_ppl(f); });

    /* flux(mom_j)^i = sqrt(g)*(
     *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
    const vec<vec<CCTK_REAL, 2>, 3> flux_moms_ppl([&](int j) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return moms_ppl(j)(f) * vtilde_ppl(f) +
               alp_sqrtg * (press_plus_pmag_ppl(f) * unit_dir_i(j) -
                            blows_ppl(j)(f) * B_over_w_lorentz_ppl(f));
      });
    });

    /* flux(tau) = sqrt(g)*(
     *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
    const vec<CCTK_REAL, 2> flux_tau_ppl([&](int f) ARITH_INLINE {
      return tau_ppl(f) * vtilde_ppl(f) +
             alp_sqrtg * (press_plus_pmag_ppl(f) * vel_ppl(f) -
                          alp_b0_ppl(f) * B_over_w_lorentz_ppl(f));
    });

    /* flux(DYe) = sqrt(g) * (D * Ye * vtilde^i) */
    const vec<CCTK_REAL, 2> DYe_ppl(
        [&](int f) ARITH_INLINE { return dens_ppl(f) * Ye_ppl(f); });
    const vec<CCTK_REAL, 2> flux_DYe_ppl(
        [&](int f) ARITH_INLINE { return DYe_ppl(f) * vtilde_ppl(f); });

    /* END REPEATED RECONSTRUCTION CODE */

    // Calc LO flux
    const CCTK_REAL fluxLOdenss = laxf_simple(dens_ppl, flux_dens_ppl);
    const CCTK_REAL fluxLODEnts = laxf_simple(DEnt_ppl, flux_DEnt_ppl);
    const CCTK_REAL fluxLODYes = laxf_simple(DYe_ppl, flux_DYe_ppl);
    const CCTK_REAL fluxLOmomxs = laxf_simple(moms_ppl(0), flux_moms_ppl(0));
    const CCTK_REAL fluxLOmomys = laxf_simple(moms_ppl(1), flux_moms_ppl(1));
    const CCTK_REAL fluxLOmomzs = laxf_simple(moms_ppl(2), flux_moms_ppl(2));
    const CCTK_REAL fluxLOtaus = laxf_simple(tau_ppl, flux_tau_ppl);

    // Get 2 * \alpha * CFL with \alpha = 3
    const CCTK_REAL a2cfl = 6 * CCTK_DELTA_TIME / p.DX[dir_i];

    // Calc dens floor
    //// At Ip
    const smat<CCTK_REAL, 3> g_avg_p([&](int i, int j) ARITH_INLINE {
      return calc_avg_v2c(gf_g(i, j), p);
    });
    const CCTK_REAL sqrtg_p = sqrt(calc_det(g_avg_p));
    const CCTK_REAL densmin_p = sqrtg_p * rho_atm(1);

    //// At Im
    const smat<CCTK_REAL, 3> g_avg_m([&](int i, int j) ARITH_INLINE {
      return calc_avg_v2c(gf_g(i, j), p, Im);
    });
    const CCTK_REAL sqrtg_m = sqrt(calc_det(g_avg_m));
    const CCTK_REAL densmin_m = sqrtg_m * rho_atm(0);

    // Calc theta from dens
    const CCTK_REAL newdens_p = dens(Ip) + a2cfl * fluxdenss(dir_i)(Ip);
    const CCTK_REAL newdens_m = dens(Im) - a2cfl * fluxdenss(dir_i)(Ip);

    const CCTK_REAL newdensLO_p = dens(Ip) + a2cfl * fluxLOdenss;
    const CCTK_REAL newdensLO_m = dens(Im) - a2cfl * fluxLOdenss;

    CCTK_REAL theta = 1.0;
    CCTK_REAL theta_m = 1.0;
    CCTK_REAL theta_p = 1.0;

    if (newdens_p < densmin_p)
      theta_p =
          min(theta,
              max(0.0, (newdensLO_p - densmin_p) /
                           (a2cfl * (fluxLOdenss - fluxdenss(dir_i)(Ip)))));

    if (newdens_m < densmin_m)
      theta_m =
          min(theta,
              max(0.0, (newdensLO_m - densmin_m) /
                           (a2cfl * (fluxdenss(dir_i)(Ip) - fluxLOdenss))));

    theta = min(theta_m, theta_p);

    // Calc theta from DYe
    // Note that this does not guarantee positivity of the primitive Ye
    if (istab) {
      const CCTK_REAL newDYe_p = DYe(Ip) + a2cfl * fluxDYes(dir_i)(Ip);
      const CCTK_REAL newDYe_m = DYe(Im) - a2cfl * fluxDYes(dir_i)(Ip);

      const CCTK_REAL newDYeLO_p = DYe(Ip) + a2cfl * fluxLODYes;
      const CCTK_REAL newDYeLO_m = DYe(Im) - a2cfl * fluxLODYes;

      const CCTK_REAL DYemin_m = densmin_m * eos_3p->rgye.min;
      const CCTK_REAL DYemin_p = densmin_p * eos_3p->rgye.min;

      if (newDYe_p < DYemin_p)
        theta_p =
            min(theta,
                max(0.0, (newDYeLO_p - DYemin_p) /
                             (a2cfl * (fluxLODYes - fluxDYes(dir_i)(Ip)))));

      if (newDYe_m < DYemin_m)
        theta_m =
            min(theta,
                max(0.0, (newDYeLO_m - DYemin_m) /
                             (a2cfl * (fluxDYes(dir_i)(Ip) - fluxLODYes))));

      theta = min(theta_m, theta_p);
    }

    // Update flux GF
    fluxdenss(dir_i)(Ip) =
        (1 - theta) * fluxLOdenss + theta * fluxdenss(dir_i)(Ip);
    fluxDEnts(dir_i)(Ip) =
        (1 - theta) * fluxLODEnts + theta * fluxDEnts(dir_i)(Ip);
    fluxDYes(dir_i)(Ip) =
        (1 - theta) * fluxLODYes + theta * fluxDYes(dir_i)(Ip);
    fluxmomxs(dir_i)(Ip) =
        (1 - theta) * fluxLOmomxs + theta * fluxmomxs(dir_i)(Ip);
    fluxmomys(dir_i)(Ip) =
        (1 - theta) * fluxLOmomys + theta * fluxmomys(dir_i)(Ip);
    fluxmomzs(dir_i)(Ip) =
        (1 - theta) * fluxLOmomzs + theta * fluxmomzs(dir_i)(Ip);
    fluxtaus(dir_i)(Ip) =
        (1 - theta) * fluxLOtaus + theta * fluxtaus(dir_i)(Ip);
    gf_theta(dir_i)(Ip) = theta;
  } else {
    gf_theta(dir_i)(p.I) = 1.0;
  }
  }
  // use_pplim=no: the PP-limiter block (and its theta write) is compiled out.
  // theta_x/y/z have no storage in that config (see schedule.ccl); theta is
  // identically 1.0, and neither consumer needs to read it -- the upwind-CT
  // drift blend below specializes on theta == 1.0 (dropping the cell-centred
  // fallback term outright) and the rhs theta_tot diagnostic uses the constant
  // directly. So there is nothing to write here.

#ifdef CCTK_DEBUG
  // Recompute the flux-CT induction quantities for diagnostics only. In a
  // production (non-CCTK_DEBUG) build these do not exist for the upwind-CT
  // instantiation (they live inside the if constexpr(!uct) block above); here
  // they are recomputed unconditionally so this NaN dump stays complete for
  // both CT schemes without carrying them in the production live-set.
  const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
      [&](int i) ARITH_INLINE { return sqrtg * Bs_rc(i); });
  const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes = calc_cross_product(
      unit_dir_i, calc_cross_product(Btildes_rc, vtildes_rc));
  bool nan_found =
      isnan(dens_rc(0)) || isnan(dens_rc(1)) || isnan(moms_rc(0)(0)) ||
      isnan(moms_rc(0)(1)) || isnan(moms_rc(1)(0)) || isnan(moms_rc(1)(1)) ||
      isnan(moms_rc(2)(0)) || isnan(moms_rc(2)(1)) || isnan(tau_rc(0)) ||
      isnan(tau_rc(1)) || isnan(Btildes_rc(0)(0)) ||
      isnan(Btildes_rc(0)(1)) || isnan(Btildes_rc(1)(0)) ||
      isnan(Btildes_rc(1)(1)) || isnan(Btildes_rc(2)(0)) ||
      isnan(Btildes_rc(2)(1)) || isnan(DYe_rc(0)) || isnan(DYe_rc(1)) ||
      isnan(flux_dens(0)) || isnan(flux_dens(1)) || isnan(flux_moms(0)(0)) ||
      isnan(flux_moms(0)(1)) || isnan(flux_moms(1)(0)) ||
      isnan(flux_moms(1)(1)) || isnan(flux_moms(2)(0)) ||
      isnan(flux_moms(2)(1)) || isnan(flux_DYe(0)) || isnan(flux_DYe(1)) ||
      isnan(flux_tau(0)) || isnan(flux_tau(1)) || isnan(flux_Btildes(0)(0)) ||
      isnan(flux_Btildes(0)(1)) || isnan(flux_Btildes(1)(0)) ||
      isnan(flux_Btildes(1)(1)) || isnan(flux_Btildes(2)(0)) ||
      isnan(flux_Btildes(2)(1)) || isnan(fluxdenss(dir_i)(p.I)) ||
      isnan(fluxmomxs(dir_i)(p.I)) || isnan(fluxmomys(dir_i)(p.I)) ||
      isnan(fluxmomzs(dir_i)(p.I)) || isnan(fluxtaus(dir_i)(p.I)) ||
      rho_rc(0) < 0.0 || rho_rc(1) < 0.0 || press_rc(0) < 0.0 ||
      press_rc(1) < 0.0;
  if constexpr (!uct) // fluxB GFs only have storage in the flux-CT config
    nan_found = nan_found || isnan(fluxB_j(dir_i)(p.I)) ||
                isnan(fluxB_k(dir_i)(p.I));
  if (nan_found) {
    printf("cctk_iteration = %i,  dir_i = %i,  ijk = %i, %i, %i, "
           "x, y, z = %16.8e, %16.8e, %16.8e.\n",
           cctk_iteration, dir_i, p.i, p.j, p.k, p.x, p.y, p.z);
    printf("  fluxdenss = %16.8e,\n", fluxdenss(dir_i)(p.I));
    printf("  fluxmoms  = %16.8e, %16.8e, %16.8e,\n", fluxmomxs(dir_i)(p.I),
           fluxmomys(dir_i)(p.I), fluxmomzs(dir_i)(p.I));
    printf("  fluxtaus  = %16.8e,\n", fluxtaus(dir_i)(p.I));
    if constexpr (!uct) // fluxB GFs only have storage in the flux-CT config
      printf("  fluxBs(j,k) = %16.8e, %16.8e\n", fluxB_j(dir_i)(p.I),
             fluxB_k(dir_i)(p.I));
    printf("  flux_denss = %16.8e, %16.8e,\n", flux_dens(0), flux_dens(1));
    printf("  flux_moms  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
           flux_moms(0)(0), flux_moms(0)(1), flux_moms(1)(0), flux_moms(1)(1),
           flux_moms(2)(0), flux_moms(2)(1));
    printf("  flux_taus  = %16.8e, %16.8e,\n", flux_tau(0), flux_tau(1));
    printf("  flux_DYes  = %16.8e, %16.8e,\n", flux_DYe(0), flux_DYe(1));
    printf("  flux_Bts   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
           flux_Btildes(0)(0), flux_Btildes(0)(1), flux_Btildes(1)(0),
           flux_Btildes(1)(1), flux_Btildes(2)(0), flux_Btildes(2)(1));
    printf("  dens_rc = %16.8e, %16.8e,\n", dens_rc(0), dens_rc(1));
    printf("  moms_rc = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
           moms_rc(0)(0), moms_rc(0)(1), moms_rc(1)(0), moms_rc(1)(1),
           moms_rc(2)(0), moms_rc(2)(1));
    printf("  tau_rc  = %16.8e, %16.8e,\n", tau_rc(0), tau_rc(1));
    printf("  DYe_rc  = %16.8e, %16.8e,\n", DYe_rc(0), DYe_rc(1));
    printf("  Bs_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
           Bs_rc(0)(0), Bs_rc(0)(1), Bs_rc(1)(0), Bs_rc(1)(1), Bs_rc(2)(0),
           Bs_rc(2)(1));
    printf("  Bts_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
           Btildes_rc(0)(0), Btildes_rc(0)(1), Btildes_rc(1)(0),
           Btildes_rc(1)(1), Btildes_rc(2)(0), Btildes_rc(2)(1));
    printf("  charmax = %16.8e, charmin = %16.8e\n", charmax, charmin);
    printf("  alp_avg = %16.8e, beta_avg = %16.8e, u_avg = %16.8e \n",
           alp_avg, beta_avg, u_avg);
    printf("  vel_rc  = %16.8e, %16.8e \n", vel_rc(0), vel_rc(1));
    printf("  rho_rc  = %16.8e, %16.8e \n", rho_rc(0), rho_rc(1));
    printf("  cs2_rc  = %16.8e, %16.8e \n", cs2_rc(0), cs2_rc(1));
    printf("  wlor_rc = %16.8e, %16.8e \n", w_lorentz_rc(0), w_lorentz_rc(1));
    printf("  h_rc    = %16.8e, %16.8e \n", h_rc(0), h_rc(1));
    printf("  bsq_rc  = %16.8e, %16.8e \n", bsq_rc(0), bsq_rc(1));
    printf("  press_rc = %16.8e, %16.8e \n", press_rc(0), press_rc(1));
    printf("  eps_rc   = %16.8e, %16.8e \n", eps_rc(0), eps_rc(1));
    printf("  rho = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
           rho(p.I - p.DI[dir_i] * 3), rho(p.I - p.DI[dir_i] * 2),
           rho(p.I - p.DI[dir_i]), rho(p.I), rho(p.I + p.DI[dir_i]),
           rho(p.I + p.DI[dir_i] * 2));
    printf("  press = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
           press(p.I - p.DI[dir_i] * 3), press(p.I - p.DI[dir_i] * 2),
           press(p.I - p.DI[dir_i]), press(p.I), press(p.I + p.DI[dir_i]),
           press(p.I + p.DI[dir_i] * 2));
    printf("  eps   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
           eps(p.I - p.DI[dir_i] * 3), eps(p.I - p.DI[dir_i] * 2),
           eps(p.I - p.DI[dir_i]), eps(p.I), eps(p.I + p.DI[dir_i]),
           eps(p.I + p.DI[dir_i] * 2));
    printf("  alp_avg, beta_avg = %16.8e, %16.8e, %16.8e, %16.8e,\n", alp_avg,
           betas_avg(0), betas_avg(1), betas_avg(2));
    printf("  g_avg = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
           g_avg(0, 0), g_avg(0, 1), g_avg(0, 2), g_avg(1, 1), g_avg(1, 2),
           g_avg(2, 2));
    printf("  sqrtg = %16.8e,\n", sqrtg);
    printf("  vlows_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
           vlows_rc(0)(0), vlows_rc(0)(1), vlows_rc(1)(0), vlows_rc(1)(1),
           vlows_rc(2)(0), vlows_rc(2)(1));
    printf("  vups_rc   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
           vels_rc(0)(0), vels_rc(0)(1), vels_rc(1)(0), vels_rc(1)(1),
           vels_rc(2)(0), vels_rc(2)(1));
    printf("  vtilde_rc = %16.8e, %16.8e.\n", vtilde_rc(0), vtilde_rc(1));
    assert(0);
  }
#endif

  /* Begin code for upwindCT */

  if constexpr (uct) { // upwind-CT only: face speeds and drift velocities
    // maxspeeds_from_lambdas gives ap = max(0, lmax) = charmax and
    // am = max(0, -lmin) = -charmin, so reuse the collapsed bounds directly.
    const CCTK_REAL ap = charmax;
    const CCTK_REAL am = -charmin;

    ap_face(dir_i)(p.I) = ap;
    am_face(dir_i)(p.I) = am;

    const CCTK_REAL vjL = vtildes_rc(dir_j)(0);
    const CCTK_REAL vjR = vtildes_rc(dir_j)(1);
    const CCTK_REAL vkL = vtildes_rc(dir_k)(0);
    const CCTK_REAL vkR = vtildes_rc(dir_k)(1);
    const CCTK_REAL vj_face = avg_upwind(vjL, vjR, ap, am);
    const CCTK_REAL vk_face = avg_upwind(vkL, vkR, ap, am);

    // theta_x/y/z are stored only when use_pplim (see schedule.ccl); with the
    // limiter off theta is identically 1.0, so use that constant and never
    // touch the (unstored) GF. For use_pplim=yes this reads the same value as
    // before -> bit-identical per config.
    //
    // The two branches are also split on pplim so that the cell-centred
    // velocity fallback is COMPILED OUT of the production kernel rather than
    // merely multiplied by zero. With theta_uct == 1.0 the blend degenerates to
    // 1.0 * v_face + 0.0 * (...), and the compiler is NOT allowed to fold the
    // 0.0 * x away: that needs both nnan (x could be Inf, and Inf * 0 -> NaN)
    // and nsz (0.0 * x is -0.0 for x < 0, so the product is not a compile-time
    // constant) -- i.e. fast-math, which this project deliberately does not
    // enable. So without this branch the compiler must emit 4 gf_vels loads per
    // direction in a configuration that provably cannot use them: 12 dead loads
    // in the production kernel.
    //
    // What this buys, measured (gfx90a, -Rpass, 2026-07-27): 12 fewer global
    // loads, and *no* change in register pressure or occupancy (256 VGPR /
    // 32 AGPR / occ 1, byte-for-byte). The register count is set by peak
    // simultaneous liveness, which lives in reconstruction; this block runs
    // after flux assembly, i.e. past the peak. Do not expect occupancy from
    // deletions down here -- see Planning/CHECKPOINT.md Lesson 9.
    //
    // Bit-identical for pplim=yes (same expression) and exactly equal for
    // pplim=no on finite data (1.0*x == x, x + 0.0 == x; the sole exception is
    // x == -0.0, which becomes +0.0, invisible to a norm-based compare).
    if constexpr (pplim) {
      const CCTK_REAL theta_uct = gf_theta(dir_i)(p.I);
      vbar_j(dir_i)(p.I) =
          theta_uct * vj_face +
          (1.0 - theta_uct) * 0.5 *
              (gf_vels(dir_j)(p.I) + gf_vels(dir_j)(p.I - p.DI[dir_i]));
      vbar_k(dir_i)(p.I) =
          theta_uct * vk_face +
          (1.0 - theta_uct) * 0.5 *
              (gf_vels(dir_k)(p.I) + gf_vels(dir_k)(p.I - p.DI[dir_i]));
    } else {
      vbar_j(dir_i)(p.I) = vj_face;
      vbar_k(dir_i)(p.I) = vk_face;
    }
  }

  /* End code for upwindCT */
}

/* The iteration box that grid.loop_mixpn_device<CI,CJ,CK>(grid.nghostzones,
 * ord) traverses (see CarpetX Loop/src/loop_device.hxx, loop_mixpn_device):
 * all points of the tile in the cell-centred (transverse) directions,
 * interior +-ord in the vertex-centred (face-normal) direction, with the
 * +-ord applied after the tile clamping. Computed with the same
 * box_int/box_all calls on the same grid object, so the fused sweep below
 * reproduces the pre-fusion per-direction domains exactly. */
template <int CI, int CJ, int CK>
static void calc_mixpn_box(const GridDescBaseDevice &grid, const int ord,
                           vect<int, dim> &imin, vect<int, dim> &imax) {
  vect<int, dim> imin_int, imax_int, imin_all, imax_all;
  grid.box_int<CI, CJ, CK>(grid.nghostzones, imin_int, imax_int);
  grid.box_all<CI, CJ, CK>(grid.nghostzones, imin_all, imax_all);
  constexpr vect<int, dim> facetype{CI, CJ, CK};
  for (int d = 0; d < dim; ++d) {
    imin[d] = facetype[d] ? imin_all[d] : imin_int[d] - ord;
    imax[d] = facetype[d] ? imax_all[d] : imax_int[d] + ord;
  }
}

// Calculate the fluxes in ALL THREE directions in one fused sweep over the
// grid: at every point of the union of the three per-direction face domains,
// each direction whose domain contains the point gets its face flux computed
// (CalcFluxAtFace<0/1/2>). This replaces the previous three per-direction
// sweeps -- the input fields are streamed through the cache once instead of
// three times, and one kernel launches instead of three. The set of faces
// computed, and the arithmetic per face, are identical to the per-direction
// sweeps (golden-master verified).
template <bool uct, bool pplim, typename EOSType>
void CalcFluxAll(CCTK_ARGUMENTS, EOSType *eos_3p, const rec_var_t rec_var,
                 const reconstruction_t reconstruction,
                 const reconstruction_t reconstruction_LO,
                 const reconstruct_params_t reconstruct_params,
                 const flux_t fluxtype) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  // Ghost zones required by the reconstruction stencil. The fused sweep
  // computes all three directions, so all of them must satisfy it (the same
  // requirement the three per-direction CalcFlux<0/1/2> calls asserted).
  const auto require_ghosts = [&](const int n) {
    for (int d = 0; d < dim; ++d)
      assert(cctk_nghostzones[d] >= n);
  };
  switch (reconstruction) {
  case reconstruction_t::Godunov:
    require_ghosts(1);
    break;
  case reconstruction_t::minmod:
    require_ghosts(2);
    break;
  case reconstruction_t::monocentral:
    require_ghosts(2);
    break;
  case reconstruction_t::ppm:
    require_ghosts(3);
    break;
  case reconstruction_t::eppm:
    require_ghosts(3);
    break;
  case reconstruction_t::wenoz:
    require_ghosts(3);
    break;
  case reconstruction_t::wenozp:
    require_ghosts(3);
    break;
  case reconstruction_t::mp5:
    require_ghosts(3);
    break;
  }

  /* grid functions for fluxes */
  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss{fxdens, fydens, fzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts{fxDEnt, fyDEnt, fzDEnt};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDYes{fxDYe, fyDYe, fzDYe};
  /* off-diagonal B-fluxes, indexed by the flux direction i: flux of B_{dir_j}
   * resp. B_{dir_k} in direction i, cf. vbar_j/vbar_k below (the diagonal
   * flux of B_{dir_i} in direction dir_i vanishes identically and has no
   * grid function) */
  const vec<GF3D2<CCTK_REAL>, dim> fluxB_j{fxBy, fyBz, fzBx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxB_k{fxBz, fyBx, fzBy};

  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_zvec{zvec_x, zvec_y, zvec_z};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_svec{svec_x, svec_y, svec_z};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_dBstags{dBx_stag, dBy_stag,
                                                    dBz_stag};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  /* grid functions for Upwind CT */
  const vec<GF3D2<CCTK_REAL>, dim> vbar_j{vbar_y_xface, vbar_z_yface,
                                          vbar_x_zface};
  const vec<GF3D2<CCTK_REAL>, dim> vbar_k{vbar_z_xface, vbar_x_yface,
                                          vbar_y_zface};
  const vec<GF3D2<CCTK_REAL>, dim> ap_face{amax_xface, amax_yface, amax_zface};
  const vec<GF3D2<CCTK_REAL>, dim> am_face{amin_xface, amin_yface, amin_zface};

  /* grid functions for PP flux limiter */
  const vec<GF3D2<CCTK_REAL>, dim> gf_theta{theta_x, theta_y, theta_z};

  // Velocity limit from Con2PrimFactory parameters
  const CCTK_REAL w_lim = sqrt(1.0 + vw_lim * vw_lim);
  const CCTK_REAL v_lim = vw_lim / w_lim;

  // Flag for tabulated EOS
  const bool istab = CCTK_EQUALS(evolution_eos, "Tabulated3d") ? true : false;

  /* Context for the per-face flux computation. Positional aggregate
   * initialization: the initializers below must stay in exactly the order of
   * the member declarations of FluxContext (each initializer is the
   * identically-named local/argument/parameter). */
  const FluxContext<EOSType> fx{
      /* grid functions for fluxes */
      fluxdenss,
      fluxDEnts,
      fluxmomxs,
      fluxmomys,
      fluxmomzs,
      fluxtaus,
      fluxDYes,
      fluxB_j,
      fluxB_k,
      /* grid functions */
      gf_vels,
      gf_zvec,
      gf_svec,
      gf_Bvecs,
      gf_dBstags,
      gf_beta,
      gf_g,
      /* grid functions for Upwind CT */
      vbar_j,
      vbar_k,
      ap_face,
      am_face,
      /* grid functions for PP flux limiter */
      gf_theta,
      /* cell-centred scalar grid functions */
      alp,
      rho,
      press,
      eps,
      entropy,
      Ye,
      temperature,
      LOflag,
      dens,
      DYe,
      /* EOS object and solver settings */
      eos_3p,
      rec_var,
      reconstruction,
      reconstruction_LO,
      reconstruct_params,
      fluxtype,
      /* flag for tabulated EOS */
      istab,
      /* velocity limit */
      v_lim,
      /* parameters */
      r_atmo,
      rho_abs_min,
      n_rho_atmo,
      recon_thresh,
      p_atmo,
      n_press_atmo,
      t_atmo,
      n_temp_atmo,
      Ye_atmo,
      atmo_tol,
      use_press_atmo,
      reconstruct_with_temperature,
      use_pplim,
      loworder_flux,
      /* flesh scalars */
      cctk_delta_time,
      cctk_timefac,
      cctk_iteration,
  };

  const int nloop = (hydro_correction_order - 2) / 2;

  /* Per-direction iteration boxes: for each direction, the
   * zero-initialization box (exactly what the pre-fusion
   * loop_all_device<face_centred...>(grid.nghostzones) init pass traversed)
   * and the flux box (exactly what loop_mixpn_device traversed; see
   * calc_mixpn_box). */
  vect<vect<int, dim>, dim> zmin, zmax;
  grid.box_all<0, 1, 1>(grid.nghostzones, zmin[0], zmax[0]);
  grid.box_all<1, 0, 1>(grid.nghostzones, zmin[1], zmax[1]);
  grid.box_all<1, 1, 0>(grid.nghostzones, zmin[2], zmax[2]);
  vect<int, dim> imin0, imax0, imin1, imax1, imin2, imax2;
  calc_mixpn_box<0, 1, 1>(grid, nloop, imin0, imax0);
  calc_mixpn_box<1, 0, 1>(grid, nloop, imin1, imax1);
  calc_mixpn_box<1, 1, 0>(grid, nloop, imin2, imax2);

  /* One fused sweep over the union of all the boxes (the shared integer
   * index space makes all three face staggerings addressable from one loop)
   * does the zeroing AND the flux computation: per point, first zero every
   * direction whose init box contains the point, then compute every
   * direction whose flux box does. This is value-identical to running the
   * init passes before the flux passes: the flux computation never READS a
   * zero-initialized value (its gf_theta and own-face flux reads are
   * preceded by same-point writes), and where boxes overlap across tiles
   * the writes are idempotent, exactly as before the fusion. Fusing saves a
   * separate full write pass over the face grid functions (the zero and the
   * flux value coalesce in cache) and a kernel launch. The PointDesc fields
   * the per-face computation uses -- p.I, p.DI, p.DX -- are
   * centering-independent; the face coordinate is NOT taken from p.X but
   * computed per direction below. */
  const vect<int, dim> fmin =
      min(min(min(zmin[0], zmin[1]), min(zmin[2], imin0)), min(imin1, imin2));
  const vect<int, dim> fmax =
      max(max(max(zmax[0], zmax[1]), max(zmax[2], imax0)), max(imax1, imax2));
  vect<int, dim> bnd_min, bnd_max;
  grid.boundary_box<0, 0, 0>(grid.nghostzones, bnd_min, bnd_max);

  /* CarpetX computes p.X = x0 + (lbnd + I - (!CI)/2) * dx for a loop of
   * centering CI (Loop/src/loop.hxx, point_desc). face_X below evaluates the
   * same expression with each direction's face centering, so it is
   * bit-identical to the p.X that the per-direction face-centred loops
   * supplied to CalcFluxAtFace before the fusion. */
  const vect<CCTK_REAL, dim> x0 = grid.x0;
  const vect<CCTK_REAL, dim> dx = grid.dx;
  const vect<int, dim> lbnd = grid.lbnd;

  grid.loop_box_device<0, 0, 0>(
      bnd_min, bnd_max, fmin, fmax, [=] CCTK_DEVICE(const PointDesc &p) {
        // Force-capture fx here, before the if constexpr blocks below: nvcc
        // rejects an extended __device__ lambda that first-captures a variable
        // inside a constexpr-if, and fx's first use would otherwise be in the
        // guarded zeroing (its unconditional use in CalcFluxAtFace comes later).
        static_cast<void>(fx);
        // initialize to zero
        for (int dir = 0; dir < dim; ++dir) {
          if (all(p.I >= zmin[dir]) && all(p.I < zmax[dir])) {
            fluxdenss(dir)(p.I) = 0;
            fluxDEnts(dir)(p.I) = 0;
            fluxmomxs(dir)(p.I) = 0;
            fluxmomys(dir)(p.I) = 0;
            fluxmomzs(dir)(p.I) = 0;
            fluxtaus(dir)(p.I) = 0;
            fluxDYes(dir)(p.I) = 0;
            // Zero the CT-scheme GFs through fx (already captured by the
            // CalcFluxAtFace calls below) rather than the local vec's: nvcc
            // forbids an extended __device__ lambda from first-capturing a
            // variable inside an if constexpr, and these are referenced only
            // in the guarded branches.
            if constexpr (!uct) { // fluxB GFs only stored in the flux-CT config
              fx.fluxB_j(dir)(p.I) = 0;
              fx.fluxB_k(dir)(p.I) = 0;
            }
            if constexpr (uct) { // face GFs only stored in the upwind-CT config
              fx.ap_face(dir)(p.I) = 0;
              fx.am_face(dir)(p.I) = 0;
              fx.vbar_j(dir)(p.I) = 0;
              fx.vbar_k(dir)(p.I) = 0;
            }
            if constexpr (pplim) // theta GFs only stored with the PP limiter
              fx.gf_theta(dir)(p.I) = 1.0;
          }
        }

        if (all(p.I >= imin0) && all(p.I < imax0)) {
          constexpr vect<bool, dim> face_centred{false, true, true};
          const vect<CCTK_REAL, dim> face_X =
              x0 + (lbnd + p.I - vect<CCTK_REAL, dim>(!face_centred) / 2) * dx;
          CalcFluxAtFace<0, uct, pplim>(fx, p, face_X);
        }
        if (all(p.I >= imin1) && all(p.I < imax1)) {
          constexpr vect<bool, dim> face_centred{true, false, true};
          const vect<CCTK_REAL, dim> face_X =
              x0 + (lbnd + p.I - vect<CCTK_REAL, dim>(!face_centred) / 2) * dx;
          CalcFluxAtFace<1, uct, pplim>(fx, p, face_X);
        }
        if (all(p.I >= imin2) && all(p.I < imax2)) {
          constexpr vect<bool, dim> face_centred{true, true, false};
          const vect<CCTK_REAL, dim> face_X =
              x0 + (lbnd + p.I - vect<CCTK_REAL, dim>(!face_centred) / 2) * dx;
          CalcFluxAtFace<2, uct, pplim>(fx, p, face_X);
        }
      });
}

extern "C" void AsterX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  eos_3param eos_3p_type;

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  rec_var_t rec_var;
  if (CCTK_EQUALS(recon_type, "v_vec")) {
    rec_var = rec_var_t::v_vec;
  } else if (CCTK_EQUALS(recon_type, "z_vec")) {
    rec_var = rec_var_t::z_vec;
  } else if (CCTK_EQUALS(recon_type, "s_vec")) {
    rec_var = rec_var_t::s_vec;
  } else {
    CCTK_ERROR("Unknown value for parameter \"recon_type\"");
  }

  // Primary reconstruction method
  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else if (CCTK_EQUALS(reconstruction_method, "eppm"))
    reconstruction = reconstruction_t::eppm;
  else if (CCTK_EQUALS(reconstruction_method, "wenoz"))
    reconstruction = reconstruction_t::wenoz;
  else if (CCTK_EQUALS(reconstruction_method, "wenozp"))
    reconstruction = reconstruction_t::wenozp;
  else if (CCTK_EQUALS(reconstruction_method, "mp5"))
    reconstruction = reconstruction_t::mp5;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  // Lower-order fallback for negative values
  reconstruction_t reconstruction_LO;
  if (CCTK_EQUALS(loworder_method, "Godunov"))
    reconstruction_LO = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(loworder_method, "minmod"))
    reconstruction_LO = reconstruction_t::minmod;
  else if (CCTK_EQUALS(loworder_method, "monocentral"))
    reconstruction_LO = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(loworder_method, "ppm"))
    reconstruction_LO = reconstruction_t::ppm;
  else
    CCTK_ERROR("Unknown value for parameter \"loworder_method\"");

  // reconstruction parameters struct
  reconstruct_params_t reconstruct_params;

  // ppm parameters
  reconstruct_params.ppm_shock_detection = ppm_shock_detection;
  reconstruct_params.ppm_zone_flattening = ppm_zone_flattening;
  reconstruct_params.poly_k = poly_k;
  reconstruct_params.poly_gamma = poly_gamma;
  reconstruct_params.ppm_eta1 = ppm_eta1;
  reconstruct_params.ppm_eta2 = ppm_eta2;
  reconstruct_params.ppm_eps = ppm_eps;
  reconstruct_params.ppm_eps_shock = ppm_eps_shock;
  reconstruct_params.ppm_small = ppm_small;
  reconstruct_params.ppm_omega1 = ppm_omega1;
  reconstruct_params.ppm_omega2 = ppm_omega2;
  reconstruct_params.enhanced_ppm_C2 = enhanced_ppm_C2;
  // wenoz parameters
  reconstruct_params.weno_eps = weno_eps;
  reconstruct_params.weno_mp = weno_mp;
  // mp5 parameters
  reconstruct_params.mp5_alpha = mp5_alpha;

  flux_t fluxtype;
  if (CCTK_EQUALS(flux_type, "LxF")) {
    fluxtype = flux_t::LxF;
  } else if (CCTK_EQUALS(flux_type, "HLLE")) {
    fluxtype = flux_t::HLLE;
  } else {
    CCTK_ERROR("Unknown value for parameter \"flux_type\"");
  }

  // The CT scheme (use_uct) is threaded as a compile-time template parameter
  // through CalcFluxAll/CalcFluxAtFace so that each configuration only compiles
  // its own half of the per-face work (flux-CT induction fluxes vs. upwind-CT
  // face speeds/drift velocities). Dispatch the runtime use_uct flag to the two
  // instantiations via a single generic lambda. Reading use_uct per call
  // respects its STEERABLE=always semantics.
  const auto run_all = [&](auto UCT, auto PPLIM) {
    constexpr bool uct = decltype(UCT)::value;
    constexpr bool pplim = decltype(PPLIM)::value;
    switch (eos_3p_type) {
    case eos_3param::IdealGas: {
      // Get local eos object
      auto eos_3p_ig = global_eos_3p_ig;

      CalcFluxAll<uct, pplim>(cctkGH, eos_3p_ig, rec_var, reconstruction,
                       reconstruction_LO, reconstruct_params, fluxtype);
      break;
    }
    case eos_3param::Hybrid: {
      // Note: the nested if conditions below could be inefficient.
      // This needs to be tested, and restructured, if required.
      if (global_eos_3p_hyb_pwpoly) {
        auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;

        CalcFluxAll<uct, pplim>(cctkGH, eos_3p_hyb, rec_var, reconstruction,
                         reconstruction_LO, reconstruct_params, fluxtype);

      } else if (global_eos_3p_hyb_poly) {
        auto eos_3p_hyb = global_eos_3p_hyb_poly;

        CalcFluxAll<uct, pplim>(cctkGH, eos_3p_hyb, rec_var, reconstruction,
                         reconstruction_LO, reconstruct_params, fluxtype);

      } else {
        CCTK_ERROR(
            "Hybrid EOS selected but no hybrid EOS object was initialized");
      }

      break;
    }
    case eos_3param::Tabulated: {
      // Get local eos object
      auto eos_3p_tab3d = global_eos_3p_tab3d;

      CalcFluxAll<uct, pplim>(cctkGH, eos_3p_tab3d, rec_var, reconstruction,
                       reconstruction_LO, reconstruct_params, fluxtype);
      break;
    }
    default:
      assert(0);
    }
  };

  // Dispatch both compile-time flags (use_uct, use_pplim) from their runtime
  // parameters. Reading them per call respects their STEERABLE=always
  // semantics. 2 (uct) x 2 (pplim) x EOS instantiations of CalcFluxAll.
  if (use_uct) {
    if (use_pplim)
      run_all(std::true_type{}, std::true_type{});
    else
      run_all(std::true_type{}, std::false_type{});
  } else {
    if (use_pplim)
      run_all(std::false_type{}, std::true_type{});
    else
      run_all(std::false_type{}, std::false_type{});
  }
}

template <int i, bool use_uct>
void CalcE_impl(CCTK_ARGUMENTS, const reconstruction_t reconstruction,
                const reconstruct_params_t reconstruct_params) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcAuxTermsForAvecPsiRHS;

  // the other two directions
  constexpr int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
  constexpr int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);

  // flux-CT: off-diagonal B-fluxes, indexed by the flux direction as in
  // CalcFlux (gf_fB_j(d) = flux of B_{dir_j(d)} in direction d, etc.; the
  // diagonal fluxes vanish identically and have no grid functions)
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fB_j{fxBy, fyBz, fzBx};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fB_k{fxBz, fyBx, fzBy};

  // upwind-CT
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> dB_stag{dBx_stag, dBy_stag, dBz_stag};
  const vec<GF3D2<const CCTK_REAL>, dim> ap_face{amax_xface, amax_yface,
                                                 amax_zface};
  const vec<GF3D2<const CCTK_REAL>, dim> am_face{amin_xface, amin_yface,
                                                 amin_zface};
  const vec<vec<GF3D2<const CCTK_REAL>, 2>, dim> vbars{
      {vbar_x_yface, vbar_x_zface},
      {vbar_y_zface, vbar_y_xface},
      {vbar_z_xface, vbar_z_yface}};
  // mapping from 3d to 2d, since we don't need iface
  constexpr int jface = 1;
  constexpr int kface = 0;

  const vec<GF3D2<CCTK_REAL>, dim> gf_E{Ex, Ey, Ez};

  // edge centered loop
  if constexpr (use_uct) { // upwind-CT
    grid.loop_int_device<i == 0, i == 1, i == 2>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // reconstruct in k-dir
          const vec<CCTK_REAL, 2> dBstag_jface_krc{
              reconstruct(dB_stag(j), p, reconstruction, k, false, false, press,
                          gf_vels(k), reconstruct_params)};
          const vec<CCTK_REAL, 2> vbar_k_jface_krc{
              reconstruct(vbars(k)(jface), p, reconstruction, k, false, false,
                          press, gf_vels(k), reconstruct_params)};
          // reconstruct in j-dir
          const vec<CCTK_REAL, 2> dBstag_kface_jrc{
              reconstruct(dB_stag(k), p, reconstruction, j, false, false, press,
                          gf_vels(j), reconstruct_params)};
          const vec<CCTK_REAL, 2> vbar_j_kface_jrc{
              reconstruct(vbars(j)(kface), p, reconstruction, j, false, false,
                          press, gf_vels(j), reconstruct_params)};

          const CCTK_REAL BjL = dBstag_jface_krc(0);
          const CCTK_REAL BjR = dBstag_jface_krc(1);
          const CCTK_REAL vkL = vbar_k_jface_krc(0);
          const CCTK_REAL vkR = vbar_k_jface_krc(1);

          const CCTK_REAL BkL = dBstag_kface_jrc(0);
          const CCTK_REAL BkR = dBstag_kface_jrc(1);
          const CCTK_REAL vjL = vbar_j_kface_jrc(0);
          const CCTK_REAL vjR = vbar_j_kface_jrc(1);

          const CCTK_REAL ap_k = ap_face(k)(p.I);
          const CCTK_REAL am_k = am_face(k)(p.I);
          const CCTK_REAL ap_j = ap_face(j)(p.I);
          const CCTK_REAL am_j = am_face(j)(p.I);

          gf_E(i)(p.I) =
              hll_upwind(BjL, BjR, vkL * BjL, vkR * BjR, ap_k, am_k) -
              hll_upwind(BkL, BkR, vjL * BkL, vjR * BkR, ap_j, am_j);
        });
  } else { // flux-CT
    grid.loop_int_device<i == 0, i == 1, i == 2>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // (i,j,k) is a cyclic permutation, so the flux of B_j in
          // direction k is direction k's dir_k-component, gf_fB_k(k), and
          // the flux of B_k in direction j is gf_fB_j(j)
          const CCTK_REAL Fjk = gf_fB_k(k)(p.I);
          const CCTK_REAL Fjk_m = gf_fB_k(k)(p.I - p.DI[j]);
          const CCTK_REAL Fkj = gf_fB_j(j)(p.I);
          const CCTK_REAL Fkj_m = gf_fB_j(j)(p.I - p.DI[k]);
          gf_E(i)(p.I) = CCTK_REAL(0.25) * ((Fjk + Fjk_m) - (Fkj + Fkj_m));
        });
  }
}

template <int i>
void CalcE(CCTK_ARGUMENTS, const bool use_uct,
           const reconstruction_t reconstruction,
           const reconstruct_params_t reconstruct_params) {
  if (use_uct) {
    CalcE_impl<i, true>(CCTK_PASS_CTOC, reconstruction, reconstruct_params);
  } else {
    CalcE_impl<i, false>(CCTK_PASS_CTOC, reconstruction, reconstruct_params);
  }
}

template <int i> void CalcFstag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcAuxTermsForAvecPsiRHS;
  DECLARE_CCTK_PARAMETERS;

  // the other two directions
  constexpr int j = (i == 0) ? 1 : ((i == 1) ? 2 : 0);
  constexpr int k = (i == 0) ? 2 : ((i == 1) ? 0 : 1);

  const vec<GF3D2<CCTK_REAL>, dim> gf_Fstag{Fx_stag, Fy_stag, Fz_stag};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Avecs{Avec_x, Avec_y, Avec_z};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_mix_device<i == 0, i == 1, i == 2>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL alp_e = calc_avg_v2e<i>(alp, p);
        const smat<CCTK_REAL, 3> g_e([&](int m, int n) ARITH_INLINE {
          return calc_avg_v2e<i>(gf_g(m, n), p);
        });
        const CCTK_REAL detg_e = calc_det(g_e);
        const CCTK_REAL sqrtg_e = sqrt(detg_e);
        const smat<CCTK_REAL, 3> ug_e = calc_inv(g_e, detg_e);

        vec<CCTK_REAL, 3> A_e;
        A_e(i) = gf_Avecs(i)(p.I);
        A_e(j) = calc_avg_e2e<i, j>(gf_Avecs(j), p);
        A_e(k) = calc_avg_e2e<i, k>(gf_Avecs(k), p);

        const vec<CCTK_REAL, 3> Aup_e = calc_contraction(ug_e, A_e);

        gf_Fstag(i)(p.I) = alp_e * sqrtg_e * Aup_e(i);
      });
}

extern "C" void AsterX_CalcAuxTermsForAvecPsiRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcAuxTermsForAvecPsiRHS;
  DECLARE_CCTK_PARAMETERS;

  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else if (CCTK_EQUALS(reconstruction_method, "eppm"))
    reconstruction = reconstruction_t::eppm;
  else if (CCTK_EQUALS(reconstruction_method, "wenoz"))
    reconstruction = reconstruction_t::wenoz;
  else if (CCTK_EQUALS(reconstruction_method, "wenozp"))
    reconstruction = reconstruction_t::wenozp;
  else if (CCTK_EQUALS(reconstruction_method, "mp5"))
    reconstruction = reconstruction_t::mp5;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  // reconstruction parameters struct
  reconstruct_params_t reconstruct_params;

  // ppm parameters
  reconstruct_params.ppm_shock_detection = ppm_shock_detection;
  reconstruct_params.ppm_zone_flattening = ppm_zone_flattening;
  reconstruct_params.poly_k = poly_k;
  reconstruct_params.poly_gamma = poly_gamma;
  reconstruct_params.ppm_eta1 = ppm_eta1;
  reconstruct_params.ppm_eta2 = ppm_eta2;
  reconstruct_params.ppm_eps = ppm_eps;
  reconstruct_params.ppm_eps_shock = ppm_eps_shock;
  reconstruct_params.ppm_small = ppm_small;
  reconstruct_params.ppm_omega1 = ppm_omega1;
  reconstruct_params.ppm_omega2 = ppm_omega2;
  reconstruct_params.enhanced_ppm_C2 = enhanced_ppm_C2;
  // wenoz parameters
  reconstruct_params.weno_eps = weno_eps;
  reconstruct_params.weno_mp = weno_mp;
  // mp5 parameters
  reconstruct_params.mp5_alpha = mp5_alpha;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_Avecs{Avec_x, Avec_y, Avec_z};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  CalcE<0>(CCTK_PASS_CTOC, use_uct, reconstruction, reconstruct_params);
  CalcE<1>(CCTK_PASS_CTOC, use_uct, reconstruction, reconstruct_params);
  CalcE<2>(CCTK_PASS_CTOC, use_uct, reconstruction, reconstruct_params);

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Fbetax(p.I) = betax(p.I) * Psi(p.I);
                                      Fbetay(p.I) = betay(p.I) * Psi(p.I);
                                      Fbetaz(p.I) = betaz(p.I) * Psi(p.I);

                                      // Init G
                                      G(p.I) = 0.0;
                                    });

  CalcFstag<0>(CCTK_PASS_CTOC);
  CalcFstag<1>(CCTK_PASS_CTOC);
  CalcFstag<2>(CCTK_PASS_CTOC);

  grid.loop_allm1_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const vec<CCTK_REAL, 3> A_vert([&](int i) ARITH_INLINE {
          return calc_avg_e2v(gf_Avecs(i), p, i);
        });
        const smat<CCTK_REAL, 3> g(
            [&](int i, int j) ARITH_INLINE { return gf_g(i, j)(p.I); });
        const vec<CCTK_REAL, 3> betas(
            [&](int i) ARITH_INLINE { return gf_beta(i)(p.I); });
        const CCTK_REAL detg = calc_det(g);
        const CCTK_REAL sqrtg = sqrt(detg);

        G(p.I) = alp(p.I) * Psi(p.I) / sqrtg - calc_contraction(betas, A_vert);
      });
}

} // namespace AsterX
