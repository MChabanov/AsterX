#ifndef ASTERX_EIGENVALUES_HXX
#define ASTERX_EIGENVALUES_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Arith;

/* Fast-magnetosonic characteristic speeds, fused form.
 *
 * Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0 is a quadratic
 * a_2 lambda^2 + a_1 lambda + a_0 = 0 with
 *
 *   a_0 = (b^2 + cs^2 h rho)(beta^2 - alp^2 u) - (cs^2-1) h rho (beta-alp v)^2 W^2
 *   a_1 = 2 beta (b^2 + cs^2 h rho) - 2 (cs^2-1) h rho (beta-alp v) W^2
 *   a_2 = b^2 + h rho (cs^2 + W^2 - cs^2 W^2)
 *
 * Introducing the total enthalpy H = rho*h + b^2 and the fast speed
 *
 *   vf2 = (b^2 + cs^2 rho h)/H = cA^2 + cs^2 (1 - cA^2),   cA^2 = b^2/H
 *
 * one has, term by term,   b^2 + cs^2 h rho = vf2 H,
 * (cs^2-1) h rho = (vf2-1) H  and  a_2 = H [vf2 + W^2 (1-vf2)], i.e. ALL THREE
 * coefficients carry the common factor H.  Since the roots are homogeneous of
 * degree zero in the a_n, H cancels identically: the quadratic is exactly the
 * pure GR-hydro acoustic quadratic with cs^2 -> vf2, and the whole magnetic
 * sector reaches the wave speeds through the single scalar vf2.
 *
 * With K = (1-vf2) W^2 and a2h = vf2 + K the discriminant collapses (the
 * beta-dependence cancels via beta - (beta - alp v) = alp v):
 *
 *   (1/4) Delta = vf2 alp^2 [u a2h - K v^2]
 *
 * and, using  beta vf2 + K(beta - alp v) = beta a2h - K alp v,
 *
 *   lambda_pm = -beta + alp [ K v +- sqrt(vf2 (u a2h - K v^2)) ] / a2h
 *
 * The radicand is strictly positive by construction:
 * u a2h - K v^2 >= u vf2 + K(u - v^2), and Cauchy-Schwarz on
 * v^{dir_i} = delta^{dir_i}_j v^j gives (v^{dir_i})^2 <= g^{dir_i dir_i} v_j v^j
 * < u.  The fmax below is therefore roundoff insurance only -- if it ever bites,
 * something upstream (EOS, velocity limiter) is wrong.  Contrast the unfused
 * form, where a_1^2 - 4 a_2 a_0 is a cancellation-prone difference of large
 * like-signed numbers all carrying the dimensional factor H and the clamp was
 * load-bearing.
 *
 * Returns the two roots (index 0 = lambda_+, 1 = lambda_-) per side.  The
 * fourfold duplication of the old 4-vector return was only ever the degenerate
 * fast pair; the caller reduces to charmax/charmin anyway. */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST vec<vec<CCTK_REAL, 2>, 2>
    eigenvalues(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> vf2,
                vec<CCTK_REAL, 2> w_lor) {
  return vec<vec<CCTK_REAL, 2>, 2>([&](int f) ARITH_INLINE {
    const CCTK_REAL K = (1 - vf2(f)) * pow2(w_lor(f));
    const CCTK_REAL a2h = vf2(f) + K;
    const CCTK_REAL rad = u_avg * a2h - K * pow2(vel(f));
    const CCTK_REAL disc = alp_avg * sqrt(fmax(0.0, vf2(f) * rad));
    const CCTK_REAL drift = alp_avg * K * vel(f);
    return vec<CCTK_REAL, 2>{-beta_avg + (drift + disc) / a2h,
                             -beta_avg + (drift - disc) / a2h};
  });
};

} // namespace AsterX

#endif // ASTERX_EIGENVALUES_HXX
