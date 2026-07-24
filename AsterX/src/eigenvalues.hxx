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

// Characteristic speeds for ONE face state.
// See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0. The two face sides
// are fully independent, so eigenvalues() below is just this called once per
// side; the flux loop (Idea 3, per-side serialization) calls it per side
// directly. Returns 4 roots (a degenerate pair) so the downstream fmax/fmin
// collapse to charmax/charmin is unchanged.
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec<CCTK_REAL, 4>
eigenvalues_oneside(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                    CCTK_REAL vel, CCTK_REAL rho, CCTK_REAL cs2,
                    CCTK_REAL w_lor, CCTK_REAL h, CCTK_REAL bsq) {
  vec<CCTK_REAL, 3> a{
      (bsq + cs2 * h * rho) * (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
          (-1 + cs2) * h * rho * pow2(beta_avg - alp_avg * vel) * pow2(w_lor),

      2 * beta_avg * (bsq + cs2 * h * rho) -
          2 * (-1 + cs2) * h * rho * (beta_avg - alp_avg * vel) * pow2(w_lor),

      bsq + h * rho * (cs2 + pow2(w_lor) - cs2 * pow2(w_lor))};

  CCTK_REAL det = pow2(a(1)) - 4 * a(2) * a(0);
  if (det < 0)
    det = 0;

  return vec<CCTK_REAL, 4>{((-a(1) + sqrt(det)) / (2 * a(2))),
                           ((-a(1) + sqrt(det)) / (2 * a(2))),
                           ((-a(1) - sqrt(det)) / (2 * a(2))),
                           ((-a(1) - sqrt(det)) / (2 * a(2)))};
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST vec<vec<CCTK_REAL, 4>, 2>
    eigenvalues(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                vec<CCTK_REAL, 2> h, vec<CCTK_REAL, 2> bsq) {
  // minus side = index 0, plus side = index 1 (the two are independent)
  vec<vec<CCTK_REAL, 4>, 2> lambda{
      eigenvalues_oneside(alp_avg, beta_avg, u_avg, vel(0), rho(0), cs2(0),
                          w_lor(0), h(0), bsq(0)),
      eigenvalues_oneside(alp_avg, beta_avg, u_avg, vel(1), rho(1), cs2(1),
                          w_lor(1), h(1), bsq(1))};
  return lambda;
};

} // namespace AsterX

#endif // ASTERX_EIGENVALUES_HXX
