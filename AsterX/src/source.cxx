#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

template <int FDORDER> void SourceTerms(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SourceTerms;
  DECLARE_CCTK_PARAMETERS;

  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_k{kxx, kxy, kxz, kyy, kyz, kzz};

  /* Loop over the entire grid (0 to n-1 cells in each direction) */
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Computing metric components at cell centers */
        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> beta_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const smat<CCTK_REAL, 3> k_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_k(i, j), p);
        });

        /* Determinant of spatial metric */
        const CCTK_REAL detg = calc_det(g_avg);
        const CCTK_REAL sqrt_detg = sqrt(detg);
        /* Upper metric */
        const smat<CCTK_REAL, 3> ug_avg = calc_inv(g_avg, detg);

        /* Derivatives of the lapse, shift and metric */
        /* calc_fd_v2c takes vertex center input, computes edge-center
         * derivatives along either dir=0, 1 or 2 using 2nd finite difference,
         * i.e, at the four edges of the cube with p.I as starting point. The
         * four edge-centered values are then interpolated to the cell-center
         * using 2nd order interpolation */
        const vec<CCTK_REAL, 3> d_alp([&](int k) ARITH_INLINE {
          return calc_fd_v2c<FDORDER>(alp, p, k);
        });
        const vec<vec<CCTK_REAL, 3>, 3> d_beta([&](int k) ARITH_INLINE {
          return vec<CCTK_REAL, 3>([&](int i) ARITH_INLINE {
            return calc_fd_v2c<FDORDER>(gf_beta(i), p, k);
          });
        });
        const vec<smat<CCTK_REAL, 3>, 3> d_g([&](int k) ARITH_INLINE {
          return smat<CCTK_REAL, 3>([&](int i, int j) ARITH_INLINE {
            return calc_fd_v2c<FDORDER>(gf_g(i, j), p, k);
          });
        });

        /* Computing v_j */
        //const vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
        //const vec<CCTK_REAL, 3> v_low = calc_contraction(g_avg, v_up);
        //const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

        /* Computing v_j */
        const vec<CCTK_REAL, 3> z_up{zvec_x(p.I), zvec_y(p.I), zvec_z(p.I)};
        const vec<CCTK_REAL, 3> z_low = calc_contraction(g_avg, z_up);
        const CCTK_REAL w_lorentz = calc_wlorentz_zvec(z_low, z_up);
        const vec<CCTK_REAL, 3> v_up{z_up/w_lorentz};
        const vec<CCTK_REAL, 3> v_low{z_low/w_lorentz};

        /* Computing [ \rho(1+\epsilon) + Pgas ]*W^2 = \rho * h * W^2 */
        const CCTK_REAL rhoenthalpyW2 =
            (rho(p.I) * (1.0 + eps(p.I)) + press(p.I)) * w_lorentz * w_lorentz;
        /* v^i - beta^i / alpha */
        const vec<CCTK_REAL, 3> velshift = v_up - beta_avg / alp_avg;

        /* cell-centered B^i and B_i */
        const vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        const vec<CCTK_REAL, 3> B_low = calc_contraction(g_avg, B_up);
        /* b^mu: b^0 and b^i */
        const CCTK_REAL bs0 =
            w_lorentz * calc_contraction(B_up, v_low) / alp_avg;
        const vec<CCTK_REAL, 3> bs =
            (B_up / w_lorentz + alp_avg * bs0 * velshift);
        /* b^2 */
        const CCTK_REAL bs2 =
            (calc_contraction(B_up, B_low) + pow2(alp_avg * bs0)) /
            pow2(w_lorentz);
        /* b_i = g_{ia} * b^a = g_{i0} * b^0 + g_{ik} * b^k */
        const vec<CCTK_REAL, 3> bs_low =
            calc_contraction(g_avg, beta_avg) * bs0 +
            calc_contraction(g_avg, bs);

        /* Computing T^{\mu\nu} */
        const CCTK_REAL t00 = ((rhoenthalpyW2 + bs2 * pow2(w_lorentz)) -
                               (press(p.I) + 0.5 * bs2)) /
                                  pow2(alp_avg) -
                              pow2(bs0);

        const vec<CCTK_REAL, 3> t0 =
            (rhoenthalpyW2 + bs2 * pow2(w_lorentz)) * velshift / alp_avg +
            (press(p.I) + 0.5 * bs2) * beta_avg / pow2(alp_avg) - bs0 * bs;

        const smat<CCTK_REAL, 3> t([&](int i, int j) ARITH_INLINE {
          return (rhoenthalpyW2 + bs2 * pow2(w_lorentz)) * velshift(i) *
                     velshift(j) +
                 (press(p.I) + 0.5 * bs2) *
                     (ug_avg(i, j) -
                      beta_avg(i) * beta_avg(j) / pow2(alp_avg)) -
                 bs(i) * bs(j);
        });

        /* Computing T^0_i = (\rho * h + b^2) * u^0 * u_i -b^0 * b_i */
        const vec<CCTK_REAL, 3> t0low =
            (rhoenthalpyW2 + bs2 * pow2(w_lorentz)) * v_low / alp_avg -
            bs_low * bs0;

        /* Contract the shift with the extrinsic curvature */
        const vec<CCTK_REAL, 3> shiftk = calc_contraction(k_avg, beta_avg);
        const CCTK_REAL shiftshiftk = calc_contraction(shiftk, beta_avg);

        /* Update term for tau */
        const CCTK_REAL tau_source =
            t00 * (shiftshiftk - calc_contraction(beta_avg, d_alp)) +
            calc_contraction(t0, -d_alp + 2.0 * shiftk) +
            calc_contraction(k_avg, t);

        /* Update term for mom */
        const vec<CCTK_REAL, 3> mom_source([&](int k) ARITH_INLINE {
          return t00 *
                     (0.5 * calc_norm(beta_avg, d_g(k)) - alp_avg * d_alp(k)) +
                 calc_contraction(d_g(k), t0, beta_avg) +
                 calc_contraction(t0low, d_beta(k)) +
                 0.5 * calc_contraction(t, d_g(k));
        });

        /* Update the RHS grid functions */
        densrhs(p.I) = 0.0;
        momxrhs(p.I) = alp_avg * sqrt_detg * mom_source(0);
        momyrhs(p.I) = alp_avg * sqrt_detg * mom_source(1);
        momzrhs(p.I) = alp_avg * sqrt_detg * mom_source(2);
        taurhs(p.I) = alp_avg * sqrt_detg * tau_source;

        // Add some debugging lines

        if (isnan(v_up(0)) || isnan(v_up(1)) || isnan(v_up(2)) ||
            isnan(eps(p.I)) || isnan(rho(p.I)) || isnan(press(p.I)) ||
            isnan(Bvecx(p.I)) || isnan(Bvecy(p.I)) || isnan(Bvecz(p.I)) ||
            rho(p.I) < 0.0 || press(p.I) < 0.0 ||
            eps(p.I) < 0.0 || isnan(momxrhs(p.I)) || isnan(momyrhs(p.I)) || isnan(momzrhs(p.I)) ||
            isnan(taurhs(p.I)) || isnan(w_lorentz)) {
      printf("cctk_iteration = %i, ijk = %i, %i, %i, "
             "x, y, z = %16.16e, %16.16e, %16.16e.\n",
             cctk_iteration, p.i, p.j, p.k, p.x, p.y, p.z);
      printf("  momxrhs  = %16.16e \n", momxrhs(p.I));
      printf("  momyrhs  = %16.16e \n", momyrhs(p.I));
      printf("  momzrhs  = %16.16e \n", momzrhs(p.I));
      printf("  taurhs  = %16.16e \n", taurhs(p.I));
      printf("  wlor  = %16.16e \n", w_lorentz);
      printf("  vupx  = %16.16e \n", v_up(0));
      printf("  vupy  = %16.16e \n", v_up(1));
      printf("  vupz  = %16.16e \n", v_up(2));
      printf("  v_lowx  = %16.16e \n", v_low(0));
      printf("  v_lowy  = %16.16e \n", v_low(1));
      printf("  v_lowz  = %16.16e \n", v_low(2));
      printf("  v*v     = %16.16e \n", calc_contraction(v_low, v_up));
      printf("  rho  = %16.16e \n", rho(p.I));
      printf("  press = %16.16e \n", press(p.I));
      printf("  eps   = %16.16e \n", eps(p.I));
      printf("  Bvecx  = %16.16e \n", Bvecx(p.I));
      printf("  Bvecy  = %16.16e \n", Bvecy(p.I));
      printf("  Bvecz  = %16.16e \n", Bvecz(p.I));
      printf("  alp_avg  = %16.16e \n", alp_avg);
      printf("  betax_avg  = %16.16e \n", beta_avg(0));
      printf("  betay_avg  = %16.16e \n", beta_avg(1));
      printf("  betaz_avg  = %16.16e \n", beta_avg(2));
      printf("  gxx_avg = %16.16e \n", g_avg(0,0));
      printf("  gxy_avg   = %16.16e \n",g_avg(0,1));
      printf("  gxz_avg  = %16.16e \n", g_avg(0,2));
      printf("  gyy_avg  = %16.16e \n", g_avg(1,1));
      printf("  gyz_avg  = %16.16e \n", g_avg(1,2));
      printf("  gzz_avg  = %16.16e \n", g_avg(2,2));
      printf("  detg  = %16.16e \n", detg);
      printf("  sqrt_detg  = %16.16e \n", sqrt_detg);     
      assert(0);
    }

      }); // end of loop over grid
}

extern "C" void AsterX_SourceTerms(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SourceTerms;
  DECLARE_CCTK_PARAMETERS;

  /* order of finite differencing */
  switch (local_spatial_order) {
  case 2:
    SourceTerms<2>(cctkGH);
    break;
  case 4:
    SourceTerms<4>(cctkGH);
    break;
  default:
    CCTK_VERROR("local_spatial_order must be set to 2 or 4.");
  }
}

} // namespace AsterX
