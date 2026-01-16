#ifndef ASTER_FD_HXX
#define ASTER_FD_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <simd.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>

namespace AsterUtils {
using namespace std;
using namespace Loop;
using namespace Arith;

// Laplace operator in direction dir
// Eq. (19) in https://arxiv.org/pdf/2310.11831
template <int dir, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
laplace_1d(const GF3D2<T> &gf, const PointDesc &p) {
  return ( gf(p.I - p.DI[dir])-T(2)*gf(p.I)+gf(p.I + p.DI[dir]) );
}

// Laplace operator in all (3) directions
// Eq. (18) in https://arxiv.org/pdf/2310.11831
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
laplace_3d(const GF3D2<T> &gf, const PointDesc &p) {
  return ( laplace_1d<0,T>(gf,p) + 
           laplace_1d<1,T>(gf,p) + 
	   laplace_1d<2,T>(gf,p) );
}

// Laplace operator in all perpendicular directions
// to dir_perp
// Eq. (20) in https://arxiv.org/pdf/2310.11831
template <int dir_perp, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
laplace_perp(const GF3D2<T> &gf, const PointDesc &p) {

  if constexpr (dir_perp == 0) {
     	return laplace_1d<1,T>(gf,p) + laplace_1d<2,T>(gf,p);
  } else if constexpr (dir_perp == 1) {
     	return laplace_1d<0,T>(gf,p) + laplace_1d<2,T>(gf,p);
  } else {
     	return laplace_1d<0,T>(gf,p) + laplace_1d<1,T>(gf,p);
  }
}

// FD2: vertex centered input, vertex centered output, oneside stencil
template <int Sign, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_v2v_oneside(const GF3D2<const T> &gf, const PointDesc &p,
                     const int dir) {
  static_assert(Sign == +1 || Sign == -1, "Sign must be +1 or -1");
  constexpr int s = Sign;

  const auto i0 = p.I;
  const auto i1 = i0 + s * p.DI[dir];
  const auto i2 = i0 + 2 * s * p.DI[dir];

  const T f0 = gf(i0);
  const T f1 = gf(i1);
  const T f2 = gf(i2);

  const T num = T(s) * (T(-3) * f0 + T(4) * f1 - f2);
  return num * (T(0.5) / p.DX[dir]);
}

// FD2: vertex centered input gridfunction along dir, cell centered output
// derivative
template <int dir, typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_forward_midpoint(const GF3D2<const T> &gf, const PointDesc &p,
                          const int order) {
  // Fill stencil
  // This assumes a maximum order/stencil size of 6
  const vect<int, dim> Im = p.I;
  array<vect<int, dim>, 6> stencil = {Im, Im, Im,
                                      Im, Im, Im}; // dummy value for init
  const int nstencil = (order / 2) - 1;
  int offset = -nstencil;
  for (int ii = 2 - nstencil; ii <= 3 + nstencil; ii++) {
    stencil[ii] = Im + offset * p.DI[dir];
    offset += 1;
  }

  return (1.0 / p.DX[dir]) *
         ((order == 2) * (gf(stencil[3]) - gf(stencil[2])) +
          (order == 4) * ((9.0 / 8.0) * (gf(stencil[3]) - gf(stencil[2])) -
                          (1.0 / 24.0) * (gf(stencil[4]) - gf(stencil[1]))) +
          (order == 6) * ((75.0 / 64.0) * (gf(stencil[3]) - gf(stencil[2])) -
                          (25.0 / 384.0) * (gf(stencil[4]) - gf(stencil[1])) +
                          (3.0 / 640.0) * (gf(stencil[5]) - gf(stencil[0]))));
}

// FD2: cell centered input gridfunction along dir, vertex centered output
// derivative
template <typename T>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
calc_fd2_backward_midpoint(const GF3D2<const T> &gf, const PointDesc &p,
                           const int dir) {
  return (gf(p.I) - gf(p.I - p.DI[dir])) / p.DX[dir];
}

// FD2: vertex centered input, cell centered output
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 2, T>
    calc_fd_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;
  const int j = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int k = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);

  const auto I = p.I;
  const auto dI = p.DI[dir];
  const auto dJ = p.DI[j];
  const auto dK = p.DI[k];
  const auto dJdK = dJ + dK;

  const auto line_fd = [&](const auto &B) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return gf(B + dI) - gf(B);
  };

  const T num =
      line_fd(I) + line_fd(I + dJ) + line_fd(I + dK) + line_fd(I + dJdK);

  const T inv4dx = T(0.25) / p.DX[dir];
  return num * inv4dx;
}

// FD4: vertex centered input, cell centered output
// Interpolation from edges to cell centers is 2nd order
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 4, T>
    calc_fd_v2c(const GF3D2<const T> &gf, const PointDesc &p, int dir) {
  T dgf1, dgf2, dgf3, dgf4;
  const int j = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
  const int k = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);

  const auto I = p.I;
  const auto dI = p.DI[dir];
  const auto dJ = p.DI[j];
  const auto dK = p.DI[k];
  const auto dJdK = dJ + dK;

  const auto line_fd4 = [&](const auto &B) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const T f_m1 = gf(B - dI);
    const T f_0 = gf(B);
    const T f_p1 = gf(B + dI);
    const T f_p2 = gf(B + dI + dI);
    return (T(1) * f_m1 + T(-27) * f_0 + T(27) * f_p1 + T(-1) * f_p2);
  };

  const T num =
      line_fd4(I) + line_fd4(I + dJ) + line_fd4(I + dK) + line_fd4(I + dJdK);

  const T inv96dx = T(1) / (T(96) * p.DX[dir]);
  return num * inv96dx;
}

// FD2: cell centered input, cell centered output
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 2, T>
    calc_fd_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  const auto I = p.I;
  const auto dI = p.DI[dir];

  const T fp = gf(I + dI);
  const T fm = gf(I - dI);

  const T inv2dx = T(0.5) / p.DX[dir];
  return (fp - fm) * inv2dx;
}

// FD4: cell centered input, cell centered output
template <int FDORDER, typename T>
CCTK_DEVICE CCTK_HOST
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::enable_if_t<FDORDER == 4, T>
    calc_fd_c2c(const GF3D2<const T> &gf, const PointDesc &p, const int dir) {
  const auto I = p.I;
  const auto dI = p.DI[dir];

  const T f_p1 = gf(I + dI);
  const T f_m1 = gf(I - dI);
  const T f_p2 = gf(I + dI + dI);
  const T f_m2 = gf(I - dI - dI);

  const T inv12dx = T(1) / (T(12) * p.DX[dir]);
  return ((T(-1) * f_p2) + (T(8) * f_p1) + (T(-8) * f_m1) + (T(1) * f_m2)) *
         inv12dx;
}

} // namespace AsterUtils

#endif // ASTER_FD_HXX
