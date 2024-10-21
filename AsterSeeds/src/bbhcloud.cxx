#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "eos.hxx"
#include "eos_idealgas.hxx"
#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

enum class bbh_id_t { Gaussian, Constant, PowerLaw };

extern "C" void BBHCloud_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_BBHCloud_Initialize;
  DECLARE_CCTK_PARAMETERS;

  bbh_id_t bbh_id;
  if (CCTK_EQUALS(initial_distribution, "gaussian")) {
    bbh_id = bbh_id_t::Gaussian;
  } else if (CCTK_EQUALS(initial_distribution, "constant")) {
    bbh_id = bbh_id_t::Constant;
  } else if (CCTK_EQUALS(initial_distribution, "powerlaw")) {
    bbh_id = bbh_id_t::PowerLaw;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_distribution\"");
  }

  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  // We'll need the metric below to compute the velocity
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          // Metric
          const smat<CCTK_REAL, 3> glo(
              [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

          // Coordinate transformation Spherical -> Cartesian
          CCTK_REAL rr   = max(sqrt(p.x*p.x + p.y*p.y + p.z*p.z),1.0e-14);
          CCTK_REAL rcyl = max(sqrt(p.x*p.x + p.y*p.y),1.0e-14);
          CCTK_REAL costheta = p.z/rr;
          CCTK_REAL sintheta = rcyl/rr;
          CCTK_REAL cosphi = p.x/rcyl;
          CCTK_REAL sinphi = p.y/rcyl;

          // Note that rad_zphi_0 = rr * z_phi_0 such that
          // zr_0, rad_ztheta_0 and rad_zphi_0 have approx. 
          // the same order of magnitude
          CCTK_REAL zx_0 = sintheta*cosphi*zr_0 + 
                           costheta*cosphi*rad_ztheta_0 - 
                           sintheta*sinphi*rad_zphi_0;
    
          CCTK_REAL zy_0 = sintheta*sinphi*zr_0 + 
                           costheta*sinphi*rad_ztheta_0 + 
                           sintheta*cosphi*rad_zphi_0;
    
          CCTK_REAL zz_0 = costheta*zr_0 - sintheta*rad_ztheta_0;
    
          vec<CCTK_REAL, 3> z_0{zx_0,zy_0,zz_0};
          const vec<CCTK_REAL, 3> zlow_0 = calc_contraction(glo, z_0);
          CCTK_REAL w_0 = calc_wlorentz_zvec(z_0,zlow_0);
    
          CCTK_REAL vx_0 = zx_0/w_0;
          CCTK_REAL vy_0 = zy_0/w_0;
          CCTK_REAL vz_0 = zz_0/w_0;

          switch (bbh_id) {
          case bbh_id_t::Gaussian : {

            rho(p.I) = rho_disk*exp(-pow(p.z,2)/disk_width);
            velx(p.I) = vx_0*exp(-pow(p.z,2)/disk_width);
            vely(p.I) = vy_0*exp(-pow(p.z,2)/disk_width);
            velz(p.I) = vz_0*exp(-pow(p.z,2)/disk_width); 
            press(p.I) = isentropic ? poly_k * pow(rho(p.I),poly_gamma) : press_disk*exp(-pow(p.z,2)/disk_width);
            eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
            break;
          };
          case bbh_id_t::Constant : {


            rho(p.I) = rho_disk;
            velx(p.I) = vx_0;
            vely(p.I) = vy_0;
            velz(p.I) = vz_0; 
            press(p.I) = isentropic ? poly_k * pow(rho(p.I),poly_gamma) : press_disk;
            eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
            break;
          };
          case bbh_id_t::PowerLaw : {

            rho(p.I) = rho_disk * pow(rr + 1e-100, -nrho);
            velx(p.I) = vx_0;
            vely(p.I) = vy_0;
            velz(p.I) = vz_0; 
            press(p.I) = isentropic ? poly_k * pow(rho(p.I),poly_gamma) : press_disk * pow(rr + 1e-100, -npress);
            eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
            break;
          };
          }

  });

  grid.loop_all<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.; });

  grid.loop_all<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = add_magnetic_fields ? B_initial*p.x : 0.; });

  grid.loop_all<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.; });

}

} // namespace AsterSeeds
