#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>

#include "setup_eos.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class ns_perturbation_t { m0, m1, m2 };

template <typename EOSType>
void AsterSeeds_PerturbNS_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_PerturbNS;
  DECLARE_CCTK_PARAMETERS;

  ns_perturbation_t perturbation_mode;

  if (CCTK_EQUALS(ns_perturbation_mode, "m0")) {
    perturbation_mode = ns_perturbation_t::m0;
    if (ns_perturbation_strength < 0.0 ||
        ns_perturbation_strength >= 1.0) {
      CCTK_VERROR(
          "AsterSeeds::ns_perturbation_strength must be in [0, 1) for m0, "
          "but is %.17g",
          static_cast<double>(ns_perturbation_strength));
    }
  } else if (CCTK_EQUALS(ns_perturbation_mode, "m1")) {
    perturbation_mode = ns_perturbation_t::m1;
  } else if (CCTK_EQUALS(ns_perturbation_mode, "m2")) {
    perturbation_mode = ns_perturbation_t::m2;
  } else {
    CCTK_VERROR("Unknown AsterSeeds::ns_perturbation_mode '%s'",
                ns_perturbation_mode);
  }

  if ((perturbation_mode == ns_perturbation_t::m1 ||
       perturbation_mode == ns_perturbation_t::m2) &&
      ns_perturbation_Req <= 0.0) {
    CCTK_VERROR("AsterSeeds::ns_perturbation_Req must be positive, but is %.17g",
                static_cast<double>(ns_perturbation_Req));
  }

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CCTK_REAL rho_atm =
            (radial_distance > r_atmo)
                ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                : rho_abs_min;
        rho_atm = max(eos_3p->rgrho.min, rho_atm);
        const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

        CCTK_REAL rhoL = rho(p.I);
        if (rhoL <= rho_atmo_cut) {
          return;
        }

        CCTK_REAL pressL = press(p.I);
        CCTK_REAL epsL = eps(p.I);
        CCTK_REAL tempL = temperature(p.I);
        CCTK_REAL entropyL = entropy(p.I);
        const CCTK_REAL YeL = Ye(p.I);

        // Perturbations from Tsokaros et al.,
        // "Dynamically stable ergostars exist: General relativistic models and
        // simulations", Phys. Rev. Lett. 123, 231103 (2019),
        // arXiv:1907.03765.
        switch (perturbation_mode) {
        case ns_perturbation_t::m0: {
          pressL *= (1.0 - ns_perturbation_strength);
          tempL = eos_3p->temp_from_rho_press_ye(rhoL, pressL, YeL);
          epsL = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          entropyL = eos_3p->kappa_from_rho_eps_ye(rhoL, epsL, YeL);
          break;
        }
        case ns_perturbation_t::m1: {
          const CCTK_REAL inv_Req = 1.0 / ns_perturbation_Req;
          const CCTK_REAL factor =
              1.0 + ns_perturbation_strength * (p.x + p.y) * inv_Req;
          rhoL = max(rhoL * factor, rho_atmo_cut);
          pressL = eos_3p->press_from_rho_temp_ye(rhoL, tempL, YeL);
          epsL = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          entropyL = eos_3p->kappa_from_rho_eps_ye(rhoL, epsL, YeL);
          break;
        }
        case ns_perturbation_t::m2: {
          const CCTK_REAL inv_Req = 1.0 / ns_perturbation_Req;
          const CCTK_REAL factor =
              1.0 + ns_perturbation_strength * (p.x * p.x - p.y * p.y) *
                        inv_Req * inv_Req;
          rhoL = max(rhoL * factor, rho_atmo_cut);
          pressL = eos_3p->press_from_rho_temp_ye(rhoL, tempL, YeL);
          epsL = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          entropyL = eos_3p->kappa_from_rho_eps_ye(rhoL, epsL, YeL);
          break;
        }
        default:
          assert(0);
        }

        rho(p.I) = rhoL;
        press(p.I) = pressL;
        eps(p.I) = epsL;
        temperature(p.I) = tempL;
        entropy(p.I) = entropyL;
      });
}

extern "C" void AsterSeeds_PerturbNS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_PerturbNS;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
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

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    // Get local eos object
    auto eos_3p_ig = global_eos_3p_ig;

    AsterSeeds_PerturbNS_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      // pwpoly cold + pwpoly-hybrid
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;

      AsterSeeds_PerturbNS_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);

    } else if (global_eos_3p_hyb_poly) {
      // poly cold + poly-hybrid
      auto eos_3p_hyb = global_eos_3p_hyb_poly;

      AsterSeeds_PerturbNS_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);

    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }

    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos object
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    AsterSeeds_PerturbNS_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace AsterSeeds
