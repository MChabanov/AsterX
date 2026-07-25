// Verify the (H, vf2) fusion against the unfused forms in fluxes.cxx /
// eigenvalues.hxx. Exact-arithmetic identity is proven; this bounds the ROUNDOFF
// deviation, which is the number needed to argue the golden rebaseline.
// Build: g++ -O2 -ffp-contract=off -o verify_fusion verify_fusion.cxx

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <algorithm>

static inline double pow2(double x) { return x * x; }

struct Roots { double lp, lm; };

// ---- OLD: GR07 Eq.(28) with b^i=0, as on opt/flux-launch-bounds -------------
static Roots eig_old(double alp, double beta, double u, double v, double rho,
                     double cs2, double W, double h, double bsq) {
  const double a0 = (bsq + cs2 * h * rho) * (pow2(beta) - pow2(alp) * u) -
                    (-1 + cs2) * h * rho * pow2(beta - alp * v) * pow2(W);
  const double a1 = 2 * beta * (bsq + cs2 * h * rho) -
                    2 * (-1 + cs2) * h * rho * (beta - alp * v) * pow2(W);
  const double a2 = bsq + h * rho * (cs2 + pow2(W) - cs2 * pow2(W));
  double det = pow2(a1) - 4 * a2 * a0;
  if (det < 0)
    det = 0;
  return {(-a1 + sqrt(det)) / (2 * a2), (-a1 - sqrt(det)) / (2 * a2)};
}

// ---- NEW: fused form ------------------------------------------------------
static Roots eig_new(double alp, double beta, double u, double v, double vf2,
                     double W) {
  const double K = (1 - vf2) * pow2(W);
  const double a2h = vf2 + K;
  const double rad = u * a2h - K * pow2(v);
  const double disc = alp * sqrt(fmax(0.0, vf2 * rad));
  const double drift = alp * K * v;
  return {-beta + (drift + disc) / a2h, -beta + (drift - disc) / a2h};
}

// ---- REFERENCE: old algebraic form in long double (64-bit mantissa) --------
// Used to decide WHICH double form is closer to truth, i.e. whether the fused
// form really is better conditioned or merely different.
struct RootsL { long double lp, lm; };
static RootsL eig_ref(long double alp, long double beta, long double u,
                      long double v, long double rho, long double cs2,
                      long double W, long double h, long double bsq) {
  const long double a0 = (bsq + cs2 * h * rho) * (beta * beta - alp * alp * u) -
                         (-1 + cs2) * h * rho * (beta - alp * v) *
                             (beta - alp * v) * W * W;
  const long double a1 = 2 * beta * (bsq + cs2 * h * rho) -
                         2 * (-1 + cs2) * h * rho * (beta - alp * v) * W * W;
  const long double a2 = bsq + h * rho * (cs2 + W * W - cs2 * W * W);
  long double det = a1 * a1 - 4 * a2 * a0;
  if (det < 0)
    det = 0;
  return {(-a1 + sqrtl(det)) / (2 * a2), (-a1 - sqrtl(det)) / (2 * a2)};
}

static uint64_t s = 0x2545F4914F6CDD1DULL;
static double urand() { // xorshift, deterministic
  s ^= s << 13; s ^= s >> 7; s ^= s << 17;
  return double(s >> 11) / double(1ULL << 53);
}
static double logu(double lo, double hi) {
  return exp(log(lo) + urand() * (log(hi) - log(lo)));
}

int main() {
  double max_lam = 0, max_Q = 0, max_tau = 0, max_rad_neg = 0;
  double max_tauQ = 0, cancel = 0;
  double err_old = 0, err_new = 0, sum_eo = 0, sum_en = 0;
  long nref = 0, new_better = 0;
  double max_lam_at[9] = {0};
  long clamps_old = 0, clamps_new = 0;
  const long N = 4000000;

  for (long n = 0; n < N; ++n) {
    // geometry
    const double alp = 0.25 + 0.75 * urand();
    const double beta = 0.6 * (urand() - 0.5);
    const double u = logu(0.3, 3.0);
    // fluid: span atmosphere -> NS core
    const double rho = logu(1e-14, 1e-2);
    const double eps = logu(1e-6, 2.0);
    const double p = logu(1e-16, 1e-3);
    const double cs2 = logu(1e-6, 0.6);
    // kinematics
    const double W = 1.0 + logu(1e-8, 6.0);
    const double v2t = 1.0 - 1.0 / pow2(W); // = v_j v^j, the FULL contraction
    // Physical bound on the single component: Cauchy-Schwarz on
    // v^{dir} = delta^{dir}_j v^j gives (v^{dir})^2 <= g^{dir dir} v_j v^j
    // = u * v2t.  Sampling outside this is unphysical (it is what made the
    // clamp fire in both forms on the first run).
    const double v = (2 * urand() - 1) * sqrt(v2t * u);
    // field: include strongly magnetically dominated (b2 >> rho h)
    const double B2 = logu(1e-20, 1e-1);
    const double alp_b0 = (2 * urand() - 1) * sqrt(B2);

    const double h = 1 + eps + p / rho;
    const double rhoh = rho + rho * eps + p;
    const double bsq = (B2 + pow2(alp_b0)) / pow2(W);
    const double H = rhoh + bsq;
    const double vf2 = (bsq + cs2 * rhoh) / H;
    if (!(vf2 >= 0.0 && vf2 < 1.0))
      continue;

    // --- eigenvalues
    const Roots o = eig_old(alp, beta, u, v, rho, cs2, W, h, bsq);
    const Roots f = eig_new(alp, beta, u, v, vf2, W);
    if (!std::isfinite(o.lp) || !std::isfinite(f.lp))
      continue;

    // did either form need its clamp?
    {
      const double a0 = (bsq + cs2 * h * rho) * (pow2(beta) - pow2(alp) * u) -
                        (-1 + cs2) * h * rho * pow2(beta - alp * v) * pow2(W);
      const double a1 = 2 * beta * (bsq + cs2 * h * rho) -
                        2 * (-1 + cs2) * h * rho * (beta - alp * v) * pow2(W);
      const double a2 = bsq + h * rho * (cs2 + pow2(W) - cs2 * pow2(W));
      if (pow2(a1) - 4 * a2 * a0 < 0) ++clamps_old;
      const double K = (1 - vf2) * pow2(W);
      const double rad = u * (vf2 + K) - K * pow2(v);
      if (vf2 * rad < 0) { ++clamps_new; max_rad_neg = fmin(max_rad_neg, rad); }
    }

    // Which double form is closer to the long-double reference?
    const RootsL r = eig_ref(alp, beta, u, v, rho, cs2, W, h, bsq);
    const double scl = fmax(fmax(fabs((double)r.lp), fabs((double)r.lm)), 1e-300);
    const double eo =
        fmax(fabs((double)(o.lp - r.lp)), fabs((double)(o.lm - r.lm))) / scl;
    const double en =
        fmax(fabs((double)(f.lp - r.lp)), fabs((double)(f.lm - r.lm))) / scl;
    err_old = fmax(err_old, eo);
    err_new = fmax(err_new, en);
    sum_eo += eo;
    sum_en += en;
    ++nref;
    if (en < eo) ++new_better;

    const double sc = fmax(fmax(fabs(o.lp), fabs(o.lm)), 1e-300);
    const double dl = fmax(fabs(o.lp - f.lp), fabs(o.lm - f.lm)) / sc;
    if (dl > max_lam) {
      max_lam = dl;
      const double v9[9] = {alp, beta, u, v, rho, cs2, W, bsq, vf2};
      for (int i = 0; i < 9; ++i) max_lam_at[i] = v9[i];
    }

    // --- Q and tau
    const double sqrtg = 0.5 + 2.0 * urand();
    const double dens = sqrtg * rho * W;
    const double Qo = dens * h * W + sqrtg * (pow2(alp_b0) + B2);
    const double Qn = sqrtg * H * pow2(W);
    max_Q = fmax(max_Q, fabs(Qo - Qn) / fmax(fabs(Qo), 1e-300));

    const double ptot = p + 0.5 * bsq;
    const double to = dens * h * W - dens + sqrtg * (B2 - ptot);
    const double tn = Qn - dens - sqrtg * (ptot + pow2(alp_b0));
    max_tau = fmax(max_tau, fabs(to - tn) / fmax(fabs(to), 1e-300));
    // tau is a difference of large like-signed numbers in BOTH forms
    // (dens_h_W - dens vs Q - dens are the same cancellation). Scale the
    // deviation by Q, the natural magnitude, to separate "my reformulation is
    // wrong" from "tau is ill-conditioned in the scheme already".
    max_tauQ = fmax(max_tauQ, fabs(to - tn) / fmax(fabs(Qn), 1e-300));
    cancel = fmax(cancel, fabs(Qn) / fmax(fabs(to), 1e-300));
  }

  printf("samples                     : %ld\n", N);
  printf("max rel dev lambda_pm       : %.3e  (%.1f eps)\n", max_lam,
         max_lam / 2.22e-16);
  printf("   at alp=%.4g beta=%.4g u=%.4g v=%.4g\n", max_lam_at[0],
         max_lam_at[1], max_lam_at[2], max_lam_at[3]);
  printf("      rho=%.4g cs2=%.4g W=%.6g bsq=%.4g vf2=%.6g\n", max_lam_at[4],
         max_lam_at[5], max_lam_at[6], max_lam_at[7], max_lam_at[8]);
  printf("max rel dev Q               : %.3e  (%.1f eps)\n", max_Q,
         max_Q / 2.22e-16);
  printf("max rel dev tau (vs tau)    : %.3e  (%.1f eps)\n", max_tau,
         max_tau / 2.22e-16);
  printf("max dev tau, scaled by Q    : %.3e  (%.1f eps)\n", max_tauQ,
         max_tauQ / 2.22e-16);
  printf("worst tau cancellation Q/tau: %.3e\n", cancel);
  printf("--- accuracy vs long-double reference (%ld samples) ---\n", nref);
  printf("  max rel err OLD form      : %.3e  (%.0f eps)\n", err_old, err_old/2.22e-16);
  printf("  max rel err NEW form      : %.3e  (%.0f eps)\n", err_new, err_new/2.22e-16);
  printf("  mean rel err OLD / NEW    : %.3e / %.3e\n", sum_eo/nref, sum_en/nref);
  printf("  NEW closer to truth in    : %.2f%% of samples\n", 100.0*new_better/nref);
  printf("clamp fired: old det<0      : %ld\n", clamps_old);
  printf("clamp fired: new vf2*rad<0  : %ld  (min rad %.3e)\n", clamps_new,
         max_rad_neg);
  return 0;
}
