# AsterX Optimization — Ideas & Future Work

Living scratchpad of assessed future work. Completed work: see
`implementation-history.md`. Current status/handoff: `CHECKPOINT.md`.

## STATUS (2026-07-24)

**The flux-kernel occupancy line of work is DONE: production ideal-gas reached
occ-2** (fusion + CT-split + eig-collapse + Ideas 1/2/3; details in `CHECKPOINT.md`
and `profiling-flux-kernel.md`). The register-reduction ideas below (CT-scheme
split = Tier-2; serialization) are implemented — kept here only as record.
Remaining genuinely-open future work: **retire flux grid functions** (blocked as a
wholesale change; scoped subset possible) and **EMF/Avec sweep fusion**
(deprioritized — cheap sector). Pending: a timed Frontier run to confirm the
occ-2 wall-clock win, then the upstream PR.

**Update 2026-07-25.** Source-level *algebraic* register removal is closed on
measurement (CHECKPOINT Lesson 8: the `(H,vf2)` fusion was a byte-for-byte
`-Rpass` null, because it deleted values the allocator was already rematerializing
and nothing downstream of the flux-assembly waist is in the AGPR overflow set).
Two occupancy levers remain: **compiler flags** (ACTIVE — own checkpoint,
`compiler-flags.md`; attractive because a flag-only winner should pass golden at
exactly 0) and **`launch_bounds`** (fallback, `launch-bounds-plan.md`).

That probe round also surfaced two items with no performance claim attached:
- **`tau` conditioning** — first section below. An accuracy defect in shipping code.
- **the `vbar` dead-code finding** — 12 provably-dead `gf_vels` loads in the UCT
  block that LLVM cannot remove (`0.0 * x` needs fast-math to fold), fixable
  bit-identically with `if constexpr (pplim)`. Recorded in `compiler-flags.md`;
  it is the only genuine REMOVAL candidate still open, and unlike the fusion it
  removes *memory loads that must stay live*, which is the category Lesson 8 says
  the allocator actually banks.

## Context

- **AsterX** is a Cactus/CarpetX thorn, built on **CarpetX** (`../CarpetX`),
  itself on **AMReX** (`../amrex`); those two are reference-only.
- Done so far (details in `implementation-history.md`): golden-master harness
  + FP-determinism CI flags; TimerReport in the test suite; ccache CI;
  **Target 1** — the flux calculation fused from 6 grid sweeps into 1 kernel,
  bit-identical (golden = exactly 0); diagonal B-fluxes retired.

## Why fusion alone shows no CPU win (Step B/C timing analysis, 2026-07-11)

Step B (one fused sweep, bit-identical) measured ≈ cost-neutral on the CI
CPU: shock tubes −0.3..−1.0 pp flux fraction, magTOV AMR +0.7..+1.4 pp,
aggregate flat. Step C2 merged the zeroing into that sweep and did help a
little relative to Step B (aggregate Fluxes −0.39 % np1 / −0.74 % np2; flux
fraction −0.28 pp np1 / −0.39 pp np2), partially recovering the magTOV AMR
overhead, but still leaving the CPU verdict cost-neutral at CI scale. Three
reasons, and what they imply:

1. **The reads are not actually shared yet.** `CalcFluxAtFace<0/1/2>` are
   three independent inlined bodies; each re-reads its own reconstruction
   stencil through the GF pointers. The per-direction stencils point along
   different axes, so the overlap is only the few central cells — and the
   compiler will not keep those in registers across three ~1000-line bodies
   with interleaved stores (aliasing between flux writes and inputs, register
   pressure). Genuine reuse = explicitly load the shared neighbourhood once
   and hand it to all three directions — that is the memory-access pass
   below. Step B builds the structure that enables it; it is not itself the
   optimization.
2. **CI test grids are cache-resident.** "Stream inputs once instead of three
   times" saves DRAM traffic only when the working set exceeds cache; on the
   small CI boxes, sweeps 2 and 3 of the old code hit L2/L3 anyway. The claim
   becomes interesting at production grid sizes.
3. **Fusion has real small-box costs.** ~3× instruction footprint (I-cache,
   unrolling), per-point box guards, and — likely dominant — ~45 interleaved
   write streams (15 face GFs × 3 staggerings) in one loop vs 15 clean
   streams per sweep; CPUs have limited write-combining fill buffers. Fits
   the data: large-box shock tubes slightly better, small-box AMR slightly
   worse. Step C2 removed one full zeroing write pass; it produced only a
   small improvement, consistent with this being a memory/cache-side cleanup
   rather than a CPU algorithmic win on small CI grids.

Expected payoffs remain: (a) **GPU** — 6→1 kernel launches and occupancy,
unmeasured in CI (cuda/rocm are build-only); (b) the
**register-level reuse pass** below, which requires the fused structure.

**(a) is now measured (2026-07-22, OLCF Frontier, 1 node/8 GPUs, production
grids).** The refactor that was cost-neutral on CI CPU is a **~27 % flux-kernel
speedup** on MI250X: `AsterX_Fluxes` −26.5 % (fixed-metric subcycling) / −27.7 %
(TOV dynamic Z4c), flux fraction −6.76/−7.74 pp, total wall time −7.5/−10.6 %.
The flux delta accounts for essentially the whole upstream Solve/rhs/total
improvement; everything else flat within noise. Table in `baseline-timings.md`
§GPU results. The reasons this was invisible on CPU (cache-resident boxes,
write-stream offset) are exactly what stops mattering on GPU, where launch
overhead and occupancy dominate.

**Occupancy follow-up (2026-07-22, full detail in `profiling-flux-kernel.md`).**
The fused flux kernel is **occupancy-bound at occ-1** for the production EOS
(tabulated3d, ideal-gas): VGPR pinned at the 256 HW max + AGPR overflow + ~4 KB
scratch (hybrid EOS are healthy, occ 7–8). occ-1 is pre-existing (pre-fusion was
also occ-1 → fusion cost no occupancy). Forcing occ 1→2 sped the flux kernel
−15.9 %, so register reduction is worth it; the peak is the MHD flux-assembly
live-set (EOS math and the PP block both excised with zero effect). Payoff (b)
below — the neighbourhood-*reuse* pass — is therefore **shelved** (it would ADD
live values to an occ-1 kernel and backfire); the lever is register *reduction*,
starting with the config-gated CT-scheme split (see Tier-2 below, now promoted).

## Two per-side flux kernels (a "true L/R split", reusing the flux GFs) — recorded 2026-07-24

Idea (user): instead of one fused kernel doing both face states, launch **two
kernels** — left (minus) then right (plus) — a true split, so each kernel holds
only ONE side ⇒ **~half the live variables**, and with **constant/no side index**
(not the rolled loop's runtime index) the arrays stay in registers ⇒ occ-2 (maybe
occ-4: half of VGPR 256 ≈ 128 → 512/128) **with LOW scratch**. The clean version of
what Ideas 2/3 attempted — avoids the register→scratch demotion. **Constraint (user):
NO new grid functions — reuse the existing flux GFs by write-then-read.**

**Why the GF reuse works — additive decomposition.** The numerical flux splits
into a left part + a right part *once the wavespeed bounds are known*:
- LxF: `F = [0.5·f_L + 0.5·c·U_L] + [0.5·f_R − 0.5·c·U_R]`, `c = max(charmax,−charmin)`.
- HLLE: same additive split with coeffs `charmax, −charmin, charmax·charmin/(cmax−cmin)`.
So **kernel A writes its side's contribution to `fluxX(dir_i)`, kernel B
read-modify-adds its side's** → one existing flux GF per conserved var, no new
storage, no atomics (sequential kernels, one write each).

**What it needs / the catches:**
- **`charmax`/`charmin` must be known to both side-kernels** (the per-side parts
  contain `c`). They're a per-face both-sides reduction, so a small eigenvalue
  pre-step is needed. For **UCT (production) they're ALREADY stored** as the
  `amax_*`/`amin_*` face GFs (`ap_face=charmax`, `am_face=−charmin`) — reuse those.
  flux-CT would need a home for them (2 GFs) — config-specific.
- **Reconstruction is both-sides + `useLO`/velocity-limit couple the sides** (OR
  over both). Each side-kernel still needs BOTH sides' reconstructed rho/entropy/
  Ye/press to make the low-order fallback decision ⇒ **each kernel reconstructs
  both sides** (using one for its flux) ⇒ **~2× reconstruction** (plus the
  eigenvalue pre-step ⇒ up to ~3×). Reconstruction is a big chunk of the kernel,
  so this is a real **compute** penalty that trades against the occupancy gain.
- Extra kernel launches (2-3 vs 1); large, config-specific refactor; bit-identity
  achievable (values unchanged, just split across kernels + additive combine) but
  a big surface (and the additive `+=` must be the exact same FP ops as the fused
  combine — under `-ffp-contract=off`, group the terms identically).

**Verdict:** the most promising *source-level* route to "occ-2 (even occ-4) without
scratch," and the GF-reuse/additive trick removes the new-storage objection. But it
trades register pressure for **redundant reconstruction compute + launches**, which
on a bandwidth/latency-sensitive kernel is an uncertain net and the biggest effort
of all options. **Try scoped `launch_bounds` first** (far cheaper, occ-2 with the
least HBM traffic); keep this as the big-swing fallback if launch_bounds disappoints,
especially if the eigenvalue pre-step + reconstruction can be shared cheaply.

## `tau` conditioning — an ACCURACY fix, not a performance one (recorded 2026-07-25)

**This is not an occupancy idea and must not be bundled with one.** Found while
measuring the `(H, vf2)` fusion probe (`probe/flux-enthalpy-fusion`, a register
null — see CHECKPOINT Lesson 8); the numerics finding is independent and survives.

**The problem, in production code today.** `fluxes.cxx:715` computes
`tau = dens_h_W - dens + sqrtg*(B2 - p_tot)`, i.e.

    tau/sqrtg = rho*W*(h*W - 1) + B^2 - p_tot

For a cold, slow fluid `h*W -> 1`, so `rho*W*(h*W-1)` is a small residual of two
large like-signed numbers. Measured over 4e6 samples spanning atmosphere→core
(`Planning/verify-fusion.cxx`): the cancellation factor `|Q/tau|` reaches **5e7**,
i.e. **`tau` loses up to ~8 decimal digits**. The atmosphere and the stellar
surface are exactly where this bites, and `tau` feeds the energy RHS.

**The fix.** With `h = 1 + eps + p/rho` and `W^2 - 1 = W^2 v^2`,

    h*W - 1 = (W - 1) + W*(eps + p/rho),    W - 1 = W^2 v^2 / (W + 1)

so

    tau/sqrtg = rho*W^3*v^2/(W+1) + W^2*(rho*eps + p) + B^2 - p_tot

No subtraction of like-signed quantities remains. Derivation and the exact-identity
check: `flux-construction.md` §10b, eqs (10.18)-(10.20).

**Two traps.**
- Build `rho*eps + p` **directly**. Writing it as `rhoh_rc - rho_rc` reintroduces
  precisely the cancellation being removed.
- `v^2` must come from the reconstruction branches of §8, where it is
  cancellation-free (`v2_rc` directly; `z^2/(1+z^2)`; likewise `s_vec`). It is
  currently a transient inside the `switch` and needs plumbing out. **Never**
  recover it as `(W^2-1)/W^2` — same cancellation again.

**Costs.** Needs `B2_rc` alive plus two new per-side values (`v^2`,
`rho*eps + p`): ~+3 doubles/side. Per Lesson 8 that is very likely occupancy-null
in both directions (all downstream of the waist, none of it in the AGPR overflow
set), so treat the register effect as noise and judge this purely on numerics.

**Not bit-identical** → cannot be golden-gated at 0. The validation case is
*conserved quantities*, which is stronger than a golden delta anyway: the gate
already evolves 100 iterations of magnetized TOV + Z4c + AMR
(`magTOV_Z4c_AMR.par:154`, `magTOV_Z4c_AMR_SC.par:157` set `cctk_itlast = 100`), so
compare rest-mass conservation and constraint norms over that run. The claim to
demonstrate is that the new form **improves** them, not merely that it differs.

**Optional companion, shares the one rebaseline.** The `vf2` eigenvalue
reparametrization (§11a) is measurably more accurate — vs a long-double reference,
max rel err 1.06e-11 → 1.36e-13 (78x), mean 4.48e-15 → 1.17e-16 (38x) — and makes
the `det < 0` clamp provably dead (0 firings in 4e6 physical samples; the radicand
`u*a2h - K*v^2` is positive by Cauchy-Schwarz). Zero occupancy benefit, so it only
makes sense folded into the same rebaseline as the `tau` fix. **Do NOT also take
the fusion's `tau` form** (`Q - dens - sqrtg*(p_tot + alp_b0^2)`): via
`b^2 W^2 = B^2 + alp_b0^2` it subtracts `alp_b0^2` back off an `H*W^2` that already
contains it, a pair the current code cancels analytically — a conditioning
*regression*. Ship `vf2` + the `tau` form above, and leave `H` out of `tau`.

## Future direction (not now)

- **Neighbourhood-reuse pass — SHELVED (2026-07-22).** Original idea: load each
  cell's neighbourhood once and reuse across directions, keeping quantities in
  registers/local arrays. Profiling killed this for the production EOS: the flux
  kernel is already occ-1 / register-bound, so *adding* live values would only
  worsen occupancy. Superseded by **register REDUCTION** (opposite direction):
  (1) config-gated CT-scheme split (Tier-2, promoted — below); (2) if needed,
  serialize the two face states in the MHD assembly to cut peak liveness ~37
  doubles for clean occ-2. Both in `profiling-flux-kernel.md` §Plan. The one
  surviving piece of the original idea: cutting the ~26–45 interleaved write
  streams — achieved as a *side effect* of the CT-scheme storage retirement.

- **Retire flux grid functions / write RHS directly? (assessed 2026-07-11)
  Blocked as a wholesale change.** The flux GFs are a multi-consumer
  interface: (1) `AsterX_RestrictFluxes` needs them as face GFs for CarpetX
  AMR flux restriction (conservation at refinement boundaries); (2) the EMF
  consumes them (flux-CT reads the B-fluxes; upwind-CT reads
  `vbar`/`amax`/`amin` face GFs); (3) with `hydro_correction_order > 2`,
  `rhs.cxx` reads a multi-face stencil per direction. Direct-to-RHS also
  means either scatter-add (GPU atomics, nondeterministic) or ~2× recompute
  of reconstruction+Riemann. If revisited: scope to unigrid + upwind-CT +
  2nd-order configs as a measured experiment.

  **Which subset retires easily (traced 2026-07-11; no test outputs flux GFs):**
  - ~~Tier 1 — diagonal B-fluxes~~ **DONE** (`1bc276a5`).
  - **Tier 2 — CT-scheme-gated retirement: PROMOTED to the active next step
    (2026-07-22).** Production **typically uses upwind-CT** (`use_uct=yes`; the
    AnalyticalSpacetime test is UCT, the TOV test was flux-CT) — so this is no
    longer conditional. `CalcFluxAtFace` currently computes BOTH schemes'
    quantities unconditionally (no `use_uct` branch); each config wastes the
    other half. Plan: templatize `CalcFluxAll<EOS, bool use_uct>` with
    `if constexpr` (mirrors `CalcE_impl<i,use_uct>`), dispatched via the existing
    runtime switch — compile-time so registers actually free. Likely
    **bit-identical per config** (dead-code removal) → golden-gate can validate.
    Two sides:
      - `use_uct=yes` (production): retire the B-fluxes (`fluxB_j/k`) — also drops
        `Es_rc`/`flux_Btildes`/`Btildes_rc` (~12–18 doubles → occupancy win toward
        occ-2). Storage retirement needs an `interface.ccl` split (off-diagonal
        B-fluxes are mixed into the `flux_x/y/z` groups with the hydro fluxes).
      - `use_uct=no` (flux-CT): retire the 12 UCT face GFs (`vbar_*`×6, `a_*`×6),
        already in 6 dedicated groups → trivial conditional STORAGE. Wins: −12
        GFs memory, ~halves face-GF write streams, modest occupancy.
    Full recipe/rationale in `profiling-flux-kernel.md` §Plan Step 1.
  - **Tier 3 — DYe (non-tab EOS) / DEnt (entropy c2p unused): NOT easy.**
    They feed evolved variables unconditionally (`DYe_rhs`, `DEntrhs` in
    `rhs.cxx`) and Ye/entropy are not identically zero in those configs —
    skipping changes results; needs a physics-level decision per config.
  - **Never: dens/mom/tau fluxes** — RHS divergence, AMR flux restriction,
    and the PP limiter (re-reads own-face `fluxdenss`/`fluxDYes`) need them.

- **Fuse the EMF / vector-potential aux sweeps (assessed 2026-07-11) —
  DEPRIORITIZED 2026-07-22.** Good technical fit, but the user reports the
  EMF/Avec sector is very cheap in the target production configs (confirmed by
  the Frontier TOV runs, where `CalcAuxTermsForAvecPsiRHS` did not even reach
  the top timers), so fusing it adds review surface for little GPU return.
  Kept below for reference; revisit only if a magnetized production profile
  shows the EMF sweeps hot.** `AsterX_CalcAuxTermsForAvecPsiRHS`
  (fluxes.cxx) runs **8 sweeps**: `CalcE<0,1,2>` (edge-centred
  `loop_int_device`; upwind-CT reconstructs `dB_stag`/`vbar` along the two
  transverse directions, flux-CT averages four B-fluxes), one vertex
  `loop_all` (Fbeta, G init), `CalcFstag<0,1,2>` (edge-centred
  `loop_mix_device`; `alp·sqrtg·A^i`), one vertex `loop_allm1` (gauge scalar
  G). The three edge centerings share the integer index space like the face
  GFs, so the Target-1 recipe transfers directly: box helpers replicating
  `loop_int`/`loop_mix`, one union sweep, per-direction guards, compile-time
  direction via a templated per-edge helper, golden gate. **Simpler than the
  flux case: none of these kernels reads `p.X`** (ReconX is index-only).
  Prerequisite: add `CCTK_HOST` to `hll_upwind` (aster_utils.hxx) — same fix
  as `maxspeeds_from_lambdas`/`avg_upwind` needed in Step A. Staging: (a)
  fuse `CalcE` 3→1 and `CalcFstag` 3→1; (b) optionally merge the vertex
  gauge loops (mind intra-routine dependencies). Payoff: top-4 kernel,
  ≈4–10 % of evolution in the MHD tests; many small kernels → launch overhead
  matters relatively more. Metric:
  `AsterX_CalcAuxTermsForAvecPsiRHS ÷ ODESolvers_Solve*`.

## Constraints & principles

- **Bit-identical results** for pure-refactor steps, gated by the golden
  check at exactly 0 (requires the FP-determinism CI flags — see CHECKPOINT).
  Performance passes that change op order are separate, fraction-validated.
- Small reviewable batches; one concern per commit; `[golden-master]` tag on
  the HEAD commit of a push that claims bit-identity.
- Verbatim-extraction discipline for refactors: anything that cannot be
  copied verbatim needs explicit sign-off before changing it.
