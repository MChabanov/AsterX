# AsterX optimization — ideas and future work

Assessed future work. Completed work: `implementation-history.md`. Status/handoff:
`CHECKPOINT.md`.

## STATUS (2026-07-27)

**The flux-kernel performance line of work is essentially done: `AsterX_Fluxes` is
−47.5 % against the original upstream code** (~−27 % fusion + config-gating, ~−28 %
occupancy via scoped `launch_bounds`), on two independent configurations, with
−13…−18 % of total wall clock. Numbers in `baseline-timings.md`; controls still owed
per `CHECKPOINT.md` §STATUS.

Every occupancy route other than `launch_bounds` is closed **on measurement**:
source-level algebraic removal (the `(H,vf2)` fusion), source-level dead-load removal
(`vbar`), and all three compiler-flag mechanisms. See `CHECKPOINT.md` Lessons 4 and 5
— removal only counts at the reconstruction pressure peak, and constraining the
allocator beat coaxing it.

**The flux kernel is no longer the dominant cost** (47.6 % → 32.7 % of Solve on TOV;
37.8 % → 24.2 % on subcycling). Further flux work now competes with
`AnalyticalSpacetimeX_SetMetric` (~29 % of Solve in the subcycling run), `Z4c_RHS`
(83 s, untouched) and `AsterX_SourceTerms` (57 s on TOV). **Pick the next target
deliberately rather than by momentum.**

**The one live item below is `tau` conditioning — an accuracy defect in shipping
code, with no performance claim attached.**

## Context

- **AsterX** is a Cactus/CarpetX thorn on **CarpetX** (`../CarpetX`), itself on
  **AMReX** (`../amrex`). AMReX is reference-only; CarpetX now carries one AsterX-driven
  change (the `min_blocks` parameter, `launch-bounds-plan.md`).
- Done: golden-master harness + FP-determinism CI flags; TimerReport in the test
  suite; ccache CI; flux calculation fused 6 sweeps → 1 kernel, bit-identical;
  diagonal B-fluxes retired; CT-scheme split; `use_pplim` compiled out; scoped
  `launch_bounds` occ-2. Provenance: `implementation-history.md`.

## `tau` conditioning — an ACCURACY fix (⭐ the live item)

**Not an occupancy idea; must not be bundled with one.** Found while measuring the
`(H,vf2)` fusion probe, which was a register null — the numerics finding is
independent and survives.

**The problem, in production code today.** `fluxes.cxx:715` computes
`tau = dens_h_W - dens + sqrtg*(B2 - p_tot)`, i.e.

    tau/sqrtg = rho*W*(h*W - 1) + B^2 - p_tot

For a cold, slow fluid `h*W → 1`, so `rho*W*(h*W-1)` is a small residual of two large
like-signed numbers. Measured over 4e6 samples spanning atmosphere→core
(`Planning/verify-fusion.cxx`): the cancellation factor `|Q/tau|` reaches **5e7**,
i.e. **`tau` loses up to ~8 decimal digits**. The atmosphere and the stellar surface
are exactly where this bites, and `tau` feeds the energy RHS.

**The fix.** With `h = 1 + eps + p/rho` and `W^2 - 1 = W^2 v^2`:

    h*W - 1 = (W - 1) + W*(eps + p/rho),    W - 1 = W^2 v^2 / (W + 1)

so

    tau/sqrtg = rho*W^3*v^2/(W+1) + W^2*(rho*eps + p) + B^2 - p_tot

No subtraction of like-signed quantities remains. Derivation and exact-identity check:
`flux-construction.md` §10b, eqs (10.18)–(10.20).

**Two traps.** Build `rho*eps + p` **directly** — writing it as `rhoh_rc - rho_rc`
reintroduces precisely the cancellation being removed. And `v^2` must come from the
reconstruction branches of §8 where it is cancellation-free (`v2_rc` directly;
`z^2/(1+z^2)`; likewise `s_vec`); it is currently a transient inside the `switch` and
needs plumbing out. **Never** recover it as `(W^2-1)/W^2` — same cancellation again.

**Costs.** Needs `B2_rc` alive plus two new per-side values (`v^2`, `rho*eps + p`):
~+3 doubles/side. Per `CHECKPOINT.md` Lesson 4 that is very likely occupancy-null in
both directions (all downstream of the waist, none of it at the pressure peak), so
treat the register effect as noise and judge this purely on numerics.

**Validation.** Not bit-identical → cannot be golden-gated at 0. The right case is
*conserved quantities*, which is stronger than a golden delta anyway: the gate already
evolves 100 iterations of magnetized TOV + Z4c + AMR (`magTOV_Z4c_AMR.par:154`,
`magTOV_Z4c_AMR_SC.par:157` set `cctk_itlast = 100`), so compare rest-mass conservation
and constraint norms over that run. The claim to demonstrate is that the new form
**improves** them, not merely that it differs.

**Optional companion, shares the one rebaseline.** The `vf2` eigenvalue
reparametrization (`flux-construction.md` §11a) is measurably more accurate — vs a
long-double reference, max rel err 1.06e-11 → 1.36e-13 (78×), mean 4.48e-15 →
1.17e-16 (38×) — and makes the `det < 0` clamp provably dead (0 firings in 4e6
physical samples; the radicand is positive by Cauchy–Schwarz). Zero occupancy benefit,
so it only makes sense folded into the same rebaseline. **Do NOT also take the
fusion's `tau` form** (`Q - dens - sqrtg*(p_tot + alp_b0^2)`): via
`b^2 W^2 = B^2 + alp_b0^2` it subtracts `alp_b0^2` back off an `H*W^2` that already
contains it, a pair the current code cancels analytically — a conditioning
*regression*. Ship `vf2` + the `tau` form above, and leave `H` out of `tau`.

### The accuracy PR — everything needed, so nothing has to be re-derived

**No work started; this is the complete shopping list.** Deliberately kept separate
from the performance PRs, because neither piece is bit-identical and both need the
same rebaseline — so they should share exactly one.

**Code already written and archived:** tag **`archive/vf2-accuracy-probe`** →
`659b48b9` (content-identical to `origin/probe/flux-enthalpy-fusion` @ `ad9b9bd6`).
Its `AsterX/src/eigenvalues.hxx` is a complete, documented `vf2` implementation:
derivation in comments, return type collapsed `vec<vec<REAL,4>,2>` → `vec<vec<REAL,2>,2>`,
clamp demoted to roundoff insurance. **Take that file. Do NOT take that branch's
`fluxes.cxx`** — its `tau` is form (B), the conditioning regression.

**Still to write:** the `tau` form (C) above, which needs `v^2` plumbed out of the §8
reconstruction `switch` (cancellation-free there) and `rho*eps + p` built directly.
Both traps are described above; both are easy to get wrong.

**Validation plan** (this is the part that makes or breaks the PR — a golden delta
alone proves nothing):
1. Golden will **not** be 0. Expect it to move; the case is that it moves *toward*
   truth. For `vf2` that is already quantified against a long-double reference.
2. The real argument is **conserved quantities** over the existing 100-iteration
   magnetized TOV + Z4c + AMR gate (`magTOV_Z4c_AMR.par:154`,
   `magTOV_Z4c_AMR_SC.par:157`): compare rest-mass conservation and constraint norms,
   and show the new forms **improve** them.
3. Rebaseline the golden `.tsv` only after (2) is convincing, and note in the PR that
   the gate exercises accumulated drift, not just one step.
4. Harness for the accuracy numbers is already committed: `Planning/verify-fusion.cxx`
   (4e6-sample sweep vs a long-double reference).

**Expect no register or timing effect** in either direction — everything here sits
downstream of the reconstruction pressure peak (`CHECKPOINT.md` Lesson 4), so judge it
purely on numerics and do not let a `-Rpass` reading influence the decision.

## Retire flux grid functions — blocked as a wholesale change

The flux GFs are a multi-consumer interface: (1) `AsterX_RestrictFluxes` needs them as
face GFs for CarpetX AMR flux restriction (conservation at refinement boundaries);
(2) the EMF consumes them (flux-CT reads the B-fluxes; upwind-CT reads
`vbar`/`amax`/`amin`); (3) with `hydro_correction_order > 2`, `rhs.cxx` reads a
multi-face stencil per direction. Direct-to-RHS also means either scatter-add (GPU
atomics, nondeterministic) or ~2× recompute of reconstruction+Riemann. If revisited:
scope to unigrid + upwind-CT + 2nd-order as a measured experiment.

Subset status (traced 2026-07-11; no test outputs flux GFs):

- ~~Tier 1, diagonal B-fluxes~~ **DONE** (`1bc276a5`).
- ~~Tier 2, CT-scheme-gated retirement~~ **DONE** (`b23fb947`/`e8b0718c`/`4a7e40ec`):
  `use_uct=yes` retires the B-fluxes, `use_uct=no` retires the 12 UCT face GFs.
  Occupancy-null but a real storage/bandwidth win. Detail:
  `profiling-flux-kernel.md` §Step 1.
- **Tier 3 — DYe (non-tab EOS) / DEnt (entropy c2p unused): NOT easy.** They feed
  evolved variables unconditionally (`DYe_rhs`, `DEntrhs` in `rhs.cxx`) and
  Ye/entropy are not identically zero in those configs, so skipping changes results.
  Needs a physics-level decision per config.
- **Never: dens/mom/tau fluxes** — RHS divergence, AMR flux restriction and the PP
  limiter (which re-reads own-face `fluxdenss`/`fluxDYes`) all need them.

## Two per-side flux kernels (a "true L/R split") — MOOT, kept for the reasoning

Idea: launch two kernels, left then right, so each holds only ONE side (~half the live
variables) with constant/no side index, giving occ-2 or even occ-4 **at low scratch** —
the clean version of what Ideas 2/3 attempted. GF reuse avoids new storage: the
numerical flux splits additively once the wavespeed bounds are known (LxF
`F = [0.5·f_L + 0.5·c·U_L] + [0.5·f_R − 0.5·c·U_R]`; HLLE the same with its
coefficients), so kernel A writes its side's contribution to `fluxX(dir_i)` and
kernel B read-modify-adds — one existing GF per conserved var, no atomics.

**Superseded:** scoped `launch_bounds` obtained occ-2 at +160 B/lane for ~90 lines in
one CarpetX header, versus this idea's large config-specific rewrite. Two catches
recorded in case it is ever revived:

- `charmax`/`charmin` must be known to both side-kernels (the per-side parts contain
  `c`), and they are a both-sides reduction → a small eigenvalue pre-step. For UCT
  they are **already stored** as the `amax_*`/`amin_*` face GFs; flux-CT would need a
  home for them.
- Reconstruction is inherently both-sides and `useLO`/velocity-limit couple the sides
  (OR over both), so **each kernel must reconstruct both sides** → ~2× reconstruction,
  ~3× with the pre-step. That trades register pressure for redundant compute, which on
  a latency-bound kernel is an uncertain net — the reason it was never the first choice.

## Deprioritized / shelved

- **Neighbourhood-reuse pass — SHELVED.** Loading each cell's neighbourhood once and
  reusing it across directions would *add* live values to an occ-1 register-bound
  kernel, so it backfires; the lever was register *reduction*, not reuse. The one
  surviving piece — cutting the ~26–45 interleaved write streams — was achieved as a
  side effect of the CT-scheme storage retirement.
- **Fuse the EMF / vector-potential aux sweeps — DEPRIORITIZED.** Good technical fit
  (`AsterX_CalcAuxTermsForAvecPsiRHS` runs 8 sweeps: `CalcE<0,1,2>`, a vertex
  `loop_all`, `CalcFstag<0,1,2>`, a vertex `loop_allm1`; the three edge centerings
  share the integer index space like the face GFs, so the fusion recipe transfers
  directly, and none of these kernels reads `p.X`). But the EMF/Avec sector is cheap
  in the target production configs — on the Frontier TOV runs `CalcAuxTermsForAvecPsiRHS`
  did not even reach the top timers — so it adds review surface for little return.
  Prerequisite if revived: add `CCTK_HOST` to `hll_upwind` (`aster_utils.hxx`), the same
  fix `maxspeeds_from_lambdas`/`avg_upwind` needed. Revisit only if a magnetized
  production profile shows the EMF sweeps hot.

## Why fusion alone showed no CPU win (historical, 2026-07-11)

Step B measured ≈ cost-neutral on CI CPU (fractions ±1 pp). Three reasons, all of
which stop applying on GPU — which is why the same code is −27 % there:

1. **The reads are not actually shared yet.** `CalcFluxAtFace<0/1/2>` are three
   independent inlined bodies, each re-reading its own stencil; the per-direction
   stencils point along different axes so the overlap is only the few central cells,
   and the compiler will not keep those in registers across three ~1000-line bodies
   with interleaved stores. Genuine reuse would need the (now shelved) neighbourhood
   pass.
2. **CI grids are cache-resident**, so "stream inputs once instead of three times"
   saves nothing — sweeps 2 and 3 of the old code hit L2/L3 anyway.
3. **Fusion has real small-box costs**: ~3× instruction footprint, per-point box
   guards, and ~45 interleaved write streams in one loop versus 15 clean streams per
   sweep. Fits the data (large-box shock tubes slightly better, small-box AMR slightly
   worse); Step C2 removed one full zeroing write pass and recovered part of it.

On GPU the launch reduction (6→1) and occupancy dominate instead, and neither is
visible in CI (cuda/rocm are build-only there).

## Constraints and principles

- **Bit-identical results** for pure-refactor steps, gated by the golden check at
  exactly 0 (requires the FP-determinism CI flags). Passes that change op order are
  separate and fraction-validated.
- Small reviewable batches, one concern per commit, `[golden-master]` on the HEAD
  commit of a push that claims bit-identity.
- Verbatim-extraction discipline for refactors: anything that cannot be copied
  verbatim needs explicit sign-off before changing it.
