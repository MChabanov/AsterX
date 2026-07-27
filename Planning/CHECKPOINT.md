# Checkpoint — AsterX flux-loop optimization (handoff)

_Last updated 2026-07-27. Read this first, then **`launch-bounds-plan.md`** (the
ACTIVE occupancy route), `compiler-flags.md` (CLOSED on measurement — kept as the
record of three probes that did not move occupancy, plus the `vbar` finding and
the fast-math analysis), `profiling-flux-kernel.md`
(all -Rpass tables + the scratch root-cause), `baseline-timings.md` (CPU/GPU
timing incl. the occ-2 sign-flip), `flux-construction.md` (the numerical pipeline;
§10a/§10b/§11a/§11b hold the fusion + `tau` findings), `implementation-history.md`
(provenance), `ideas.md` (future work). `ideas-2-3-serialization.md` = the
(executed) serialization design note._

## Mission

Optimize `AsterX/src/fluxes.cxx` (the fused flux kernel `CalcFluxAll`), starting
from bit-identity-gated refactors (golden-master = exactly 0) and now
register/occupancy work measured on Frontier/MI250X. CI CPU is cost-neutral by
design; the payoff is GPU.

## STATUS: occ-2 is reachable but carries a scratch tradeoff — one route left (`launch_bounds`)

The production ideal-gas flux kernel was occupancy-bound at **occ-1**. A
source-only serialization stack (Ideas 1/2/3) reached **occ-2**, but the way it
got there (rolled loops) traded on-chip registers for **off-chip scratch**, and
the timing verdict is **grid-size-dependent**:

- **TOV Z4c, flux-CT, larger grid (ideal-gas): occ-2 WON, AsterX_Fluxes −4.8%**
  (163.3 vs 171.6 s), −1.34 pp of Solve, Z4c untouched (flat).
- **fixed-metric subcycling, UCT, smaller grid (ideal-gas): occ-2 REGRESSED +8.9%**
  (137.2 vs 125.9 s), +1.82 pp.

So occ-2 helps where there's enough memory latency to hide (large grid) and hurts
where the extra scratch bandwidth dominates (small grid). And the *genuine* occ-2
(−4.8%) is far short of the *forced* occ-2 (−15.9%, earlier global-launch_bounds
test on the same TOV run) — that gap is the scratch. **=> the lever is to reach
occ-2 WITHOUT the serialization scratch.**

### Where that stands as of 2026-07-27 — two of three routes are now closed

Three routes to occ-2-at-low-scratch, in the order they were tried. **Only #3 is
still open**, and both closures are measurements, not guesses:

1. **Source removal — CLOSED on measurement, twice.** (a) The `(H,vf2)` algebraic
   fusion (`flux-construction.md` §10a/§11a) was a byte-for-byte `-Rpass` null —
   Lesson 8: it deleted values the allocator was already rematerializing. Branch
   `probe/flux-enthalpy-fusion` kept as the record. (b) The `vbar` fix
   (2026-07-27) removed 12 provably-dead *memory loads*, the category Lesson 8
   said the allocator banks, and was **also byte-for-byte null** — Lesson 9: the
   register count is peak simultaneous liveness, and the UCT epilogue is past the
   peak, so position beats category. **No removal candidates remain**: everything
   at the peak is reconstruction-side, which is separately blocked.
2. **Compiler flags — CLOSED on measurement, 2026-07-27.** Own doc:
   **`compiler-flags.md`** (§SWEEP RESULTS + §VERDICT). All three mechanisms run
   against the target `uct1` (needs AGPR 32→0):
   scheduling `--misched=gcn-max-occupancy` = **NULL** (byte-for-byte identical in
   all 11 kernels — because that strategy is already AMDGPU's *default* machine
   scheduler, so the flag re-selected what was running);
   allocation `--enable-deferred-spilling` = **AGPR 32→30**, occ-1, and it made
   `CalcE uct1` SGPR spills worse (20/22/24 → 30/26/33) → reject;
   relaxed FP `-ffinite-math-only` = **AGPR 32→34**, occ-1 → negative.
   Spread across three independent mechanisms: ±2 registers against a 32-register
   gap. Probe 2 doubled as the positive control (counters moved), so probe 1's null
   is a real codegen null, not a broken build knob.
   Two durable findings: the 32 AGPRs are **load-bearing live values, not allocator
   sloppiness** (the flag that re-ran allocation found 2; the flag that relaxed
   semantics found −2) — the compiler-side confirmation of Lesson 8; and
   **TU-scoped flags cannot be aimed at the flux kernel alone** — all 11 kernels
   live in `fluxes.cxx`, so collateral is unavoidable, which is a *structural*
   advantage for `launch_bounds`, not merely a strength-of-knob one.
   ⚠ Probe 3 did not actually test the `vbar` finding: folding `0.0 * x` needs
   `nsz` (`-fno-signed-zeros`), which `-ffinite-math-only` does not set. Resolve by
   writing the fix and measuring it, not by another flag build.
3. **`launch_bounds` plumbing — ⭐ BUILT AND MEASURED 2026-07-27: occ-2 REACHED
   AT LOW SCRATCH.** `launch-bounds-plan.md` §TRIAL RESULT owns the numbers.
   Production `uct1` went **256/32/3448/occ-1 → 128/128/3608/occ-2** with both
   spill counters still 0, and **all 9 EMF kernels byte-for-byte unchanged** (they
   keep `MB=0` and the ordinary `ParallelFor` path). Scratch cost +160 B/lane
   against the serialization's +1536 for the same occupancy — the Lesson-3(b)
   hypothesis, confirmed. Implementation: CarpetX branch
   `opt/loop-box-device-min-blocks` @ `19267243` (pushed, `MChabanov/CarpetX`),
   AsterX side `c96df5d5` (**breaks CI by design**, needs the forked CarpetX).
   The two plan corrections held (the replication calls
   `amrex::detail::call_f_intvect_handler`; the 7th template argument does break
   stock-CarpetX builds).
   **Open: TIMING at both grid sizes — the only remaining question.** Note the
   likely sign flip vs the serialization: TOV-large is flux-CT (`uct0`, which the
   reduced build says was already occ-2, so possibly nothing to gain and +48
   B/lane to lose), while UCT-small is the kernel that actually gained a wave.

Two by-products worth more than the occupancy hunt so far, both in
`flux-construction.md`: **`tau` loses up to ~8 digits to cancellation in
production today** (measured `|Q/tau| ~ 5e7`, §10b has a cancellation-free
rewrite), and the **`vbar` dead-code finding** (`compiler-flags.md`) — 12
provably-dead `gf_vels` loads that LLVM cannot remove without fast-math, fixable
bit-identically with `if constexpr (pplim)`.

## Branches (on `origin` = MChabanov/AsterX)

- **`opt/flux-eig-collapse` (HEAD `210a013c`)** — fusion + CT-split + eig-collapse
  + Idea 1 (`use_pplim` template) + Idea 2 (flux-by-flux) + Idea 3 (face-state
  serialize). **occ-2 via serialization, all golden PASS at 0**, but scratch
  4984 B/lane; large-grid win / small-grid regress. Kept as the serialization
  record and the large-grid data point.
- **`opt/flux-launch-bounds` (HEAD `defd9f60`) — the ACTIVE branch.** Reset to the
  Idea-1 state (fusion + CT-split + eig-collapse + `use_pplim` template + theta
  storage gating), i.e. **before Ideas 2/3**: occ-1 but **low scratch 3448**,
  unrolled/constant-index code. This is the clean base for the `launch_bounds`
  experiment. (`eigenvalues.hxx` here is the both-sides-only original; no
  `eigenvalues_oneside` — that only exists on `opt/flux-eig-collapse`.)
- **`probe/flux-enthalpy-fusion` (HEAD `659b48b9`) — MEASURED NULL, keep as record.**
  Branched off `opt/flux-launch-bounds`. Algebraic removal via the `(H, vf2)`
  waist (`flux-construction.md` §10a/§11a/§11b): total enthalpy `H = rho*h + b^2`
  and fast speed `vf2 = cA^2 + cs^2(1-cA^2)` absorb the whole magnetic sector;
  `cs2_rc`, `h_rc`, `dens_h_W_rc`, `dens_h_W_plus_sqrtg_W2b2_rc` deleted,
  `B2_rc`/`bsq_rc` demoted to transients, eigenvalue signature 15→9 doubles.
  **-Rpass: 256/40/3448/occ-1 — IDENTICAL to the Idea-1 baseline in all four
  counters.** Not bit-identical, so never golden-gated. Do not retry this shape;
  see Lesson 8. Its *accuracy* findings are real and live on in `ideas.md`.
- **`4f902ca1` "Only temporary trick for testing"** — sits on top of `defd9f60` on
  `opt/flux-launch-bounds`. Adds the INERT `ASTERX_PROBE_IDEALGAS_ONLY` guards to
  `fluxes.cxx`; with `-DASTERX_PROBE_IDEALGAS_ONLY` (set only in the uncommitted
  `make.code.deps`) the CalcFluxAll instantiation matrix drops 16 → 2 (idealgas,
  `pplim=0`, both CT schemes) for ~8x faster probe rebuilds. **NOT codegen-neutral:
  AGPR reads 32 in the reduced build vs 40 in the full one** — see
  `compiler-flags.md`. Aborts at runtime on pplim/hybrid/tabulated. Revert (or
  leave inert) before upstreaming.

Commit stack on `opt/flux-eig-collapse`: `b23fb947`/`e8b0718c`/`4a7e40ec` (CT
split compute/hoist/storage), `549a28ca` (eig-collapse), `defd9f60` (Idea 1),
`e8613114` (Idea 2), `210a013c` (Idea 3). **Idea 4** (hoist eig+UCT ahead of the
assembly) was tried, golden-PASS but occupancy-null (pure reorder), and REVERTED.

## The -Rpass progression (production `CalcFluxAll<uct=1,pplim=0,idealgas>`)

| stage | VGPR | AGPR | scratch B/lane | occ |
|---|---:|---:|---:|---:|
| eig-collapse (pre-Idea1) | 256 | ~64 | ~3960 | 1 |
| Idea 1 (template use_pplim) | 256 | 40 | 3448 | 1 |
| Idea 2 (flux-by-flux `for(j)`) | 254 | 1 | 3816 | 1 |
| Idea 3 (face `for(f)`) | 253 | 0 | 4984 | 2 |
| **(H,vf2) fusion, off Idea 1** | **256** | **40** | **3448** | **1** |

Reduced build (`-DASTERX_PROBE_IDEALGAS_ONLY`, AGPR 32 not 40 — a *different*
compilation, comparable only within this block):

| stage | VGPR | AGPR | scratch B/lane | occ |
|---|---:|---:|---:|---:|
| Idea-1 baseline | 256 | 32 | 3448 | 1 |
| `--misched=gcn-max-occupancy` | 256 | 32 | 3448 | 1 |
| `--enable-deferred-spilling` | 256 | 30 | 3448 | 1 |
| `-ffinite-math-only` | 256 | 34 | 3448 | 1 |
| **`vbar` fix (12 dead loads removed)** | **256** | **32** | **3448** | **1** |
| **`launch_bounds(256,2)` — MB=2** | **128** | **128** | **3608** | **2** ⭐ |

Four independent attempts to *coax* the allocator — three compiler-side, one
source-side — all landed within ±2 of a 32-register gap. **Constraining it
instead worked on the first try**: total pressure 288 → 256 (exactly the occ-2
budget), occ 1 → 2, and scratch up only 160 B/lane versus the serialization's
+1536. The 128/128 split with both spill counters at 0 means the allocator used
AGPRs as *on-chip* spill space rather than scratch memory — the best available
outcome. Detail and caveats: `launch-bounds-plan.md` §TRIAL RESULT.

The fusion row is the Idea-1 row **byte for byte** — a true null, not a small
move (2026-07-25, `probe/flux-enthalpy-fusion`). Staleness ruled out: `-Rpass`
remarks are emitted at compile time, so getting output at all proves
`fluxes.cxx` recompiled.

flux-CT+noPP idealgas also hit occ-2 (251/0/4984). PP-on idealgas stays occ-1
(28-29 AGPR; not production). tabulated3d stays occ-1 (AGPR ~172-244,
EOS-table-bound — occ-2 unreachable there by any flux change). hybrid 54/0/occ8.

## KEY LESSONS (the accumulated knowledge — most important for a new agent)

1. **The register allocator (cce/gfx90a) banks REMOVAL and STRUCTURAL change, not
   REORDERING.** REMOVAL (Idea 1, PP compiled out) → AGPR 64→40 AND scratch
   3976→3448 (genuine, both down). Rolled loops (Ideas 2/3) → AGPR 40→0. Pure
   reorders (CT-hoist, eig-collapse, Idea-4) → occupancy-null. Don't waste effort
   on reorders / hand-rematerialization. **⚠ See Lesson 8 — "removal" as stated
   here is too coarse and cost one wasted probe.**
2. **occ-2 gauge = AGPR → exactly 0** (VGPR welded at the 256 ceiling; reducible
   overflow is AGPR). occ = floor(512/(VGPR+AGPR)), granule-rounded (VGPR granule
   8). At AGPR 0, VGPR 253→256 → 512/256 = 2 waves. One residual AGPR granule
   blocks it (we sat at 254/1/occ-1 one granule short until Idea 3).
3. **⚠ THE CENTRAL INSIGHT — serialization RELOCATES, it doesn't REMOVE.** The
   rolled loops (`#pragma unroll 1`) make the loop indices `f`/`j` RUNTIME
   variables. **GPU registers cannot be dynamically indexed**, so any array read
   with a runtime index — the reconstructed vectors (`rho_rc`, `vels_rc`, `Bs_rc`,
   `vlows_rc`, …) and the cut-set accumulators (`moms_rc`, `flux_moms`, `dens_rc`,
   …) — is demoted from registers to **scratch (off-chip HBM)**. So AGPR↓ and
   scratch↑ are the SAME data moving register→scratch, not shrinking. That is why
   scratch rose 3448→3816→4984 (Idea 3 = +1168, the bulk) while AGPR fell. The
   *reconstruction arrays* (demoted by the runtime-`f` reads) dominate the scratch;
   the cut-set is a fraction. Unrolling would put them back in registers (constant
   index) but then both sides coexist → back to occ-1. **The seesaw:**
   unrolled = registers/occ-1/low-scratch ; rolled = scratch/occ-2/high-scratch —
   the data has to live somewhere. The only way to get occ-2 with LOW scratch is
   (a) genuine removal (exhausted) or (b) `launch_bounds` (force occ-2 on the
   unrolled kernel, letting the allocator spill SELECTIVELY — hot data stays in
   registers, only cold values spill). **"(exhausted)" was an assumption when
   written; it was retested 2026-07-25 (the `(H,vf2)` fusion) and is now CLOSED on
   measurement — with a REASON, in Lesson 8.**

   **⚠ CORRECTION 2026-07-25 — the baseline 3448 is `alloca`, NOT spill, and
   scratch does NOT cap occupancy.** The `-Rpass` block reports spill counts and
   `ScratchSize` as SEPARATE fields, and the baseline flux kernel reads
   **scratch 3448 with `VGPRs Spill: 0` AND `SGPRs Spill: 0`**. If VGPRs were
   being spilled to scratch, `VGPRs Spill` would be non-zero. So those 3448
   B/lane (~431 doubles/lane) are private-memory ALLOCATIONS, almost certainly
   the ReconX stencil arrays — runtime-indexed loops over the 6-point stencil
   defeat SROA and force them to private memory. Two things follow:

   - **The seesaw is real but SMALLER than "the data has to live somewhere"
     implies.** A large part of the reconstruction state was ALREADY in scratch
     before Ideas 2/3; the rolled loops added ~1536 B/lane *on top of* an
     existing 3448 rather than moving the bulk of it there. Read the
     3448→3816→4984 progression as an increment on a pre-existing floor.
   - **Occupancy and scratch are INDEPENDENT levers.** AMD occupancy is set by
     VGPR+AGPR (and LDS), not private memory. Reaching occ-2 is purely the
     "clear the 256-register cliff" problem; scratch is a latency/bandwidth cost
     that is what made the small grid regress (+8.9%). So "beat 4984" is a
     TIMING criterion, not an occupancy one — do not conflate them.

   Comparison point (reduced build, 2026-07-25): `CalcE uct1` sits at 72 B/lane
   scratch WITH 20-24 SGPR spills at occ 5 — i.e. that kernel's scratch really is
   spill. The flux kernel's is not.
4. **Production is upwind-CT (`use_uct=yes`), `use_pplim=no`, ideal-gas +
   tabulated3d.** Default `use_pplim=no` and no test par overrides it, so the
   golden gate exercises the production `pplim=false` kernel; `theta_x/y/z` +
   `theta_tot` are not in the golden .tsv set (made the theta storage-gating safe).
5. **Bit-identity of serialization** holds under `-ffp-contract=off` (rolled loops
   reorder independent sub-computations, never reassociate within an expression;
   per-side contractions reuse the SAME scalar `calc_contraction` overloads in
   `AsterUtils/src/aster_utils.hxx`; eigenvalue sides independent). But the
   `CCTK_DEBUG` NaN-dump on `opt/flux-eig-collapse` is disabled (`#if 0`) because
   it references now-serialized quantities — rewrite before upstreaming.
6. **NVIDIA H100/GH100 has no AGPR** (unified 64K-reg file, 255/thread, spill to
   cached local mem). AGPR-drain framing is AMD-only; on Hopper measure
   `registers/thread` + spill bytes via `--resource-usage`. Detail: `GPUHardwareDict.md`.
7. **My hand double-counts were unreliable** (predicted Idea-2 40→16, got 40→1;
   predicted Idea-3 net ~0-6, got occ-2). Trust `-Rpass` + timing, not arithmetic.
8. **⚠ REMOVAL ONLY COUNTS IF THE DATA IS RESIDENT IN THE OVERFLOW SET.**
   **⚠ SUPERSEDED IN PART BY LESSON 9 — read that first.** This lesson's
   load-vs-rematerializable-arithmetic distinction turned out not to be the
   operative one; *position relative to the pressure peak* is. Lesson 8's
   conclusions about the fusion remain correct, but its prediction that removing
   non-rematerializable memory loads would work was tested and failed. Sharpens
   Lesson 1, learned the expensive way from the null in the `-Rpass` table above.
   What the allocator banks is removal of **long-lived values that are expensive to
   REMATERIALIZE**. Removing a source-level value that is short-lived, or cheaply
   recomputable from operands that stay live anyway, is indistinguishable from a
   reorder → exactly zero. The fusion deleted `h_rc` (`1+eps+p/rho`: two flops
   from live inputs), `B2_rc`/`bsq_rc` (contractions of live `Bs_rc`/`Blows_rc`/
   `alp_b0_rc`/`W`), `dens_h_W_rc` and `dens_h_W_plus_sqrtg_W2b2_rc` (products of
   live values), and the inlined eigenvalue transients `a_m/a_p/det`. The
   allocator was **already rematerializing every one of those**; only `cs2_rc` (an
   EOS call) was expensive, and that is 2 doubles. A count of NAMES IN THE SOURCE
   is not a count of VALUES IN REGISTERS.
   **Corollary — where the 40 AGPRs actually are.** Everything the fusion touched
   is downstream of reconstruction. Per Lesson 3 the overflow set is the
   *reconstruction* state (S2–S6, longest and most-interfering live ranges); Idea 1
   moved AGPR 64→40 because compiling out `pplim` deleted a duplicated copy of the
   whole S5–S9 chain *including its own reconstruction state*, not because it
   deleted a dozen doubles of algebra. **Nothing downstream of the §10a waist is in
   the overflow set, so no amount of algebraic fusion there will move AGPR.**
   Route (a) of Lesson 3 is now genuinely CLOSED: the only removable data that
   matters is reconstruction-side, and that path is separately blocked (both-sides
   `reconstruct()`, `useLO` OR-over-both-sides, ReconX risk — see "Ideas considered
   and set aside"). Compiler flags (`compiler-flags.md`) and `launch_bounds`
   (route b) are the remaining levers.
   **⚠ Caveat on the corollary, 2026-07-25.** Per the Lesson-3 correction the
   baseline 3448 B/lane is `alloca`, not spill — so the ReconX stencil arrays are
   ALREADY in private memory, not in the AGPR overflow. That means
   reconstruction-side *data removal* would cut scratch traffic (the small-grid
   regression lever) more than it cuts registers (the occupancy lever). Those may
   be two separate wins needing two separate changes. Before spending effort on
   either, dump the ISA and look: `v_accvgpr_read/write_b32` shows what is truly
   AGPR-resident, `scratch_load/store` shows what is in private memory. Stop
   inferring from counters — that is what produced the two wrong predictions.

9. **⚠ REGISTERS ARE SET BY PEAK SIMULTANEOUS LIVENESS — POSITION BEATS
   CATEGORY.** The final sharpening of Lessons 1/8, measured 2026-07-27 by the
   `vbar` fix. That change deleted **12 provably-dead global loads** from the
   production kernel (`if constexpr (pplim)` around the UCT drift blend, so the
   cell-centred fallback is compiled out instead of multiplied by zero) and the
   `-Rpass` result was **byte-for-byte identical in all 11 kernels** —
   256/32/3448/occ-1 unchanged.
   Lesson 8 predicted this would work: unlike the `(H,vf2)` fusion's
   rematerializable algebra, these were *memory loads that cannot be
   rematerialized*, "exactly the category the allocator banks". **That prediction
   was wrong, and the load-vs-arithmetic distinction is not the operative one.**
   AGPR count = the maximum number of values live *at one point*. The pressure
   peak is the reconstruction state (S2–S6). The UCT epilogue runs *after* flux
   assembly, i.e. after the peak, so values born and consumed there occupy
   registers that are already free — removing them lowers total work, never the
   peak, and therefore never the register count.
   **Operational rule: before proposing any removal, ask WHERE it lives relative
   to the reconstruction peak, not what kind of value it is.** Everything
   downstream of the flux-assembly waist — algebra (Lesson 8) *and* memory loads
   (this lesson) — is register-null. That is now measured twice from two
   different directions.
   **Consequence: route (a) of Lesson 3, "genuine removal", is CLOSED with no
   candidates remaining.** The `vbar` fix was the last one. Every removable thing
   at the peak is reconstruction-side, and that path is separately blocked
   (both-sides `reconstruct()`, `useLO` OR-over-both-sides, ReconX risk — see
   "Ideas considered and set aside"). Combined with the compiler-flag route
   closing the same day, **`launch_bounds` is the only remaining occupancy lever.**
   Corollary for the *scratch* problem (the small-grid regression, a separate
   criterion per the Lesson-3 correction): reconstruction-side data removal would
   cut the 3448 B/lane of `alloca`, and that is a different change from anything
   tried so far. Dump the ISA before guessing — `v_accvgpr_read/write_b32` shows
   what is truly AGPR-resident, `scratch_load/store` what is in private memory.

## Ideas considered and set aside (do not re-derive from scratch)

- **Per-side (incremental) Riemann combine** — accumulate `fsum=flux0+flux1`,
  `vacc=var1-var0` per side, apply `c`/`0.5` post-loop. Bit-identical for LxF
  (NOT the naive `0.5*flux0+0.5*flux1` — that reassociates); HLLE saves less
  (per-side coeffs `charmax/-charmin` don't separate). Trims the both-sides
  cut-set ~7-14 doubles (~15% of the scratch increase) but NOT the dominant
  reconstruction-array scratch → partial, not worth the redo.
- **Grade the atmosphere once at the face** (instead of per cell-center) — a
  physically-defensible SCHEME change (non-bit-identical → needs maintainer
  sign-off + physics validation, not golden=0); modest compute win (halves the
  atmosphere `pow`s). Out of the bit-identical track.
- **Full per-side kernel rewrite touching reconstruction** — blocked by (a) the
  `useLO`/velocity-limit fallback being an OR over BOTH sides (reconstruction is
  inherently both-sides; `reconstruct()` returns both), (b) constant-index
  two-block unroll relies on the allocator NOT interleaving (likely occ-null, à la
  Idea 4) unless forced by `noinline` (a dead probe), + huge bit-identity/ReconX
  risk. Not recommended.

## FP-determinism prerequisite (do not undo)

CI cpu config compiles C++ with `-ffp-contract=off`, no `-funsafe-math` (`scripts/
actions-cpu-real64.cfg`); golden `.tsv` regenerated under those flags. Changing
them invalidates golden + baseline.

## How to operate

- **Golden check:** `[golden-master]` in the HEAD commit message of the push. Pass
  = `2235 ok, 0 fail`, worst `|abs|=|rel|=0`, np1+np2, cpu/rocm/cuda. Watch
  `gh run watch <id> --exit-status`. All Ideas 1/2/3 passed.
- **-Rpass (user drives on Frontier):** `fluxes.cxx.o: CXXFLAGS +=
  -Rpass-analysis=kernel-resource-usage` in `AsterX/src/make.code.deps` (local,
  not committed). `c++filt`; `Lb1E`=true/`Lb0E`=false → `CalcFluxAll<uct,pplim,EOS>`.
  frontier.cfg DEBUG=no (production-path numbers).
- **Timing (user drives):** side-by-side final TimerReport, same iteration; judge
  by `AsterX_Fluxes / ODESolvers::Solve*` fraction (`_Subcycling` for subcycling).

## NEXT STEPS for the new agent

1. **DONE — compiler-flag sweep, all 3 probes, route CLOSED and DISMANTLED.**
   Results and the verdict table are in `compiler-flags.md` (§SWEEP RESULTS,
   §VERDICT); summary in the STATUS section above. The `PROBE_FLAGS` switchboard
   and the tiered candidate inventory were **deleted** from
   `AsterX/src/make.code.deps` and from that doc — deliberately, so nobody
   re-walks a list measurement has already priced. `make.code.deps` keeps only
   `-Rpass-analysis=kernel-resource-usage` and `-DASTERX_PROBE_IDEALGAS_ONLY`,
   which are still what you need to measure anything else. Reduced-build
   baseline: `uct1` 256/32/3448/occ-1, needing AGPR → 0; `uct0` 252/2/3448/occ-2
   is the canary.
2. **DONE — the `vbar` fix: written and measured 2026-07-27. Register NULL.**
   `if constexpr (pplim)` around the UCT drift-velocity blend in `fluxes.cxx`, so
   the cell-centred fallback is compiled out of the production kernel instead of
   multiplied by zero. Removed 12 provably-dead `gf_vels` loads and **moved not
   one counter** (all 11 kernels identical) → **Lesson 9**, and route (a) "genuine
   removal" is now closed with no candidates left.
   **Keep the change**, with the claim restated: bit-identical, deletes 12 real
   global loads (traffic, which `-Rpass` does not measure), removes a
   multiply-by-zero. **No occupancy claim.** Still owed for it: a full build (the
   reduced build never instantiates the `pplim=true` branch, which is where
   bit-identity has to hold exactly) and a golden gate, expected exactly 0. Do not
   spend a dedicated timing run on it — fold it into the next one.
3. **DONE — the `launch_bounds` plumbing: BUILT, occ-2 AT LOW SCRATCH** (see the
   STATUS section and `launch-bounds-plan.md` §TRIAL RESULT).
   **NEXT ACTION IS TIMING**, at both grid sizes, vs the Idea-1 baseline, recorded
   in `baseline-timings.md`. Reference points: forced occ-2 via a *global*
   `__launch_bounds__` measured −15.9% on TOV (with Z4c collateral this scoped
   version avoids); occ-2 via serialization only −4.8% because of scratch. This
   sits at 3608 B/lane, near the low-scratch end. Also owed: a **full-build**
   re-measurement (record `uct0`'s full-build baseline too — it was never taken),
   then the golden gate, which is expected at exactly 0 since `launch_bounds` is a
   codegen directive.
   Historical plan detail below:
   **THE ORIGINAL PLAN — the `launch_bounds` plumbing, now the primary route**
   (`launch-bounds-plan.md` Changes 1+2) — add a `min_blocks` template param to
   CarpetX `loop_box_device` in the user's fork at `../CarpetX`, routing to AMReX's
   existing 3-arg `launch_global<NT,MB>`. Two corrections recorded in that doc: the
   replication must call `amrex::detail::call_f_intvect_handler` (CarpetX passes an
   `(i,j,k)` lambda, not the 1D form the plan sketched), and passing a 7th template
   argument **breaks the CI build against stock CarpetX** — decide the CI story up
   front. The plan's "optional zero-code pre-check" `--amdgpu-waves-per-eu` **does
   not exist** in this toolchain. Note the sweep also produced a new argument *for*
   this route: it is per-launch-site, so unlike a TU flag it cannot damage
   `CalcE`/`CalcFstag`/`CalcAux`.
4. **Re-measure any winner in a FULL build** (comment out
   `-DASTERX_PROBE_IDEALGAS_ONLY`) before believing it — the reduced build is a
   different compilation (AGPR 32 vs 40).
5. Golden gate, then timing at BOTH grid sizes (UCT-small + TOV-large) vs the
   Idea-1 baseline. Record in `baseline-timings.md`. The reduced build IS usable for
   the timing runs (both are idealgas, `use_pplim=no`, and both CT schemes are
   compiled) — but remove the `-D` before any test-suite or golden run.
   **Exit condition:** if `launch_bounds` reaches occ-2 but the small grid still
   regresses, ship Idea-1 alone (occ-1, low scratch, fusion + `use_pplim` removal
   already banked) and stop chasing occ-2.
6. **PR cleanup** (whichever route ships): revert or neutralise `4f902ca1`;
   revert `2fd4595e` (`make.code.deps` — still hardwires `-Rpass` +
   `-DASTERX_PROBE_IDEALGAS_ONLY`, so the binary aborts on pplim/hybrid/tabulated);
   rewrite the `#if 0` CCTK_DEBUG block (only on the serialization branch); reword
   `[golden-master]`/`PROBE` commit messages for upstream.
7. **Independent of all occupancy work — the `tau` conditioning fix.** The fusion
   probe was a register null but turned up a real numerical finding: `tau` loses up
   to ~8 decimal digits to cancellation *in production today* (measured
   `|Q/tau| ~ 5e7`), because `rho*W*(h*W-1)` is a small residual of two large
   like-signed numbers whenever the fluid is cold and slow. A cancellation-free
   rewrite exists. **This is an accuracy improvement to shipping code with no
   occupancy claim attached** — see `ideas.md` §"tau conditioning" and
   `flux-construction.md` §10b. Do not bundle it with performance work.

## Superseded/historical

`Planning/Screenshot*.png` (the "is occ-2 worth it" deliberation) — resolved; occ-2
was reached, and the real question turned out to be scratch, not reachability.
