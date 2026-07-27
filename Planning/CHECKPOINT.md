# Checkpoint — AsterX flux-loop optimization (handoff)

_Last updated 2026-07-27. Read this first. Then: `baseline-timings.md` (all
timings), `launch-bounds-plan.md` (the winning occupancy change),
`flux-construction.md` (numerical pipeline; §10b = the open `tau` accuracy fix),
`profiling-flux-kernel.md` (profiling methodology + EOS occupancy table),
`ideas.md` (open future work), `implementation-history.md` (provenance),
`compiler-flags.md` + `ideas-2-3-serialization.md` (records of closed routes),
`GPUHardwareDict.md` (terminology)._

## Mission

Optimize `AsterX/src/fluxes.cxx` (the fused flux kernel `CalcFluxAll`), starting
from bit-identity-gated refactors (golden-master = exactly 0) and then
register/occupancy work measured on Frontier/MI250X. CI CPU is cost-neutral by
design; the payoff is GPU.

## STATUS: ⭐ SOLVED (pending controls) — **−47.5 % on the flux kernel**

A `min_blocks` template parameter on CarpetX `loop_box_device` plus `MB=2` at the
one fused flux site forces **occ-2 at scratch 3608 B/lane** and, combined with the
earlier fusion work, cuts `AsterX_Fluxes` in half versus the original upstream
code. Two independent configurations, near-identical decomposition:

| config | `AsterX_Fluxes` | fusion + Idea 1 | + occupancy | cumulative | CCTK total |
|---|---|---:|---:|---:|---:|
| UCT-small (subcycling) | 157.3 → 82.0 s | −26.5 % | −29.1 % | **−47.9 %** | **−13.4 %** |
| TOV-large (Z4c, flux-CT) | 237.4 → 124.8 s | −27.7 % | −27.3 % | **−47.4 %** | **−17.6 %** |

The two halves are near-equal, the two configs agree to 0.5 pp, `Z4c_RHS` is flat,
and in every run the flux delta accounts for the whole `Solve::rhs` delta. Against
the *same* TOV baseline the serialization route managed only −4.8 % — same
occupancy, 5.7× less benefit, the difference being purely where the register
overflow lives. **The flux kernel is no longer the dominant cost** (47.6 % → 32.7 %
of Solve on TOV, 37.8 % → 24.2 % on subcycling), which is worth weighing before
any further flux work.

**Owed before this is quotable / PR-ready:**

1. **The `MB=0` control** — same reduced build, same par file. Both timing runs are
   reduced-build (`-DASTERX_PROBE_IDEALGAS_ONLY`) against a full-build baseline, so
   the headline is (instantiation cut + `MB=2`). The cut alone is probably worth
   little (at `MB=0` the reduced build was still 256/32/3448/occ-1, changing neither
   occupancy nor scratch) but that is inference. One run settles both configs.
2. **Full-build re-measurement**, including `uct0`'s full-build baseline, never taken.
3. **Golden gate**, expected exactly 0 (`launch_bounds` is a codegen directive). Needs
   the probe commits reverted or the CarpetX branch upstreamed — golden CI builds
   against stock CarpetX.

## Branches — PR-BOUND STACK (built 2026-07-27)

Two clean branches, both containing **AsterX files only** (verified: no `Planning/`,
no `make.code.deps` — the probe machinery is deliberately excluded entirely):

| branch | tip | contents |
|---|---|---|
| **`opt/flux-registers`** | `1ed53981` | the bit-identical stack: `b23fb947` → `e8b0718c` → `4a7e40ec` (CT-scheme split: compute / hoist / storage) → `549a28ca` (eigenvalue collapse) → `defd9f60` (Idea 1, `use_pplim` compiled out) → `1ed53981` (`vbar` fix). **Fast-forwarded from `4a7e40ec` to `defd9f60` on 2026-07-27** — it previously held only the CT-scheme split, i.e. it was two commits short of the state every timing number was measured against. |
| **`opt/flux-min-blocks`** | `cd6b3b7a` | `opt/flux-registers` + `MB=2` at the flux site. ⚠ **Requires the CarpetX branch**; cannot build against stock CarpetX. |

`opt/flux-registers` should pass golden at exactly 0 end to end (`1ed53981` is the only
row not yet gated). `opt/flux-min-blocks` needs CarpetX
`opt/loop-box-device-min-blocks` @ `19267243` present first.

**How the two branches differ, and which to use when.** They share the base
`defd9f60`, both carry the `vbar` fix and `MB=2`, and the difference is **purely
additive** — 111 insertions, 0 deletions: `opt/flux-launch-bounds` adds 51 lines of
`#ifdef ASTERX_PROBE_IDEALGAS_ONLY` guard in `fluxes.cxx`, the 60-line
`make.code.deps`, and the 11 Planning files. The guard is `#ifdef`-inert, so **without
the `-D` both branches preprocess to the same translation unit**. But
`opt/flux-launch-bounds` ships a `make.code.deps` that *does* define it, so as checked
out it builds the reduced 2-instantiation binary that **aborts at runtime** on
`use_pplim=yes`, hybrid or tabulated.

- **`opt/flux-launch-bounds`** → `-Rpass` sweeps and the timing runs (both are
  idealgas/`use_pplim=no`, both CT schemes compiled). **Run the owed `MB=0` control
  here** — flip `MB` back to 0 at the call site, cheap rebuild, same par file.
- **`opt/flux-min-blocks`** → full 16-instantiation build: the **full-build
  re-measurement**, the test suite, golden. Same code, without the instantiation cut,
  which is exactly the confound the control removes.

The two histories are disjoint after `defd9f60` — do not merge or rebase one onto the
other. Further bit-identical work goes on `opt/flux-registers` and is replayed onto both.

**Stale but still present** (deletable whenever): `opt/flux-loops` — its tip
`c2614613` is **already upstream**, merged into `EinsteinToolkit/AsterX` as PR #146
("Opt/flux loops: Update CI test suite", merged 2026-07-13), verified as an ancestor of
`etk/dev`; and `opt/flux-loops-impl`, an ancestor of `opt/flux-registers`. Both exist
local and on `origin`. The branch names still appear in these docs as *labels for code
states* (e.g. the GPU fusion comparison) — that usage stays correct.

`549a28ca`'s subject says "PROBE" — **decided 2026-07-27: leave it as is**, no rebase.

## Archived branches — retrieve from tags, not from branches

Both dead-end branches were **deleted and replaced by annotated tags**, local and on
`origin`. Their commits are unreachable from any branch, so **`git fetch --tags` is
required** — plain `git fetch` only follows tags that point into branch history.

| tag | commit | what it holds |
|---|---|---|
| **`archive/vf2-accuracy-probe`** | `659b48b9` | the `(H,vf2)` fusion probe (was `probe/flux-enthalpy-fusion`; the remote tip `ad9b9bd6` had an identical tree). **`eigenvalues.hxx` here is the asset for the accuracy PR** — see §Accuracy PR in `ideas.md`. Register motivation dead; its `tau` is form (B), a regression — don't take it. |
| **`archive/serialization-occ2`** | `210a013c` | Ideas 2+3, the rolled-loop serialization (was `opt/flux-eig-collapse`). The only source-only occ-2, golden PASS at 0, but scratch 4984 → −4.8 % / +8.9 %. Kept as the counterfactual the `launch_bounds` result is measured against, plus `eigenvalues_oneside` and the rolled-loop bit-identity argument. |

Each tag message records what to take and what to avoid. To restore either as a
working branch:

```bash
git fetch --tags origin
git branch <name> archive/serialization-occ2     # or: git checkout -b <name> <tag>
git show archive/vf2-accuracy-probe:AsterX/src/eigenvalues.hxx   # single file, no checkout
```

## Working branch (probe machinery + these docs)

**AsterX** (`origin` = MChabanov/AsterX), branch **`opt/flux-launch-bounds`**:

| commit | content |
|---|---|
| `b23fb947`/`e8b0718c`/`4a7e40ec` | CT-scheme split (compute / hoist / storage gating) |
| `549a28ca` | eigenvalue collapse to `charmax`/`charmin` |
| `defd9f60` | **Idea 1** — template `use_pplim`, PP compiled out. The occ-1 low-scratch (3448) reference build for every timing comparison. Upstream-bound work ends here. |
| `4f902ca1` | ⚠ probe guard `ASTERX_PROBE_IDEALGAS_ONLY` (inert without the `-D`) |
| `2fd4595e` | ⚠ `make.code.deps` — hardwires `-Rpass` + the `-D`; binary **aborts** on pplim/hybrid/tabulated |
| `c96df5d5` | ⚠ `MB=2` at the flux site — **breaks CI by design**, needs the forked CarpetX |
| + `vbar` fix | `if constexpr (pplim)` around the UCT drift blend; register-null but keep (12 dead loads) |

This branch is the **working/measurement** branch: it keeps the probe guard, the
`make.code.deps` with `-Rpass`, and these Planning docs. **None of that is upstream
bound** — the PR-bound content was extracted to the two branches above by path
(`AsterX/src/fluxes.cxx` hunks only), not by cherry-picking commits, because
`3a6039c8`/`86bb2c48` mix the `vbar` fix with `make.code.deps` and doc edits.
**These docs are deliberately kept local and are not committed onto the PR branches.**

**CarpetX** (`origin` = MChabanov/CarpetX): **`opt/loop-box-device-min-blocks`** @
`19267243`, off `dev` @ `55e7434e`, pushed. One file, `Loop/src/loop_device.hxx`.

**Dead-end branches are gone — see §Archived branches above** for the two tags that
replaced them (`archive/serialization-occ2`, `archive/vf2-accuracy-probe`). Also
deleted: `backup/flux-registers-prelinear`, a pre-rebase copy whose three commits were
tree-identical to `b23fb947`/`e8b0718c`/`4a7e40ec`, so nothing was archived from it.
**Idea 4** (hoist eig+UCT ahead of the assembly) was golden-PASS but occupancy-null
and was REVERTED — it exists only in this history, not as a branch.

## The `-Rpass` progression (production `CalcFluxAll<uct=1,pplim=0,idealgas>`)

Full build (all 16 instantiations):

| stage | VGPR | AGPR | scratch B/lane | occ |
|---|---:|---:|---:|---:|
| eig-collapse (pre-Idea 1) | 256 | ~64 | ~3960 | 1 |
| Idea 1 (`defd9f60`) | 256 | 40 | 3448 | 1 |
| Idea 2 (flux-by-flux `for(j)`) | 254 | 1 | 3816 | 1 |
| Idea 3 (face `for(f)`) | 253 | 0 | 4984 | **2** |
| `(H,vf2)` fusion off Idea 1 | 256 | 40 | 3448 | 1 (null) |

Reduced build (`-DASTERX_PROBE_IDEALGAS_ONLY` — a *different compilation*, AGPR 32
not 40; comparable only within this block):

| stage | VGPR | AGPR | scratch | occ |
|---|---:|---:|---:|---:|
| Idea-1 baseline | 256 | 32 | 3448 | 1 |
| `--misched=gcn-max-occupancy` | 256 | 32 | 3448 | 1 |
| `--enable-deferred-spilling` | 256 | 30 | 3448 | 1 |
| `-ffinite-math-only` | 256 | 34 | 3448 | 1 |
| `vbar` fix (12 dead loads removed) | 256 | 32 | 3448 | 1 |
| **`launch_bounds(256,2)`** | **128** | **128** | **3608** | **2** ⭐ |

Other configs: flux-CT+noPP idealgas also reached occ-2 under serialization
(251/0/4984). PP-on idealgas stays occ-1 (28–29 AGPR; not production).
**tabulated3d stays occ-1 (AGPR ~172–244), EOS-table-bound — occ-2 is unreachable
there by any flux change.** hybrid is healthy at 54/0/occ-8.

## KEY LESSONS

1. **The allocator banks REMOVAL and STRUCTURAL change, never REORDERING.**
   Removal (Idea 1, PP compiled out) → AGPR 64→40 *and* scratch 3976→3448, both
   down. Rolled loops (Ideas 2/3) → AGPR 40→0. Pure reorders (CT-hoist,
   eig-collapse, Idea 4) → occupancy-null: the allocator re-derives its own
   schedule. Don't spend effort on reorders or hand-rematerialization. See
   Lesson 4 for what "removal" actually requires.

2. **Occupancy arithmetic.** `occ = floor(512/(VGPR+AGPR))`, granule-rounded (VGPR
   granule 8), cap 8. **At occ-1 the gauge was AGPR → exactly 0** — VGPR welds at
   the 256 ceiling while demand exceeds it, so reducible overflow lands in AGPR,
   and one residual AGPR granule blocks the crossing (254/1 was still occ-1).
   **⚠ That gauge does NOT survive `launch_bounds`.** Under `MB=2` the kernel sits
   at VGPR 128 / AGPR 128 = 256 total = the occ-2 budget, with both spill counters
   0: the allocator used the AGPR half of the unified gfx90a file as *on-chip*
   spill space. An even 128/128 split is a healthy allocation, not overflow. Read
   AGPR as "overflow to drain" only when VGPR is pinned at 256.

3. **Serialization RELOCATES, it doesn't REMOVE — and scratch is independent of
   occupancy.** Rolled loops (`#pragma unroll 1`) make the indices `f`/`j` runtime
   variables; GPU registers cannot be dynamically indexed, so every runtime-indexed
   array (the reconstructed vectors, the cut-set accumulators) is demoted to
   scratch. AGPR↓ and scratch↑ are the same data moving off-chip, which is why
   scratch rose 3448→3816→4984. The seesaw: unrolled = registers/occ-1/low-scratch;
   rolled = scratch/occ-2/high-scratch.
   Two refinements, both measured:
   - **The baseline 3448 B/lane is `alloca`, NOT spill** — `-Rpass` reports
     `ScratchSize` and spill counts separately, and the baseline reads 3448 with
     *both* spill counters 0. Those ~431 doubles/lane are private-memory
     allocations, almost certainly the ReconX stencil arrays whose runtime-indexed
     6-point loops defeat SROA. So a large part of the reconstruction state was
     *already* in scratch before Ideas 2/3; read 3448→4984 as an increment on a
     pre-existing floor. (Contrast `CalcE uct1`: 72 B/lane *with* 20–24 SGPR
     spills — that kernel's scratch really is spill.)
   - **Occupancy and scratch are INDEPENDENT levers.** Occupancy is VGPR+AGPR (and
     LDS) only. "Reach occ-2" and "beat 4984 B/lane" are separate criteria — the
     first is occupancy, the second is the timing/bandwidth cost that made the
     small grid regress. Never conflate them.
   **This lesson is now confirmed end-to-end:** occ-2 at 3608 B/lane timed −27.3 %
   on TOV where occ-2 at 4984 timed −4.8 %. Occupancy was never the problem; the
   scratch was.

4. **⚠ REMOVAL ONLY COUNTS AT THE PRESSURE PEAK — position beats category.**
   Register count is the maximum number of values live *at one point*. The peak is
   the **reconstruction** state (S2–S6, longest and most-interfering live ranges).
   Anything born and consumed downstream of the flux-assembly waist occupies
   registers that are already free, so removing it lowers total work but never the
   peak. Measured twice, from two different directions:
   - the `(H,vf2)` algebraic fusion deleted `h_rc`, `B2_rc`/`bsq_rc`,
     `dens_h_W_rc`, the inlined `a_m`/`a_p`/`det` transients — all of which the
     allocator was **already rematerializing** — and was byte-for-byte null;
   - the `vbar` fix deleted **12 provably-dead global loads**, i.e. values that
     *cannot* be rematerialized, and was **also** byte-for-byte null.
   So the load-vs-arithmetic distinction is irrelevant; only position is. Idea 1
   worked because compiling out `pplim` deleted a duplicated copy of the whole
   S5–S9 chain *including its own reconstruction state* — not because it deleted a
   dozen doubles of algebra. A count of names in the source is not a count of
   values in registers.
   **Consequence: source-level removal is CLOSED with no candidates left.**
   Everything at the peak is reconstruction-side, and that path is separately
   blocked (see "Ideas considered and set aside"). Before proposing any removal,
   ask *where* it lives relative to the peak.
   Corollary for the *scratch* problem: reconstruction-side data removal would cut
   the 3448 B/lane of `alloca` more than it cuts registers — a different change
   from anything tried. Dump the ISA before guessing:
   `v_accvgpr_read/write_b32` shows what is truly AGPR-resident, `scratch_load/store`
   what is in private memory.

5. **⭐ CONSTRAINING THE ALLOCATOR BEAT COAXING IT, decisively.** Four independent
   attempts to persuade the allocator to find 32 registers — three compiler-side
   (scheduling, allocation, relaxed FP) and one source-side (`vbar`) — all landed
   within ±2 registers of the 32-register gap. A hard `__launch_bounds__(256,2)`
   found all 32 on the first try. The registers were not recoverable by *asking*;
   the peak was compressible once a budget was *imposed*. Two supporting findings:
   - **The 32 AGPRs were load-bearing live values, not allocator sloppiness** — the
     one flag that genuinely re-ran allocation found 2 registers; the one that
     relaxed semantics found −2.
   - **TU-scoped flags cannot be aimed at one kernel.** `CalcFluxAll`,
     `CalcE_impl`, `CalcFstag` and `CalcAux` all live in `fluxes.cxx`, so every
     `-mllvm` flag hits all 11 kernels — and `--enable-deferred-spilling`
     measurably damaged `CalcE` (SGPR spills 20/22/24 → 30/26/33) while helping the
     target by 2. `launch_bounds` is per-launch-site by construction: measured
     zero collateral, all 9 EMF kernels byte-for-byte unchanged, `Z4c_RHS` flat in
     timing. That is structural, not merely a stronger knob.

6. **⚠ A REDUCED BUILD IS A DIFFERENT COMPILATION — never carry occupancy across.**
   `-DASTERX_PROBE_IDEALGAS_ONLY` cuts 16 instantiations → 2 and moves AGPR 40→32
   with no flags, because inlining of the shared ReconX/EOS callees changes (LLVM
   scores call sites partly on user count). This cost one wrong prediction: `uct0`
   reads 252/2/occ-2 in the reduced build (two registers under the cliff, i.e.
   threshold luck), which was extrapolated to a full-build baseline where it had
   never been measured — predicting TOV-large would be neutral, where it measured
   −27.3 %. Reduced-build rows are comparable to each other only.

7. **Trust `-Rpass` and timing, not arithmetic.** Hand double-counts were wrong in
   both directions: predicted Idea-2 AGPR 40→16, got 40→1; predicted Idea-3 net
   ~0–6, got occ-2; predicted the fusion at 24–28 registers, got 0.

8. **Production is upwind-CT (`use_uct=yes`), `use_pplim=no`, ideal-gas +
   tabulated3d.** No test par overrides `use_pplim`, so the golden gate exercises
   the production `pplim=false` kernel, and `theta_x/y/z` + `theta_tot` are not in
   the golden `.tsv` set (which is what made theta storage-gating safe). Measured
   aside: **UCT is the expensive CT config by ~34 registers** (`uct=1` totals 288
   against `uct=0` at 254 in the reduced build).

9. **Bit-identity scoping.** The serialization is bit-identical under
   `-ffp-contract=off`: rolled loops reorder independent sub-computations without
   reassociating within an expression, per-side contractions reuse the same scalar
   `calc_contraction` overloads (`AsterUtils/src/aster_utils.hxx`), and the
   eigenvalue sides were already independent. **But `-ffp-contract` is unset on
   Frontier (so FMA contraction is ON) and explicitly `off` in the CI CPU config —
   golden = 0 on CI implies nothing about bit-identity of the Frontier binary.**
   CI validates, Frontier measures; never conflate the two.

10. **NVIDIA H100/GH100 has no AGPR** (unified 64K register file, 255/thread, spills
    to cached local memory). The AGPR framing is AMD-only; on Hopper measure
    `registers/thread` + spill bytes via `--resource-usage`. Detail:
    `GPUHardwareDict.md`.

## Ideas considered and set aside (do not re-derive)

- **Per-side (incremental) Riemann combine** — accumulate `fsum`/`vacc` per side,
  apply `c`/`0.5` post-loop. Bit-identical for LxF (but *not* the naive
  `0.5*flux0+0.5*flux1`, which reassociates); HLLE saves less since the per-side
  coefficients don't separate. Trims the both-sides cut-set ~7–14 doubles, i.e.
  ~15 % of the scratch increase, but not the dominant reconstruction-array
  scratch → not worth the redo.
- **Grade the atmosphere once at the face** — a defensible *scheme* change
  (non-bit-identical, needs maintainer sign-off + physics validation); modest
  compute win (halves the atmosphere `pow`s). Outside the bit-identical track.
- **Full per-side kernel rewrite touching reconstruction** — blocked by the
  `useLO`/velocity-limit fallback being an OR over BOTH sides (`reconstruct()`
  returns both), plus a constant-index two-block unroll relying on the allocator
  not interleaving (likely occ-null à la Idea 4), plus large bit-identity/ReconX
  risk. The GF-reuse variant is assessed in `ideas.md`; `launch_bounds` obtained
  the same occupancy for a fraction of the effort, so this is now moot.

## FP-determinism prerequisite (do not undo)

The CI cpu config compiles C++ with `-ffp-contract=off` and no `-funsafe-math`
(`scripts/actions-cpu-real64.cfg`); the golden `.tsv` were regenerated under those
flags. Changing them invalidates golden and the baseline.

## How to operate

- **Golden check:** put `[golden-master]` in the HEAD commit message of the push.
  Pass = `2235 ok, 0 fail`, worst `|abs|=|rel|=0`, np1+np2, cpu/rocm/cuda. Watch
  with `gh run watch <id> --exit-status`. All of Ideas 1/2/3 passed.
- **`-Rpass` (user drives on Frontier):**
  `fluxes.cxx.o: CXXFLAGS += -Rpass-analysis=kernel-resource-usage` in
  `AsterX/src/make.code.deps`. Demangle with `c++filt`; `Lb1E`=true / `Lb0E`=false
  → `CalcFluxAll<uct,pplim,EOS>`. frontier.cfg has DEBUG=no, so the numbers are
  production-path. Remarks are emitted at compile time, so getting output at all
  proves the TU recompiled. Parsing gotchas and the extraction one-liner:
  `compiler-flags.md`.
- **Timing (user drives):** side-by-side final TimerReport, same iteration; judge by
  `AsterX_Fluxes` and its fraction of `ODESolvers::Solve*` (`_Subcycling` for the
  subcycling test). Confirm comparability by checking that unrelated timers are flat.

## NEXT STEPS

1. **The `MB=0` control**, then the full-build re-measurement, then golden — the
   three items in §STATUS. Nothing else should be quoted until the control is in.
2. **PR structuring — decide the split before opening anything.**
   `opt/flux-registers` against `etk/dev` is **18 commits, 2038 files, ~51.7k
   insertions**, because Target 1 (the fusion) was never upstreamed and the bulk of
   those files are the regenerated golden `.tsv` from `cc435a06`. Only the
   *infrastructure* went up (PR #146). So this is one very large PR as it stands.
   Natural seam: **(i)** fusion + FP-determinism flags + golden rebaseline
   (`9fcb41d6`…`1ce741de`, `b9c200ef`, `cc435a06` — the 2038-file part), then
   **(ii)** the register work (CT split, eig-collapse, Idea 1, `vbar`), which is a
   handful of source files. Reviewers can actually read (ii); (i) is mostly
   regenerated data with a one-paragraph rationale.
   `549a28ca`'s subject still says "PROBE" — **decided 2026-07-27: leave as is.**
3. **CarpetX `min_blocks` — routing decided 2026-07-27: propose it to `lwJi/CarpetX`
   first** (the `liwei` remote), *if* the change is wanted at all; not directly to
   `EinsteinToolkit/CarpetX`. Note for whoever does it: the branch
   `opt/loop-box-device-min-blocks` was cut from the local fork's `dev`, which is
   hundreds of commits ahead of `EinsteinToolkit/CarpetX:main` (fork-only TimerReport
   and ccache work, plus lwJi's agent_scripts/subcycling commits), and upstream has no
   `dev` branch. So a PR must be a **cherry-pick of `19267243` alone** onto a branch cut
   from whichever base the receiving repo wants — it is one self-contained file, so this
   is trivial. `opt/flux-min-blocks` cannot merge anywhere until that parameter exists in
   the CarpetX that AsterX CI resolves.
4. **Independent of all performance work — the `tau` conditioning fix.** `tau`
   loses up to ~8 decimal digits to cancellation *in production today* (measured
   `|Q/tau| ~ 5e7`), because `rho*W*(h*W-1)` is a small residual of two large
   like-signed numbers whenever the fluid is cold and slow. A cancellation-free
   rewrite exists: `flux-construction.md` §10b, `ideas.md` §tau. Not
   bit-identical, so it needs a conserved-quantity argument rather than golden = 0.
   **Do not bundle it with performance work.**
5. Optional, now that the flux kernel is down to a quarter/third of Solve: pick the
   next target deliberately rather than by momentum. `AnalyticalSpacetimeX_SetMetric`
   (99 s, ~29 % of Solve in the subcycling run), `Z4c_RHS` (83 s, untouched) and
   `AsterX_SourceTerms` (57 s on TOV) are all now larger than further flux gains.
