# Flux-kernel GPU occupancy — investigation, results

Frontier / MI250X (gfx90a, CDNA2), ROCm, cce/20. Concerns the fused flux kernel
`CalcFluxAll<uct,pplim,EOSType>` in `AsterX/src/fluxes.cxx`. Terms: `GPUHardwareDict.md`.

## Bottom line

- **occ-2 ACHIEVED (2026-07-24) for the production ideal-gas kernel** via a
  source-only register-reduction stack (no `launch_bounds`). Progression in
  §Serialization result. Live summary + lessons: `CHECKPOINT.md`.
- History: the fused kernel was occ-1 for the production EOS (ideal-gas,
  tabulated3d) — VGPR pinned 256 + AGPR overflow + ~4 KB scratch; hybrid healthy
  (54/occ7-8). occ-1 pre-existing (pre-fusion also occ-1) -> fusion cost no
  occupancy. Global `__launch_bounds__(MT,2)` forcing occ 1->2 gave the flux
  kernel **-15.9%** (the motivation).
- **tabulated3d is EOS-table-bound (AGPR ~172-244) -> occ-2 unreachable there;
  occupancy work is ideal-gas-only.**
- **⚠ UPDATE 2026-07-25 — the "~4 KB scratch" is `alloca`, NOT spill, and scratch
  does NOT cap occupancy.** The remark block reports spill counts and `ScratchSize`
  as separate fields, and the baseline kernel reads scratch 3448 with *both*
  `VGPRs Spill: 0` and `SGPRs Spill: 0`. So those ~431 doubles/lane are private
  memory ALLOCATIONS (almost certainly ReconX stencil arrays, whose runtime-indexed
  loops defeat SROA), not spilled registers. Occupancy is VGPR+AGPR (and LDS) only;
  scratch is a pure latency/bandwidth cost and is what made the small grid regress.
  Treat "reach occ-2" and "beat 4984 B/lane" as two INDEPENDENT criteria. Detail +
  the corrected Lesson 3: `CHECKPOINT.md`.
- **Two further measured facts (2026-07-25, reduced build — see
  `compiler-flags.md`).** (i) **UCT is the expensive CT config by ~34 registers**:
  `uct=1` totals 288 (256 VGPR + 32 AGPR, occ-1) against `uct=0` at 254 (252+2,
  occ-2). (ii) Cutting the instantiation matrix 16→2 moved AGPR 40→32 with no
  flags, i.e. **inlining of shared ReconX/EOS callees measurably drives register
  pressure** — so reduced-build numbers are not comparable to the tables below.
- **Active route is now compiler flags** (`compiler-flags.md`), with
  `launch_bounds` (`launch-bounds-plan.md`) as fallback; source-level algebraic
  removal is closed on measurement (`CHECKPOINT.md` Lesson 8).

## Production configuration (user)

Production typically **upwind-CT (`use_uct=yes`)** with **`use_pplim=no`**; EOS
ideal-gas + tabulated3d (both were occ-1). The -27% fusion win was measured on
both a UCT run (fixed-metric subcycling) and a flux-CT run (TOV Z4c).

## How to profile (methodology, reusable)

- **Tier 0 — static -Rpass (workhorse).** `fluxes.cxx.o: CXXFLAGS +=
  -Rpass-analysis=kernel-resource-usage` in `AsterX/src/make.code.deps` (local,
  not committed). Rebuild; read `VGPRs`/`AGPRs`/`ScratchSize`/`Occupancy`/spills
  per kernel. `c++filt` the mangled names; template args `Lb1E`=true/`Lb0E`=false
  -> `CalcFluxAll<uct,pplim,EOS>`. frontier.cfg DEBUG=no -> production-path numbers.
- **Compile-only source bisection** — `#if`-guard a region, rebuild, watch the
  remark. Localizes register pressure without a GPU run (how EOS/PP were exonerated).
- **Tier 1 rocprof** (`--stats --hip-trace`): launch count + per-launch time.
- **Tier 2 omniperf/rocprof-compute** (`-k <substr>`): roofline, achieved
  occupancy + limiter, cache/VALU/MemUnit — to disambiguate occupancy- vs
  bandwidth-bound.
- **NVIDIA (H100/GH100):** no -Rpass equivalent for AGPR; use `--resource-usage`
  (`-Xptxas -v`) -> `registers/thread` + `stack frame`/spill bytes, and `ncu` for
  achieved occupancy. See `GPUHardwareDict.md` §H100.

## Findings (evidence)

Baseline register/occupancy by EOS, fused `CalcFluxAll` (pre-serialization):

| EOS | VGPR | AGPR | Scratch B/lane | Occ |
|---|---:|---:|---:|---:|
| ideal-gas | 256 | 74 | 3976 | 1 |
| hybrid<poly>/<piecewise_poly> | 54 | 0 | 72 | 7-8 |
| tabulated3d | 256 | 243 | 4024 | 1 |

EMF/aux kernels (CalcE_impl, CalcFstag, CalcAux) all occ 6-8 -> EMF sector not
worth fusing.

**Exonerations (source bisection):** stubbing ideal-gas EOS math (`csnd`/`eps`/
`press`/`kappa`) -> peak byte-identical (not EOS math). Excising the whole PP
block -> peak VGPR unchanged at 256, only AGPR overflow dropped 74->46 (so PP
feeds the spillover, not the 256 core). The 256 core peak is the **main MHD
flux-assembly live-set** (many simultaneously-live `vec<vec<2>,3>` quantities).

**Why idealgas(256) vs hybrid(54) on identical MHD source:** register allocation
is a whole-function heuristic; ideal-gas's transparent trivial EOS lets the
scheduler stretch live ranges (higher peak), hybrid's opaque `ColdEOS*` calls act
as scheduling barriers (shorter live ranges). `CalcFluxAtFace` is
`ALWAYS_INLINE`, so the ~1000-line worker flat-inlines into one region.

**Forced occ-2 (global `__launch_bounds__(MT,2)`, TOV flux-CT):** AsterX_Fluxes
-15.9% (171.6->144.3s) BUT Z4c_RHS +48%, total +7.9% worse -> the directive must
be **scoped to the flux kernel** (global regressed the other fat kernels).

## Plan — EXECUTED

The go-forward plan (config-gated CT split; then serialize) is done. Sequence and
per-step -Rpass are in §Step 1 result and §Serialization result. **Not levers**
(ruled out by evidence, do not re-probe): trimming EOS math (exonerated);
un-fusing (pre-fusion was also occ-1); the neighbourhood-reuse pass (adds live
values -> backfires on occ-1); pure REORDERS and hand-rematerialization (the
allocator re-derives its own schedule — see CHECKPOINT Lesson 1; Idea-4 hoist was
bit-identical + golden but occupancy-null, reverted). **Fallback if ever needed:**
scoped `launch_bounds` (thread a `min_blocks` template param through CarpetX
`loop_box_device` -> the 3-arg `launch_global<NT,MB>`, applied only at the flux
site) — forces occ-2 by spilling; not needed now that occ-2 is reached genuinely.

## Step 1 — RESULT (CT-scheme split, commits b23fb947 / e8b0718c / 4a7e40ec)

Templatized `CalcFluxAll<uct,EOS>` + `if constexpr(uct)` guard (UCT face-speed/
drift block) / `if constexpr(!uct)` (flux-CT B-fluxes), dispatched via a generic
lambda over `std::true_type/false_type`. Template param named `uct` (NOT
`use_uct` — `DECLARE_CCTK_PARAMETERS` in `CalcFluxAll` shadows it). Layer 2 gated
STORAGE (`vbar_*`/`a_*` iff use_uct, `fluxBs_*` iff !use_uct) via duplicated
schedule blocks + interface.ccl B-flux split + sync.cxx restriction. nvcc gotcha:
an extended `__device__` lambda can't first-capture a var inside `if constexpr`
-> the fused-sweep zeroing needed `static_cast<void>(fx)` up top (nvcc-only).

**Both layers golden PASS at 0 (cpu/rocm/cuda). Occupancy: NO WIN** — production
idealgas stayed VGPR 256 / occ 1 (AGPR only 74->68->68, scratch 3976->3528). The
B-flux/UCT-face quantities are NOT what pins the 256 peak. Layer 1 is a
bit-identical config-gating refactor + enabler for layer 2 (storage/bandwidth:
-6 GFs, ~halved face-write streams, -11% scratch); the occupancy win needed the
serialization below.

## Serialization result (Ideas 1/2/3, 2026-07-24) — occ-2 achieved

Production kernel `CalcFluxAll<uct=1, pplim=0, idealgas>`, Frontier -Rpass:

| stage (commit) | VGPR | AGPR | scratch | occ |
|---|---:|---:|---:|---:|
| eig-collapse (549a28ca) | 256 | ~64 | ~3960 | 1 |
| Idea 1 template use_pplim (defd9f60) | 256 | 40 | 3448 | 1 |
| Idea 2 flux-by-flux momentum (e8613114) | 254 | 1 | 3816 | 1 |
| **Idea 3 face-state serialize (210a013c)** | **253** | **0** | 4984 | **2** |

- **Idea 1** compiled out the PP `_ppl` live-set for use_pplim=no (matches the
  historical PP-excision 74->46): AGPR 64->40. flux-CT+PP idealgas 64->39 too.
- **Idea 2** rolled `for(j)` over the momentum component (blows/moms/flux_moms one
  at a time): AGPR 40->1, VGPR 256->254. One granule from occ-2 (VGPR 254 rounds
  to 256; residual 1 AGPR rounds up -> total >256 -> occ 1).
- **Idea 3** rolled `for(f)` over the two face sides (group-B intermediates as
  per-side scalars; `eigenvalues_oneside` per side; momentum `for(j)` nested;
  combine post-loop): drained the last AGPR 1->0 -> VGPR 253 alone = 512/256 = 2
  waves = **occ 2**. flux-CT+noPP idealgas also crossed (251/0/occ2).
- **Scratch rose 3448->4984 B/lane** (overflow moved AGPR->off-chip HBM). At occ-1
  a slight loss; at occ-2 the 2nd wave hides it. Net win pending a timed run.
- PP-on idealgas (`<1,1>`=28 AGPR, `<0,1>`=29) stays occ-1 (not production).
  tabulated `<1,0>` 256/172/occ1 (down from 244 but EOS-bound). hybrid 54/0/occ8.

**Lesson confirmed:** the allocator banks REMOVAL (Idea 1) and STRUCTURAL rolled
loops (Ideas 2/3), not reorders (CHECKPOINT Lesson 1). Every step golden PASS at 0.

**⚠ Why scratch ROSE (3448→3816→4984) while AGPR fell — the central caveat.** The
rolled loops make the loop indices `f`/`j` RUNTIME variables. GPU registers can't
be dynamically indexed, so every array read with a runtime index — the
reconstructed vectors (`rho_rc`/`vels_rc`/`Bs_rc`/`vlows_rc`/…) and the cut-set
accumulators (`moms_rc`/`flux_moms`/`dens_rc`/…) — is demoted from registers to
scratch (off-chip HBM). So AGPR↓ and scratch↑ are the SAME data relocating, not
shrinking. Idea 3's +1168 B is dominated by the reconstruction arrays (demoted by
the runtime-`f` reads); the cut-set is a fraction. Idea 1 (removal) reduced BOTH
registers and scratch — the "good kind"; Ideas 2/3 (serialization) relocate.

**Timing outcome (Frontier, ideal-gas): occ-2 is grid-size-dependent** —
TOV/flux-CT/large **−4.8%** flux; fixed-metric/UCT/small **+8.9%** (scratch
bandwidth dominates at small grid). Genuine occ-2 (−4.8%) ≪ forced occ-2 (−15.9%,
same TOV run) — the gap is scratch. Full numbers: `baseline-timings.md`. ⇒ the
active lever is scoped `launch_bounds` on the low-scratch Idea-1 kernel
(`launch-bounds-plan.md`); a two-per-side-kernel split reusing the flux GFs is a
recorded big-swing alternative (`ideas.md`).

## Frontier practicalities

- Profile single-rank/single-GCD (`ROCR_VISIBLE_DEVICES=0`), short `cctk_itlast`
  but enough that every FMR level's flux launch fires.
- `-ffp-contract=off` is a CI-cpu-only concern; use normal Frontier opt flags so
  the profile reflects production codegen. frontier.cfg: DEBUG=no, OPTIMISE=yes -O3.
- The AMReX launch_bounds edit (fallback) is header-only (template
  `launch_global` instantiated in `fluxes.cxx`) -> rebuild AsterX only.
