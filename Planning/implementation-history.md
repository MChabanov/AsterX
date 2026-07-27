# Implementation history — AsterX flux-loop optimization

Compact provenance record. Remotes: `origin` = MChabanov/AsterX (CI runs here),
`etk` = EinsteinToolkit/AsterX. Full design discussions lived in per-target plan files
(deleted after completion; the git history of `Planning/` is not kept — this file is
the record). Current branch/commit map: `CHECKPOINT.md` §Branches.

## Supporting infrastructure (branch `opt/flux-loops`, 2026-07-10)

- **Golden-master harness:** `scripts/compare_tsv.py` — stdlib-only comparator of
  produced vs committed test `.tsv`, keyed by grid locator (handles 1p/2p row
  reordering), `|Δ| ≤ atol + rtol·|b|`, default 1e-13. Wired into CI as an opt-in
  step. Committed data regenerated from the CI container (`bea0b913`), since data from
  other machines drifts ~1e-9. Decisions: reuse the existing 37 test pars and their
  committed `.tsv` as the golden master (no new tests); TSV precision fixed at 16
  digits (hardcoded in CarpetX `io_tsv.cxx`, not tunable) — accepted.
- **Per-routine timing:** all 18 evolution test pars activate `TimerReport`
  (`out_every=50`, top 50); `scripts/print_timers.py` prints the hottest AsterX
  routines per test in the CI log; full timer files ride in the `asterx-test-output`
  artifact. (`cctk_timer_output` alone was insufficient — it hides `AsterX_Fluxes`
  inside ODESolvers.)
- **ccache CI:** `scripts/build.sh` + `actions/cache` on `.ccache`; warm cuda/rocm
  builds ~25 min → ~6–7 min.
- Commits: `bae19fe5`, `bc162f74`, `eadd1fa7`, `f882c088`, `5830bc95`, `bea0b913`,
  `f4a469f7`, `c2614613`.

## Target 1 — fused single-sweep flux calculation (`opt/flux-loops-impl`, 2026-07-11)

Hard requirement: bit-identical, each step gated by golden diff = 0. Fidelity spec for
the extraction: verbatim copy-paste, only sanctioned substitutions (`p.X` → explicit
`face_X`; captures → `FluxContext` members).

- `0930f4c3` **face_X made explicit** in the flux loop (used by the atmosphere-radius
  code, the one place needing the *face* coordinate).
- `b6739ec9` **`FluxContext<EOSType>`** aggregate: everything the per-face code needs,
  members named as the locals they carry; positional aggregate init.
- `9fcb41d6` **Step A extraction**: loop body (1002 lines) moved verbatim
  (mechanically verified, only a uniform 2-space dedent) into
  `CalcFluxAtFace<dir_i, EOSType>(fx, p, face_X)`; `dir_i` stays compile-time.
- `e25ce7ce` cuda fix: `maxspeeds_from_lambdas`/`avg_upwind` in AsterUtils were the
  only helpers without `CCTK_HOST`; nvcc rejects device-only calls from the host path
  of the H+D helper. (`hll_upwind` is still device-only — it will need the same if the
  EMF sweeps are ever fused.)
- `1bc276a5` **diagonal B-fluxes retired**: fxBx/fyBy/fzBz were identically zero
  (`(dir_i != i) * flux`) and unread. Replaced variable-indexed `fluxB{x,y,z}s` vecs by
  direction-indexed off-diagonal pairs `fluxB_j = {fxBy, fyBz, fzBx}`,
  `fluxB_k = {fxBz, fyBx, fzBy}` (the `vbar_j/k` idiom); CalcE mapping with cyclic
  `(i,j,k)`: `gf_fBs(j)(k) = gf_fB_k(k)`, `gf_fBs(k)(j) = gf_fB_j(j)`. Drops 3/30 flux
  GFs and 1/10 Riemann evals per face.
- `0d1f41eb` CI: `[golden-master]` in the HEAD commit message of a push triggers the
  golden step (besides dispatch input and PR label).

### FP-determinism lesson (critical)

The first golden run on the provably-verbatim Step A **FAILED**: 123/2235 files,
ulp-level early deviations growing to ~1e-6 rel (shock fronts amplify last-bit
perturbations; limiter branches flip on thresholds). Cause: the CI cpu cfg compiled C++
with `-funsafe-math-optimizations -ffp-contract=fast`, under which the compiler may
legally change FP results whenever code is restructured — so bit-identity across ANY
refactor is unattainable.

- `b9c200ef` `scripts/actions-cpu-real64.cfg`: CXX flags → `-ffp-contract=off`, no
  unsafe-math (C/Fortran and cuda/rocm cfgs unchanged; CI-only config).
- `cc435a06` golden `.tsv` regenerated from PRE-refactor code under the new flags (run
  29158250436 on temp branch `tmp-golden-baseline`; 1-proc == 2-proc byte-identical;
  2030 files replaced, 205 unchanged). That run is also the canonical performance
  baseline (`baseline-timings.md`).
- Validation run 29159401173 then reported **2235 ok, 0 fail, worst |abs|=|rel|=0** for
  Step A + the B-flux retirement. (gh note: the keyring token expired mid-session once;
  `gh auth login` fixes it.)

### Steps B and C — the fusion

- `1b6ca2ee` **Step B**: three per-direction sweeps → one fused `CalcFluxAll<EOSType>`.
  Per-direction domains reproduced exactly by `calc_mixpn_box<CI,CJ,CK>` (the same
  `box_int`/`box_all` ± `nloop`-after-tile-clamping computation as CarpetX's
  `loop_mixpn_device`, on the same grid object); fused loop = `loop_box_device<0,0,0>`
  over the per-dim union with per-direction guards. `face_X` computed per direction
  with CarpetX's own `point_desc` formula `x0 + (lbnd + I - (!CI)/2) * dx` →
  bit-identical to the old face-centred `p.X`. Verified prerequisites: the flux body
  and ReconX use only the centering-independent `p.I`/`p.DI`/`p.DX`; no
  cross-face/cross-direction data flow. **Golden: exactly 0.**
- `6728188b` **Step C1**: three zeroing kernels → one (runtime dir loop, guards from
  the same `box_all` boxes). Launches 4 → 2.
- `1ce741de` **Step C2**: zeroing merged into the main sweep → **1 kernel launch per
  `AsterX_Fluxes` call**, eliminating a separate full write pass over the face GFs.
  Value-identity: the flux code never reads an init-written value (`gf_theta` and
  own-face flux reads are preceded by same-point writes in both PP-limiter branches);
  cross-tile overlaps idempotent. **Golden: exactly 0.**

CPU verdict: cost-neutral by design; GPU verdict: **−27 %**. Both in
`baseline-timings.md`.

## Register reduction to occ-2 (branch `opt/flux-eig-collapse`, 2026-07-23/24)

Goal: drain the ideal-gas AGPR overflow to reach occ-2. All bit-identical (golden PASS
at 0, cpu/rocm/cuda). `-Rpass` tables: `profiling-flux-kernel.md`.

- `b23fb947`/`e8b0718c`/`4a7e40ec` **CT-scheme split** (compute + hoist + storage
  gating), `CalcFluxAll<uct,EOS>`. Golden PASS; occupancy-null (AGPR 74→68) but a
  storage/bandwidth win and the enabler for what follows.
- `549a28ca` **eigenvalue collapse** — `eigenvalues()` is degenerate (4 stored = 2
  distinct per side) and every consumer reduces to `{charmax=max(0,·),
  charmin=min(0,·)}`. Collapsed to those two scalars, threaded through
  calcflux/laxf/UCT. Bit-identical (max/min/fabs exact). ~2–4 AGPR.
- `defd9f60` **Idea 1** — template `use_pplim` (like `uct`): `if constexpr(pplim)`
  around the whole PP-limiter block, so for `use_pplim=no` (production) the `_ppl` live
  set compiles out. AGPR 64→40 *and* scratch 3976→3448. Also storage-gated
  `theta_x/y/z` on `use_pplim` (their other consumers — the UCT drift blend and the
  `rhs.cxx` `theta_tot` diagnostic — are guarded, since theta ≡ 1 when pplim is off).
  Golden-safe: all tests default `use_pplim=no`, and theta is not in the golden `.tsv`
  set. Files: `fluxes.cxx`, `rhs.cxx` (`AsterX_RHS_impl<pplim>`), `schedule.ccl`.
  **This is the reference build for every subsequent timing comparison.**
- **Idea 4** (hoist eig+UCT ahead of the assembly) — bit-identical, golden PASS, but
  **occupancy-null** (AGPR 40→40, a pure reorder) → REVERTED. The probe that nailed
  `CHECKPOINT.md` Lesson 1.
- `e8613114` **Idea 2** — momentum flux-by-flux, rolled `#pragma unroll 1 for(j)`.
  AGPR 40→1, VGPR 256→254.
- `210a013c` **Idea 3** — serialize the two face states, rolled `for(f)` with the
  momentum `for(j)` nested; `eigenvalues_oneside` split out in `eigenvalues.hxx`.
  AGPR 1→0 → **occ 2**. The CCTK_DEBUG NaN-dump is disabled (`#if 0`) because it
  references now-serialized quantities — rewrite before upstreaming this branch.
  Tradeoff: scratch 3448→4984 (AGPR→HBM), which is what made the timing
  grid-size-dependent. Design note: `ideas-2-3-serialization.md`.

## Closed probes (2026-07-25/27) — measured negative results

Not shipped; recorded so the shapes are not retried. Reasons: `CHECKPOINT.md`
Lessons 4 and 5.

- **`probe/flux-enthalpy-fusion` @ `659b48b9`** — the `(H,vf2)` algebraic fusion
  (`flux-construction.md` §10a/§11a): total enthalpy and fast speed absorb the whole
  magnetic sector, deleting `cs2_rc`/`h_rc`/`dens_h_W_rc`/`dens_h_W_plus_…` and
  demoting `B2_rc`/`bsq_rc` to transients, eigenvalue signature 15→9 doubles.
  **`-Rpass`: byte-for-byte identical (256/40/3448/occ-1).** Not bit-identical, so
  never golden-gated. Its *accuracy* findings are real and live on in `ideas.md`.
- **Compiler-flag sweep** — three mechanisms (scheduling, allocation, relaxed FP), all
  within ±2 registers of a 32-register gap; `--enable-deferred-spilling` also damaged
  `CalcE`. Full tables: `compiler-flags.md`. The `PROBE_FLAGS` switchboard and the
  candidate inventory were deleted rather than archived.
- **The `vbar` fix** — `if constexpr (pplim)` around the UCT drift blend, removing 12
  provably-dead `gf_vels` loads that LLVM cannot fold without `nnan`+`nsz`.
  **Register-null**, but KEPT: bit-identical, deletes 12 real global loads, removes a
  multiply-by-zero. No occupancy claim.

## Scoped `launch_bounds` — occ-2 at low scratch (2026-07-27)

The change that worked. Full detail: `launch-bounds-plan.md`.

- **CarpetX** `opt/loop-box-device-min-blocks` @ `19267243` (pushed): `loop_box_device`
  gains `int MB = 0` as a 7th template parameter; `MB == 0` keeps the existing
  `ParallelFor<NT>` path byte-for-byte, `MB != 0` replicates AMReX's
  `ParallelFor<MT>(BoxND,L)` body and launches through the already-existing 3-arg
  `launch_global<NT,MB>`. SYCL falls back to `ParallelFor`.
- **AsterX** `c96df5d5`: `MB=2` at the one fused flux site → `__launch_bounds__(256,2)`
  on `CalcFluxAll` alone. **⚠ Breaks CI against stock CarpetX by design** (the 7th
  template argument binds to `typename F`); revert or upstream CarpetX before golden.
- Result: `uct1` 256/32/3448/occ-1 → **128/128/3608/occ-2**, spills 0/0, all 9 EMF
  kernels in the same TU byte-for-byte unchanged. Timed **−29.1 %** (UCT-small) and
  **−27.3 %** (TOV-large) on `AsterX_Fluxes`, `Z4c_RHS` flat; cumulative **−47.5 %**
  against the pre-fusion original. Controls still owed (`CHECKPOINT.md` §STATUS).
