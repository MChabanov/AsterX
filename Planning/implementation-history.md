# Implementation history — AsterX flux-loop optimization

Compact provenance record. Working branch: `opt/flux-loops-impl` (stacked on
`opt/flux-loops`, based on `main` @ `bc6f2d2f`). Remotes: `origin` =
MChabanov/AsterX (CI runs here), `etk` = EinsteinToolkit/AsterX. Full design
discussions lived in per-target plan files (deleted after completion; git
history of `Planning/` is NOT kept — this file is the record).

## Supporting infrastructure (branch `opt/flux-loops`, 2026-07-10)

- **Golden-master harness (Phase 0)**: `scripts/compare_tsv.py` — stdlib-only
  comparator of produced vs committed test `.tsv`, keyed by grid locator
  (handles 1p/2p row reordering), `|Δ| ≤ atol + rtol·|b|`, default 1e-13.
  Wired into CI as an opt-in step. Committed data regenerated from the CI
  container (`bea0b913`) since data from other machines drifts ~1e-9.
  Decisions: reuse the existing 37 test pars and their committed `.tsv` as
  the golden master (no new tests); TSV precision fixed at 16 digits
  (hardcoded in CarpetX `io_tsv.cxx`, not tunable) — accepted.
- **Per-routine timing (Target 2)**: all 18 evolution test pars activate
  `TimerReport` (`out_every=50`, readable output, top 50);
  `scripts/print_timers.py` prints the hottest AsterX routines per test in
  the CI log; full timer files ride in the `asterx-test-output` artifact.
  (`cctk_timer_output` alone was insufficient — hides AsterX_Fluxes inside
  ODESolvers.)
- **ccache CI (Target 3)**: `scripts/build.sh` + `actions/cache` on
  `.ccache`; warm cuda/rocm builds ~25 min → ~6–7 min.
- Commits: `bae19fe5`, `bc162f74`, `eadd1fa7`, `f882c088`, `5830bc95`,
  `bea0b913`, `f4a469f7`, `c2614613`.

## Target 1 — fused single-sweep flux calculation (branch `opt/flux-loops-impl`, 2026-07-11)

Hard requirement: bit-identical results, each step gated by golden diff = 0.
Fidelity spec for the extraction: verbatim copy-paste, only sanctioned
substitutions (`p.X` → explicit `face_X`; captures → `FluxContext` members).

- `0930f4c3` **face_X made explicit** in the flux loop (used by the
  atmosphere-radius code, the one place needing the *face* coordinate).
- `b6739ec9` **FluxContext<EOSType>** aggregate: everything the per-face code
  needs, members named as the locals they carry; positional aggregate init.
- `9fcb41d6` **Step A extraction**: loop body (1002 lines) moved verbatim
  (mechanically verified, only a uniform 2-space dedent) into
  `CalcFluxAtFace<dir_i, EOSType>(fx, p, face_X)`; dir_i stays compile-time.
- `e25ce7ce` cuda fix: `maxspeeds_from_lambdas`/`avg_upwind` in AsterUtils
  were the only helpers without `CCTK_HOST`; nvcc rejects device-only calls
  from the host path of the H+D helper. (`hll_upwind` is still device-only —
  will need the same when the EMF sweeps are fused.)
- `1bc276a5` **diagonal B-fluxes retired**: fxBx/fyBy/fzBz were identically
  zero ((dir_i != i) * flux) and unread (RHS consumes only hydro fluxes; B
  evolves via Avec; flux-CT CalcE reads only off-diagonals). Replaced
  variable-indexed fluxB{x,y,z}s vecs by direction-indexed off-diagonal pairs
  `fluxB_j = {fxBy, fyBz, fzBx}`, `fluxB_k = {fxBz, fyBx, fzBy}` (vbar_j/k
  idiom); CalcE mapping (cyclic (i,j,k)): gf_fBs(j)(k) = gf_fB_k(k),
  gf_fBs(k)(j) = gf_fB_j(j). Drops 3/30 flux GFs and 1/10 Riemann evals/face.
- `0d1f41eb` CI: `[golden-master]` in the HEAD commit message of a push
  triggers the golden step (besides dispatch input and PR label).

### FP-determinism lesson (critical)

The first golden run on the provably-verbatim Step A FAILED: 123/2235 files,
ulp-level early deviations growing to ~1e-6 rel (shock fronts amplify
last-bit perturbations; limiter branches flip on thresholds). Cause: the CI
cpu cfg compiled C++ with `-funsafe-math-optimizations -ffp-contract=fast` —
under those flags the compiler may legally change FP results whenever code is
restructured, so bit-identity across ANY refactor is unattainable.

- `b9c200ef` `scripts/actions-cpu-real64.cfg`: CXX flags → `-ffp-contract=off`,
  no unsafe-math (C/Fortran and cuda/rocm cfgs unchanged; CI-only config).
- `cc435a06` golden `.tsv` regenerated from PRE-refactor code under the new
  flags (run 29158250436 on temp branch `tmp-golden-baseline`; 1-proc ==
  2-proc byte-identical; 2030 files replaced, 205 unchanged). This run is
  also the canonical performance baseline (`baseline-timings.md`).
- Validation run 29159401173 then reported **2235 ok, 0 fail, worst
  |abs|=|rel|=0** for Step A + the B-flux retirement. gh auth note: the gh
  keyring token expired mid-session once; `gh auth login` fixes it.

### Steps B and C — the fusion

- `1b6ca2ee` **Step B**: three per-direction sweeps → one fused
  `CalcFluxAll<EOSType>`. Per-direction domains reproduced exactly by
  `calc_mixpn_box<CI,CJ,CK>` (same `box_int`/`box_all` ± `nloop`-after-tile-
  clamping computation as CarpetX's `loop_mixpn_device`, on the same grid
  object); fused loop = `loop_box_device<0,0,0>` over the per-dim union with
  per-direction guards. `face_X` computed per direction with CarpetX's own
  point_desc formula `x0 + (lbnd + I - (!CI)/2) * dx` → bit-identical to the
  old face-centred `p.X`. Verified prerequisites: flux body + ReconX use only
  the centering-independent `p.I/p.DI/p.DX`; no cross-face/cross-direction
  data flow. **Golden: exactly 0** (run 29161158922).
- `6728188b` **Step C1**: three zeroing kernels → one (runtime dir loop,
  guards from the same `box_all` boxes). Launches 4 → 2.
- `1ce741de` **Step C2**: zeroing merged into the main sweep → **1 kernel
  launch per AsterX_Fluxes call**; eliminates the separate full write pass
  over the face GFs (zero + value coalesce in cache). Value-identity: the
  flux code never reads an init-written value (gf_theta and own-face flux
  reads are preceded by same-point writes in both PP-limiter branches);
  cross-tile overlaps idempotent. **Golden: exactly 0** (run 29162200669;
  both np1 and np2: `2235 ok, 0 fail, 0 missing`, worst
  `|abs|=|rel|=0.000e+00`).

### Timing verdicts (CPU, CI scale — see baseline-timings.md)

Step A + B-flux retirement: no regression (fractions ±0.9 pp; total Fluxes
−1..−2 %). Step B fusion: ≈ cost-neutral (shock tubes slightly better, small
AMR boxes slightly worse); analysis of why in `ideas.md` §Why fusion alone
shows no CPU win. Step C fused zeroing slightly improved over Step B:
aggregate Fluxes −0.39 % (np1) / −0.74 % (np2), flux fraction −0.28 pp
(np1) / −0.39 pp (np2), and the Step B magTOV AMR overhead was partially
recovered.

### GPU validation + occupancy investigation (2026-07-22, OLCF Frontier)

Not committed code — measurement/analysis; recorded in `baseline-timings.md`
§GPU and `profiling-flux-kernel.md`. Fusion validated on GPU: `AsterX_Fluxes`
**−26.5 %** (fixed-metric, UCT) / **−27.7 %** (TOV, flux-CT). Occupancy
profiling: the fused flux kernel is occ-1 / register-bound for the production
EOS (tabulated3d, ideal-gas), occ-1 pre-existing (fusion cost no occupancy);
forcing occ 1→2 gave −15.9 %. Peak is the MHD flux-assembly live-set (EOS math
and PP block exonerated). **Next work item** (see CHECKPOINT §Next steps): the
config-gated CT-scheme split (`CalcFluxAll<EOS, use_uct>`) — production is UCT,
so retire B-fluxes (the promoted `ideas.md` Tier-2); flux-CT retires the 12 UCT
face GFs. The neighbourhood-reuse pass is shelved (would backfire on occ-1).

## Register reduction to occ-2 (branch `opt/flux-eig-collapse`, 2026-07-23/24)

Goal: drain the ideal-gas AGPR overflow to reach occ-2. All bit-identical
(golden PASS at 0, cpu/rocm/cuda). Full -Rpass tables: `profiling-flux-kernel.md`
§Serialization result. Result: **production `<uct=1,pplim=0,idealgas>` occ 1->2**.

- `b23fb947`/`e8b0718c`/`4a7e40ec` **CT-scheme split** (compute + hoist + storage
  gating), `CalcFluxAll<uct,EOS>`. Golden PASS; occupancy-null (AGPR 74->68) but a
  storage/bandwidth win + the enabler for what follows.
- `549a28ca` **eigenvalue collapse** — `eigenvalues()` degenerate (4 stored = 2
  distinct/side); every consumer reduces to `{charmax=max(0,·), charmin=min(0,·)}`.
  Collapsed to those 2 scalars, threaded through calcflux/laxf/UCT. Bit-identical
  (max/min/fabs exact). ~2-4 AGPR.
- `defd9f60` **Idea 1** — template `use_pplim` (like `uct`): `if constexpr(pplim)`
  around the whole PP-limiter block -> for use_pplim=no (production) the `_ppl`
  live-set compiles out. AGPR 64->40. Also storage-gated `theta_x/y/z` on
  use_pplim (only PP-associated GFs; also read by the UCT drift blend + the rhs
  `theta_tot` diagnostic, both guarded since theta==1 when pplim off). Golden-safe:
  all tests default use_pplim=no; theta not in the golden .tsv set. Files:
  fluxes.cxx, rhs.cxx (AsterX_RHS_impl<pplim>), schedule.ccl.
- **Idea 4 (hoist eig+UCT block ahead of the assembly)** — bit-identical + golden
  PASS but **occupancy-null** (production AGPR 40->40, only reorder) -> REVERTED.
  The probe that nailed Lesson 1 (allocator ignores reorders).
- `e8613114` **Idea 2** — momentum flux-by-flux, rolled `#pragma unroll 1 for(j)`.
  AGPR 40->1, VGPR 256->254. (CCTK_DEBUG recomputes the full moms/flux_moms.)
- `210a013c` **Idea 3** — serialize the two face states, rolled `for(f)` (group-B
  as per-side scalars) with the momentum `for(j)` nested; `eigenvalues_oneside`
  split in `eigenvalues.hxx`. AGPR 1->0 -> **occ 2**. CCTK_DEBUG NaN-dump disabled
  (`#if 0`, references serialized quantities; rewrite before upstream).
- Tradeoff: scratch 3448->4984 B/lane (AGPR->HBM); hidden at occ-2. Timed Frontier
  run pending to confirm the wall-clock win (~-16% expected).
