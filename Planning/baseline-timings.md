# AsterX flux timings — baseline and all measured results

Canonical performance record for the flux-loop work, taken under the deterministic
CI FP flags (`-ffp-contract=off`, no unsafe-math) that all runs use. Supersedes the
earlier baseline taken under the old value-changing flags
(`implementation-history.md`).

## Provenance and methodology

- **CPU baseline run:** 29158250436 (branch `tmp-golden-baseline`, deleted after
  use), code `c2614613` (pre-refactor) + `b9c200ef` (FP flags only), 2026-07-11,
  container `einsteintoolkit/carpetx:cpu-real64`. np1 = 1 rank/2 threads, np2 = 2
  ranks/1 thread. The same run's output regenerated the committed golden `.tsv`
  (`cc435a06`).
- **⚠ CI wall-times carry ~±30 % run-to-run noise** (shared runners). Judge CPU
  changes by the **flux fraction** = `AsterX_Fluxes ÷ ODESolvers_Solve*` (stable to
  ±1–3 pp; `_Subcycling` for subcycling tests), or by comparing runs whose total
  `Solve*` happens to agree. To reproduce: download the `asterx-test-output`
  artifact and run `scripts/print_timers.py`.
- **GPU runs (Frontier)** are side-by-side final TimerReport blocks at the same
  iteration; absolute `AsterX_Fluxes` is meaningful there, unlike CI. Confirm
  comparability by checking that unrelated timers are flat.

## CPU baseline per test (run 29158250436)

| Test | it | np1 flux % | np2 flux % | np1 Fluxes [s] | np2 Fluxes [s] |
|------|---:|----:|----:|----:|----:|
| Alfven_wave | 801 | 31.8 | 31.2 | 4.2439 | 3.3059 |
| Balsara1_shocktube_xdir_godunov | 80 | 25.1 | 30.0 | 0.9207 | 1.1156 |
| Balsara1_shocktube_xdir_monocentral | 80 | 23.6 | 29.3 | 1.0295 | 1.1765 |
| Balsara1_shocktube_xdir_mp5 | 80 | 29.3 | 30.3 | 1.6296 | 2.0548 |
| Balsara1_shocktube_xdir_ppm | 80 | 16.9 | 19.1 | 1.8439 | 2.3492 |
| Balsara1_shocktube_xdir_ppm_PalenzuelaC2P | 80 | 29.8 | 30.0 | 1.8171 | 2.4153 |
| Balsara1_shocktube_xdir_wenoz | 80 | 31.7 | 32.5 | 1.8281 | 2.3381 |
| Balsara1_shocktube_ydir_ppm | 80 | 31.9 | 32.4 | 1.8571 | 2.2559 |
| Balsara1_shocktube_zdir_ppm | 80 | 32.0 | 32.5 | 1.8512 | 2.2524 |
| Balsara2_shocktube_xdir_ppm | 80 | 38.8 | 39.4 | 1.8310 | 2.2442 |
| Balsara3_shocktube_xdir_minmod_LxF | 129 | 34.6 | 35.5 | 2.4948 | 2.7720 |
| Balsara4_shocktube_xdir_ppm | 80 | 11.3 | 14.2 | 1.8199 | 2.2079 |
| Balsara5_shocktube_xdir_ppm | 110 | 39.1 | 39.5 | 2.5145 | 3.0702 |
| Cylindrical_blast | 96 | 40.3 | 40.2 | 9.7776 | 13.0685 |
| Magnetic_rotor | 72 | 46.6 | 44.9 | 9.6093 | 11.8189 |
| Sound_wave | 10 | 32.1 | 30.8 | 0.4474 | 0.3168 |
| magTOV_Z4c_AMR | 100 | 22.1 | 23.3 | 21.8045 | 21.8673 |
| magTOV_Z4c_AMR_SC | 100 | 23.3 | 23.9 | 14.4883 | 12.6992 |

Aggregates: np1 81.8 s Fluxes of 303.1 s Solve (27.0 %); np2 89.3 / 310.0 (28.8 %).

## CPU verdict on the fusion: cost-neutral, by design

Compressed from three CI runs (#286 Step A + B-flux retirement, #287 Step B fused
sweep, #288 Step C fused zeroing):

| Run | Step | np1 Fluxes / Solve / % | np2 Fluxes / Solve / % |
|-----|------|---|---|
| 29159401173 | A + B-flux retirement | 80.80 / 297.88 / 27.13 | 87.50 / 305.60 / 28.63 |
| 29161158922 | B fused sweep | 89.05 / 319.18 / 27.90 | 97.15 / 333.49 / 29.13 |
| 29162200669 | C fused zeroing | 88.71 / 321.19 / 27.62 | 96.43 / 335.47 / 28.74 |

Fractions moved ±0.9 pp at most; absolute times track runner speed, not code. Step C
vs B: Fluxes −0.39 % (np1) / −0.74 % (np2), fraction −0.28 / −0.39 pp, which
partially recovered Step B's small-box magTOV AMR overhead (`magTOV_Z4c_AMR` np1
22.30 → 23.50 → 22.91 %; `_SC` 23.39 → 24.71 → 24.33 %). Step C golden PASS at
exactly 0 for np1 and np2. **Why CPU shows nothing** (cache-resident CI boxes,
write-stream costs offsetting launch savings): `ideas.md`.

## GPU — the fusion (2026-07-22, OLCF Frontier, 1 node / 8 GPUs, production grids)

`opt/flux-loops` (pre-fusion, tip `c2614613`) vs `opt/flux-loops-impl` (fused,
Steps A–C). The baseline branch lacks only the CI-only FP flags, so apples-to-apples.

| Test (it) | CT | Fluxes base | fused | Δ | frac base → fused | CCTK total Δ |
|-----------|----|----:|----:|----:|---|----:|
| Fixed-metric subcycling (3488) | UCT | 169.790 s | 124.794 | **−26.5 %** | 37.74 → 30.98 % (−6.76 pp) | −7.5 % |
| TOV dynamic Z4c (128) | flux-CT | 237.357 s | 171.622 | **−27.7 %** | 47.63 → 39.89 % (−7.74 pp) | −10.6 % |

**~27 % off the flux kernel, stable across problem types.** In both runs the flux
delta accounts for essentially the entire upstream improvement — rhs and Solve
deltas track it to within 1–2 s, and `Z4c_RHS`, `SourceTerms`, `Tmunu`, `Con2Prim`,
`OutputSilo/Norms` are flat within noise. Subcycling: Fluxes −45.0 s, rhs −45.6,
Solve_Subcycling −47.1, total −46.0. TOV: Fluxes −65.7 s, rhs −65.6, Solve −68.2,
total −69.8. Caveats: Frontier runs, so no run IDs; grid sizes and rank decomposition
not recorded; two problem types measured.

## GPU — occupancy

### Global `__launch_bounds__(MT,2)` probe (the motivation)

Forcing occ 1→2 globally on the TOV flux-CT run cut `AsterX_Fluxes` **−15.9 %**
(171.622 → 144.283 s) — but the *global* directive regressed `Z4c_RHS` +48 % and
others, for a net total **+7.9 % worse**. Confirms the flux kernel is
occupancy-bound and that the directive must be **scoped**. ⚠ The code state under
this probe is not recorded (probably pre-Idea-1), so treat −15.9 % as indicative
rather than a matched measurement.

### occ-2 via serialization (Ideas 1/2/3, `210a013c`) — grid-size-dependent

occ-2 without `launch_bounds` (VGPR 253 / AGPR 0), but scratch 3448 → 4984 B/lane:

| run | grid | CT | Fluxes occ-2 | baseline | Δ | Δ fraction |
|---|---|---|---:|---:|---:|---:|
| TOV Z4c | larger | flux-CT | 163.3 s | 171.6 s | **−4.8 %** | −1.34 pp |
| fixed-metric subcycling | smaller | UCT | 137.2 s | 125.9 s | **+8.9 %** | +1.82 pp |

**Sign flip:** occ-2 won where there was enough memory latency to hide the +45 %
scratch and lost where scratch bandwidth dominated. `Z4c_RHS` flat (83.27 vs 83.31),
so the serialization is scoped by construction. Confound: the two runs differ in grid
*and* CT scheme *and* subcycling. The genuine occ-2 (−4.8 %) fell far short of the
forced occ-2 (−15.9 %) — the gap is the scratch, which set up the whole
`launch_bounds` route.

### ⭐ occ-2 via scoped `launch_bounds(256,2)` (2026-07-27) — BOTH grids win

CarpetX `opt/loop-box-device-min-blocks` @ `19267243` + AsterX `c96df5d5`. `-Rpass`:
`uct1` 128 VGPR / 128 AGPR / 3608 B/lane / occ 2, spills 0/0
(`launch-bounds-plan.md` §TRIAL RESULT).

**Baseline for both runs is the same occ-1 full build at `defd9f60`** (`uct1` =
256/40/3448/occ-1) — *not* an occ-2 build, and unrelated to the serialization branch.
Same reference for both grid sizes, which is what makes the two Δ's comparable.

| timer | UCT-small `MB=2` → base | Δ | TOV-large `MB=2` → base | Δ |
|---|---|---:|---|---:|
| **AsterX_Fluxes** | 81.958 → 115.639 s | **−29.1 %** | 124.848 → 171.622 s | **−27.3 %** |
| Solve::rhs | 137.648 → 171.656 | −19.8 % | 295.379 → 342.543 | −13.8 % |
| ODESolvers::Solve | 339.148 → 373.177 | −9.1 % | 382.350 → 430.279 | −11.1 % |
| CCTK total | 489.082 → 521.988 | −6.3 % | 542.565 → 588.706 | −7.8 % |
| Fluxes / Solve | 24.2 → 31.0 % | −6.8 pp | 32.7 → 39.9 % | −7.2 pp |

**Isolated to the flux kernel in both, and the bookkeeping closes.** UCT-small:
Fluxes −33.68 s, rhs −34.01, Solve −34.03. TOV-large: Fluxes −46.77 s, rhs −47.16,
Solve −47.93 (= rhs + poststep −0.76). Unrelated timers flat inside ±1 %: UCT-small
`SetMetric` 99.244/99.502, `SourceTerms` 30.471/30.441, `Con2Prim` 17.818/17.866,
`CalcAux` 10.374/10.443; TOV-large **`Z4c_RHS` 83.189/83.313**, `SourceTerms`
56.667/56.613, `Tmunu` 32.670/32.644, `OutputSilo` 88.218/88.233. That flatness also
evidences the runs are comparable in size and iteration count, and confirms at the
timing level the zero collateral seen in `-Rpass`. (`OutputNorms`, `OutputGH`,
`Initialise`, `Sync` swing a few percent — I/O and startup variance.)

**⭐ The two occ-2 routes, head-to-head on the same TOV baseline:**

| route to occ-2 on TOV-large | scratch B/lane | Fluxes | Δ vs `defd9f60` |
|---|---:|---:|---:|
| serialization (rolled loops, `210a013c`) | 4984 | 163.3 s | −4.8 % |
| **scoped `launch_bounds(256,2)`** | **3608** | **124.8 s** | **−27.3 %** |
| (global `__launch_bounds__`, indicative) | — | ~144 s | −15.9 % |

**Same occupancy, 5.7× the benefit, and the only difference is where the register
overflow lives.** The cleanest confirmation of `CHECKPOINT.md` Lesson 3: occ-2 was
never the problem, the scratch was. It also beats the global probe while leaving
`Z4c_RHS` untouched, which that probe did not.

**The grid-size sign flip is gone**: −29.1 % small vs −27.3 % large, within 2 pp.
The serialization's dramatic flip (−4.8 % vs +8.9 %) was a scratch-bandwidth
artifact; at +160 B/lane instead of +1536 the grid dependence largely disappears.

**⚠ CONFOUND, applies to both runs.** The `MB=2` runs are the **reduced** build
(`-DASTERX_PROBE_IDEALGAS_ONLY`) while the baseline is a **full** build, contrary to
this file's own methodology note, so the headline is (instantiation cut + `MB=2`).
The cut alone is probably worth little — at `MB=0` the reduced build was still
256/32/3448/occ-1, changing neither occupancy nor scratch, and only 8 AGPRs separate
it from the full build's 40, both occ-1 — but that is inference. **The control is one
cheap run: same reduced build, `MB=0`, same par file.** It covers both configs. Do it
before quoting these numbers. Also owed: a full-build re-measurement including
`uct0`'s never-recorded baseline, and the golden gate (expected exactly 0).

## ⭐⭐ Cumulative: pre-fusion original → fused + occ-2

Both configurations against the **original upstream code, before the fused flux
loop**, with `defd9f60` as the midpoint. Same iteration in every column (TimerReport
headers agree: it 3232 / t = 484.8 for UCT-small, it 128 / t = 5.01953 for TOV):

| config | Fluxes original | → fusion + Idea 1 | → + `launch_bounds` occ-2 | **cumulative** |
|---|---:|---:|---:|---:|
| fixed-metric subcycling / UCT / smaller | 157.335 s | 115.639 (−26.5 %) | **81.958 (−29.1 %)** | **−47.9 %** |
| TOV Z4c / flux-CT / larger | 237.357 s | 171.622 (−27.7 %) | **124.848 (−27.3 %)** | **−47.4 %** |

| config | CCTK total | Δ | ODESolvers::Solve | Δ | Fluxes / Solve |
|---|---:|---:|---:|---:|---|
| UCT-small | 489.082 vs 564.423 | **−13.4 %** | 339.148 vs 416.709 | −18.6 % | 24.2 vs 37.8 % (−13.6 pp) |
| TOV-large | 542.565 vs 658.478 | **−17.6 %** | 382.350 vs 498.430 | −23.3 % | 32.7 vs 47.6 % (−15.0 pp) |

**Two independent configurations, essentially the same decomposition:** ~−27 % from
the source work, ~−28 % from occupancy, compounding to **−47.5 %** in both. The
halves being near-equal in both configs, and the two configs agreeing to 0.5 pp, is
strong evidence neither number is an artifact. Non-flux timers are flat against the
original too (UCT-small `SetMetric` 99.244/99.642, `SourceTerms` 30.471/30.680,
`Con2Prim` 22.968/23.061; TOV-large `Z4c_RHS` 83.189/83.229, `OutputSilo`
88.218/88.124, `Tmunu` 32.670/32.636), so the whole gain is the flux kernel.

**The flux kernel is no longer the dominant cost** — 47.6 % → 32.7 % of Solve on TOV,
37.8 % → 24.2 % on subcycling. Weigh that before further flux-kernel work.
