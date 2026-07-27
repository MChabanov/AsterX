# AsterX flux baseline timings (bitwise-stable FP flags)

Canonical **pre-refactor** performance baseline for the flux-loop work, taken
with the deterministic CI FP flags (`-ffp-contract=off`, no unsafe-math in
`CXX_OPTIMISE_FLAGS`) that all current and future runs use. This supersedes
the earlier baseline taken under the old value-changing flags (see
`implementation-history.md`).

## Provenance

- **Run:** 29158250436 (branch `tmp-golden-baseline`, deleted after use)
- **Code:** `c2614613` (pre-refactor flux code) + `b9c200ef` (FP flags only)
- **Date:** 2026-07-11, container `einsteintoolkit/carpetx:cpu-real64`
- **Ranks:** np1 = 1 rank / 2 threads, np2 = 2 ranks / 1 thread
- The same run's output regenerated the committed golden `.tsv` (`cc435a06`).

## ⚠ Methodology (unchanged)

Absolute CI wall-times carry **~±30 % run-to-run noise** (shared runners).
Judge changes by the **flux fraction** = `AsterX_Fluxes ÷ ODESolvers_Solve*`
(stable to ±1–3 pp; use `ODESolvers_Solve_Subcycling` for subcycling tests),
or by comparing two runs whose total `ODESolvers_Solve*` happens to agree. To
reproduce a table like the one below: download the `asterx-test-output`
artifact of a run and parse the final TimerReport block per test
(`scripts/print_timers.py` prints the raw numbers).

## Baseline flux fraction and absolute Fluxes time (run 29158250436)

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

Aggregates: np1 total Fluxes 81.8 s of Solve 303.1 s (27.0 %); np2 total
Fluxes 89.3 s of Solve 310.0 s (28.8 %).

## Measurements against this baseline so far

- **Step A + diagonal-B retirement** (run 29159401173): fractions within
  ±0.9 pp everywhere; total Fluxes −1.2 % (np1) / −2.0 % (np2) — no
  regression, consistent with dropping the wasted diagonal-B Riemann solve.
- **Step B fused sweep** (run 29161158922; runner ~5–8 % slower overall, use
  fractions): shock tubes −0.3..−1.0 pp, magTOV AMR +0.7..+1.4 pp, aggregate
  ≈ flat. Analysis in `ideas.md` §Why fusion alone shows no CPU win.
- **Step C fused zeroing** (run 29162200669): **golden-master PASS exactly
  0** for both np1 and np2 (`2235 ok, 0 fail, 0 missing`; worst
  `|abs| = |rel| = 0.000e+00`). Compared with Step B, aggregate Fluxes
  improved slightly (np1 −0.39 %, np2 −0.74 %) and the flux fraction dropped
  by 0.28 pp (np1) / 0.39 pp (np2). It partially recovered the Step B magTOV
  AMR overhead, especially for np2, but remains within CI noise.

## Recent optimization-run aggregates

These are final TimerReport-block sums over the 18 evolution tests, using
`AsterX_Fluxes ÷ ODESolvers_Solve*`.

| Run | Step | np1 Fluxes [s] | np1 Solve [s] | np1 flux % | np2 Fluxes [s] | np2 Solve [s] | np2 flux % |
|-----|------|----:|----:|----:|----:|----:|----:|
| 29159401173 (#286) | Step A + B-flux retirement | 80.8030 | 297.8766 | 27.13 | 87.5047 | 305.6041 | 28.63 |
| 29161158922 (#287) | Step B fused sweep | 89.0543 | 319.1754 | 27.90 | 97.1508 | 333.4862 | 29.13 |
| 29162200669 (#288) | Step C fused zeroing | 88.7110 | 321.1875 | 27.62 | 96.4272 | 335.4664 | 28.74 |

Step C vs Step B: aggregate Fluxes −0.39 % (np1) / −0.74 % (np2), and flux
fraction −0.28 pp (np1) / −0.39 pp (np2). Step C vs Step A is still higher
in absolute Fluxes time because runs #287/#288 landed on slower CI runners;
use the fractions for the signal.

## magTOV AMR flux-fraction check

The key regression watchpoint from Step B was the small-box AMR cases. Step C
reduced those fractions relative to Step B:

| Test | #286 np1 % | #287 np1 % | #288 np1 % | #286 np2 % | #287 np2 % | #288 np2 % |
|------|----:|----:|----:|----:|----:|----:|
| magTOV_Z4c_AMR | 22.30 | 23.50 | 22.91 | 23.40 | 24.05 | 23.43 |
| magTOV_Z4c_AMR_SC | 23.39 | 24.71 | 24.33 | 23.95 | 25.05 | 24.12 |

## GPU results — OLCF Frontier (the headline result), 2026-07-22

The CPU verdict above is cost-neutral **by design** (see `ideas.md` §Why
fusion alone shows no CPU win): CI grids are cache-resident and small-box
write-stream costs offset the launch savings. `CHECKPOINT §Next steps #1`
deferred the real test — **6→1 kernel launches + occupancy** — to GPU
hardware at production grid size. This is that measurement.

- **Machine:** OLCF Frontier, **1 node = 8 GPUs** (4× MI250X, 8 GCDs), ROCm
  build. Production-scale grids (not the CI test pars).
- **Comparison:** `opt/flux-loops` (pre-fusion baseline, tip `c2614613` — has
  the TimerReport infra so the timer blocks are byte-for-byte comparable) vs
  `opt/flux-loops-impl` (fused single sweep, Steps A–C). The baseline branch
  lacks only the FP-determinism flags (`b9c200ef`), which are cpu-cfg-only and
  do not affect ROCm timings. Apples-to-apples.
- **Method:** side-by-side final TimerReport block, same iteration in each
  pair. Metric = flux fraction `AsterX_Fluxes ÷ ODESolvers::Solve*`
  (`_Subcycling` for the subcycling test), plus absolute `AsterX_Fluxes` — on
  real GPU the launch reduction shows in absolute time, unlike CPU.

| Test (it) | CT | Fluxes base [s] | Fluxes fused [s] | Fluxes Δ | frac base | frac fused | frac Δ | CCTK total Δ |
|-----------|----|----:|----:|----:|----:|----:|----:|----:|
| Fixed-metric subcycling (3488) | UCT | 169.790 | 124.794 | **−26.5 %** | 37.74 % | 30.98 % | **−6.76 pp** | −7.5 % |
| TOV dynamic Z4c (128) | flux-CT | 237.357 | 171.622 | **−27.7 %** | 47.63 % | 39.89 % | **−7.74 pp** | −10.6 % |

(CT scheme per user, 2026-07-22: the fixed-metric/AnalyticalSpacetime test runs
upwind-CT — production's usual mode; the TOV test was run flux-CT. The −27 %
fusion win holds for both.)

**Verdict: fusion is a large GPU win, ~27 % off the flux kernel, stable across
problem types.** In both runs the flux delta accounts for essentially the
*entire* improvement upstream — the ODESolvers rhs and Solve deltas track the
flux delta to within ~1–2 s, and everything else (Z4c_RHS, SourceTerms, Tmunu,
Con2Prim, OutputSilo/Norms) is flat within noise — so this is unambiguously the
flux fusion, not runner variance or a shift elsewhere:

- **Fixed-metric subcycling** (it 3488): Fluxes 169.790→124.794 s (−45.0 s);
  rhs 230.860→185.245 (−45.6 s); Solve_Subcycling 449.842→402.750 (−47.1 s);
  CCTK total 608.995→563.045 (−46.0 s).
- **TOV dynamic Z4c** (it 128): Fluxes 237.357→171.622 s (−65.7 s); rhs
  408.159→342.543 (−65.6 s); Solve 498.430→430.279 (−68.2 s); CCTK total
  658.478→588.706 (−69.8 s). Here flux is the *dominant* RHS cost (47.6 % of
  Solve at baseline), even under a full spacetime evolution, so the fusion
  yields the larger flux-fraction drop and the bigger total win (−10.6 %).

This validates the core hypothesis of the whole flux-loop effort: the refactor
that was cost-neutral on CI CPU (±0.3 pp) is worth **~27 % of the flux kernel /
6.8–7.7 pp of flux fraction** on Frontier GPUs at production scale. Caveats:
run-IDs are not captured (Frontier runs, not CI); exact grid sizes and rank
decomposition not recorded here; two problem types measured — consistent, but
add MHD / AMR cases if a third run is available.

### Occupancy follow-up (headroom beyond fusion)

Forcing the flux kernel occ 1→2 (global `__launch_bounds__(MT,2)`) on the TOV
flux-CT run cut **AsterX_Fluxes −15.9 %** (171.622→144.283 s) — but the *global*
directive regressed Z4c_RHS (+48 %) and others → net total +7.9 % worse, so the
win must be **scoped to the flux kernel**. Confirms the flux kernel is
occupancy-bound; projected scoped gain ~−4.6 % total on top of fusion. Full
investigation, EOS-by-EOS occupancy table, and the register-reduction plan:
`profiling-flux-kernel.md`.

### occ-2 (serialization, Ideas 1/2/3, `210a013c`) — timing: GRID-SIZE-DEPENDENT (2026-07-24)

The register-reduction stack reached **occ-2 without `launch_bounds`** (`-Rpass`:
VGPR 253, AGPR 0, occ 2; but scratch rose 3448→4984 B/lane as overflow moved
register→HBM — see `CHECKPOINT.md` §KEY LESSONS 3). Two timed Frontier runs
(occ-2 build vs fused baseline), both ideal-gas:

| run | grid | CT | AsterX_Fluxes occ-2 | baseline | Δ Fluxes | Δ fraction |
|---|---|---|---:|---:|---:|---:|
| TOV Z4c | larger | flux-CT | 163.3 s | 171.6 s | **−4.8 %** | −1.34 pp |
| fixed-metric subcycling | smaller | UCT | 137.2 s | 125.9 s | **+8.9 %** | +1.82 pp |

**Sign flip:** occ-2 WON on the larger grid (enough memory latency to hide the
+45% scratch) and REGRESSED on the smaller grid (scratch bandwidth dominates).
Z4c_RHS was flat in the TOV run (83.27 vs 83.31) — the serialization is scoped to
the flux kernel by construction, so no collateral (unlike the earlier *global*
`__launch_bounds__` test). Caveat/confound: the two runs differ in grid AND CT
scheme AND subcycling, but both kernels reached occ-2 with near-identical
registers/scratch, so grid/latency regime is the likely driver.

Crucially, the *genuine* occ-2 (−4.8%) is far short of the *forced* occ-2
(−15.9%, global-launch_bounds on the same TOV run) — the gap is the scratch. So
the lever is to reach occ-2 without the wholesale-demotion scratch.

**Update 2026-07-25 (SUPERSEDED — see the `launch_bounds` section below) — the
lever is now compiler flags, not `launch_bounds`.** See `compiler-flags.md`;
`launch-bounds-plan.md` is the fallback. Whichever route wins, its timings (both
grid sizes, vs the Idea-1 baseline) go here next. Three things to carry into that
measurement:

- **The reduced probe build IS usable for timing** (idealgas, `use_pplim=no`, both
  CT schemes compiled) — but it is a *different compilation* from a full build
  (AGPR 32 vs 40), so time the FULL build, not the probe build.
- **occ-2 and scratch are independent** (`CHECKPOINT.md` Lesson 3 correction): the
  baseline 3448 B/lane is `alloca`, not spill, and does not affect occupancy. So
  "beat 4984" is the criterion for the small-grid regression, and "clear the
  256-register cliff" is the criterion for occupancy. A change can win one and lose
  the other.
- **The two runs below are confounded** (grid AND CT scheme AND subcycling all
  differ, per the caveat above). Now that `uct=0` is measured at occ-2 and `uct=1`
  at occ-1 in the same build, that confound is sharper than it looked: decide
  whether a matched pair is needed before drawing a grid-size conclusion.

### ⭐ occ-2 via scoped `launch_bounds(256,2)` — fixed-metric subcycling / UCT / smaller grid (2026-07-27)

CarpetX `opt/loop-box-device-min-blocks` @ `19267243` + AsterX `c96df5d5`
(`MB=2` at the fused flux site). `-Rpass`: `uct1` 128 VGPR / 128 AGPR / 3608
B/lane / **occ 2**, spills 0/0 (`launch-bounds-plan.md` §TRIAL RESULT).
Side-by-side TimerReport, same par file, same iteration count:

| timer | occ-2 (`MB=2`) | baseline | Δ |
|---|---:|---:|---:|
| **AsterX_Fluxes** | **81.958 s** | **115.639 s** | **−29.1 %** |
| ODESolvers::Solve::rhs | 137.648 | 171.656 | −19.8 % |
| ODESolvers::Solve | 339.148 | 373.177 | −9.1 % |
| CCTK total time | 489.082 | 521.988 | −6.3 % |
| Fluxes as fraction of Solve | 24.2 % | 31.0 % | **−6.8 pp** |

**The win is isolated to the flux kernel, and the bookkeeping closes.** Fluxes
−33.68 s, `Solve::rhs` −34.01 s, `Solve` −34.03 s: the entire RHS improvement is
the flux kernel and nothing else. Every unrelated timer is flat inside ±1 % —
`SetMetric` 99.244 vs 99.502, `SourceTerms` 30.471 vs 30.441, `Con2Prim` 17.818 vs
17.866, `CalcAux` 10.374 vs 10.443, `OutputNorms` 89.578 vs 88.926, `AsterX_RHS`
14.685 vs 14.967. That flatness is also the evidence the two runs are comparable
in size and iteration count, and it independently confirms the scoping seen in
`-Rpass`: no collateral anywhere (`Initialise` 5.8 vs 4.1 is startup, ignore).

**This is the configuration the serialization LOST on.** Same fixed-metric
subcycling / UCT / smaller grid where Ideas 2/3 measured **+8.9 %**; scoped
`launch_bounds` measures **−29.1 %**. Sign flipped and multiplied. It also beats
the old *global* forced-`launch_bounds` reference (−15.9 %, on TOV) in relative
terms, without that test's Z4c collateral. Mechanism, consistent with everything
in `CHECKPOINT.md` Lesson 3: this kernel is latency-bound, a second resident wave
buys latency hiding, and the +160 B/lane of extra scratch is small enough not to
eat it — unlike the serialization's +1536.

**⚠ CONFOUND — two changes at once; the headline is not yet a clean measurement.**
The occ-2 run is the **reduced** build (`-DASTERX_PROBE_IDEALGAS_ONLY`, 2
instantiations) while the baseline is a **full** build (the `use_pplim`
compile-time commit). The methodology note above says to time the full build, and
this run does not. So −29.1 % is (instantiation cut + `MB=2`) versus neither.

How much could the cut alone be worth? Probably little: at `MB=0` the reduced
build was still 256/**32**/3448/**occ-1**, i.e. the cut changed neither occupancy
nor scratch, and only 8 AGPRs separate it from the full build's 40 — both occ-1.
So the occupancy doubling is the plausible cause of essentially all of it. But
"plausible" is not "measured", and the reduced build IS a different compilation
(different inlining of the shared ReconX/EOS callees).

**The control that settles it is cheap: same reduced build, flip `MB` to 0, rerun
the same par file.** That isolates `launch_bounds` from the instantiation cut in
one run, and it doubles as the sanity check that `MB=0` reproduces occ-1. Do this
before quoting −29 % anywhere.

Also still owed: **TOV-large / flux-CT** (where `uct0` was already occ-2 in the
reduced build, so expect neutral-to-slightly-negative — a flat result there is not
failure), a **full-build** re-measurement, and the **golden gate** (expected
exactly 0; `launch_bounds` is a codegen directive).

Cross-run caveat: the small-grid baseline in the serialization table above reads
125.9 s for `AsterX_Fluxes` where this baseline reads 115.639 s — different
baseline builds/runs, so compare the two experiments by **fraction, not seconds**.
