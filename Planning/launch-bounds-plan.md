# Scoped launch_bounds plan — reach occ-2 WITHOUT the serialization scratch

**Status (2026-07-27): PROMOTED BACK TO PRIMARY.** The compiler-flag route is
closed on measurement — all three probes ran and `uct1` stayed at occ-1
(scheduling 0 registers, allocation −2 with `CalcE` collateral, relaxed FP +2;
full tables in `compiler-flags.md` §SWEEP RESULTS). Two results from that sweep
bear directly on this plan:

- **The 32 AGPRs are load-bearing live values, not allocator sloppiness** — so the
  premise here is the right one: don't ask the allocator to find registers that
  aren't there, *force* the occupancy target and let it spill selectively.
- **A TU-scoped flag cannot be aimed at the flux kernel** — `CalcFluxAll`,
  `CalcE_impl`, `CalcFstag` and `CalcAux` all live in `fluxes.cxx`, and probe 2
  measurably damaged `CalcE` (SGPR spills 20/22/24 → 30/26/33) while helping the
  target by 2. `launch_bounds` is per-launch-site by construction, so it is not
  merely a stronger knob than a flag — it is the *only* correctly scoped one.

Do the `vbar` source fix first (`compiler-flags.md`, ~10 lines, expected golden 0);
it is independent of this plumbing and the sweep never validly tested it.

**Status (2026-07-25, superseded): DEMOTED TO FALLBACK.** The active route is
compiler flags — see `compiler-flags.md` — because a flag-only winner needs no
CarpetX change, no CI-build workaround, and should pass golden at exactly 0. Come
back here only if those probes all leave `uct1` at occ-1.

Branch: `opt/flux-launch-bounds` (at `defd9f60`, the Idea-1 low-scratch/occ-1
state; `4f902ca1` on top adds the inert probe guard). User has a CarpetX fork at
`../CarpetX` for the trial. Rationale + why this beats the serialization:
`CHECKPOINT.md` §KEY LESSONS 3 (the register↔scratch seesaw).

**PR hygiene:** the upstream-bound work is `b23fb947`..`defd9f60` only; `4f902ca1`
(probe guard), `25f08a66` (Planning docs) and `2fd4595e` (`make.code.deps`, which
hardwires `-Rpass` + `-DASTERX_PROBE_IDEALGAS_ONLY` and makes the binary abort on
pplim/hybrid/tabulated) sit on top and must be branched off or reverted before
submitting — and must never reach a golden or test-suite run.

### ⚠ THREE CORRECTIONS, verified against the sources 2026-07-25

1. **The pseudocode below is wrong about the lambda form.** CarpetX passes an
   `(i,j,k)` lambda to `ParallelFor` (`loop_device.hxx:83`), not the 1D `icell`
   form. AMReX's body (`AMReX_GpuLaunchFunctsG.H:1011-1030`) is 1D internally and
   bridges via `detail::call_f_intvect_handler`, so the replication **must call
   `amrex::detail::call_f_intvect_handler(KERNEL, iv, Gpu::Handler(...))`** rather
   than invoking the lambda directly. As written it will not compile.
2. **The AsterX-side change BREAKS THE CI BUILD.** Stock CarpetX's signature is
   `<int CI, int CJ, int CK, int VS, int N, int NT, typename F>`, so a 7th template
   argument binds `2` to `typename F` — a hard compile error, not a graceful
   fallback. The "Validate" section below claims CI can confirm the MB path
   compiles; it cannot. **Decide up front** how to keep CI green: keep `MB=2`
   uncommitted for the trial, upstream the CarpetX change first, or feature-detect.
3. **The "optional zero-code pre-check" does not exist.**
   `-mllvm --amdgpu-waves-per-eu=2` is rejected outright ("Unknown command line
   argument"); it is an IR function attribute, not a cl::opt. Full inventory of what
   *does* exist in this toolchain: `compiler-flags.md`.

Also note `AMREX_LAUNCH_KERNEL` has **no** min-blocks variant
(`AMReX_GpuLaunch.H:37-45`), only `MT` and `NOBOUND` — confirming the body must be
replicated rather than a macro reused. `launch_global<max_threads, min_blocks, L>`
itself is confirmed present in `AMReX_GpuLaunchGlobal.H`.

## IMPLEMENTATION STATUS (2026-07-27)

**Change 1 is COMMITTED AND PUSHED** — CarpetX branch
**`opt/loop-box-device-min-blocks`** @ `19267243`, off `dev` @ `55e7434e`, on
`origin` = `MChabanov/CarpetX`. `Loop/src/loop_device.hxx`, +85/−5, single file:

- `loop_box_device` gains `int MB = 0` as the 7th template parameter, before
  `typename F`.
- The device lambda is hoisted to a named `kernel` (was passed inline), then
  `if constexpr (MB == 0)` → the untouched `amrex::ParallelFor<NT>(box, kernel)`
  path; `else` → the replicated launch.
- The replication mirrors `ParallelFor<MT>(BoxND, L)` from
  `AMReX_GpuLaunchFunctsG.H` (~line 1008): `amrex::BoxIndexer indexer(box)`,
  `amrex::Gpu::makeNExecutionConfigs<NT>(box)`, the grid-stride
  `icell = NT*blockIdx.x + threadIdx.x + start_idx` guard, and
  **`amrex::detail::call_f_intvect_handler(kernel, iv, Gpu::Handler(...))`** —
  correction #1 above, applied.
- Launch: `amrex::launch_global<NT, MB><<<...>>>` under CUDA;
  `hipLaunchKernelGGL(HIP_KERNEL_NAME(amrex::launch_global<NT, MB>), ...)` under
  HIP, where `HIP_KERNEL_NAME` protects the comma from the macro's argument
  splitting. `AMREX_GPU_ERROR_CHECK()` after the loop.
- **SYCL/other backends fall back to `ParallelFor<NT>`** — there is no
  `launch_global` there. MB is a codegen hint, so dropping it costs performance,
  never correctness.
- `static_assert(MB >= 0)`; the `#ifndef AMREX_USE_GPU` CPU path and the
  `gpu_sync_after_every_kernel` tail are untouched.

**Caller audit (settles the blast radius).** Across every sibling thorn there is
exactly **one** `loop_box_device` call outside `loop_device.hxx` itself:
`AsterX/src/fluxes.cxx:1534`, the fused flux site, and it passes `<0, 0, 0>` —
three explicit arguments, everything else defaulted. The 11 internal callers
(`loop_all_device`, `loop_int_device`, `loop_mixpn_device`, …) all pass
`<CI, CJ, CK, VS, N, NT>` and deliberately do **not** forward MB. So a 7th
parameter with a default is safe everywhere, and the flux kernel is reachable
directly — no wrapper needs plumbing, contrary to what a `loop_mixpn_device`
call site would have required.

**NOT compile-tested — NOTHING has been compiled.** No ROCm and no Cactus build
on the dev machine: CarpetX's `agent_scripts/build.sh` aborts with
`missing required env var: CACTUSX` (and, note for next time, still exits 0 —
do not read a zero exit from it as a pass). Only a brace/paren balance check and
a read-through were done. **The Frontier build is the first compile of this code**,
so expect ordinary first-compile breakage there, most likely candidates:
`HIP_KERNEL_NAME` not being defined (fallback: the triple-chevron form, which
`CC -x hip` accepts), or `amrex::min` needing an explicit
`#include <AMReX_Algorithm.H>`.

**Change 2 (AsterX side) IS APPLIED — `fluxes.cxx` now passes `MB=2`.**
`grid.loop_box_device<0, 0, 0, /*VS*/ 1, /*N*/ 1, AMREX_GPU_MAX_THREADS, /*MB*/ 2>`
at the fused flux site, under a loud comment block. Deliberate, at the user's
instruction, for the Frontier trial:

- **CI WILL FAIL on every job** until the CarpetX branch is upstreamed or this
  line is reverted — against stock CarpetX the `2` binds to `typename F`, a hard
  compile error (correction #2). Accepted knowingly; this is not a regression to
  investigate.
- `AMREX_GPU_MAX_THREADS` is 256 on Frontier and 0 on CPU builds, and the
  `static_assert(NT > 0)` sits inside the GPU branch, so writing the default out
  explicitly reproduces the previous behaviour exactly on both.
- Revert to `grid.loop_box_device<0, 0, 0>(` before any golden or test-suite run,
  along with `4f902ca1` and `2fd4595e`.

**To build the trial on Frontier:** point the Cactus CarpetX checkout at
`opt/loop-box-device-min-blocks` (fetch + checkout in
`$CACTUS/repos/CarpetX`, or `origin` if the arrangement differs), rebuild, and
read `-Rpass` for `CalcFluxAll<uct=1,pplim=0,idealgas>`. **Success = occ 2 with
scratch ~3448** (not 4984). A useful intermediate step if the build misbehaves:
temporarily set `MB` back to 0 at the call site — that must reproduce
256/32/3448/occ-1 exactly, which isolates "the replication is wrong" from
"forcing occ-2 does not help".

## ⭐ TRIAL RESULT (2026-07-27) — occ-2 REACHED AT LOW SCRATCH. Criterion MET.

First compile of the plumbing, reduced build (`-DASTERX_PROBE_IDEALGAS_ONLY`),
`MB=2` at the flux site. **It built, and it did what the hypothesis said it would.**

| kernel | SGPR | VGPR | AGPR | VGPR+AGPR | scratch | spill v/s | occ |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLUX uct1 pp0 idealgas | 100 | **128** | **128** | **256** | **3608** | 0/0 | **2** |
| FLUX uct0 pp0 idealgas | 100 | 128 | 128 | 256 | 3496 | 0/0 | **2** |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 78 | 72 | 0/20,22,24 | 5 |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | — | 0 | 0/0 | 8 |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | — | 0 | 0/0 | 8 |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | — | 0 | 0/0 | 8 |

Against the reduced baseline (`uct1` = 100/256/32/3448/occ-1):

- **occ 1 → 2.** The target is met.
- **scratch 3448 → 3608, i.e. +160 B/lane.** The serialization's route to the same
  occupancy cost **+1536** (4984). So this reaches occ-2 for **~10% of the scratch
  the rolled loops cost**. That was the entire bet in Lesson 3(b), and it paid.
- **Total register pressure 288 → 256**, exactly the occ-2 budget (512/2). The
  allocator found the 32 registers that three compiler flags could not, because
  it was *told* to rather than asked.
- **The split is the interesting part: VGPR 256→128, AGPR 32→128, both spill
  counters still 0.** The allocator did not spill to memory; it used the AGPR half
  of the unified gfx90a register file as **on-chip** spill space. That is the best
  available outcome — `v_accvgpr_read/write` is register-speed, unlike
  `scratch_load/store`. It also retires the old framing of AGPR as "the overflow
  set to be drained to zero" (Lesson 2's gauge): at occ-2 the budget is 256 total
  and an even 128/128 split is a *healthy* allocation, not overflow.
- **Zero collateral: all 9 EMF kernels are byte-for-byte unchanged.** Their
  mangled names carry `Li0E` for MB and still route through
  `launch_global<256>` (the 2-arg overload, line 21) via `ParallelFor`, while the
  two flux kernels now use `launch_global<256, 2>` (the 3-arg overload, line 25).
  This is the structural advantage over a TU-scoped `-mllvm` flag, demonstrated
  rather than argued: probe 2 damaged `CalcE`, this does not touch it.

**Two caveats before this is called a win.**

1. **This is a reduced build.** The full build is a different compilation (AGPR 40
   vs 32 with no flags). `launch_bounds` forces the budget either way so occ-2
   should hold, but the scratch number may differ — re-measure.
2. **`uct0` was ALREADY occ-2 in the reduced build** (252/2/3448) and is now
   128/128/3496. Same occupancy, +48 B/lane scratch, more AGPR traffic — so for
   flux-CT this change is plausibly **neutral-to-slightly-negative**. That matters
   for how the timing runs will read: **the large-grid TOV run is flux-CT
   (`use_uct=no`)**, i.e. the kernel that may have had nothing to gain, while the
   *small*-grid subcycling run is UCT, the kernel that just gained a wave. Expect
   the sign of the effect to be the OPPOSITE of the serialization's (which won on
   TOV-large and lost on UCT-small). Do not read a flat TOV result as failure.
   ⚠ `uct0`'s full-build baseline was never recorded — only `uct1`'s
   (256/40/3448/occ-1). Get it while re-measuring, or this comparison stays
   guesswork.

**Next: timing, which is now the only open question.** Both grid sizes, against
the Idea-1 baseline, per `baseline-timings.md`. Reference points: forced occ-2 via
a global `__launch_bounds__` measured **−15.9%** on the TOV run (with Z4c
collateral, which this scoped version avoids); genuine occ-2 via serialization
managed only −4.8% because of the scratch. This sits at 3608 B/lane, much nearer
the low-scratch end, so the honest expectation is somewhere between −5% and −16%
on whichever kernel actually gained.

## Goal

Force the fused flux kernel to occ-2 on the **unrolled, constant-index, low-scratch**
Idea-1 kernel (scratch 3448 B/lane), letting the allocator spill **selectively**
(cold values only, hot reconstructed vectors stay in registers) — instead of the
serialization's wholesale register→scratch demotion (scratch 4984). Target the
same occ-2 the global-`__launch_bounds__(MT,2)` test reached at −15.9%, but SCOPED
to the flux kernel (no Z4c collateral).

## Mechanism (verified in the code)

- AsterX flux site: `grid.loop_box_device<0,0,0>(bnd_min,bnd_max,fmin,fmax,lambda)`.
- CarpetX `loop_box_device<CI,CJ,CK,VS,N,NT,F>` (`../CarpetX/Loop/src/loop_device.hxx`,
  ~line 42) → `amrex::ParallelFor<NT>(box, lambda)` (~line 84).
- AMReX `ParallelFor<MT>(box,f)` → `AMREX_LAUNCH_KERNEL(MT, nblocks, MT, …)` →
  **`launch_global<MT>`** = `__launch_bounds__(MT)` (max threads only, NO min-blocks).
- AMReX ALREADY has the knob: `AMReX_GpuLaunchGlobal.H` defines a 3-arg
  `template<int max_threads,int min_blocks,class L> __launch_bounds__(max_threads,
  min_blocks) launch_global(L)`. It's just not exposed through `ParallelFor`.
- `MB=2` is what the earlier global `__launch_bounds__(MT,2)` used to reach occ-2.

## Changes (3, backward-compatible)

**1. CarpetX `Loop/src/loop_device.hxx` — add `int MB = 0` to `loop_box_device`:**
```cpp
template <int CI, int CJ, int CK, int VS = 1, int N = 1,
          int NT = AMREX_GPU_MAX_THREADS, int MB = 0, typename F>   // +MB, default 0
void loop_box_device(...) const {
  ... build amrex::Box box ...
  if constexpr (MB == 0) {
    amrex::ParallelFor<NT>(box, KERNEL);       // unchanged path
  } else {
    // mirror amrex::ParallelFor<NT>(box,·) (CUDA/HIP body, AMReX_GpuLaunchFunctsG.H
    // ~lines 1008-1030) but launch via launch_global<NT,MB>:
    //   const auto& nec = amrex::Gpu::makeNExecutionConfigs<NT>(box);
    //   for (auto const& ec : nec) {
    //     const auto start_idx = std::uint64_t(ec.start_idx);
    //     const amrex::BoxIndexerND<3> indexer(box);
    //     auto KERNEL2 = [=] AMREX_GPU_DEVICE(){ grid-stride: icell = NT*blockIdx.x
    //        + threadIdx.x + start_idx; if (icell<indexer.numPts()) f(point_desc(...)); };
    // #if defined(AMREX_USE_CUDA)
    //     amrex::launch_global<NT,MB><<<ec.nblocks, NT, 0, amrex::Gpu::gpuStream()>>>(KERNEL2);
    // #elif defined(AMREX_USE_HIP)
    //     hipLaunchKernelGGL(HIP_KERNEL_NAME(amrex::launch_global<NT,MB,decltype(KERNEL2)>),
    //                        ec.nblocks, NT, 0, amrex::Gpu::gpuStream(), KERNEL2);
    // #endif
    //   }
    //   AMREX_GPU_ERROR_CHECK();
  }
  ... existing gpu_sync_after_every_kernel tail ...
}
```
`MB=0` default → every wrapper caller (`loop_int_device`/`loop_all_device`/… which
pass ≤ `NT`) is unchanged. KERNEL is the same wrapped device lambda the current
code builds (the `point_desc` grid-stride body). Keep the `#ifndef AMREX_USE_GPU`
CPU path unchanged (MB irrelevant there).

**2. AsterX `fluxes.cxx` (on `opt/flux-launch-bounds`) — pass MB=2 at the one flux site:**
```cpp
grid.loop_box_device<0, 0, 0, /*VS*/1, /*N*/1, AMREX_GPU_MAX_THREADS, /*MB*/2>(
    bnd_min, bnd_max, fmin, fmax, [=] CCTK_DEVICE(const PointDesc &p){ ... });
```
Applies `__launch_bounds__(256, 2)` to ONLY the fused flux kernel.

**3. No AMReX edit** — reuse AMReX's existing `launch_global<NT,MB>`; the ~15 new
lines live in the CarpetX fork.

## Validate

- Golden gate (`[golden-master]`) — must stay 0 (launch_bounds is a codegen
  directive; results unchanged). Note: golden CI uses stock CarpetX, not the fork,
  so the CarpetX change must be exercised via the user's Frontier build; the CI
  golden confirms the AsterX-side call compiles/runs identically when MB path is
  present.
- Frontier `-Rpass` on `<uct=1,pplim=0,idealgas>`: expect **occ 2** with scratch
  ~**3448** (NOT 4984). If scratch stays low and occ=2 → the hypothesis holds.
- Timing at BOTH grid sizes (UCT-small subcycling + TOV-large flux-CT) vs the
  Idea-1 baseline. Success = flux-kernel speedup (toward −15.9%) AND no small-grid
  regression. Record in `baseline-timings.md`.

## Risks / open decisions (raise with user)

- **Replicate vs upstream:** the MB!=0 branch duplicates ~12 lines of AMReX's
  `ParallelFor<MT>(box)` internals (`makeNExecutionConfigs`, grid-stride, error
  check). Cleaner long-term = an `amrex::ParallelFor<MT,MB>` overload upstreamed,
  but that touches AMReX (pinned release, not forked). For the trial: CarpetX-replicate.
- **MB value:** 2 (occ-2 target; matches the −15.9% test). Confirm occ via -Rpass.
- **Scoping:** guaranteed — only the flux site passes MB≠0.
- ~~**Optional zero-code pre-check (Frontier, cce):** `fluxes.cxx.o: CXXFLAGS +=
  -mllvm --amdgpu-waves-per-eu=2`~~ — **DEAD, 2026-07-25.** That cl::opt does not
  exist in this toolchain (hard error: "Unknown command line argument"); it is an
  IR function attribute only. Neither do `--amdgpu-spill-vgpr-to-agpr` or
  `--amdgpu-num-vgpr`. The TU-scoped-flag idea itself is sound and is now its own
  checkpoint — **`compiler-flags.md`** — using the options that DO exist
  (`--misched=gcn-max-occupancy`, `--enable-deferred-spilling`,
  `--inline-threshold=N`). Run those before building any of the plumbing below.

## Fallbacks

If launch_bounds does NOT win (scratch-from-selective-spill still eats it, or the
small grid still regresses): the serialization branch (`opt/flux-eig-collapse`) is
the large-grid-only win; or ship Idea-1 alone (occ-1, low scratch, the -27% fusion
+ Idea-1 removal still stand) and stop chasing occ-2.
