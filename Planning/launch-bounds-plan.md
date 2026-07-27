# Scoped `launch_bounds` — occ-2 without the serialization scratch. ⭐ IMPLEMENTED, WORKS.

_The occupancy result. Reaches occ-2 on the unrolled low-scratch kernel by forcing
the target rather than coaxing the allocator (`CHECKPOINT.md` Lesson 5), and times
−27…−29 % on the flux kernel with zero collateral. Timings live in
`baseline-timings.md`; this doc owns the mechanism, the implementation and the
outstanding controls._

## Goal and mechanism

Force occ-2 on the **unrolled, constant-index, low-scratch** Idea-1 kernel (scratch
3448 B/lane), letting the allocator spill **selectively**, instead of the
serialization's wholesale register→scratch demotion (4984). Scoped to the flux
kernel, unlike the earlier *global* `__launch_bounds__` probe that bought −15.9 % on
TOV at the cost of regressing `Z4c_RHS` +48 %.

The knob already existed in AMReX and was simply not reachable:

- AsterX flux site → CarpetX `loop_box_device<CI,CJ,CK,VS,N,NT,F>` →
  `amrex::ParallelFor<NT>(box, lambda)` → `AMREX_LAUNCH_KERNEL(MT,…)` →
  **`launch_global<MT>`** = `__launch_bounds__(MT)`, max threads only.
- `AMREX_LAUNCH_KERNEL` has **no** min-blocks variant (`AMReX_GpuLaunch.H:37-45`),
  only `MT` and `NOBOUND` — so the body must be replicated, not a macro reused.
- But `AMReX_GpuLaunchGlobal.H` already defines the 3-arg
  `template<int max_threads,int min_blocks,class L> __launch_bounds__(max_threads,
  min_blocks) launch_global(L)`. It is just not exposed through `ParallelFor`.

## Implementation (both changes done)

**CarpetX** — branch `opt/loop-box-device-min-blocks` @ `19267243`, off `dev` @
`55e7434e`, pushed to `MChabanov/CarpetX`. One file, `Loop/src/loop_device.hxx`,
+85/−5:

- `loop_box_device` gains `int MB = 0` as the 7th template parameter, before
  `typename F`.
- The device lambda is hoisted to a named `kernel` (was passed inline), then
  `if constexpr (MB == 0)` → the untouched `amrex::ParallelFor<NT>(box, kernel)`
  path; `else` → a replicated launch.
- The replication mirrors `ParallelFor<MT>(BoxND, L)` from
  `AMReX_GpuLaunchFunctsG.H` (~line 1008): `amrex::BoxIndexer indexer(box)`,
  `amrex::Gpu::makeNExecutionConfigs<NT>(box)`, the grid-stride
  `icell = NT*blockIdx.x + threadIdx.x + start_idx` guard, and
  **`amrex::detail::call_f_intvect_handler(kernel, iv, Gpu::Handler(...))`** —
  necessary because CarpetX passes an `(i,j,k)` lambda while AMReX's body is 1D
  internally; calling the lambda directly does not compile.
- Launch: `amrex::launch_global<NT, MB><<<...>>>` under CUDA;
  `hipLaunchKernelGGL(HIP_KERNEL_NAME(amrex::launch_global<NT, MB>), ...)` under
  HIP, where `HIP_KERNEL_NAME` protects the comma from the macro's argument
  splitting. `AMREX_GPU_ERROR_CHECK()` after the loop.
- **SYCL and other backends fall back to `ParallelFor<NT>`** — no `launch_global`
  exists there, and MB is a codegen hint, so dropping it costs performance, never
  correctness.
- `static_assert(MB >= 0)`; the `#ifndef AMREX_USE_GPU` CPU path and the
  `gpu_sync_after_every_kernel` tail untouched.

**AsterX** — `c96df5d5` on `opt/flux-launch-bounds`:

```cpp
grid.loop_box_device<0, 0, 0, /*VS*/ 1, /*N*/ 1, AMREX_GPU_MAX_THREADS, /*MB*/ 2>(
    bnd_min, bnd_max, fmin, fmax, [=] CCTK_DEVICE(const PointDesc &p){ ... });
```

- **⚠ This breaks CI on every job** until the CarpetX branch is upstreamed or the
  line is reverted: stock CarpetX's signature is
  `<int CI, int CJ, int CK, int VS, int N, int NT, typename F>`, so the `2` binds to
  `typename F` — a hard compile error, not a graceful fallback. Accepted knowingly
  for the trial; **decide the CI story before the PR** (upstream CarpetX first, or
  feature-detect).
- `AMREX_GPU_MAX_THREADS` is 256 on Frontier and 0 on CPU builds, and the
  `static_assert(NT > 0)` sits inside the GPU branch, so spelling the default out
  explicitly reproduces the previous behaviour exactly on both.

**Caller audit (settles the blast radius).** Across every sibling thorn there is
exactly **one** `loop_box_device` call outside `loop_device.hxx` itself — the fused
flux site, passing `<0,0,0>`. The 11 internal wrappers (`loop_all_device`,
`loop_int_device`, `loop_mixpn_device`, …) pass `<CI,CJ,CK,VS,N,NT>` and deliberately
do not forward MB. So a 7th parameter with a default is safe everywhere, and the flux
kernel is reachable **directly** — had it gone through `loop_mixpn_device`, MB would
have needed plumbing through that wrapper too.

**Never compile-tested off Frontier.** No ROCm or Cactus build on the dev machine
(CarpetX's `agent_scripts/build.sh` aborts with `missing required env var: CACTUSX`
— and still exits 0, so don't read a zero exit from it as a pass). The Frontier
build was the first compile and it succeeded.

## ⭐ TRIAL RESULT — occ-2 at low scratch. Criterion MET.

Reduced build (`-DASTERX_PROBE_IDEALGAS_ONLY`), `MB=2` at the flux site:

| kernel | SGPR | VGPR | AGPR | VGPR+AGPR | scratch | spill v/s | occ |
|---|---:|---:|---:|---:|---:|---:|---:|
| FLUX uct1 pp0 idealgas | 100 | **128** | **128** | **256** | **3608** | 0/0 | **2** |
| FLUX uct0 pp0 idealgas | 100 | 128 | 128 | 256 | 3496 | 0/0 | **2** |
| all 9 EMF kernels | — | — | — | — | — | — | **unchanged** |

Against the reduced baseline (`uct1` = 100/256/32/3448/occ-1):

- **occ 1 → 2**, and **scratch 3448 → 3608, i.e. +160 B/lane** where the
  serialization cost **+1536** for the same occupancy: ~10 % of the price. That was
  the entire bet, and it paid — see the head-to-head in `baseline-timings.md`.
- **Total pressure 288 → 256**, exactly the occ-2 budget. The allocator found the 32
  registers three compiler flags could not, because it was *told* to rather than
  asked (`CHECKPOINT.md` Lesson 5).
- **The 128/128 split with both spill counters at 0** means the allocator used the
  AGPR half of the unified gfx90a file as **on-chip** spill space, not scratch
  memory — the best available outcome, and it retires the old "drain AGPR to 0"
  gauge (Lesson 2).
- **Zero collateral, demonstrated twice.** The 9 EMF kernels are byte-for-byte
  unchanged — their mangled names carry `Li0E` for MB and still route through
  `launch_global<256>` via `ParallelFor`, while only the two flux kernels use
  `launch_global<256, 2>` — and `Z4c_RHS` is flat in the TOV timing. Contrast the
  TU-scoped flag, which damaged `CalcE`.

**Outstanding controls** (also in `CHECKPOINT.md` §STATUS): the `MB=0` control in the
same reduced build; a full-build re-measurement, since the full build is a different
compilation (AGPR 40 vs 32) and `uct0`'s full-build baseline was never taken; then
golden, expected exactly 0. A useful diagnostic if a future build misbehaves: set
`MB` back to 0 at the call site — that must reproduce 256/32/3448/occ-1 exactly,
which isolates "the replication is wrong" from "forcing occ-2 does not help".

**One prediction from this doc was wrong, and it is worth remembering why.** It
predicted TOV-large would be neutral because `uct0` "was already occ-2" — taking
`uct0`'s occupancy from the *reduced* build (252/2, threshold luck at two registers
under the cliff) and applying it to a full-build baseline where it had never been
measured. TOV-large measured −27.3 %. `CHECKPOINT.md` Lesson 6.

## Open design question for upstreaming

The `MB != 0` branch duplicates ~12 lines of AMReX's `ParallelFor<MT>(box)`
internals (`makeNExecutionConfigs`, grid-stride, error check), so it must be kept in
sync if those internals change. The cleaner long-term shape is an
`amrex::ParallelFor<MT,MB>` overload upstreamed into AMReX — but that touches a
pinned release rather than a fork. CarpetX-side replication was the right call for
the trial; for upstream, raise both options with the maintainers.
