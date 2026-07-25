# Scoped launch_bounds plan — reach occ-2 WITHOUT the serialization scratch

**Status (2026-07-25): DEMOTED TO FALLBACK.** The active route is compiler flags —
see `compiler-flags.md` — because a flag-only winner needs no CarpetX change, no
CI-build workaround, and should pass golden at exactly 0. Come back here only if
those probes all leave `uct1` at occ-1.

Branch: `opt/flux-launch-bounds` (at `defd9f60`, the Idea-1 low-scratch/occ-1
state; `4f902ca1` on top adds the inert probe guard). User has a CarpetX fork at
`../CarpetX` for the trial. Rationale + why this beats the serialization:
`CHECKPOINT.md` §KEY LESSONS 3 (the register↔scratch seesaw).

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
