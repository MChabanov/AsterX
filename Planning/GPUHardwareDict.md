# GPU hardware dictionary — MI250X / gfx90a (CDNA2)

Reference for the terms in the `-Rpass-analysis=kernel-resource-usage` output
and the flux-kernel profiling (`profiling-flux-kernel.md`). All numbers are for
AMD MI250X (Frontier), architecture gfx90a / CDNA2.

## Hardware hierarchy

```
GPU (MI250X = 2 GCDs) → each GCD has ~110 Compute Units (CU)
   CU → 4 SIMD units
      SIMD → runs up to 8 wavefronts at once   ← "occupancy"
         wavefront (wave) → 64 lanes (work-items) in lockstep
```

## Execution units

- **wave / wavefront** — a bundle of **64 work-items (lanes)** executing the
  same instruction in lockstep. AMD is 64-wide (NVIDIA "warps" are 32-wide).
  In the flux loop, one lane ≈ one grid cell.
- **SIMD** — the physical execution unit that runs waves. Each CU has 4 of
  them; each SIMD owns the register files below and keeps several waves
  resident at once.
- **waves/SIMD (occupancy, "occ")** — **how many wavefronts are resident on a
  SIMD simultaneously**, 1–8 on gfx90a. The single most important number. GPUs
  hide memory latency by *oversubscription*: when one wave stalls on an HBM
  load, the scheduler instantly runs another resident wave. **Occ 8 = up to 8
  waves hiding behind each other; occ 1 = when the wave stalls, the SIMD
  idles.** Occupancy is set by whichever resource runs out first — VGPRs,
  SGPRs, LDS, or the hardware cap of 8.

## Registers (the scarce resources that cap occupancy)

Two axes: **per-lane vs uniform**, and **on-chip vs off-chip**.

| Term | What it holds | Key fact |
|---|---|---|
| **VGPR** (Vector GPR) | **Per-lane** values — each of the 64 lanes has its own copy (per-cell prims, fluxes, coords) | The scarce one. A SIMD has a **512-VGPR pool** shared across resident waves, so `waves ≈ 512 ÷ VGPRs-per-wave` (capped at 8). Max **256** addressable by one wave. |
| **SGPR** (Scalar GPR) | **Uniform** values — one copy shared by all 64 lanes (loop bounds, base pointers, grid metadata, constants) | Cheap and plentiful (~800/SIMD); rarely the limiter. A small "SGPR Spill" is minor. |
| **AGPR** (Accumulation VGPR) | A **second** vector register file, originally for matrix (MFMA) accumulation | On gfx90a the compiler repurposes it as **on-chip overflow** when 256 arch-VGPRs aren't enough. arch+acc *both* draw down the 512 pool. **Read it in context** (see below): with VGPR pinned at 256 it is an overflow signal, but under `__launch_bounds__` an even split such as 128/128 is a *healthy* allocation — on-chip spill space in place of off-chip scratch. |
| **Scratch (B/lane)** | Per-lane private memory living in **off-chip HBM**, not registers | ⚠ **Not necessarily spill** — see the `alloca` note below. When it *is* spill, it is the expensive kind: ~4000 B/lane means every work-item shuttling 4 KB to/from HBM mid-kernel. |

## Worked examples (from the flux-kernel profile)

- **Healthy** (`eos_3p_hybrid`): VGPR 54 → 512÷54 ≈ 9 → capped near 8 → **occ 7**;
  no AGPR; scratch 72 B (ABI baseline, ignore).
- **Register/scratch-bound** (`eos_3p_tabulated3d`): VGPR 256 (maxed) + AGPR 243
  ≈ 499 of the 512 pool → 512÷499 ≈ 1 → **occ 1**, and it *still* overflowed
  ~4 KB to scratch. Worst case: too few waves to hide latency **and** constant
  HBM scratch traffic with nothing to hide it behind — two penalties stacked.

## Notes

- **"Reduce register pressure"** = shrink the kernel's live-value footprint
  until VGPR drops enough that another wave fits **and** scratch spilling stops
  — turning an occ-1 double-penalty into latency-hidden throughput.
- **Granularity:** registers are allocated in fixed granules (~multiples of 8
  VGPR), so usage rounds up and occupancy comes back in steps, not per-register.
- **⚠ `ScratchSize` and the spill counters are SEPARATE fields, and scratch is
  often not spill at all.** The flux kernel reads 3448 B/lane with *both*
  `VGPRs Spill: 0` and `SGPRs Spill: 0` — those are private-memory `alloca`s
  (runtime-indexed stencil arrays that defeat SROA), not spilled registers.
  Contrast `CalcE`, at 72 B/lane *with* 20–24 SGPR spills, where the scratch
  really is spill. **Scratch does not cap occupancy either way** — occupancy is
  VGPR+AGPR (and LDS) only, so "reach occ-2" and "cut scratch" are independent
  goals. Measured consequence: two builds at the same occ-2 differing only in
  scratch (3608 vs 4984 B/lane) timed −27 % vs −5 %.
- **LDS (Local Data Share)** — on-chip scratchpad shared by a workgroup (not
  used by these kernels, LDS Size 0); can also be an occupancy limiter when a
  kernel allocates a lot of it.

## occ-2 arithmetic (gfx90a, the flux kernel)

occ ≈ `floor(512 / (VGPR + AGPR))` with **granule rounding** (VGPR granule 8), cap
8. occ-2 needs granule-rounded `VGPR + AGPR ≤ 256`.

**Two regimes, and the gauge differs between them:**

- *Unconstrained.* VGPR welds at the 256 ceiling while demand exceeds it, so the
  reducible overflow lands in AGPR and the gauge is **AGPR → exactly 0**. With
  AGPR 0, VGPR 253 rounds to 256 → 512/256 = 2 waves. A single residual AGPR
  rounds up to a granule → total > 256 → back to occ-1 (which is why 254 VGPR /
  1 AGPR was still occ-1, one granule short).
- *Under `__launch_bounds__(256, 2)`.* The budget is imposed, so the allocator
  distributes 256 total however it likes — measured **128 VGPR / 128 AGPR with
  both spill counters 0**, i.e. it used the AGPR half as *on-chip* spill space
  rather than off-chip scratch. Here a large AGPR count is the good outcome, not
  overflow. **Do not read "AGPR → 0" as the gauge once occupancy is forced.**

## NVIDIA H100 / GH100 (sm_90) contrast — TACC Vista

Different shape; the AMD AGPR framing does **not** transfer.

| | MI250X gfx90a | H100 GH100 sm_90 |
|---|---|---|
| warp/wave | 64 lanes | 32 threads |
| register file | 512 vec/SIMD (256 arch VGPR + 256 AGPR) | **65536 32-bit/SM, one unified file** |
| per-thread cap | 256 VGPR | **255 regs/thread** |
| overflow | spills to **AGPR** then scratch | **no AGPR** → spills to local memory (L1/L2-cached, DRAM) |
| occupancy unit | waves/SIMD (1-8) | warps/SM (≤64) |
| gauge to watch | **AGPR → 0** | **registers/thread + spill bytes** |

H100 occupancy ≈ `min(64, 256 / ceil(regs_per_thread/8))` warps (256-reg warp
granule): 32 regs→64 warps(100%), 64→32(50%), 128→16(25%), 255→8(12.5%). Measure
with **`--resource-usage`** (≡ `-Xptxas -v`, nvcc) or `-Rpass-analysis=...` if
clang-CUDA → per-kernel `registers`, `stack frame`, `spill stores/loads`; achieved
occupancy via `ncu`. Force occupancy with `__launch_bounds__(maxT, minBlocks)` or
`-maxrregcount`. The live-set reductions (Ideas 1/2/3) should also cut regs/thread
+ spills on Hopper, but "AGPR→0 crosses a granule" is AMD-only — expect a smoother
regs/thread→warps curve.
