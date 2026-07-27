# Flux-kernel GPU occupancy — methodology and evidence

Frontier / MI250X (gfx90a, CDNA2), ROCm, cce/20. Concerns `CalcFluxAll<uct,pplim,EOS>`
in `AsterX/src/fluxes.cxx`. Terminology: `GPUHardwareDict.md`. Status, lessons and
the final result: `CHECKPOINT.md`. Timings: `baseline-timings.md`.

## Where this ended up

occ-2 is reached and timed (**−27…−29 % on the flux kernel**) via scoped
`launch_bounds` — `launch-bounds-plan.md`. The serialization stack also reached occ-2
but at +1536 B/lane of scratch, which cost most of the benefit. Two constraints from
this investigation still bound any future work:

- **tabulated3d is EOS-table-bound (AGPR ~172–244) → occ-2 is unreachable there by
  any flux change.** Occupancy work is ideal-gas-only.
- **occ-1 was pre-existing**: the pre-fusion code was also occ-1, so the fusion cost
  no occupancy.

## How to profile (reusable)

- **Tier 0 — static `-Rpass` (the workhorse).**
  `fluxes.cxx.o: CXXFLAGS += -Rpass-analysis=kernel-resource-usage` in
  `AsterX/src/make.code.deps`; rebuild and read
  `VGPRs`/`AGPRs`/`ScratchSize`/`Occupancy`/spills per kernel. Demangle with
  `c++filt`; `Lb1E`=true / `Lb0E`=false → `CalcFluxAll<uct,pplim,EOS>`. frontier.cfg
  DEBUG=no → production-path numbers. Full recipe, mangled-name patterns and two
  parsing gotchas: `compiler-flags.md`.
- **Compile-only source bisection** — `#if`-guard a region, rebuild, watch the
  remark. Localizes register pressure with no GPU run; this is how the EOS math and
  the PP block were exonerated.
- **Tier 1 `rocprof`** (`--stats --hip-trace`): launch count + per-launch time.
- **Tier 2 `omniperf`/`rocprof-compute`** (`-k <substr>`): roofline, achieved
  occupancy + limiter, cache/VALU/MemUnit — to separate occupancy- from
  bandwidth-bound.
- **ISA inspection** (not yet done, and the recommended next step before any further
  register guessing): `v_accvgpr_read/write_b32` shows what is truly AGPR-resident,
  `scratch_load/store` what is in private memory.
- **NVIDIA (H100/GH100):** no `-Rpass` AGPR equivalent; use `--resource-usage`
  (`-Xptxas -v`) → `registers/thread` + stack-frame/spill bytes, and `ncu` for
  achieved occupancy. `GPUHardwareDict.md` §H100.

## Evidence

Register/occupancy by EOS, fused `CalcFluxAll`, pre-serialization:

| EOS | VGPR | AGPR | scratch B/lane | occ |
|---|---:|---:|---:|---:|
| ideal-gas | 256 | 74 | 3976 | 1 |
| hybrid `<poly>`/`<piecewise_poly>` | 54 | 0 | 72 | 7–8 |
| tabulated3d | 256 | 243 | 4024 | 1 |

EMF/aux kernels (`CalcE_impl`, `CalcFstag`, `CalcAux`) are all occ 6–8, so the EMF
sector is not worth fusing for occupancy.

**Exonerations (source bisection).** Stubbing the ideal-gas EOS math
(`csnd`/`eps`/`press`/`kappa`) left the peak byte-identical → not EOS math. Excising
the whole PP block left peak VGPR unchanged at 256 and dropped only AGPR overflow
74→46 → PP feeds the spillover, not the 256 core. **The 256 core peak is the main MHD
flux-assembly live set** (many simultaneously-live `vec<vec<2>,3>` quantities), and per
`CHECKPOINT.md` Lesson 3 the reconstruction stencils dominate what is resident.

**Why idealgas (256) vs hybrid (54) on identical MHD source.** Register allocation is
a whole-function heuristic: ideal-gas's transparent trivial EOS lets the scheduler
stretch live ranges (higher peak), while hybrid's opaque `ColdEOS*` calls act as
scheduling barriers (shorter live ranges). `CalcFluxAtFace` is `ALWAYS_INLINE`, so the
~1000-line worker flat-inlines into one region.

**Not levers — ruled out by evidence, do not re-probe:** trimming EOS math
(exonerated); un-fusing (pre-fusion was also occ-1); the neighbourhood-reuse pass
(*adds* live values → backfires on an occ-1 kernel); pure reorders and hand
rematerialization (`CHECKPOINT.md` Lesson 1 — Idea 4 was bit-identical and golden-PASS
but occupancy-null, reverted).

## Step 1 result — CT-scheme split (`b23fb947` / `e8b0718c` / `4a7e40ec`)

Templatized `CalcFluxAll<uct,EOS>` with `if constexpr(uct)` around the UCT
face-speed/drift block and `if constexpr(!uct)` around the flux-CT B-fluxes,
dispatched via a generic lambda over `std::true_type`/`false_type`. The template
parameter must be named `uct`, not `use_uct` — `DECLARE_CCTK_PARAMETERS` in
`CalcFluxAll` shadows the latter. Layer 2 gated *storage* (`vbar_*`/`a_*` iff
use_uct, `fluxBs_*` iff !use_uct) via duplicated schedule blocks, an interface.ccl
B-flux split and a sync.cxx restriction. **nvcc gotcha:** an extended `__device__`
lambda cannot first-capture a variable inside `if constexpr`, so the fused-sweep
zeroing needs `static_cast<void>(fx)` up top.

**Both layers golden PASS at 0; occupancy NO WIN** — production idealgas stayed VGPR
256 / occ 1 (AGPR 74→68, scratch 3976→3528). The B-flux/UCT-face quantities are not
what pins the 256 peak. Layer 1 is a bit-identical config-gating refactor and the
enabler for layer 2's storage/bandwidth win (−6 GFs, ~halved face-write streams,
−11 % scratch).

## Serialization result (Ideas 1/2/3) — occ-2 achieved, at a price

Production `CalcFluxAll<uct=1,pplim=0,idealgas>`:

| stage (commit) | VGPR | AGPR | scratch | occ |
|---|---:|---:|---:|---:|
| eig-collapse (`549a28ca`) | 256 | ~64 | ~3960 | 1 |
| Idea 1, template `use_pplim` (`defd9f60`) | 256 | 40 | 3448 | 1 |
| Idea 2, flux-by-flux momentum (`e8613114`) | 254 | 1 | 3816 | 1 |
| **Idea 3, face-state serialize (`210a013c`)** | **253** | **0** | 4984 | **2** |

- **Idea 1** compiled out the PP `_ppl` live set for `use_pplim=no`, matching the
  historical PP-excision (74→46): AGPR 64→40, *and* scratch down — genuine removal.
- **Idea 2** rolled `for(j)` over the momentum component: AGPR 40→1, VGPR 256→254.
  One granule short of occ-2 (254 rounds to 256, residual 1 AGPR rounds up).
- **Idea 3** rolled `for(f)` over the two face sides with the momentum `for(j)`
  nested and `eigenvalues_oneside` per side: drained the last AGPR → occ-2.
  flux-CT+noPP idealgas also crossed (251/0). Design note:
  `ideas-2-3-serialization.md`.
- **But scratch rose 3448→4984** because the rolled loops make `f`/`j` runtime
  indices and GPU registers cannot be dynamically indexed — relocation, not removal
  (`CHECKPOINT.md` Lesson 3). Idea 3's +1168 B is dominated by the reconstruction
  arrays; the cut-set is a fraction.
- PP-on idealgas stays occ-1 (`<1,1>` 28 AGPR, `<0,1>` 29; not production).
  tabulated `<1,0>` 256/172/occ-1, down from 244 but still EOS-bound. hybrid 54/0/occ-8.

Every step golden PASS at exactly 0. Timing verdict (grid-size-dependent, −4.8 % large
/ +8.9 % small) and the comparison against scoped `launch_bounds`:
`baseline-timings.md`.

## Frontier practicalities

- Profile single-rank/single-GCD (`ROCR_VISIBLE_DEVICES=0`), short `cctk_itlast` but
  long enough that every FMR level's flux launch fires.
- `-ffp-contract=off` is a CI-cpu-only concern; use normal Frontier opt flags so the
  profile reflects production codegen. frontier.cfg: DEBUG=no, OPTIMISE=yes, -O3.
- The `launch_bounds` change is header-only on the CarpetX side (template
  `launch_global` instantiated in `fluxes.cxx`) → rebuilding AsterX is enough.
