# Compiler-flag route to occ-2 — CLOSED. Record of a spent route.

_Opened 2026-07-25, closed on measurement 2026-07-27: all three mechanisms probed,
none moved occupancy. Kept for the negative results (expensive to re-derive), the
`-Rpass` recipe (still the way to measure anything), and the `vbar` finding. The
forward-looking parts — the `PROBE_FLAGS` switchboard and the tiered `-mllvm`
candidate inventory — were **deleted, not archived**, so nobody re-walks a list
that measurement has priced. The occupancy result is in `launch-bounds-plan.md`._

## Why it was tried, and the one property it had

Source-level removal was closed (`CHECKPOINT.md` Lesson 4) and `launch_bounds`
needed a CarpetX fork change with an unsolved CI problem. Flags had a property
neither had: **scheduling and register allocation do not alter arithmetic**, so a
flag-only winner should pass golden at exactly 0 — no rebaseline, no physics
validation. That is why it went first.

## Baselines — two sets, NOT comparable

`CalcFluxAll<uct=1,pplim=0,idealgas>`, full build (16 instantiations):
**256 VGPR / 40 AGPR / 3448 B/lane / occ 1**, spills 0/0.

Reduced build (`-DASTERX_PROBE_IDEALGAS_ONLY`, 2 instantiations):

| kernel | SGPR | VGPR | AGPR | scratch | spill v/s | occ |
|---|---:|---:|---:|---:|---:|---:|
| FLUX uct1 pp0 idealgas | 100 | 256 | 32 | 3448 | 0/0 | **1** |
| FLUX uct0 pp0 idealgas | 100 | 252 | 2 | 3448 | 0/0 | 2 |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 72 | 0/20,22,24 | 5 |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | 0 | 0/0 | 8 |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | 0 | 0/0 | 8 |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | 0 | 0/0 | 8 |

**AGPR is 40 full / 32 reduced with no flags** — the instantiation cut changes
inlining of the shared ReconX/EOS callees, so this is a different compilation
(`CHECKPOINT.md` Lesson 6). Calibration: `uct1` needed **32 registers** to clear the
cliff — not noise distance. `uct0` at 254 is two registers under the cliff, i.e.
threshold luck, but a useful **canary**: a flag that helps `uct1` while pushing
`uct0` back to occ-1 is doing real harm.

## SWEEP RESULTS — three mechanisms, all ±2 registers

All rows reduced build. **In all three probes every one of the 9 EMF kernels was
byte-for-byte identical to the baseline except where noted**, so only the target
rows are tabulated:

| probe | mechanism | `uct1` SGPR/VGPR/AGPR/scratch/occ | verdict |
|---|---|---|---|
| 1 `--misched=gcn-max-occupancy` | scheduling | 100/256/**32**/3448/**1** | **NULL**, byte-identical everywhere |
| 2 `--enable-deferred-spilling` | allocation | 102/256/**30**/3448/**1** | moves 2, **rejected** |
| 3 `-ffinite-math-only` | relaxed FP | 98/256/**34**/3448/**1** | **negative** |

- **Probe 1 was probably a null BY CONSTRUCTION.** In LLVM's
  `AMDGPUTargetMachine.cpp`, `GCNTargetMachine::createMachineScheduler` already
  returns `createGCNMaxOccupancyMachineScheduler` — the same function the
  `"gcn-max-occupancy"` registry entry names — and the generic pass only overrides
  the target's choice when `-misched` is given. So the flag re-selected the strategy
  already running, which explains a *perfectly* identical result across kernels with
  nothing in common. It does not retire scheduling as a mechanism; it mostly re-ran
  the baseline. (Corollary: AMDGPU already runs an unclustered-high-register-pressure
  reduction stage, so the scheduler is *already* fighting this.)
- **Probe 2 rejected despite moving in the right direction.** 32 → 30 against a
  requirement of 32 → 0 is 6 % of the distance with occupancy unchanged, and it
  *costs*: +2 SGPR on both flux kernels and **`CalcE uct1` SGPR spills 20/22/24 →
  30/26/33**, a real traffic increase in a kernel that runs three times per RHS.
  Probe 2 also served as the **positive control**: counters moved, so `-mllvm` flags
  demonstrably reach the compiler and probe 1's null is a genuine codegen null.
- **Probe 3 came out worse and its premise was flawed.** It was meant to price the
  `vbar` finding, but folding `0.0 * x` to a constant needs `nnan` **and `nsz`**
  (without no-signed-zeros the result is `±0.0` depending on the sign of `x`, so it
  is not a compile-time constant). `-ffinite-math-only` supplies only the first, so
  the dead loads most likely stayed live and the +2 AGPR is unrelated perturbation.
  What it does establish: **there is no free lunch in relaxed FP for this kernel's
  register pressure.**

## VERDICT

| mechanism | best effect on the target (needs AGPR 32→0) |
|---|---|
| scheduling | 0 registers (already the default strategy) |
| allocation | −2 registers, plus +4…+10 SGPR spills in `CalcE uct1` |
| relaxed FP | +2 registers (worse), and never shippable anyway |

Nothing here covers a 32-register gap; the spread across three independent
mechanisms is ±2. Two structural findings, both now folded into `CHECKPOINT.md`
Lesson 5: **the 32 AGPRs were load-bearing live values, not allocator sloppiness**,
and **TU-scoped flags cannot be aimed at one kernel** (all 11 kernels share
`fluxes.cxx`; probe 2's `CalcE` damage proves the collateral is real). The
remaining candidate list was deleted rather than deferred: every entry was a
variation on one of the two mechanisms already measured at ±2. **If a future agent
wants occupancy, the instrument is `launch_bounds`, not another flag.**

## Ruled out as `-mllvm` options (re-deriving each costs a failed build)

These exist only as IR **function attributes**, unreachable without attaching an
attribute to the kernel — which *is* the CarpetX plumbing, since the kernel is
`amrex::launch_global<>` inside AMReX headers:

- `--amdgpu-waves-per-eu` — hard error, *"Unknown command line argument"*. **This
  was the original plan's "zero-code pre-check"; it does not exist.**
- `--amdgpu-spill-vgpr-to-agpr`, `--amdgpu-num-vgpr` / `--amdgpu-num-sgpr`
- `--vgpr-regalloc=pbqp` — only `basic`/`greedy`/`fast` offered for the split
  AMDGPU allocators
- `--amdgpu-schedule-relaxed-occupancy` — *relaxes* occupancy, wrong direction

## Fast math: why not just enable it

- **The main perf component is already on.** `-ffp-contract` is unset in
  `frontier.cfg`, so HIP device code gets clang's `fast` default and FMA
  contraction is ON. What's left buys much less and carries nearly all the risk.
- **`-ffinite-math-only` disables safety machinery.** The `CCTK_DEBUG` NaN dump uses
  `isnan()`, which may fold to `false`; it also undermines `if (det_m < 0) det_m = 0`
  in `eigenvalues.hxx` and the atmosphere/`useLO` edge comparisons. That is a
  category difference from "last-bits shift".
- **`-fassociative-math` hits exactly the two expressions measured as fragile:**
  `tau` (cancels by up to `|Q/tau| ~ 5e7`) and the eigenvalue discriminant
  `a1^2 − 4a2a0` (~6-digit cancellation between two `O(H^2)` terms).
  See `flux-construction.md` §10b / §11a.
- **And it is measurably worthless here:** probe 3 moved the target from AGPR 32 to
  **34**. The case against shipping fast-math is empirical as well as semantic.

## The `vbar` finding — implemented, register-null, still worth keeping

**Found** by comparing `uct1` (288 registers) with `uct0` (254). In the
`if constexpr (uct)` block with `pplim=false` (production), `theta_uct` is the
literal `1.0`, so the drift blend degenerates to
`1.0*vj_face + (1.0-1.0)*0.5*(gf_vels(...) + gf_vels(...))`. `(1.0-1.0)` folds to
`0.0`, but **`0.0 * x` cannot be folded away** without `nnan` *and* `nsz` (see probe
3 above) — so the compiler is *required* to keep those loads: 4 per direction,
**12 total, provably dead** in production. Fix (in `fluxes.cxx` since 2026-07-27):

```cpp
if constexpr (pplim) {
  const CCTK_REAL theta_uct = gf_theta(dir_i)(p.I);
  vbar_j(dir_i)(p.I) = theta_uct * vj_face + (1.0 - theta_uct) * 0.5 * (...);
  vbar_k(dir_i)(p.I) = theta_uct * vk_face + (1.0 - theta_uct) * 0.5 * (...);
} else {
  vbar_j(dir_i)(p.I) = vj_face;
  vbar_k(dir_i)(p.I) = vk_face;
}
```

**MEASURED: a register null** — all 11 kernels byte-for-byte identical
(`uct1` still 100/256/32/3448/occ-1). The prediction that these would count because
they are non-rematerializable *memory loads* was wrong; the reason is
`CHECKPOINT.md` **Lesson 4** (position relative to the pressure peak is what
matters, and the UCT epilogue is past it).

**Keep the change, with the claim restated:** bit-identical (`1.0*x == x`;
`x + 0.0 == x` except `x == -0.0` → `+0.0`, invisible to a norm-based compare, so
expect golden = 0), it deletes 12 real global loads (traffic, which `-Rpass` does
not measure), and it removes a multiply-by-zero that only looked meaningful. **No
occupancy claim.** Affects `uct1` only. Don't spend a dedicated timing run on it.

## Reading `-Rpass` output — the recipe (still current)

`AsterX/src/make.code.deps` (local, uncommitted as far as upstream is concerned)
carries the two surviving lines:

```make
fluxes.cxx.o: CXXFLAGS += -Rpass-analysis=kernel-resource-usage
fluxes.cxx.o: CXXFLAGS += -DASTERX_PROBE_IDEALGAS_ONLY
```

```bash
export CACTUS=/ccs/home/mchabanov/EinsteinToolkit/Cactus
export CONFIG=cray20-adios-register
cd $CACTUS
rm -f configs/$CONFIG/build/AsterX/fluxes.cxx.o
gmake $CONFIG 2>&1 | tee ~/probe.log
```

Remarks are emitted during compilation, before linking — safe to `Ctrl-C` once they
scroll past, and getting output at all proves the TU recompiled.

`-DASTERX_PROBE_IDEALGAS_ONLY` (guard `4f902ca1`) cuts 16 instantiations → 2:
idealgas, `pplim=0`, **both** CT schemes — both kept deliberately, since both timing
runs are idealgas/`use_pplim=no` and TOV-large is flux-CT, so the fast build serves
the timing runs too. It **aborts at runtime** on `use_pplim=yes`, hybrid or
tabulated: never run the test suite or golden against it.

Kernel identification (template order `<use_uct, use_pplim, EOSType>`):

```
CalcFluxAllILb1ELb0E...idealgas   -> the production kernel
CalcFluxAllILb0ELb0E...idealgas   -> flux-CT
CalcE_implILi[0-2]ELb[01]E        -> EMF, per direction, per CT scheme
```

```bash
grep -A 12 'Function Name.*CalcFluxAllILb1ELb0E.*idealgas' ~/probe.log \
  | grep -E 'VGPRs:|AGPRs:|ScratchSize|Occupancy|Spill' \
  | sed 's/.*remark: *//; s/ \[-Rpass.*//'
```

Two parsing gotchas that cost real time: the last whitespace field is the
`[-Rpass-...]` tag, not the value (strip it first, then take what follows the LAST
`": "`, since `ScratchSize [bytes/lane]:` and `Occupancy [waves/SIMD]:` carry their
own brackets); and awk has no block scope, so a helper function's temporaries are
globals unless declared as extra parameters — using `s` inside one silently
clobbered the ScratchSize global and printed `scratch 0` for every kernel.
