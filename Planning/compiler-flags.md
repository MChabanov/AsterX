# Compiler-flag route to occ-2 — CLOSED, kept as the record

_Opened 2026-07-25, **closed on measurement 2026-07-27**. Sibling of
`CHECKPOINT.md`; read that first for the overall mission. This doc owned the
attempt to reach occ-2 on the flux kernel via **compiler flags** rather than
source changes; all three mechanisms were measured and none moved occupancy
(§SWEEP RESULTS, §VERDICT). The live occupancy route is now
`launch-bounds-plan.md`._

_The forward-looking parts of this doc — the `PROBE_FLAGS` switchboard in
`make.code.deps` and the tiered `-mllvm` candidate inventory — were **deleted**,
not archived, so nobody re-walks a list that measurement has priced. What remains
is evidence: the two baselines, the three probe tables, the verdict, the
`-Rpass` reading recipe (still needed for any future measurement), the `vbar`
finding, and the fast-math analysis._

## Why this route exists

`CHECKPOINT.md` Lesson 3 named two ways to occ-2 with low scratch: (a) genuine
removal, (b) `launch_bounds`. Route (a) closed on measurement (Lesson 8: the
`(H,vf2)` fusion was a byte-for-byte `-Rpass` null, because nothing downstream of
the flux-assembly waist is in the AGPR overflow set). Route (b) needs a CarpetX
fork change and has an unsolved CI-build problem. **Compiler flags are a third
route**, and they have a property neither of the others has:

> **Scheduling and register allocation do not alter arithmetic.** With no
> `-ffast-math` anywhere (confirmed, see below), a flag-only winner should pass
> the golden gate at **exactly 0** — no rebaseline, no physics validation, no
> maintainer sign-off on changed numbers. That is a materially better outcome
> than either source-change route.

## The two targets, and they are INDEPENDENT

Per the Lesson-3 correction: AMD occupancy is set by VGPR+AGPR (and LDS), **not**
by scratch.

- **Occupancy**: clear the 256-register cliff. `occ = floor(512/(VGPR+AGPR))`,
  VGPR granule 8. VGPR is welded at the 256 ceiling, so the gauge is **AGPR → 0**.
- **Scratch**: a latency/bandwidth cost only. It is what made the small grid
  regress +8.9% under Ideas 2/3 (4984 B/lane). "Beat 4984" is a **timing**
  criterion, not an occupancy one. Do not conflate them.

Baseline scratch of 3448 B/lane is **`alloca`, not spill** (spill counters both
read 0) — ~431 doubles/lane of private memory, almost certainly ReconX stencil
arrays whose runtime-indexed loops defeat SROA.

## Baselines — TWO sets, NOT comparable

**Full build** (all 16 instantiations), `CalcFluxAll<uct=1,pplim=0,idealgas>`:

| VGPR | AGPR | scratch | spill v/s | occ |
|---:|---:|---:|---:|---:|
| 256 | 40 | 3448 | 0/0 | 1 |

**Reduced build** (`-DASTERX_PROBE_IDEALGAS_ONLY`, 2 instantiations), 2026-07-25:

| kernel | SGPR | VGPR | AGPR | scratch | spill v/s | occ |
|---|---:|---:|---:|---:|---:|---:|
| FLUX uct1 pp0 idealgas | 100 | 256 | 32 | 3448 | 0/0 | **1** |
| FLUX uct0 pp0 idealgas | 100 | 252 | 2 | 3448 | 0/0 | **2** |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 72 | 0/20,22,24 | 5 |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | 0 | 0/0 | 8 |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | 0 | 0/0 | 8 |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | 0 | 0/0 | 8 |

**⚠ AGPR is 40 in the full build and 32 in the reduced one, with no flags.**
Cutting the instantiation matrix changes inlining of the shared ReconX/EOS
callees (LLVM scores call sites partly on user count), so the reduced build is a
*different compilation*. Sweep rows are comparable to each other **only**; any
winner must be re-measured in a full build before it is believed.

Two calibrations that follow:

- **`uct1` needs 32 registers to clear the cliff.** That is not noise distance.
- **`uct0` was occ-1 before the cut and is now occ-2 at 254** — two registers
  under the cliff, i.e. threshold luck, not a mechanism. It is however a useful
  **canary**: a flag that helps `uct1` while pushing `uct0` back to occ-1 is
  doing real harm.

## Mechanism — how to measure (kept: still the recipe for any `-Rpass` run)

`AsterX/src/make.code.deps` (LOCAL, uncommitted) carries the two surviving lines:

```make
fluxes.cxx.o: CXXFLAGS += -Rpass-analysis=kernel-resource-usage
fluxes.cxx.o: CXXFLAGS += -DASTERX_PROBE_IDEALGAS_ONLY
```

```bash
export CACTUS=/ccs/home/mchabanov/EinsteinToolkit/Cactus
export REPO=$CACTUS/repos/AsterX
export CONFIG=cray20-adios-register

cd $CACTUS
rm -f configs/$CONFIG/build/AsterX/fluxes.cxx.o
gmake $CONFIG 2>&1 | tee ~/probe.log
```

Remarks are emitted during compilation, before linking — safe to `Ctrl-C` once
they scroll past.

**The `PROBE_FLAGS` switchboard was removed 2026-07-27** along with the `-mllvm`
candidate inventory, when the route closed (§VERDICT). Do not re-add it.

`-DASTERX_PROBE_IDEALGAS_ONLY` (guard committed as `4f902ca1`, "Only temporary
trick for testing") cuts 16 instantiations → 2: idealgas, `pplim=0`, **both** CT
schemes. Both are kept deliberately, because both timing runs are idealgas with
`use_pplim=no` and the large-grid TOV run is flux-CT — so the fast build is
usable for the timing runs too. It **aborts at runtime** on `use_pplim=yes`,
hybrid, or tabulated: never run the test suite or a golden check against it.

## Reading `-Rpass` output — two gotchas that cost real time

Remark lines look like:

```
/path:21:1: remark:     ScratchSize [bytes/lane]: 3448 [-Rpass-analysis=kernel-resource-usage]
```

1. **The last whitespace field is the `[-Rpass-...]` tag, not the value.** Strip
   the trailing bracketed tag first, then take what follows the LAST `": "` (two
   lines carry their own brackets: `ScratchSize [bytes/lane]:` and
   `Occupancy [waves/SIMD]:`).
2. **awk has no block scope.** A helper function's temporaries are globals unless
   declared as extra parameters. Using `s` inside such a helper silently clobbered
   the global holding ScratchSize and printed `scratch 0` for every kernel.

Kernel identification from the mangled name (template order
`<use_uct, use_pplim, EOSType>`, `Lb1E`=true, `Lb0E`=false):

```
CalcFluxAllILb1ELb0E...idealgas   -> the production kernel
CalcFluxAllILb0ELb0E...idealgas   -> flux-CT
CalcE_implILi[0-2]ELb[01]E        -> EMF, per direction, per CT scheme
```

Quick extraction without a helper script:

```bash
grep -A 12 'Function Name.*CalcFluxAllILb1ELb0E.*idealgas' ~/probe.log \
  | grep -E 'VGPRs:|AGPRs:|ScratchSize|Occupancy|Spill' \
  | sed 's/.*remark: *//; s/ \[-Rpass.*//'
```

## Flag inventory — DELETED 2026-07-27

The tiered candidate list (Tier 1 schedulers, Tier 2 allocator options, Tier 2b
inlining, Tier 3 scheduler weights) and the `llc --help-hidden` enumeration
recipe were removed when the route closed on measurement. They described work
that will not be done: every remaining entry was a variation on one of the two
mechanisms already measured at ±2 registers against a 32-register gap. See
§SWEEP RESULTS for what was run and §VERDICT for why the rest is not worth a
build.

One negative result from the inventory is kept, because re-deriving it costs a
failed build each time. These are IR **function attributes**, not `-mllvm`
cl::opts, so they are unreachable without attaching an attribute to the kernel
(which *is* the CarpetX plumbing, since the kernel is `amrex::launch_global<>`
inside AMReX headers):

- `--amdgpu-waves-per-eu` — hard error: *"Unknown command line argument"*,
  suggested `--amdgpu-dce-in-ra`. **This was the plan's "zero-code pre-check"; it
  does not exist.**
- `--amdgpu-spill-vgpr-to-agpr`
- `--amdgpu-num-vgpr` / `--amdgpu-num-sgpr`
- `--vgpr-regalloc=pbqp` — only `basic`/`greedy`/`fast` are offered for the split
  AMDGPU allocators.
- `--amdgpu-schedule-relaxed-occupancy` — *relaxes* occupancy targets, wrong
  direction.

**Diagnostics — never ship.**

- `-ffinite-math-only` / `-ffast-math` — bounds what dead-code elimination and
  reassociation could buy. See the fast-math section below.
- `--stress-regalloc=216` — hard-caps *all* regclasses; an LLVM testing option,
  expect poor code. Only answers "is occ-2 reachable at all".

## Fast math: why not just enable it

Asked and answered 2026-07-25. **Not in `frontier.cfg`**, but worth one TU-scoped
diagnostic build.

- **You already have the main perf component.** `-ffp-contract` is unset in
  `frontier.cfg`, so HIP device code gets clang's `fast` default: FMA contraction
  is ON. What's left (`-ffinite-math-only`, `-fassociative-math`,
  `-freciprocal-math`, `-fno-signed-zeros`) buys much less and carries nearly all
  the risk.
- **`-ffinite-math-only` disables safety machinery.** The `CCTK_DEBUG` NaN dump
  uses `isnan()`, which the compiler may then fold to `false` — NaN detection
  silently stops working. It also undermines `if (det_m < 0) det_m = 0` in
  `eigenvalues.hxx` and the atmosphere/`useLO` edge comparisons. This is a
  category difference from "last-bits shift".
- **`-fassociative-math` hits exactly the two expressions measured as fragile:**
  `tau` (cancels by up to `|Q/tau| ~ 5e7`) and the eigenvalue discriminant
  `a1^2 - 4*a2*a0` (a ~6-digit cancellation between two `O(H^2)` terms). See
  `flux-construction.md` §10b / §11a.
- **The one register item it would buy is obtainable exactly.** See the `vbar`
  finding below.
- **MEASURED 2026-07-27: the register upside is not merely small, it is negative.**
  `-ffinite-math-only` moved the target from AGPR 32 to **34** (probe 3). Whatever
  relaxed FP deletes here, the allocator does not turn into occupancy. The case
  against shipping fast-math is now empirical as well as semantic.

Also note for the record: `-ffp-contract` is unset on Frontier but explicitly
`off` in the CI CPU config, so **golden = 0 on CI implies nothing about
bit-identity of the Frontier binary**. `CHECKPOINT.md` Lesson 5's bit-identity
claim is CI-scoped only. Fine for the workflow (CI validates, Frontier measures)
as long as the two are never conflated.

## The `vbar` finding — a real removal, independent of all flags

**STATUS: IMPLEMENTED in `fluxes.cxx` 2026-07-27 (uncommitted at time of writing);
`-Rpass` and golden not yet run.** This survives the closure of the flag route —
it is a source change, and the flag sweep never validly tested it (probe 3 needed
`-fno-signed-zeros`, see above). Measure it directly.

Found by comparing `uct1` (288 registers) against `uct0` (254): **UCT is the
expensive config by ~34 registers**, which was not previously documented.

In the `if constexpr (uct)` block, with `pplim=false` (production) `theta_uct` is
the literal `1.0`, so:

```cpp
vbar_j(dir_i)(p.I) = 1.0 * vj_face
                   + (1.0 - 1.0) * 0.5 * (gf_vels(dir_j)(p.I) + gf_vels(dir_j)(p.I - p.DI[dir_i]));
```

`(1.0 - 1.0)` constant-folds to `0.0`, but **`0.0 * x` cannot be folded away**
without fast-math. Two conditions are needed and clang gives neither by default:
`nnan` (because `x` could be ±Inf and `Inf * 0 → NaN`) **and `nsz`** (because
`0.0 * x` is `-0.0` for `x < 0`, so the product is not a compile-time constant).
`-ffinite-math-only` supplies only the first — which is why probe 3 was not a
valid test of this. So the compiler is *required* to keep those `gf_vels` loads:
4 per direction, **12 total, provably dead** in the production config.

These are **memory loads that must stay live** across the UCT epilogue — exactly
the category Lesson 8 says the allocator banks, unlike the rematerializable
arithmetic the `(H,vf2)` fusion deleted. Fix:

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

Effectively bit-identical: `1.0*x == x` exactly, and `x + 0.0 == x` except for
`x == -0.0` (which yields `+0.0` before, `-0.0` after) — invisible to a norm-based
golden compare. **Expect golden = 0.** Affects `uct1` only; `uct0` does not
compile that block.

## SWEEP RESULTS

All rows are the **reduced build** (`-DASTERX_PROBE_IDEALGAS_ONLY`), so they are
comparable to each other and to the reduced baseline above — never to the full
build.

### Probe 1 — `-mllvm --misched=gcn-max-occupancy` — **NULL** (2026-07-27)

| kernel | SGPR | VGPR | AGPR | scratch | spill v/s | occ | vs baseline |
|---|---:|---:|---:|---:|---:|---:|---|
| FLUX uct1 pp0 idealgas | 100 | 256 | 32 | 3448 | 0/0 | **1** | identical |
| FLUX uct0 pp0 idealgas | 100 | 252 | 2 | 3448 | 0/0 | 2 | identical |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 72 | 0/20,22,24 | 5 | identical |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | 0 | 0/0 | 8 | identical |

**Byte-for-byte identical to the no-flag reduced baseline in every counter of
every one of the 11 kernels** — not a small move. Target unchanged at
`uct1` 256/32/3448/occ-1; canary `uct0` unmoved at occ-2. Recompilation is proven
(remarks are emitted at compile time, so getting output at all proves
`fluxes.cxx` rebuilt).

**⚠ Probably a null BY CONSTRUCTION — this flag likely re-selects the default.**
In LLVM's `AMDGPUTargetMachine.cpp`, `GCNTargetMachine::createMachineScheduler`
already returns `createGCNMaxOccupancyMachineScheduler`, and the same function is
what the `"gcn-max-occupancy"` `MachineSchedRegistry` entry names. The generic
`MachineScheduler` pass only overrides the target's choice when `-misched` is
given — so passing `-misched=gcn-max-occupancy` asks for the strategy that was
already running. That is consistent with a *perfectly* identical result across
kernels that have nothing in common, and it is a better explanation than "the
flag never arrived". Corollary: this probe carries **much less information than
the STATUS table assumed** — it does not retire the scheduling mechanism, it
mostly re-ran the baseline. Tier 3's premise ("the scheduler is already fighting
this") is now direct evidence rather than inference.

**Positive control — SATISFIED by probe 2 (below).** The worry was that
"default-already" and "flag-never-propagated" predict the same identical output,
and only one is benign. Probe 2 used the identical mechanism (in-file
`PROBE_FLAGS ?=` in `make.code.deps`) and **did change counters**, so `-mllvm`
arguments demonstrably reach the compiler. Probe 1's null is therefore a real
codegen null, and the default-scheduler explanation stands. No separate control
build was needed, and none is owed now that the route is closed.

### Probe 2 — `-mllvm --enable-deferred-spilling` — **MOVES, BUT LOSES** (2026-07-27)

| kernel | SGPR | VGPR | AGPR | scratch | spill v/s | occ | vs baseline |
|---|---:|---:|---:|---:|---:|---:|---|
| FLUX uct1 pp0 idealgas | 102 | 256 | **30** | 3448 | 0/0 | **1** | SGPR +2, **AGPR −2**, occ unchanged |
| FLUX uct0 pp0 idealgas | 102 | 251 | 2 | 3448 | 0/0 | 2 | SGPR +2, VGPR −1, occ unchanged |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 72 | 0/**30,26,33** | 5 | **SGPR spills +10/+4/+9 — WORSE** |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | 0 | 0/0 | 8 | identical |

**Verdict: reject.** The direction is right — this is the first flag to touch the
target's AGPR at all — but **32 → 30 against a requirement of 32 → 0** is 6% of
the distance, and occupancy is unchanged on both flux kernels. Meanwhile it
*costs*: +2 SGPR on both flux kernels and **+4 to +10 SGPR spills in `CalcE uct1`**,
which is a real traffic increase in a kernel that runs three times per RHS at
occ-5. Nothing to gain, something to lose.

Two things worth more than the verdict:

- **The allocator-side lever is alive but weak.** Unlike the scheduler (probe 1,
  already at its occupancy-max default), deferred spilling genuinely changed the
  allocation and still only found 2 registers. Combined with Lesson 8's "nothing
  cheap is left to remove", this is evidence that the flux kernel's 32 AGPRs are
  *load-bearing live values*, not allocator sloppiness. That predicts the
  remaining Tier-2 flags (`--regalloc-eviction-max-interference-cutoff`,
  `--vgpr-regalloc=basic`, `--split-spill-mode`) will behave the same way: a
  couple of registers, not 32. **Nothing in the flag route looks capable of
  covering a 32-register gap.**
- **⚠ TU-scoped flags CANNOT be aimed at the flux kernel alone.** `CalcFluxAll`,
  `CalcE_impl`, `CalcFstag` and `CalcAuxTermsForAvecPsiRHS` all live in
  `fluxes.cxx`, so every `-mllvm` flag hits all 11 kernels — and probe 2 shows
  that collateral is not hypothetical. This is a *structural* disadvantage of the
  flag route versus `launch_bounds`, which is per-launch-site by construction and
  therefore scoped exactly to the one kernel. Judge any future flag winner on
  **all 11 rows**, not just the target row.

### Probe 3 — `-ffinite-math-only` — **NEGATIVE** (2026-07-27)

| kernel | SGPR | VGPR | AGPR | scratch | spill v/s | occ | vs baseline |
|---|---:|---:|---:|---:|---:|---:|---|
| FLUX uct1 pp0 idealgas | 98 | 256 | **34** | 3448 | 0/0 | **1** | SGPR −2, **AGPR +2 (WORSE)** |
| FLUX uct0 pp0 idealgas | 98 | 252 | 2 | 3448 | 0/0 | 2 | SGPR −2, else identical |
| EMF CalcE uct1 (×3) | 106 | 77 | 1 | 72 | 0/**24,22,20** | 5 | same multiset {20,22,24}, permuted across directions |
| EMF CalcE uct0 (×3) | 36 | 16–19 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcFstag (×3) | 72 | 52–56 | 0 | 0 | 0/0 | 8 | identical |
| EMF CalcAux (×2) | 58, 92 | 15, 38 | 0 | 0 | 0/0 | 8 | identical |

Relaxing FP semantics made the target **worse** (AGPR 32 → 34), traded 2 SGPRs
for it, and left occupancy at 1. The `CalcE uct1` spill counts are the same three
values reassigned to different directions — allocation-order noise, not a trend.

**⚠ This probe probably did NOT test what it was designed to test.** The STATUS
table claimed it would price the `vbar` finding. Re-deriving the fold: turning
`0.0 * x` into a constant needs `nnan` (so `Inf * 0 → NaN` can't happen) **and
`nsz`**, because without no-signed-zeros the result is `+0.0` or `-0.0` depending
on the sign of `x` and is therefore not a compile-time constant.
`-ffinite-math-only` sets `nnan`/`ninf` but **not `nsz`** — that is
`-fno-signed-zeros`. So the dead `gf_vels` loads were most likely still required
to stay live, and the +2 AGPR is unrelated perturbation from other finite-math
simplifications.

**Do not spend another build closing this gap.** Two ways forward, and the second
is strictly better:

1. Re-run with `-ffinite-math-only -fno-signed-zeros` (or `-ffast-math`) to make
   the diagnostic actually valid. Costs a build, answers only a question about a
   flag we will never ship.
2. **Just write the `vbar` fix and measure it.** It is ~10 lines of
   `if constexpr (pplim)` (recipe above), effectively bit-identical, expected
   golden = 0, and it is *the actual candidate* — the flag was only ever a proxy
   for it. A direct `-Rpass` on the real change is a definitive answer where the
   flag is an inference. **Recommended.**

What probe 3 *does* establish, independent of the `vbar` question: **there is no
free lunch in relaxed FP for this kernel's register pressure.** Whatever
`-ffinite-math-only` deletes, the allocator does not convert into occupancy — it
came out 2 AGPRs behind. Combined with the accuracy hazards already documented
(`isnan()` folding, `tau` cancellation, the eigenvalue discriminant), the case for
shipping any fast-math flag here is now empirically dead as well as
theoretically bad.

## STATUS / next actions

Three builds queued, deliberately three *different mechanisms* so a flat result
on all three is itself conclusive:

| # | `PROBE_FLAGS` | tests | result |
|---|---|---|---|
| 1 | `-mllvm --misched=gcn-max-occupancy` | scheduling | **NULL** — but ≈ the default; see above |
| 2 | `-mllvm --enable-deferred-spilling` | allocation | **AGPR 32→30, occ-1. Reject** (costs CalcE spills) |
| 3 | `-ffinite-math-only` | dead code (diagnostic) | **AGPR 32→34, occ-1. Negative** (and premise flawed) |

## VERDICT: the compiler-flag route is CLOSED (2026-07-27)

All three mechanisms are spent, and the score across 11 kernels is:

| mechanism | best effect on the target (`uct1`, needs AGPR 32→0) |
|---|---|
| scheduling | 0 registers (already the default strategy) |
| allocation | −2 registers, plus +4…+10 SGPR spills in `CalcE uct1` |
| relaxed FP | **+2 registers** (worse), and never shippable anyway |

Nothing here can cover a 32-register gap; the spread across three independent
mechanisms is ±2. Two structural findings explain why, and both should be treated
as settled:

- **The 32 AGPRs are load-bearing live values, not allocator sloppiness.** The one
  flag that genuinely re-ran the allocation found 2 registers, and the one that
  relaxed the *semantics* found −2. This is the same conclusion Lesson 8 reached
  from the source side, now confirmed from the compiler side. Register pressure
  here is a property of the algorithm's live ranges, not of a bad heuristic.
- **TU-scoped flags cannot be aimed at the flux kernel.** All 11 kernels compile
  in `fluxes.cxx`, so every `-mllvm` flag hits `CalcE`/`CalcFstag`/`CalcAux` too,
  and probe 2 shows that collateral is real. `launch_bounds` is per-launch-site by
  construction — a structural advantage, not just a stronger knob.

The remaining candidate list has been **deleted, not deferred** (2026-07-27), and
`PROBE_FLAGS` is gone from `make.code.deps`. Every entry was a variation on one of
the two mechanisms already measured at ±2, so working down the list buys builds,
not registers. Treat this route as exhausted: if a future agent wants occupancy,
the instrument is `launch_bounds`, not another flag.

## NEXT (in order)

1. **Write the `vbar` fix and `-Rpass` it directly.** ~10 lines of
   `if constexpr (pplim)`, recipe above. This is the last genuine removal
   candidate and the only outstanding question the flag sweep failed to answer
   (probe 3's premise was flawed — it needed `-fno-signed-zeros`). Effectively
   bit-identical → golden expected exactly 0, so it is worth doing on its own
   merits (12 dead `gf_vels` loads) even at AGPR-neutral.
2. **`launch-bounds-plan.md` Changes 1+2** — now the primary route, not the
   fallback. The two corrections recorded there still apply (the replication must
   call `amrex::detail::call_f_intvect_handler`; a 7th template argument breaks
   the CI build against stock CarpetX, so decide the CI story up front).
3. If `launch_bounds` reaches occ-2 at scratch ≈ 3448: golden gate, then timing at
   **both** grid sizes vs the Idea-1 baseline. If it reaches occ-2 but the small
   grid still regresses, the honest outcome is to ship Idea-1 alone (occ-1, low
   scratch, the fusion + `use_pplim` removal already banked) and stop chasing
   occ-2.

Build 3 is the highest-information diagnostic: a large AGPR drop says the dead
`vbar` FP is the bulk of the problem and the bit-identical `if constexpr` fix
captures it without shipping fast-math; a flat result retires both fast-math and
`vbar` as *occupancy* levers (the `vbar` fix may still be worth it for the removed
memory traffic).

Then, in order: the `vbar` source fix (independent of flags); full-build
re-measurement of any winner; golden gate (expect exactly 0 for flag-only
winners); timing at both grid sizes.

**If nothing reaches occ-2:** fall back to `launch-bounds-plan.md` Changes 1+2,
with the two corrections recorded there (the replication must call
`amrex::detail::call_f_intvect_handler`; and passing a 7th template argument
breaks the CI build against stock CarpetX).

## Deleted helpers (2026-07-25)

`Planning/rpass-table.sh` and `Planning/sweep-regalloc-flags.sh` were removed as
unnecessary once the sweep went manual. Everything they encoded — the baselines,
the flag tiers, the mangled-name patterns, the two `-Rpass` parsing gotchas, and
the extraction one-liner — is preserved above.
