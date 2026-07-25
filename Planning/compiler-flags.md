# Compiler-flag route to occ-2 — checkpoint

_Opened 2026-07-25. Sibling of `CHECKPOINT.md`; read that first for the overall
mission. This doc owns everything about reaching occ-2 on the flux kernel via
**compiler flags** rather than source changes. `launch-bounds-plan.md` is the
fallback if this route closes._

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

## Mechanism — how to run a probe

`AsterX/src/make.code.deps` (LOCAL, uncommitted) carries:

```make
fluxes.cxx.o: CXXFLAGS += -Rpass-analysis=kernel-resource-usage
fluxes.cxx.o: CXXFLAGS += -DASTERX_PROBE_IDEALGAS_ONLY
fluxes.cxx.o: CXXFLAGS += $(PROBE_FLAGS)
```

so candidates go on the gmake command line, no file edit per probe:

```bash
export CACTUS=/ccs/home/mchabanov/EinsteinToolkit/Cactus
export REPO=$CACTUS/repos/AsterX
export CONFIG=cray20-adios-register

cd $CACTUS
rm -f configs/$CONFIG/build/AsterX/fluxes.cxx.o
gmake $CONFIG PROBE_FLAGS="-mllvm --misched=gcn-max-occupancy" 2>&1 | tee ~/probe.log
```

Remarks are emitted during compilation, before linking — safe to `Ctrl-C` once
they scroll past.

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

## Flag inventory (Frontier, cce/20.0.0 + rocm/6.4.2, gfx90a)

Enumerate with:

```bash
$ROCM_PATH/llvm/bin/llc -march=amdgcn -mcpu=gfx90a --help-hidden 2>&1 \
  | grep -iE '^ +--[^ =]*(vgpr|agpr|occup|wave|spill|regalloc|pressure|sched)' | sort -u
# allowed values for =<value> options (there is no `=help`):
$ROCM_PATH/llvm/bin/llc -march=amdgcn -mcpu=gfx90a --help-hidden 2>&1 \
  | grep -A 10 -E '^ +--(vgpr-regalloc|regalloc|split-spill-mode|misched)='
```

### RULED OUT — do not retry

These exist only as IR **function attributes**, not `-mllvm` cl::opts, so they are
unreachable without attaching an attribute to the kernel (= the CarpetX plumbing,
since the kernel is `amrex::launch_global<>` inside AMReX headers):

- `--amdgpu-waves-per-eu` — hard error: *"Unknown command line argument"*,
  suggested `--amdgpu-dce-in-ra`. **This was the plan's "zero-code pre-check"; it
  does not exist.**
- `--amdgpu-spill-vgpr-to-agpr`
- `--amdgpu-num-vgpr` / `--amdgpu-num-sgpr`
- `--vgpr-regalloc=pbqp` — only `basic`/`greedy`/`fast` are offered for the split
  AMDGPU allocators. `--regalloc=pbqp` exists globally but AMDGPU overrides with
  separate `sgpr-`/`vgpr-`/`wwm-regalloc`, so it is likely a silent no-op.
- `--amdgpu-schedule-relaxed-occupancy` — *relaxes* occupancy targets, wrong
  direction.

### AVAILABLE — tiered by prior

**Tier 1 — purpose-built occupancy/register schedulers.** A different scheduler,
not a weight on the default one.

- `--misched=gcn-max-occupancy` — *"Run GCN scheduler to maximize occupancy"*,
  non-experimental. **Highest prior.**
- `--misched=gcn-iterative-minreg` — *"minimal register usage"* (experimental)
- `--misched=gcn-iterative-max-occupancy-experimental`

Caveat on the experimental two: slow to compile, and they can hit their stated
metric while producing worse code — an occ-2 from those needs the timing run more
than the others, not less.

**Tier 2 — allocator-side.** Attacks "the allocator gave up finding a coloring".

- `--enable-deferred-spilling` — *"defer the actual code insertion to the end of
  the allocation. That way the allocator might still find a suitable coloring…
  because of other evicted variables"*. Directly the failure mode.
- `--regalloc-eviction-max-interference-cutoff=100000` — the allocator bails after
  N interferences (default ~10); *"To disable, pass a very large number"*.
- `--split-spill-mode=size|speed`
- `--vgpr-regalloc=basic`

**Tier 2b — inlining.** Direct evidence it matters: the instantiation cut, which
changes *only* inlining of shared callees, moved `uct0` from occ-1 to occ-2.

- `--inline-threshold=100 | 50 | 25`

⚠ Watch **scratch**, not just AGPR: AMDGPU calls need stack frames, so de-inlining
can convert registers into private memory — the Ideas-2/3 seesaw in a new costume.

**Tier 3 — default-scheduler tuning.** Low prior: the flag list shows AMDGPU
*already* runs an unclustered-high-register-pressure reduction stage
(`--amdgpu-disable-unclustered-high-rp-reschedule` exists to turn it off), so the
scheduler is already fighting this.

- `--amdgpu-schedule-metric-bias=100` — *"Set it to 100 to chase the occupancy
  only"*
- `--amdgpu-opt-vgpr-liverange` — *"VGPR liverange optimizations for if-else
  structure"*; the kernel is branch-heavy (`useLO`, `resetL/R`, the `rec_var`
  switch, per-direction box guards)
- `--sink-insts-to-avoid-spills` — sinks "into cycles"; the kernel body is not a
  loop, so probably inert

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

Also note for the record: `-ffp-contract` is unset on Frontier but explicitly
`off` in the CI CPU config, so **golden = 0 on CI implies nothing about
bit-identity of the Frontier binary**. `CHECKPOINT.md` Lesson 5's bit-identity
claim is CI-scoped only. Fine for the workflow (CI validates, Frontier measures)
as long as the two are never conflated.

## The `vbar` finding — a real removal, independent of all flags

Found by comparing `uct1` (288 registers) against `uct0` (254): **UCT is the
expensive config by ~34 registers**, which was not previously documented.

In the `if constexpr (uct)` block, with `pplim=false` (production) `theta_uct` is
the literal `1.0`, so:

```cpp
vbar_j(dir_i)(p.I) = 1.0 * vj_face
                   + (1.0 - 1.0) * 0.5 * (gf_vels(dir_j)(p.I) + gf_vels(dir_j)(p.I - p.DI[dir_i]));
```

`(1.0 - 1.0)` constant-folds to `0.0`, but **`0.0 * x` cannot be folded away**
without fast-math, because `x` could be NaN or ±Inf. So the compiler is *required*
to keep those `gf_vels` loads: 4 per direction, **12 total, provably dead** in the
production config.

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

## STATUS / next actions

Three builds queued, deliberately three *different mechanisms* so a flat result
on all three is itself conclusive:

| # | `PROBE_FLAGS` | tests |
|---|---|---|
| 1 | `-mllvm --misched=gcn-max-occupancy` | scheduling |
| 2 | `-mllvm --enable-deferred-spilling` | allocation |
| 3 | `-ffinite-math-only` | dead code (diagnostic) |

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
