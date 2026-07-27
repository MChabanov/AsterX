# Ideas 2 + 3 — serialize the MHD flux assembly — DONE, superseded

Design note for the flux-by-flux (Idea 2) and face-state (Idea 3) serialization of
`CalcFluxAtFace`. Both implemented and golden PASS at 0, and they did deliver occ-2 for
the production ideal-gas kernel — but at +1536 B/lane of scratch, which cost most of the
benefit (−4.8 % on TOV-large, **+8.9 %** on UCT-small). **Scoped `launch_bounds` later
obtained the same occupancy at +160 B/lane and −27 %**, so this shape is superseded;
tag `archive/serialization-occ2` @ `210a013c` is kept as the record. Retained here: the
structure, the bit-identity argument, and one open cleanup item.

## What was built

- **Idea 2 (`e8613114`)** — serialize the momentum-component axis `j`: a rolled
  `#pragma unroll 1 for(j=0..2)` builds `blows_j`/`moms_j`/`flux_moms_j` one component
  at a time instead of the full `vec<vec<2>,3>`. AGPR 40→1.
- **Idea 3 (`210a013c`)** — serialize the two face states `f`: the whole
  conserved-var/flux assembly becomes a rolled `for(f=0..1)`, with the group-B
  intermediates as per-side scalars and the momentum `for(j)` nested inside. AGPR 1→0,
  **occ 2**.

```
// pre-loop (side-independent + dir_i components): beta_avg, vel_rc, B_rc,
//   vtilde_rc, unit_dir_i, alp_sqrtg, u_avg
// cut-set accumulators (survive the loop for the combine):
//   dens_rc,DEnt_rc,tau_rc,DYe_rc + flux_* (vec<2>); moms_rc,flux_moms
//   (vec<vec<2>,3>); charmax,charmin
#pragma unroll 1
for (f=0..1) {                       // Idea 3
  Bs_f, vlows_f = per-side slices    // reindex, no arithmetic
  group-B scalars: alp_b0_f, Blows_f(vec<3>), B2_f, bsq_f, cs2_f, h_f,
                   dens_h_W_f, dens_h_W_plus_f, press_plus_pmag_f, B_over_w_f
  write dens_rc(f),DEnt_rc(f),tau_rc(f),DYe_rc(f), flux_dens/DEnt/tau/DYe(f)
  #pragma unroll 1
  for (j=0..2) { blows_fj, moms_fj -> moms_rc(j)(f); flux_moms(j)(f) }  // Idea 2
  lambda_f = eigenvalues_oneside(...); fold into charmax/charmin
}
// post-loop combine (needs both sides + charmax/charmin):
//   calcflux/laxf for dens/DEnt/tau/DYe, then rolled for(j) for the 3 moms
```

## Why it is bit-identical

- Rolled loops change only the order and lifetime of *independent* sub-computations,
  never reassociate within an expression → exact under `-ffp-contract=off`.
- Per-side contractions call the SAME scalar `calc_contraction` overloads
  (`AsterUtils/src/aster_utils.hxx`, sequential `sum<3>`) on extracted single-side
  `vec<3>` slices → identical summation order.
- `eigenvalues.hxx`: the two sides were already independent; `eigenvalues_oneside` is
  the verbatim per-side body and `eigenvalues()` now calls it twice. `fmax`/`fmin` over
  both sides equals the old reduction (max/min are exact, hence order-independent).
- **Loops MUST stay rolled** (`#pragma unroll 1`; nvcc/clang/cce honor it, GCC ignores
  it harmlessly) or `-O3` re-merges the live ranges and the saving is lost.

## Two things this taught

- **Ideas 2 and 3 pay off only together.** Idea 3 forces `eigenvalues` per side, so
  `charmax`/`charmin` are ready only after the loop, so the HLL and momentum combines
  must be post-loop, so `moms_rc`/`flux_moms` are held as `vec<vec<2>,3>` for the
  combine — giving back part of Idea 2's saving. The hand estimate of the net was ~0–6
  doubles (ambiguous); empirically it crossed the last granule to occ-2.
- **The rolled indices are exactly why scratch rose**: runtime `f`/`j` cannot index
  registers, so the reconstruction arrays were demoted to private memory
  (`CHECKPOINT.md` Lesson 3). Relocation, not removal.

## Open item, if this branch is ever revived

The `#if 0` CCTK_DEBUG NaN-dump in `fluxes.cxx` references now-serialized quantities
(`cs2_rc`/`h_rc`/`bsq_rc`/`alp_b0_rc`/`Blows_rc`/`dens_h_W_plus`/`press_plus_pmag`/
`B_over_w`). Rewrite it — recompute those as full `vec<2>` locally, or print the
per-side scalars — before upstreaming, so `-DCCTK_DEBUG` builds compile again.
