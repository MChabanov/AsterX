#!/usr/bin/env python3
"""
make_beta_slice_from_compose_h5.py

Extract a beta-equilibrium slice from an already-created full 3D CompOSE HDF5
table by solving

    mu_l(T, n_b, Y_q) = 0

along the Y_q direction at each (T, n_b).

For the HS(DD2) HDF5 table generated with the selection

    index_thermo = 1, 2, 3, 4, 5, 7, 12

CompOSE thermodynamic index 5 is the lepton chemical potential mu_l. For
neutrino-less beta equilibrium, mu_l = 0 is the condition used here.

Example usage:

    python make_beta_slice_from_compose_h5.py eos.table.HSDD2.h5 \
      --out-prefix HSDD2_beta_allT

Lowest-temperature/cold slice only:

    python make_beta_slice_from_compose_h5.py eos.table.HSDD2.h5 \
      --cold --out-prefix HSDD2_beta_cold

A specific temperature index:

    python make_beta_slice_from_compose_h5.py eos.table.HSDD2.h5 \
      --temperature-index 0 --out-prefix HSDD2_beta_T0

Outputs:

    <out-prefix>.dat
        ASCII table with columns:
        T, nb, Yq_beta, status,
        all selected thermodynamic quantities,
        optional composition pair quantities,
        optional average-nucleus quantities.

    <out-prefix>.h5
        Compact postprocessed HDF5 file with /Beta_slice/table.
        This is NOT the same 3D CompOSE layout as the input file.

status column:
    0 = mu_l=0 crossing was bracketed and linearly interpolated.
    1 = no crossing was bracketed; nearest |mu_l| point was used.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np


THERMO_NAMES = {
    1: "press_MeV_fm^-3",
    2: "entropy",
    3: "mu_b_minus_m_n_MeV",
    4: "mu_q_MeV",
    5: "mu_l_MeV",
    7: "eps_scaled_E_over_mn_minus_1",
    12: "cs2",
    15: "Gamma",
    21: "E_MeV",
}


def map_axes_by_unique_size(shape: Tuple[int, ...], sizes: Dict[str, int]) -> Dict[str, int]:
    """
    Map raw HDF5 dataset axes to nb/t/yq/q by matching dimension sizes.

    For HS(DD2), the dimensions are usually unique:
        nb=326, t=81, yq=60, q=7 or 3 or 1.
    """
    axes: Dict[str, int] = {}
    used = set()

    for name, size in sizes.items():
        matches = [i for i, s in enumerate(shape) if s == size and i not in used]
        if len(matches) != 1:
            raise RuntimeError(
                f"Could not uniquely identify axis {name} of size {size}. "
                f"Dataset shape is {shape}; matches={matches}. "
                "If two dimensions have the same size, edit the transpose logic manually."
            )
        axes[name] = matches[0]
        used.add(matches[0])

    return axes


def to_nb_t_yq_q(raw: np.ndarray, nrho: int, ntemp: int, nye: int, nq: int) -> np.ndarray:
    axes = map_axes_by_unique_size(raw.shape, {"nb": nrho, "t": ntemp, "yq": nye, "q": nq})
    return np.transpose(raw, (axes["nb"], axes["t"], axes["yq"], axes["q"]))


def interp_along_yq_to_mu_target(
    yq: np.ndarray,
    values: np.ndarray,
    mu_l: np.ndarray,
    target_mu_l: float,
) -> Tuple[np.ndarray, float, int]:
    """
    Interpolate all values to mu_l = target_mu_l along Y_q.

    values has shape (nye, ncols).
    """
    f = mu_l - target_mu_l
    finite = np.isfinite(yq) & np.isfinite(f) & np.all(np.isfinite(values), axis=1)

    if np.count_nonzero(finite) < 2:
        return np.full(values.shape[1], np.nan), np.nan, 1

    y = yq[finite]
    ff = f[finite]
    vv = values[finite, :]

    # Exact/nearly exact point.
    iz = int(np.argmin(np.abs(ff)))
    if abs(ff[iz]) < 1.0e-12:
        return vv[iz].copy(), float(y[iz]), 0

    # Locate sign changes.
    changes = np.where(ff[:-1] * ff[1:] <= 0.0)[0]

    if len(changes) == 0:
        # No bracket; use nearest |mu_l|.
        return vv[iz].copy(), float(y[iz]), 1

    # If multiple crossings exist, choose the local bracket closest to zero.
    i0 = int(min(changes, key=lambda i: min(abs(ff[i]), abs(ff[i + 1]))))
    i1 = i0 + 1

    f0 = ff[i0]
    f1 = ff[i1]
    if f1 == f0:
        w = 0.0
    else:
        w = -f0 / (f1 - f0)

    w = float(np.clip(w, 0.0, 1.0))
    y_beta = y[i0] + w * (y[i1] - y[i0])
    out = vv[i0, :] + w * (vv[i1, :] - vv[i0, :])

    return out, float(y_beta), 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_h5", help="Full 3D CompOSE HDF5 table")
    parser.add_argument("--out-prefix", default="beta_slice", help="Output prefix")
    parser.add_argument("--target-mul", type=float, default=0.0, help="Target mu_l in MeV")
    parser.add_argument("--cold", action="store_true", help="Only extract the lowest-temperature slice")
    parser.add_argument("--temperature-index", type=int, default=None, help="Extract one temperature index")
    args = parser.parse_args()

    input_h5 = Path(args.input_h5)
    out_prefix = Path(args.out_prefix)

    with h5py.File(input_h5, "r") as f:
        nb = np.asarray(f["/Parameters/nb"][:], dtype=float)
        temp = np.asarray(f["/Parameters/t"][:], dtype=float)
        yq = np.asarray(f["/Parameters/yq"][:], dtype=float)

        index_thermo = np.asarray(f["/Thermo_qty/index_thermo"][:], dtype=int)
        thermo_raw = np.asarray(f["/Thermo_qty/thermo"][:], dtype=float)

        nrho = len(nb)
        ntemp = len(temp)
        nye = len(yq)
        nthermo = len(index_thermo)

        thermo = to_nb_t_yq_q(thermo_raw, nrho, ntemp, nye, nthermo)

        if 5 not in index_thermo:
            raise RuntimeError("Thermodynamic index 5, mu_l, is not present.")
        i_mul = int(np.where(index_thermo == 5)[0][0])

        comp_pair_indices = np.array([], dtype=int)
        comp_pairs = None
        if "/Composition_pairs" in f:
            comp_pair_indices = np.asarray(f["/Composition_pairs/index_yi"][:], dtype=int)
            yi_raw = np.asarray(f["/Composition_pairs/yi"][:], dtype=float)
            comp_pairs = to_nb_t_yq_q(yi_raw, nrho, ntemp, nye, len(comp_pair_indices))

        av_indices = np.array([], dtype=int)
        av_data = None
        av_group_name = None
        if "/Composition_quadrupels" in f:
            av_group_name = "/Composition_quadrupels"
        elif "/Composition_quadruples" in f:
            av_group_name = "/Composition_quadruples"

        if av_group_name is not None:
            g = f[av_group_name]
            av_indices = np.asarray(g["index_av"][:], dtype=int)
            navsets = len(av_indices)

            yav = to_nb_t_yq_q(np.asarray(g["yav"][:], dtype=float), nrho, ntemp, nye, navsets)
            aav = to_nb_t_yq_q(np.asarray(g["aav"][:], dtype=float), nrho, ntemp, nye, navsets)
            zav = to_nb_t_yq_q(np.asarray(g["zav"][:], dtype=float), nrho, ntemp, nye, navsets)
            nav = to_nb_t_yq_q(np.asarray(g["nav"][:], dtype=float), nrho, ntemp, nye, navsets)

            # Interleave by set: Yav, Aav, Zav, Nav for each set.
            av_parts = []
            for iset in range(navsets):
                av_parts.extend(
                    [
                        yav[:, :, :, iset:iset+1],
                        aav[:, :, :, iset:iset+1],
                        zav[:, :, :, iset:iset+1],
                        nav[:, :, :, iset:iset+1],
                    ]
                )
            av_data = np.concatenate(av_parts, axis=3)

    if args.temperature_index is not None:
        t_indices = [args.temperature_index]
    elif args.cold:
        t_indices = [0]
    else:
        t_indices = list(range(ntemp))

    for jt in t_indices:
        if jt < 0 or jt >= ntemp:
            raise RuntimeError(f"Temperature index {jt} is outside [0, {ntemp - 1}]")

    thermo_col_names = [THERMO_NAMES.get(int(idx), f"thermo_{int(idx)}") for idx in index_thermo]
    comp_col_names: List[str] = [f"Y_particle_{int(idx)}" for idx in comp_pair_indices]

    av_col_names: List[str] = []
    for idx in av_indices:
        av_col_names.extend(
            [
                f"Yav_set_{int(idx)}",
                f"Aav_set_{int(idx)}",
                f"Zav_set_{int(idx)}",
                f"Nav_set_{int(idx)}",
            ]
        )

    columns = ["T_MeV", "nb_fm^-3", "Yq_beta", "status"]
    columns += thermo_col_names + comp_col_names + av_col_names

    rows = []
    n_cross = 0
    n_fallback = 0

    for jt in t_indices:
        for ir in range(nrho):
            pieces = [thermo[ir, jt, :, :]]

            if comp_pairs is not None:
                pieces.append(comp_pairs[ir, jt, :, :])

            if av_data is not None:
                pieces.append(av_data[ir, jt, :, :])

            values = np.concatenate(pieces, axis=1)
            mu_l = thermo[ir, jt, :, i_mul]

            out_values, yq_beta, status = interp_along_yq_to_mu_target(
                yq, values, mu_l, args.target_mul
            )

            if status == 0:
                n_cross += 1
            else:
                n_fallback += 1

            row = [temp[jt], nb[ir], yq_beta, float(status)]
            row.extend(out_values.tolist())
            rows.append(row)

    table = np.asarray(rows, dtype=float)

    dat_path = out_prefix.with_suffix(".dat")
    h5_path = out_prefix.with_suffix(".h5")

    np.savetxt(dat_path, table, header=" ".join(columns), fmt="%.16e")

    with h5py.File(h5_path, "w") as fout:
        fout.attrs["source_file"] = str(input_h5)
        fout.attrs["beta_condition"] = "mu_l = target_mu_l"
        fout.attrs["target_mu_l_MeV"] = args.target_mul
        fout.attrs["status_meaning"] = "0=crossing found; 1=no crossing, nearest |mu_l| used"

        gp = fout.create_group("Parameters")
        gp.create_dataset("nb", data=nb)
        gp.create_dataset("t", data=temp[t_indices])
        gp.create_dataset("columns", data=np.asarray(columns, dtype="S"))

        gb = fout.create_group("Beta_slice")
        gb.create_dataset("table", data=table)
        gb.create_dataset("index_thermo", data=index_thermo)
        if comp_pair_indices.size:
            gb.create_dataset("index_yi", data=comp_pair_indices)
        if av_indices.size:
            gb.create_dataset("index_av", data=av_indices)

    print(f"Wrote {dat_path}")
    print(f"Wrote {h5_path}")
    print(f"Rows: {len(rows)}")
    print(f"Crossings found: {n_cross}")
    print(f"Fallback/min-|mu_l| points: {n_fallback}")
    if n_fallback:
        print("WARNING: Some points did not bracket mu_l=0 over the available Yq range.")


if __name__ == "__main__":
    main()
