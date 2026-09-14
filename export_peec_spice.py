#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

try:
    from scipy.io import mmread
except ImportError as exc:
    raise SystemExit(
        "SciPy is required.\n"
        "Install it in the environment with:\n"
        "    python -m pip install scipy numpy\n"
    ) from exc


def load_mtx(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(path)

    value = mmread(path)

    if hasattr(value, "toarray"):
        value = value.toarray()

    return np.asarray(value, dtype=float)


def spice_number(value: float) -> str:
    """
    Use scientific notation explicitly.
    This avoids ambiguity of SPICE suffixes (m = milli, meg = mega).
    """
    return f"{value:.16e}"


def branch_nodes_from_A(A: np.ndarray) -> list[tuple[int, int]]:
    """
    PEEC convention:
        A[e,a] = -1
        A[e,b] = +1

    Returns:
        (a,b) for every edge e, corresponding to global orientation a -> b.
    """
    ne, nv = A.shape
    result: list[tuple[int, int]] = []

    for e in range(ne):
        row = A[e]

        neg = np.flatnonzero(np.isclose(row, -1.0, atol=1e-12))
        pos = np.flatnonzero(np.isclose(row, +1.0, atol=1e-12))

        if len(neg) != 1 or len(pos) != 1:
            raise RuntimeError(
                f"A row {e}: expected one -1 and one +1, "
                f"got {len(neg)} negative and {len(pos)} positive entries."
            )

        result.append((int(neg[0]), int(pos[0])))

    return result


def validate_inputs(
    L: np.ndarray,
    P: np.ndarray,
    A: np.ndarray,
    R: np.ndarray,
) -> None:
    if L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise RuntimeError(f"L must be square, got {L.shape}")

    if P.ndim != 2 or P.shape[0] != P.shape[1]:
        raise RuntimeError(f"P must be square, got {P.shape}")

    ne = L.shape[0]
    nv = P.shape[0]

    if A.shape != (ne, nv):
        raise RuntimeError(
            f"A shape {A.shape}, expected {(ne, nv)}"
        )

    if R.shape != (ne, ne):
        raise RuntimeError(
            f"R shape {R.shape}, expected {(ne, ne)}"
        )

    for name, matrix in (
        ("L", L),
        ("P", P),
        ("A", A),
        ("R", R),
    ):
        if not np.all(np.isfinite(matrix)):
            raise RuntimeError(f"{name} contains NaN/Inf")

    l_sym = float(np.max(np.abs(L - L.T)))
    p_sym = float(np.max(np.abs(P - P.T)))

    if l_sym > 1e-8 * max(1e-300, float(np.max(np.abs(L)))):
        raise RuntimeError(f"L is not symmetric enough: {l_sym:.3e}")

    if p_sym > 1e-8 * max(1e-300, float(np.max(np.abs(P)))):
        raise RuntimeError(f"P is not symmetric enough: {p_sym:.3e}")


def make_capacitance_network(
    P: np.ndarray,
    relative_cutoff: float,
) -> tuple[
    np.ndarray,
    list[tuple[int, int, float]],
    list[tuple[int, float]],
    float,
]:
    """
    Convert coefficient-of-potential matrix

        V = P Q

    to Maxwell capacitance matrix

        Q = C V,  C = P^{-1}.

    Standard capacitor realization:
        C_ij(branch) = -C[i,j]  for i != j
        C_i0         = sum_j C[i,j]

    Signed branch values are retained: positive definiteness belongs to the
    complete nodal matrix C, not to each branch in this algebraic realization.
    max_negative_branch reports the largest negative branch magnitude.
    """
    if not np.isfinite(relative_cutoff) or relative_cutoff < 0:
        raise ValueError("Capacitance cutoff must be finite and nonnegative")
    if not np.all(np.isfinite(P)) or not np.allclose(P, P.T, rtol=1e-12, atol=0):
        raise ValueError("Potential matrix must be finite and symmetric")
    try:
        np.linalg.cholesky(P)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Potential matrix is not positive definite; no passive capacitance model") from exc
    C = np.linalg.solve(P, np.eye(P.shape[0]))
    C = 0.5 * (C + C.T)

    scale = max(
        float(np.max(np.abs(C))),
        1e-300,
    )
    cutoff = relative_cutoff * scale

    nv = C.shape[0]

    mutual_caps: list[tuple[int, int, float]] = []
    ground_caps: list[tuple[int, float]] = []

    max_negative_branch = 0.0

    for i in range(nv):
        for j in range(i + 1, nv):
            c = -float(C[i, j])

            if c < 0.0:
                max_negative_branch = max(
                    max_negative_branch,
                    -c,
                )

            if abs(c) > cutoff:
                mutual_caps.append(
                    (i, j, c)
                )

    row_sum = np.sum(C, axis=1)

    for i, c0 in enumerate(row_sum):
        c0 = float(c0)

        if abs(c0) > cutoff:
            ground_caps.append(
                (i, c0)
            )

    exported = np.zeros_like(C)
    for i, j, c in mutual_caps:
        exported[i, i] += c
        exported[j, j] += c
        exported[i, j] -= c
        exported[j, i] -= c
    for i, c in ground_caps:
        exported[i, i] += c
    try:
        np.linalg.cholesky(exported)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Capacitance cutoff destroys positive definiteness; reduce --cap-cutoff") from exc
    return C, mutual_caps, ground_caps, max_negative_branch


def write_netlist(
    output: Path,
    L: np.ndarray,
    P: np.ndarray,
    A: np.ndarray,
    R: np.ndarray,
    strike_node: int | None,
    return_node: int | None,
    shunt_a_node: int | None,
    shunt_b_node: int | None,
    shunt_R: float,
    K_amp: float,
    alpha: float,
    beta: float,
    delay: float,
    tstep: float,
    tstop: float,
    r_floor: float,
    cap_cutoff: float,
    k_cutoff: float,
) -> dict[str, float | int]:
    ne = L.shape[0]
    nv = P.shape[0]

    edge_nodes = branch_nodes_from_A(A)

    diag_L = np.diag(L)

    if np.any(diag_L <= 0.0):
        bad = int(np.flatnonzero(diag_L <= 0.0)[0])
        raise RuntimeError(
            f"L[{bad},{bad}] must be positive, got {diag_L[bad]:.6e}"
        )

    # R in the present PEEC implementation is diagonal.
    offdiag_R = R - np.diag(np.diag(R))
    max_offdiag_R = float(np.max(np.abs(offdiag_R)))

    if max_offdiag_R > 1e-12 * max(
        1.0,
        float(np.max(np.abs(R))),
    ):
        raise RuntimeError(
            "The SPICE exporter currently supports diagonal R only. "
            f"max off-diagonal R = {max_offdiag_R:.6e}"
        )

    diag_R = np.diag(R)

    C, mutual_caps, ground_caps, max_negative_cap = (
        make_capacitance_network(
            P,
            cap_cutoff,
        )
    )

    lines: list[str] = []

    lines.append("* ============================================================")
    lines.append("* PEEC -> SPICE netlist")
    lines.append("* Generated from C-solver Matrix Market matrices")
    lines.append("*")
    lines.append("* PEEC equations represented:")
    lines.append("*   L dI/dt + R I + A V = 0")
    lines.append("*   Q = C V,  C = inv(P)")
    lines.append("*")
    lines.append("* Node PEEC j is SPICE node N{j}.")
    lines.append("* SPICE node 0 represents electrostatic infinity/reference.")
    lines.append("* ============================================================")
    lines.append("")

    # --------------------------------------------------------
    # PEEC edge branches
    #
    # Each edge gets a private internal node XEe so R and L are
    # in series between PEEC end nodes.
    # --------------------------------------------------------
    lines.append("* --- PEEC edge branches: R_e + L_ee ---")

    for e, (a, b) in enumerate(edge_nodes):
        r = max(float(diag_R[e]), r_floor)
        l = float(diag_L[e])

        lines.append(
            f"RPEEC{e} N{a} XE{e} {spice_number(r)}"
        )
        lines.append(
            f"LPEEC{e} XE{e} N{b} {spice_number(l)}"
        )

    lines.append("")

    # --------------------------------------------------------
    # Mutual inductive coupling
    # --------------------------------------------------------
    lines.append("* --- Mutual partial inductances ---")

    n_k = 0
    max_abs_k = 0.0

    for i in range(ne):
        li = float(diag_L[i])

        for j in range(i + 1, ne):
            lj = float(diag_L[j])

            kij = float(
                L[i, j] / np.sqrt(li * lj)
            )

            max_abs_k = max(
                max_abs_k,
                abs(kij),
            )

            if abs(kij) <= k_cutoff:
                continue

            # Numerical integration can occasionally produce |k| just over 1.
            if abs(kij) >= 1.0:
                if abs(kij) <= 1.0 + 1e-9:
                    kij = np.copysign(
                        1.0 - 1e-12,
                        kij,
                    )
                else:
                    raise RuntimeError(
                        f"Nonphysical coupling |k| >= 1 for edges "
                        f"{i},{j}: k={kij:.12e}"
                    )

            lines.append(
                f"KPEEC{i}_{j} LPEEC{i} LPEEC{j} "
                f"{spice_number(kij)}"
            )
            n_k += 1

    lines.append("")

    # --------------------------------------------------------
    # Electrostatic capacitance network
    # --------------------------------------------------------
    lines.append("* --- Electrostatic network C = inv(P) ---")
    lines.append("* Signed capacitors preserve the complete positive-definite nodal matrix.")
    lines.append("* Individual negative branches are algebraic elements, not physical parts.")

    for index, (i, j, value) in enumerate(mutual_caps):
        lines.append(
            f"CPEECM{index} N{i} N{j} {spice_number(value)}"
        )

    for index, (i, value) in enumerate(ground_caps):
        lines.append(
            f"CPEECG{index} N{i} 0 {spice_number(value)}"
        )

    lines.append("")

    # --------------------------------------------------------
    # Lightning source
    # --------------------------------------------------------
    if (strike_node is None) != (return_node is None):
        raise RuntimeError(
            "Specify both --strike-node and --return-node, or neither."
        )

    if strike_node is not None:
        if not (0 <= strike_node < nv):
            raise RuntimeError("--strike-node out of range")

        if not (0 <= return_node < nv):
            raise RuntimeError("--return-node out of range")

        lines.append("* --- Lightning current source ---")
        lines.append(
            "* SPICE source current is defined from first node to second."
        )
        lines.append(
            "* RETURN -> STRIKE injects +I_L into the PEEC strike node."
        )

        expression = (
            f"(time<{spice_number(delay)} ? 0 : "
            f"{spice_number(K_amp)}*("
            f"exp(-{spice_number(alpha)}*(time-{spice_number(delay)}))"
            f"-exp(-{spice_number(beta)}*(time-{spice_number(delay)}))))"
        )

        lines.append(
            f"BFLASH N{return_node} FLASH_DRIVE I={expression}"
        )
        lines.append(f"VFLASH FLASH_DRIVE N{strike_node} 0")
        lines.append("")

    # --------------------------------------------------------
    # Internal shunt
    # --------------------------------------------------------
    if (shunt_a_node is None) != (shunt_b_node is None):
        raise RuntimeError(
            "Specify both --shunt-a-node and --shunt-b-node, or neither."
        )

    if shunt_a_node is not None:
        if not (0 <= shunt_a_node < nv):
            raise RuntimeError("--shunt-a-node out of range")

        if not (0 <= shunt_b_node < nv):
            raise RuntimeError("--shunt-b-node out of range")

        if shunt_R <= 0.0:
            raise RuntimeError("--shunt-R must be > 0")

        lines.append("* --- Internal resistive shunt ---")
        lines.append(
            f"RSHUNT N{shunt_a_node} N{shunt_b_node} "
            f"{spice_number(shunt_R)}"
        )
        lines.append("")

    # --------------------------------------------------------
    # Transient analysis
    # --------------------------------------------------------
    lines.append("* --- Transient analysis ---")
    lines.append(
        f".tran {spice_number(tstep)} {spice_number(tstop)} 0 {spice_number(tstep)} uic"
    )
    lines.append("")

    lines.append(".control")
    lines.append("set wr_singlescale")
    lines.append("set wr_vecnames")
    lines.append("run")

    if shunt_a_node is not None:
        lines.append(f"let i_shunt = v(N{shunt_a_node},N{shunt_b_node})/{spice_number(shunt_R)}")
        lines.append(
            "wrdata peec_spice_shunt.dat i_shunt "
            f"v(N{shunt_a_node},N{shunt_b_node})"
        )
    elif strike_node is not None:
        lines.append(
            "wrdata peec_spice_source.dat "
            "i(VFLASH)"
        )

    lines.append(".endc")
    lines.append("")
    lines.append(".end")
    lines.append("")

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    return {
        "n_nodes": nv,
        "n_edges": ne,
        "n_mutual_inductors": n_k,
        "n_mutual_capacitors": len(mutual_caps),
        "n_ground_capacitors": len(ground_caps),
        "max_abs_k": max_abs_k,
        "max_negative_cap": max_negative_cap,
        "capacitance_condition": float(np.linalg.cond(P)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export a PEEC model saved by the C solver in Matrix Market "
            "format to an ngspice-compatible RLC netlist."
        )
    )

    parser.add_argument(
        "--matrix-dir",
        required=True,
        help="Directory containing L.mtx, A.mtx, P.mtx and R.mtx.",
    )

    parser.add_argument(
        "--output",
        default="results/spice/peec_model.cir",
    )

    parser.add_argument("--strike-node", type=int, default=None)
    parser.add_argument("--return-node", type=int, default=None)

    parser.add_argument("--shunt-a-node", type=int, default=None)
    parser.add_argument("--shunt-b-node", type=int, default=None)
    parser.add_argument("--shunt-R", type=float, default=1000.0)

    parser.add_argument("--K", type=float, default=102900.0)
    parser.add_argument("--alpha", type=float, default=1500.0)
    parser.add_argument("--beta", type=float, default=1.0e6)
    parser.add_argument("--delay", type=float, default=0.0)

    parser.add_argument("--dt", type=float, default=1.0e-7)
    parser.add_argument("--t-end", type=float, default=1.0e-4)

    parser.add_argument(
        "--r-floor",
        type=float,
        default=1.0e-12,
        help=(
            "Small series resistance [ohm] used when PEEC R_e is zero. "
            "Helps avoid singular ideal-inductor loops in SPICE."
        ),
    )

    parser.add_argument(
        "--cap-cutoff",
        type=float,
        default=0.0,
        help=(
            "Relative cutoff for capacitances with respect to max |C|. "
            "Use 0 for exact export."
        ),
    )

    parser.add_argument(
        "--k-cutoff",
        type=float,
        default=0.0,
        help=(
            "Absolute cutoff for mutual coupling coefficient |k|. "
            "Use 0 for exact export."
        ),
    )

    parser.add_argument(
        "--max-edges",
        type=int,
        default=200,
        help="Safety limit for direct SPICE export.",
    )

    parser.add_argument(
        "--max-nodes",
        type=int,
        default=200,
        help="Safety limit for direct SPICE export.",
    )

    parser.add_argument(
        "--force-large",
        action="store_true",
        help="Allow export beyond the safety limits.",
    )

    args = parser.parse_args()

    matrix_dir = Path(args.matrix_dir)

    print()
    print("=" * 72)
    print("PEEC -> SPICE EXPORT")
    print("=" * 72)
    print(f"matrix dir : {matrix_dir}")

    L = load_mtx(
        matrix_dir / "L.mtx"
    )
    A = load_mtx(
        matrix_dir / "A.mtx"
    )
    P = load_mtx(
        matrix_dir / "P.mtx"
    )
    R = load_mtx(
        matrix_dir / "R.mtx"
    )

    validate_inputs(
        L,
        P,
        A,
        R,
    )

    ne = L.shape[0]
    nv = P.shape[0]

    print(f"Ne         : {ne}")
    print(f"Nv         : {nv}")

    if not args.force_large:
        if ne > args.max_edges or nv > args.max_nodes:
            n_k = ne * (ne - 1) // 2
            n_c = nv * (nv - 1) // 2

            raise SystemExit(
                "\nDirect SPICE export stopped by safety limit.\n"
                f"Model has Ne={ne}, Nv={nv}.\n"
                f"Potential mutual-inductor pairs: {n_k:,}\n"
                f"Potential mutual-capacitor pairs : {n_c:,}\n\n"
                "Use a coarser/smaller validation mesh first, or pass "
                "--force-large deliberately."
            )

    info = write_netlist(
        output=Path(args.output),
        L=L,
        P=P,
        A=A,
        R=R,
        strike_node=args.strike_node,
        return_node=args.return_node,
        shunt_a_node=args.shunt_a_node,
        shunt_b_node=args.shunt_b_node,
        shunt_R=args.shunt_R,
        K_amp=args.K,
        alpha=args.alpha,
        beta=args.beta,
        delay=args.delay,
        tstep=args.dt,
        tstop=args.t_end,
        r_floor=args.r_floor,
        cap_cutoff=args.cap_cutoff,
        k_cutoff=args.k_cutoff,
    )

    print("-" * 72)
    print(f"mutual K    : {info['n_mutual_inductors']}")
    print(f"mutual C    : {info['n_mutual_capacitors']}")
    print(f"ground C    : {info['n_ground_capacitors']}")
    print(f"max |k|     : {info['max_abs_k']:.9e}")
    print(f"cond(P)     : {info['capacitance_condition']:.9e}")
    print(
        "max negative capacitor (retained) : "
        f"{info['max_negative_cap']:.9e} F"
    )
    print(f"netlist     : {args.output}")
    print("=" * 72)
    print()
    print("Run with:")
    print(f"    ngspice {args.output}")


if __name__ == "__main__":
    main()
