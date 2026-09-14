from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def lightning_current(
    t: np.ndarray,
    K: float,
    alpha: float,
    beta: float,
    delay: float,
) -> np.ndarray:
    tau = t - delay

    current = np.zeros_like(
        t,
        dtype=float,
    )

    mask = tau >= 0.0

    current[mask] = (
        K
        * (
            np.exp(
                -alpha * tau[mask]
            )
            -
            np.exp(
                -beta * tau[mask]
            )
        )
    )

    return current


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Publication-style plot of the PEEC lightning "
            "double-exponential current waveform."
        )
    )

    parser.add_argument(
        "--K",
        type=float,
        default=102900.0,
        help="Amplitude coefficient K [A].",
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=1500.0,
        help="Slow exponential coefficient alpha [1/s].",
    )

    parser.add_argument(
        "--beta",
        type=float,
        default=1.0e6,
        help="Fast exponential coefficient beta [1/s].",
    )

    parser.add_argument(
        "--delay",
        type=float,
        default=0.0,
        help="Time delay [s].",
    )

    parser.add_argument(
        "--t-end",
        type=float,
        default=1.0e-4,
        help="Final time [s].",
    )

    parser.add_argument(
        "--samples",
        type=int,
        default=2000,
        help="Number of plotted samples.",
    )

    parser.add_argument(
        "--output",
        default="results/lightning_double_exponential.pdf",
        help="Output figure file.",
    )

    parser.add_argument(
        "--png",
        default=None,
        help="Optional additional PNG output.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=400,
    )

    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title.",
    )

    args = parser.parse_args()

    if args.samples < 2:
        raise ValueError(
            "--samples must be >= 2."
        )

    if args.t_end <= 0.0:
        raise ValueError(
            "--t-end must be > 0."
        )

    if args.alpha <= 0.0:
        raise ValueError(
            "--alpha must be > 0."
        )

    if args.beta <= 0.0:
        raise ValueError(
            "--beta must be > 0."
        )

    t = np.linspace(
        0.0,
        args.t_end,
        args.samples,
    )

    current = lightning_current(
        t,
        args.K,
        args.alpha,
        args.beta,
        args.delay,
    )

    time_us = 1.0e6 * t
    current_kA = 1.0e-3 * current

    peak_index = int(
        np.argmax(
            current
        )
    )

    peak_time = float(
        t[peak_index]
    )

    peak_current = float(
        current[peak_index]
    )

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "mathtext.fontset": "dejavuserif",
            "figure.dpi": 150,
            "savefig.dpi": args.dpi,
            "axes.linewidth": 0.8,
        }
    )

    fig, ax = plt.subplots(
        figsize=(6.8, 4.2),
        constrained_layout=True,
    )

    ax.plot(
        time_us,
        current_kA,
        linewidth=1.8,
        label=r"$I_L(t)$",
    )

    ax.set_xlabel(
        r"Time, $t$ [$\mu$s]"
    )

    ax.set_ylabel(
        r"Lightning current, $I_L$ [kA]"
    )

    if args.title:
        ax.set_title(
            args.title
        )

    ax.grid(
        True,
        linewidth=0.5,
        alpha=0.3,
    )

    ax.tick_params(
        direction="in",
        top=True,
        right=True,
    )

    ax.legend(
        frameon=False,
        loc="best",
    )

    output = Path(
        args.output
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig.savefig(
        output,
        bbox_inches="tight",
    )

    if args.png is not None:
        png_path = Path(
            args.png
        )

        png_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig.savefig(
            png_path,
            bbox_inches="tight",
        )

    print()
    print("=" * 72)
    print("LIGHTNING DOUBLE-EXPONENTIAL WAVEFORM")
    print("=" * 72)
    print(
        f"K                  : {args.K:.9e} A"
    )
    print(
        f"alpha              : {args.alpha:.9e} 1/s"
    )
    print(
        f"beta               : {args.beta:.9e} 1/s"
    )
    print(
        f"delay              : {args.delay:.9e} s"
    )
    print(
        f"sampled peak time  : {peak_time * 1e6:.9f} us"
    )
    print(
        f"sampled peak       : {peak_current:.9e} A"
    )
    print(
        f"figure             : {output}"
    )

    if args.png is not None:
        print(
            f"PNG                : {args.png}"
        )

    print("=" * 72)

    plt.show()


if __name__ == "__main__":
    main()
