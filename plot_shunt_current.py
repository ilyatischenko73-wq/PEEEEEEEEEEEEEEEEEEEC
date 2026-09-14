from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def frame_number(path: Path) -> int:
    match = re.search(r"_(\d+)\.dat$", path.name)

    if match is None:
        return -1

    return int(match.group(1))


def read_shunt_file(path: Path) -> tuple[float, float]:
    time_value = None
    current_value = None

    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            parts = raw_line.split()

            if len(parts) < 2:
                continue

            if parts[0] == "time":
                time_value = float(parts[1])

            elif parts[0] == "shunt_current_A":
                current_value = float(parts[1])

    if time_value is None:
        raise RuntimeError(
            f"{path}: field 'time' not found."
        )

    if current_value is None:
        raise RuntimeError(
            f"{path}: field 'shunt_current_A' not found."
        )

    return time_value, current_value


def load_series(directory: Path, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    files = sorted(
        directory.glob(f"{prefix}_shunt_*.dat"),
        key=frame_number,
    )

    if not files:
        raise FileNotFoundError(
            f"No files found: {directory}/{prefix}_shunt_*.dat"
        )

    times = []
    currents = []

    for path in files:
        t, current = read_shunt_file(path)

        times.append(t)
        currents.append(current)

    times = np.asarray(times, dtype=float)
    currents = np.asarray(currents, dtype=float)

    order = np.argsort(times)

    return times[order], currents[order]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Publication-style plot of PEEC shunt current versus time."
        )
    )

    parser.add_argument(
        "--directory",
        default="results/TEST",
        help="Directory containing lightning_stage2_shunt_*.dat",
    )

    parser.add_argument(
        "--prefix",
        default="lightning_stage2",
        help="VTK/shunt output prefix.",
    )

    parser.add_argument(
        "--output",
        default="results/TEST/I_shunt_vs_time.pdf",
        help="Output figure file (.pdf, .png, .svg, ...).",
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

    directory = Path(args.directory)
    output = Path(args.output)

    time_s, current_a = load_series(
        directory,
        args.prefix,
    )

    time_us = 1.0e6 * time_s

    peak_index = int(
        np.argmax(
            np.abs(current_a)
        )
    )

    peak_time_us = float(
        time_us[peak_index]
    )

    peak_current_a = float(
        current_a[peak_index]
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

    current_abs_a = np.abs(
        current_a
    )

    ax.plot(
        time_us,
        current_abs_a,
        linewidth=1.8,
        label=r"$|I_{\mathrm{shunt}}(t)|$",
    )

    ax.set_xlabel(
        r"Time, $t$ [$\mu$s]"
    )

    ax.set_ylabel(
        r"Shunt current magnitude, $|I_{\mathrm{shunt}}|$ [A]"
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
    print("PEEC SHUNT CURRENT PLOT")
    print("=" * 72)
    print(
        f"directory          : {directory}"
    )
    print(
        f"samples            : {len(time_s)}"
    )
    print(
        f"time interval      : {time_us[0]:.6f} ... "
        f"{time_us[-1]:.6f} us"
    )
    print(
        f"peak |I_shunt|     : {abs(peak_current_a):.9e} A"
    )
    print(
        f"I_shunt at peak    : {peak_current_a:.9e} A"
    )
    print(
        f"peak time          : {peak_time_us:.9f} us"
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
