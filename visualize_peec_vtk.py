from __future__ import annotations

import argparse
import re
import time
from pathlib import Path

import numpy as np
import pyvista as pv


# ============================================================
# HELPERS
# ============================================================

def cell_centers(poly: pv.PolyData) -> np.ndarray:
    return np.asarray(
        poly.cell_centers().points,
        dtype=float,
    )


def model_scale(poly: pv.PolyData) -> float:
    bounds = np.asarray(
        poly.bounds,
        dtype=float,
    )

    extent = np.array(
        [
            bounds[1] - bounds[0],
            bounds[3] - bounds[2],
            bounds[5] - bounds[4],
        ],
        dtype=float,
    )

    scale = float(
        np.max(extent)
    )

    if scale <= 0.0:
        scale = 1.0

    return scale


def safe_unit_vectors(
    vectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    vectors = np.asarray(
        vectors,
        dtype=float,
    )

    magnitudes = np.linalg.norm(
        vectors,
        axis=1,
    )

    unit = np.zeros_like(
        vectors
    )

    mask = magnitudes > 0.0

    unit[mask] = (
        vectors[mask]
        /
        magnitudes[mask, None]
    )

    return unit, magnitudes


def build_arrow_glyphs(
    edge_poly: pv.PolyData,
    vectors: np.ndarray,
    scalar_values: np.ndarray,
    threshold_fraction: float,
    arrow_scale: float,
) -> pv.PolyData | None:
    """
    Build arrows at edge-cell centers.

    Important:
        arrow size is normalized by the maximum vector magnitude,
        so very small physical values do not make invisible glyphs.
    """

    vectors = np.asarray(
        vectors,
        dtype=float,
    )

    scalar_values = np.asarray(
        scalar_values,
        dtype=float,
    )

    if vectors.shape != (edge_poly.n_cells, 3):
        raise ValueError(
            f"vector shape {vectors.shape} != "
            f"({edge_poly.n_cells}, 3)"
        )

    if scalar_values.shape != (edge_poly.n_cells,):
        raise ValueError(
            f"scalar shape {scalar_values.shape} != "
            f"({edge_poly.n_cells},)"
        )

    directions, magnitudes = safe_unit_vectors(
        vectors
    )

    maximum = float(
        np.max(magnitudes)
    )

    if maximum <= 0.0:
        return None

    threshold = (
        threshold_fraction
        *
        maximum
    )

    mask = (
        magnitudes >= threshold
    )

    if not np.any(mask):
        return None

    centers = cell_centers(
        edge_poly
    )

    selected_centers = centers[mask]
    selected_directions = directions[mask]
    selected_magnitudes = magnitudes[mask]
    selected_scalars = scalar_values[mask]

    relative = (
        selected_magnitudes
        /
        maximum
    )

    scale_length = (
        arrow_scale
        *
        model_scale(edge_poly)
    )

    vectors_for_glyph = (
        selected_directions
        *
        relative[:, None]
        *
        scale_length
    )

    cloud = pv.PolyData(
        selected_centers
    )

    cloud.point_data[
        "current_vector"
    ] = vectors_for_glyph

    cloud.point_data[
        "value"
    ] = selected_scalars

    return cloud.glyph(
        orient="current_vector",
        scale="current_vector",
        factor=1.0,
        geom=pv.Arrow(),
    )


def add_base_surface(
    plotter: pv.Plotter,
    surface_file: Path | None,
) -> None:
    if surface_file is None:
        return

    if not surface_file.exists():
        return

    surface = pv.read(
        surface_file
    )

    plotter.add_mesh(
        surface,
        color="lightgray",
        opacity=0.16,
        show_edges=False,
        pickable=False,
    )


def set_camera(
    plotter: pv.Plotter,
    dataset: pv.DataSet,
) -> None:
    bounds = np.asarray(
        dataset.bounds,
        dtype=float,
    )

    z_extent = (
        bounds[5]
        -
        bounds[4]
    )

    total_extent = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )

    if z_extent < 1.0e-10 * total_extent:
        plotter.view_xy()
    else:
        plotter.view_isometric()

    plotter.reset_camera()


# ============================================================
# FIELD DETECTION
# ============================================================

def dataset_kind(
    edge_poly: pv.PolyData,
) -> str:
    names = set(
        edge_poly.cell_data.keys()
    )

    if (
        "current_A" in names
        and
        "J_A_m" in names
    ):
        return "transient"

    if (
        "current_re_A" in names
        and
        "current_im_A" in names
        and
        "J_re_A_m" in names
        and
        "J_im_A_m" in names
    ):
        return "harmonic"

    raise RuntimeError(
        "Unknown PEEC VTK format.\n"
        f"Cell arrays: {sorted(names)}"
    )


# ============================================================
# SHUNT METADATA
# ============================================================

def matching_shunt_file(
    edge_file: Path,
) -> Path | None:
    """
    lightning_stage2_edges_000123.vtk
        ->
    lightning_stage2_shunt_000123.dat
    """

    name = edge_file.name

    if "_edges_" not in name:
        return None

    candidate = edge_file.with_name(
        name.replace(
            "_edges_",
            "_shunt_",
        ).replace(
            ".vtk",
            ".dat",
        )
    )

    return (
        candidate
        if candidate.exists()
        else None
    )


def read_shunt_metadata(
    edge_file: Path,
) -> dict | None:
    path = matching_shunt_file(
        edge_file
    )

    if path is None:
        return None

    data: dict[str, object] = {
        "path": path,
    }

    with path.open(
        "r",
        encoding="utf-8",
    ) as file:
        for raw_line in file:
            parts = raw_line.split()

            if not parts:
                continue

            key = parts[0]

            if key in (
                "node_a",
                "node_b",
            ):
                data[key] = int(
                    parts[1]
                )

            elif key in (
                "xyz_a",
                "xyz_b",
            ):
                data[key] = np.asarray(
                    [
                        float(parts[1]),
                        float(parts[2]),
                        float(parts[3]),
                    ],
                    dtype=float,
                )

            elif key in (
                "time",
                "shunt_current_A",
                "shunt_resistance_ohm",
            ):
                data[key] = float(
                    parts[1]
                )

    required = {
        "node_a",
        "node_b",
        "xyz_a",
        "xyz_b",
        "shunt_current_A",
    }

    if not required.issubset(
        data.keys()
    ):
        return None

    return data


def add_shunt_overlay(
    plotter: pv.Plotter,
    edge_file: Path,
    actor_prefix: str = "shunt",
) -> dict | None:
    data = read_shunt_metadata(
        edge_file
    )

    if data is None:
        return None

    xyz_a = np.asarray(
        data["xyz_a"],
        dtype=float,
    )

    xyz_b = np.asarray(
        data["xyz_b"],
        dtype=float,
    )

    current = float(
        data["shunt_current_A"]
    )

    node_a = int(
        data["node_a"]
    )

    node_b = int(
        data["node_b"]
    )

    # A terminal.
    point_a = pv.PolyData(
        np.asarray(
            [xyz_a],
            dtype=float,
        )
    )

    plotter.add_mesh(
        point_a,
        color="green",
        point_size=22,
        render_points_as_spheres=True,
        name=f"{actor_prefix}_a",
    )

    plotter.add_point_labels(
        np.asarray(
            [xyz_a],
            dtype=float,
        ),
        [
            f"SHUNT A\nnode {node_a}"
        ],
        point_size=0,
        font_size=11,
        always_visible=True,
        shape_color="green",
        text_color="black",
        name=f"{actor_prefix}_a_label",
    )

    # B terminal.
    point_b = pv.PolyData(
        np.asarray(
            [xyz_b],
            dtype=float,
        )
    )

    plotter.add_mesh(
        point_b,
        color="orange",
        point_size=22,
        render_points_as_spheres=True,
        name=f"{actor_prefix}_b",
    )

    plotter.add_point_labels(
        np.asarray(
            [xyz_b],
            dtype=float,
        ),
        [
            f"SHUNT B\nnode {node_b}"
        ],
        point_size=0,
        font_size=11,
        always_visible=True,
        shape_color="orange",
        text_color="black",
        name=f"{actor_prefix}_b_label",
    )

    # Draw the physical shunt branch A -> B.
    line = pv.Line(
        xyz_a,
        xyz_b,
    )

    plotter.add_mesh(
        line,
        color="yellow",
        line_width=6,
        name=f"{actor_prefix}_branch",
    )

    # Current label at the branch midpoint.
    midpoint = (
        0.5
        * (
            xyz_a
            +
            xyz_b
        )
    )

    plotter.add_point_labels(
        np.asarray(
            [midpoint],
            dtype=float,
        ),
        [
            (
                f"I_shunt = {current:.6e} A\n"
                f"positive A -> B"
            )
        ],
        point_size=0,
        font_size=11,
        always_visible=True,
        shape_color="white",
        text_color="black",
        name=f"{actor_prefix}_current_label",
    )

    return data


# ============================================================
# TRANSIENT
# ============================================================

def transient_state(
    edge_poly: pv.PolyData,
) -> tuple[np.ndarray, np.ndarray, float]:
    currents = np.asarray(
        edge_poly.cell_data["current_A"],
        dtype=float,
    )

    vectors = np.asarray(
        edge_poly.cell_data["J_A_m"],
        dtype=float,
    )

    if "time" in edge_poly.cell_data:
        time_array = np.asarray(
            edge_poly.cell_data["time"],
            dtype=float,
        )

        time_value = float(
            time_array[0]
        )
    else:
        time_value = 0.0

    return currents, vectors, time_value


def show_transient_static(
    edge_file: Path,
    surface_file: Path | None,
    threshold: float,
    arrow_scale: float,
    show_arrows: bool,
) -> None:
    edges = pv.read(
        edge_file
    )

    currents, vectors, time_value = (
        transient_state(
            edges
        )
    )

    maximum = float(
        np.max(
            np.abs(currents)
        )
    )

    max_edge = int(
        np.argmax(
            np.abs(currents)
        )
    )

    print()
    print("=" * 72)
    print("PEEC TRANSIENT CURRENT")
    print("=" * 72)
    print(f"file          = {edge_file}")
    print(f"time          = {time_value:.9e} s")
    print(f"time          = {time_value * 1e6:.6f} us")
    print(f"max |I|       = {maximum:.9e} A")
    print(f"max edge      = {max_edge}")
    print("=" * 72)

    plotter = pv.Plotter(
        title="PEEC transient current"
    )

    add_base_surface(
        plotter,
        surface_file,
    )

    edge_view = edges.copy()

    edge_view.cell_data[
        "I_abs_A"
    ] = np.abs(
        currents
    )

    plotter.add_mesh(
        edge_view,
        scalars="I_abs_A",
        preference="cell",
        line_width=5,
        cmap="turbo",
        scalar_bar_args={
            "title": "|I_i| [A]",
        },
    )

    if show_arrows:
        glyphs = build_arrow_glyphs(
            edges,
            vectors,
            np.abs(currents),
            threshold,
            arrow_scale,
        )

        if glyphs is not None:
            plotter.add_mesh(
                glyphs,
                scalars="value",
                cmap="turbo",
                show_scalar_bar=False,
            )

    shunt_data = add_shunt_overlay(
        plotter,
        edge_file,
    )

    if shunt_data is not None:
        print(
            "shunt current = "
            f"{float(shunt_data['shunt_current_A']):.9e} A"
        )

    plotter.add_text(
        (
            f"t = {time_value * 1e6:.3f} us\n"
            f"max |I| = {maximum:.4e} A\n"
            f"edge = {max_edge}"
        ),
        position="upper_left",
        font_size=11,
    )

    plotter.show_axes()

    set_camera(
        plotter,
        edges,
    )

    plotter.show()


# ============================================================
# HARMONIC
# ============================================================

def harmonic_state(
    edge_poly: pv.PolyData,
    phase_deg: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    """
    e^{j omega t} convention:

        I(t) = Re{ I_hat exp(j phase) }

             = I_re cos(phase)
               - I_im sin(phase)

    Same formula for J.
    """

    theta = np.deg2rad(
        phase_deg
    )

    c = np.cos(
        theta
    )

    s = np.sin(
        theta
    )

    I_re = np.asarray(
        edge_poly.cell_data[
            "current_re_A"
        ],
        dtype=float,
    )

    I_im = np.asarray(
        edge_poly.cell_data[
            "current_im_A"
        ],
        dtype=float,
    )

    J_re = np.asarray(
        edge_poly.cell_data[
            "J_re_A_m"
        ],
        dtype=float,
    )

    J_im = np.asarray(
        edge_poly.cell_data[
            "J_im_A_m"
        ],
        dtype=float,
    )

    currents = (
        I_re * c
        -
        I_im * s
    )

    vectors = (
        J_re * c
        -
        J_im * s
    )

    return currents, vectors


def show_harmonic_static(
    edge_file: Path,
    surface_file: Path | None,
    phase_deg: float,
    threshold: float,
    arrow_scale: float,
    show_arrows: bool,
) -> None:
    edges = pv.read(
        edge_file
    )

    currents, vectors = harmonic_state(
        edges,
        phase_deg,
    )

    I_abs_phasor = np.asarray(
        edges.cell_data[
            "current_abs_A"
        ],
        dtype=float,
    )

    maximum = float(
        np.max(
            I_abs_phasor
        )
    )

    instantaneous_max = float(
        np.max(
            np.abs(currents)
        )
    )

    plotter = pv.Plotter(
        title="PEEC harmonic current"
    )

    add_base_surface(
        plotter,
        surface_file,
    )

    edge_view = edges.copy()

    edge_view.cell_data[
        "I_instant_A"
    ] = currents

    plotter.add_mesh(
        edge_view,
        scalars="I_instant_A",
        preference="cell",
        line_width=5,
        cmap="coolwarm",
        clim=[
            -maximum,
            maximum,
        ],
        scalar_bar_args={
            "title": "I_i(phase) [A]",
        },
    )

    if show_arrows:
        glyphs = build_arrow_glyphs(
            edges,
            vectors,
            np.abs(currents),
            threshold,
            arrow_scale,
        )

        if glyphs is not None:
            plotter.add_mesh(
                glyphs,
                scalars="value",
                cmap="turbo",
                show_scalar_bar=False,
            )

    plotter.add_text(
        (
            f"phase = {phase_deg:.1f} deg\n"
            f"max |I_hat| = {maximum:.4e} A\n"
            f"max |I(phase)| = {instantaneous_max:.4e} A"
        ),
        position="upper_left",
        font_size=11,
    )

    plotter.show_axes()

    set_camera(
        plotter,
        edges,
    )

    plotter.show()


def animate_harmonic(
    edge_file: Path,
    surface_file: Path | None,
    frames: int,
    cycles: int,
    frame_delay: float,
    threshold: float,
    arrow_scale: float,
) -> None:
    edges = pv.read(
        edge_file
    )

    maximum = float(
        np.max(
            np.asarray(
                edges.cell_data[
                    "current_abs_A"
                ],
                dtype=float,
            )
        )
    )

    if maximum <= 0.0:
        raise RuntimeError(
            "All harmonic currents are zero."
        )

    currents, vectors = harmonic_state(
        edges,
        0.0,
    )

    edge_view = edges.copy()

    edge_view.cell_data[
        "I_instant_A"
    ] = currents

    plotter = pv.Plotter(
        title="PEEC harmonic current animation"
    )

    add_base_surface(
        plotter,
        surface_file,
    )

    plotter.add_mesh(
        edge_view,
        scalars="I_instant_A",
        preference="cell",
        line_width=5,
        cmap="coolwarm",
        clim=[
            -maximum,
            maximum,
        ],
        scalar_bar_args={
            "title": "I_i(t) [A]",
        },
    )

    arrow_name = "current_arrows"
    text_name = "current_info"

    plotter.add_text(
        "phase = 0.0 deg",
        position="upper_left",
        font_size=11,
        name=text_name,
    )

    plotter.show_axes()

    set_camera(
        plotter,
        edges,
    )

    plotter.show(
        auto_close=False,
        interactive_update=True,
    )

    total_frames = (
        frames
        *
        cycles
    )

    for frame_index in range(
        total_frames
    ):
        phase_deg = (
            360.0
            *
            (frame_index % frames)
            /
            frames
        )

        currents, vectors = harmonic_state(
            edges,
            phase_deg,
        )

        edge_view.cell_data[
            "I_instant_A"
        ] = currents

        edge_view.Modified()

        try:
            plotter.remove_actor(
                arrow_name
            )
        except Exception:
            pass

        glyphs = build_arrow_glyphs(
            edges,
            vectors,
            np.abs(currents),
            threshold,
            arrow_scale,
        )

        if glyphs is not None:
            plotter.add_mesh(
                glyphs,
                scalars="value",
                cmap="turbo",
                show_scalar_bar=False,
                name=arrow_name,
            )

        try:
            plotter.remove_actor(
                text_name
            )
        except Exception:
            pass

        plotter.add_text(
            (
                f"phase = {phase_deg:.1f} deg\n"
                f"max |I(t)| = "
                f"{np.max(np.abs(currents)):.4e} A"
            ),
            position="upper_left",
            font_size=11,
            name=text_name,
        )

        try:
            plotter.render()
            plotter.update()
        except Exception:
            break

        time.sleep(
            frame_delay
        )

    try:
        plotter.interactor.start()
    except Exception:
        pass

    try:
        plotter.close()
    except Exception:
        pass


# ============================================================
# TRANSIENT SERIES
# ============================================================

def extract_frame_number(
    path: Path,
) -> int:
    match = re.search(
        r"_(\d+)\.vtk$",
        path.name,
    )

    if match is None:
        return -1

    return int(
        match.group(1)
    )


def find_transient_files(
    directory: Path,
    prefix: str,
) -> list[Path]:
    files = list(
        directory.glob(
            f"{prefix}_edges_*.vtk"
        )
    )

    files.sort(
        key=extract_frame_number
    )

    return files


def matching_surface_file(
    edge_file: Path,
) -> Path:
    return edge_file.with_name(
        edge_file.name.replace(
            "_edges_",
            "_surface_",
        )
    )


def animate_transient(
    directory: Path,
    prefix: str,
    step: int,
    frame_delay: float,
    threshold: float,
    arrow_scale: float,
) -> None:
    files = find_transient_files(
        directory,
        prefix,
    )

    if not files:
        raise FileNotFoundError(
            f"No files: "
            f"{directory}/{prefix}_edges_*.vtk"
        )

    first_edges = pv.read(
        files[0]
    )

    currents, vectors, time_value = transient_state(
        first_edges
    )

    global_max = 0.0

    for path in files[::step]:
        data = pv.read(
            path
        )

        current, _, _ = transient_state(
            data
        )

        global_max = max(
            global_max,
            float(
                np.max(
                    np.abs(current)
                )
            ),
        )

    if global_max <= 0.0:
        global_max = 1.0

    edge_view = first_edges.copy()

    edge_view.cell_data[
        "I_abs_A"
    ] = np.abs(
        currents
    )

    plotter = pv.Plotter(
        title=f"PEEC transient: {prefix}"
    )

    first_surface = matching_surface_file(
        files[0]
    )

    add_base_surface(
        plotter,
        first_surface,
    )

    plotter.add_mesh(
        edge_view,
        scalars="I_abs_A",
        preference="cell",
        line_width=5,
        cmap="turbo",
        clim=[
            0.0,
            global_max,
        ],
        scalar_bar_args={
            "title": "|I_i| [A]",
        },
    )

    arrow_name = "current_arrows"
    text_name = "current_info"

    plotter.add_text(
        (
            f"t = {time_value * 1e6:.3f} us\n"
            f"max |I| = {np.max(np.abs(currents)):.4e} A"
        ),
        position="upper_left",
        font_size=11,
        name=text_name,
    )

    plotter.show_axes()

    set_camera(
        plotter,
        first_edges,
    )

    plotter.show(
        auto_close=False,
        interactive_update=True,
    )

    for frame_index, path in enumerate(
        files[::step]
    ):
        data = pv.read(
            path
        )

        currents, vectors, time_value = transient_state(
            data
        )

        edge_view.cell_data[
            "I_abs_A"
        ] = np.abs(
            currents
        )

        edge_view.Modified()

        try:
            plotter.remove_actor(
                arrow_name
            )
        except Exception:
            pass

        glyphs = build_arrow_glyphs(
            data,
            vectors,
            np.abs(currents),
            threshold,
            arrow_scale,
        )

        if glyphs is not None:
            plotter.add_mesh(
                glyphs,
                scalars="value",
                cmap="turbo",
                show_scalar_bar=False,
                name=arrow_name,
            )

        try:
            plotter.remove_actor(
                text_name
            )
        except Exception:
            pass

        plotter.add_text(
            (
                f"frame = {frame_index}\n"
                f"t = {time_value * 1e6:.3f} us\n"
                f"max |I| = "
                f"{np.max(np.abs(currents)):.4e} A"
            ),
            position="upper_left",
            font_size=11,
            name=text_name,
        )

        for suffix in (
            "_a",
            "_a_label",
            "_b",
            "_b_label",
            "_branch",
            "_current_label",
        ):
            try:
                plotter.remove_actor(
                    "shunt" + suffix
                )
            except Exception:
                pass

        add_shunt_overlay(
            plotter,
            path,
        )

        try:
            plotter.render()
            plotter.update()
        except Exception:
            break

        time.sleep(
            frame_delay
        )

    try:
        plotter.interactor.start()
    except Exception:
        pass

    try:
        plotter.close()
    except Exception:
        pass


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "PyVista visualization for PEEC VTK output "
            "written by the C solver."
        )
    )

    parser.add_argument(
        "--edges",
        default=None,
        help=(
            "Single *_edges.vtk or *_edges_XXXXXX.vtk file."
        ),
    )

    parser.add_argument(
        "--surface",
        default=None,
        help=(
            "Optional matching surface VTK file."
        ),
    )

    parser.add_argument(
        "--directory",
        default=None,
        help=(
            "Directory containing transient VTK frames."
        ),
    )

    parser.add_argument(
        "--prefix",
        default=None,
        help=(
            "Transient prefix, e.g. lightning_stage1 "
            "or lightning_stage2."
        ),
    )

    parser.add_argument(
        "--phase",
        type=float,
        default=90.0,
        help=(
            "Harmonic physical phase [deg]. "
            "Used with e^(j omega t)."
        ),
    )

    parser.add_argument(
        "--animate",
        action="store_true",
    )

    parser.add_argument(
        "--arrows",
        action="store_true",
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.10,
        help=(
            "Show arrows for vectors >= threshold*max."
        ),
    )

    parser.add_argument(
        "--arrow-scale",
        type=float,
        default=0.08,
        help=(
            "Maximum arrow length as fraction of model size."
        ),
    )

    parser.add_argument(
        "--frames",
        type=int,
        default=60,
        help="Harmonic animation frames per cycle.",
    )

    parser.add_argument(
        "--cycles",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Use every N-th transient VTK frame.",
    )

    parser.add_argument(
        "--frame-delay",
        type=float,
        default=0.03,
    )

    args = parser.parse_args()

    if not (
        0.0 <= args.threshold <= 1.0
    ):
        raise ValueError(
            "--threshold must be between 0 and 1."
        )

    if args.arrow_scale <= 0.0:
        raise ValueError(
            "--arrow-scale must be > 0."
        )

    if args.step < 1:
        raise ValueError(
            "--step must be >= 1."
        )

    # --------------------------------------------------------
    # Transient directory animation.
    # --------------------------------------------------------
    if args.directory is not None:
        if args.prefix is None:
            raise ValueError(
                "--directory requires --prefix."
            )

        animate_transient(
            directory=Path(args.directory),
            prefix=args.prefix,
            step=args.step,
            frame_delay=args.frame_delay,
            threshold=args.threshold,
            arrow_scale=args.arrow_scale,
        )

        return

    # --------------------------------------------------------
    # Single VTK file.
    # --------------------------------------------------------
    if args.edges is None:
        raise ValueError(
            "Use either --edges FILE or "
            "--directory DIR --prefix PREFIX."
        )

    edge_file = Path(
        args.edges
    )

    if not edge_file.exists():
        raise FileNotFoundError(
            edge_file
        )

    if args.surface is not None:
        surface_file = Path(
            args.surface
        )
    else:
        if "_edges_" in edge_file.name:
            candidate = matching_surface_file(
                edge_file
            )
        else:
            candidate = edge_file.with_name(
                edge_file.name.replace(
                    "_edges.vtk",
                    "_surface.vtk",
                )
            )

        surface_file = (
            candidate
            if candidate.exists()
            else None
        )

    edges = pv.read(
        edge_file
    )

    kind = dataset_kind(
        edges
    )

    print()
    print("=" * 72)
    print("PEEC PYVISTA VIEWER")
    print("=" * 72)
    print(f"file   = {edge_file}")
    print(f"type   = {kind}")
    print(
        "arrays = "
        f"{list(edges.cell_data.keys())}"
    )
    print("=" * 72)

    if kind == "transient":
        show_transient_static(
            edge_file=edge_file,
            surface_file=surface_file,
            threshold=args.threshold,
            arrow_scale=args.arrow_scale,
            show_arrows=args.arrows,
        )

        return

    if args.animate:
        animate_harmonic(
            edge_file=edge_file,
            surface_file=surface_file,
            frames=args.frames,
            cycles=args.cycles,
            frame_delay=args.frame_delay,
            threshold=args.threshold,
            arrow_scale=args.arrow_scale,
        )
    else:
        show_harmonic_static(
            edge_file=edge_file,
            surface_file=surface_file,
            phase_deg=args.phase,
            threshold=args.threshold,
            arrow_scale=args.arrow_scale,
            show_arrows=args.arrows,
        )


if __name__ == "__main__":
    main()
