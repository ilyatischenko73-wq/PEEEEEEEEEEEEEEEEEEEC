from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyvista as pv


# ============================================================
# DATA
# ============================================================

@dataclass
class SolverMesh:
    """
    Mesh data in EXACTLY the same local-node numbering used by
    src/mesh/mesh_gmsh.c.

    xyz[i] is the coordinate of PEEC local node i.

    gmsh_tags[i] is the original Gmsh node tag corresponding
    to PEEC local node i.

    quads stores PEEC local node indices.
    """
    xyz: np.ndarray
    gmsh_tags: np.ndarray
    quads: np.ndarray


# ============================================================
# GMSH 4.x ASCII PARSER
# ============================================================

def _find_line(lines: list[str], marker: str) -> int:
    for i, line in enumerate(lines):
        if line.strip() == marker:
            return i

    raise RuntimeError(
        f"Section {marker} not found."
    )


def read_solver_mesh(mesh_file: Path) -> SolverMesh:
    """
    Reproduce mesh_gmsh.c exactly.

    Node numbering in mesh_gmsh.c:

        local_node_index = 0

        for every $Nodes entity block, in file order:
            read all node tags of the block
            read coordinates in the same order

            tag_map[gmsh_tag] = local_node_index
            local_node_index += 1

    Therefore the PEEC node index is NOT assumed to be:
        gmsh_tag - 1
    and is NOT taken from the Gmsh Python API order.

    It is reconstructed directly from the ASCII .msh file.
    """

    if not mesh_file.exists():
        raise FileNotFoundError(
            mesh_file
        )

    lines = mesh_file.read_text(
        encoding="utf-8",
    ).splitlines()

    # --------------------------------------------------------
    # Mesh format
    # --------------------------------------------------------
    mesh_format_pos = _find_line(
        lines,
        "$MeshFormat",
    )

    fmt = lines[
        mesh_format_pos + 1
    ].split()

    if len(fmt) < 3:
        raise RuntimeError(
            "Invalid $MeshFormat line."
        )

    version = float(
        fmt[0]
    )

    file_type = int(
        fmt[1]
    )

    if not (
        4.0 <= version < 5.0
    ):
        raise RuntimeError(
            f"Expected Gmsh 4.x, got {version}."
        )

    if file_type != 0:
        raise RuntimeError(
            "Binary Gmsh is not supported by the current C solver."
        )

    # --------------------------------------------------------
    # Nodes
    # --------------------------------------------------------
    node_pos = _find_line(
        lines,
        "$Nodes",
    )

    header = lines[
        node_pos + 1
    ].split()

    if len(header) != 4:
        raise RuntimeError(
            "Invalid $Nodes header."
        )

    n_blocks = int(
        header[0]
    )

    n_nodes = int(
        header[1]
    )

    cursor = (
        node_pos + 2
    )

    xyz = np.empty(
        (n_nodes, 3),
        dtype=float,
    )

    gmsh_tags = np.empty(
        n_nodes,
        dtype=np.int64,
    )

    tag_to_local: dict[int, int] = {}

    local_node_index = 0

    for _block in range(
        n_blocks
    ):
        block_header = lines[
            cursor
        ].split()

        cursor += 1

        if len(block_header) != 4:
            raise RuntimeError(
                "Invalid node block header."
            )

        entity_dim = int(
            block_header[0]
        )

        entity_tag = int(
            block_header[1]
        )

        parametric = int(
            block_header[2]
        )

        n_block_nodes = int(
            block_header[3]
        )

        del entity_dim, entity_tag

        if parametric != 0:
            raise RuntimeError(
                "Parametric Gmsh nodes are not supported "
                "by mesh_gmsh.c."
            )

        tags: list[int] = []

        for _ in range(
            n_block_nodes
        ):
            values = lines[
                cursor
            ].split()

            cursor += 1

            if not values:
                raise RuntimeError(
                    "Invalid Gmsh node tag."
                )

            tags.append(
                int(values[0])
            )

        for i in range(
            n_block_nodes
        ):
            values = lines[
                cursor
            ].split()

            cursor += 1

            if len(values) < 3:
                raise RuntimeError(
                    "Invalid Gmsh node coordinates."
                )

            if local_node_index >= n_nodes:
                raise RuntimeError(
                    "Number of nodes exceeds $Nodes header."
                )

            x = float(
                values[0]
            )

            y = float(
                values[1]
            )

            z = float(
                values[2]
            )

            gmsh_tag = int(
                tags[i]
            )

            xyz[
                local_node_index
            ] = (
                x,
                y,
                z,
            )

            gmsh_tags[
                local_node_index
            ] = gmsh_tag

            tag_to_local[
                gmsh_tag
            ] = local_node_index

            local_node_index += 1

    if local_node_index != n_nodes:
        raise RuntimeError(
            f"Expected {n_nodes} nodes, "
            f"read {local_node_index}."
        )

    if lines[
        cursor
    ].strip() != "$EndNodes":
        raise RuntimeError(
            "Expected $EndNodes."
        )

    # --------------------------------------------------------
    # Elements
    # --------------------------------------------------------
    elem_pos = _find_line(
        lines,
        "$Elements",
    )

    header = lines[
        elem_pos + 1
    ].split()

    if len(header) != 4:
        raise RuntimeError(
            "Invalid $Elements header."
        )

    n_element_blocks = int(
        header[0]
    )

    cursor = (
        elem_pos + 2
    )

    quads: list[
        tuple[int, int, int, int]
    ] = []

    for _block in range(
        n_element_blocks
    ):
        block_header = lines[
            cursor
        ].split()

        cursor += 1

        if len(block_header) != 4:
            raise RuntimeError(
                "Invalid element block header."
            )

        entity_dim = int(
            block_header[0]
        )

        entity_tag = int(
            block_header[1]
        )

        element_type = int(
            block_header[2]
        )

        n_block_elements = int(
            block_header[3]
        )

        del entity_dim, entity_tag

        for _ in range(
            n_block_elements
        ):
            values = lines[
                cursor
            ].split()

            cursor += 1

            if element_type != 3:
                continue

            if len(values) != 5:
                raise RuntimeError(
                    "Invalid Gmsh quadrangle."
                )

            node_tags = [
                int(values[1]),
                int(values[2]),
                int(values[3]),
                int(values[4]),
            ]

            try:
                quad = tuple(
                    tag_to_local[tag]
                    for tag in node_tags
                )
            except KeyError as exc:
                raise RuntimeError(
                    f"Unknown Gmsh node tag {exc.args[0]}."
                ) from exc

            quads.append(
                quad
            )

    if not quads:
        raise RuntimeError(
            "Mesh has no 4-node quadrangles."
        )

    return SolverMesh(
        xyz=xyz,
        gmsh_tags=gmsh_tags,
        quads=np.asarray(
            quads,
            dtype=np.int64,
        ),
    )


# ============================================================
# PYVISTA SURFACE
# ============================================================

def build_surface(
    mesh: SolverMesh,
) -> pv.PolyData:
    faces = np.empty(
        (len(mesh.quads), 5),
        dtype=np.int64,
    )

    faces[:, 0] = 4
    faces[:, 1:] = mesh.quads

    surface = pv.PolyData(
        mesh.xyz,
        faces.ravel(),
    )

    surface.point_data[
        "peec_solver_index"
    ] = np.arange(
        len(mesh.xyz),
        dtype=np.int64,
    )

    surface.point_data[
        "gmsh_node_tag"
    ] = mesh.gmsh_tags

    return surface


# ============================================================
# VISUAL HELPERS
# ============================================================

def set_camera(
    plotter: pv.Plotter,
    dataset: pv.DataSet,
) -> None:
    bounds = np.asarray(
        dataset.bounds,
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

    total = max(
        float(np.max(extent)),
        1.0,
    )

    if extent[2] < 1.0e-10 * total:
        plotter.view_xy()
    else:
        plotter.view_isometric()

    plotter.reset_camera()


# ============================================================
# NODE PICKER
# ============================================================

def run_node_picker(
    mesh_file: Path,
    point_size: float,
    surface_opacity: float,
) -> None:
    mesh = read_solver_mesh(
        mesh_file
    )

    surface = build_surface(
        mesh
    )

    selected: dict[
        str,
        tuple[int, int, np.ndarray] | None
    ] = {
        "strike": None,
        "return": None,
        "shunt_a": None,
        "shunt_b": None,
    }

    labels = {
        "strike": "STRIKE",
        "return": "RETURN",
        "shunt_a": "SHUNT A",
        "shunt_b": "SHUNT B",
    }

    colors = {
        "strike": "red",
        "return": "blue",
        "shunt_a": "green",
        "shunt_b": "orange",
    }

    actor_names = {
        key: f"selected_{key}"
        for key in selected
    }

    pick_mode = "strike"

    plotter = pv.Plotter(
        title=(
            "PEEC exact C-node picker: "
            f"{mesh_file.name}"
        )
    )

    plotter.add_mesh(
        surface,
        color="lightgray",
        opacity=surface_opacity,
        show_edges=True,
        edge_color="gray",
        pickable=True,
    )

    node_cloud = pv.PolyData(
        mesh.xyz
    )

    node_cloud.point_data[
        "peec_solver_index"
    ] = np.arange(
        len(mesh.xyz),
        dtype=np.int64,
    )

    plotter.add_mesh(
        node_cloud,
        color="white",
        point_size=point_size,
        render_points_as_spheres=True,
        pickable=False,
    )

    status_name = "picker_status"

    def nearest_node(
        picked_point,
    ) -> tuple[int, int, np.ndarray]:
        point = np.asarray(
            picked_point,
            dtype=float,
        )

        delta = (
            mesh.xyz
            -
            point[None, :]
        )

        d2 = np.einsum(
            "ij,ij->i",
            delta,
            delta,
        )

        solver_index = int(
            np.argmin(
                d2
            )
        )

        gmsh_tag = int(
            mesh.gmsh_tags[
                solver_index
            ]
        )

        xyz = mesh.xyz[
            solver_index
        ].copy()

        return (
            solver_index,
            gmsh_tag,
            xyz,
        )

    def status_text() -> str:
        lines = [
            "PEEC EXACT C-NODE PICKER",
            "",
            f"mode: {labels[pick_mode]}",
            "",
            "S : strike",
            "R : return",
            "A : shunt A",
            "B : shunt B",
            "P : print CLI",
            "C : clear",
            "",
            "Left click: nearest PEEC node",
            "",
        ]

        for key in (
            "strike",
            "return",
            "shunt_a",
            "shunt_b",
        ):
            value = selected[
                key
            ]

            if value is None:
                lines.append(
                    f"{labels[key]:8s}: ---"
                )
                continue

            solver_index, gmsh_tag, xyz = value

            lines.append(
                (
                    f"{labels[key]:8s}: "
                    f"PEEC={solver_index}, "
                    f"Gmsh={gmsh_tag}"
                )
            )

            lines.append(
                (
                    " " * 10
                    +
                    f"xyz=({xyz[0]:.6g}, "
                    f"{xyz[1]:.6g}, "
                    f"{xyz[2]:.6g})"
                )
            )

        return "\n".join(
            lines
        )

    def update_status() -> None:
        try:
            plotter.remove_actor(
                status_name
            )
        except Exception:
            pass

        plotter.add_text(
            status_text(),
            position="upper_left",
            font_size=10,
            name=status_name,
        )

    def redraw_marker(
        key: str,
    ) -> None:
        actor = actor_names[
            key
        ]

        for name in (
            actor,
            actor + "_label",
        ):
            try:
                plotter.remove_actor(
                    name
                )
            except Exception:
                pass

        value = selected[
            key
        ]

        if value is None:
            return

        solver_index, gmsh_tag, xyz = value

        marker = pv.PolyData(
            np.asarray(
                [xyz],
                dtype=float,
            )
        )

        plotter.add_mesh(
            marker,
            color=colors[key],
            point_size=22,
            render_points_as_spheres=True,
            name=actor,
        )

        plotter.add_point_labels(
            np.asarray(
                [xyz],
                dtype=float,
            ),
            [
                (
                    f"{labels[key]}\n"
                    f"PEEC {solver_index}\n"
                    f"Gmsh {gmsh_tag}"
                )
            ],
            point_size=0,
            font_size=10,
            always_visible=True,
            shape_color=colors[key],
            text_color="black",
            name=actor + "_label",
        )

    def set_mode(
        key: str,
    ) -> None:
        nonlocal pick_mode

        pick_mode = key

        print()
        print(
            f"Pick mode: {labels[key]}"
        )

        update_status()
        plotter.render()

    def on_pick(
        picked_point,
    ) -> None:
        if picked_point is None:
            return

        solver_index, gmsh_tag, xyz = nearest_node(
            picked_point
        )

        selected[
            pick_mode
        ] = (
            solver_index,
            gmsh_tag,
            xyz,
        )

        print()
        print("=" * 72)
        print(
            f"{labels[pick_mode]} SELECTED"
        )
        print("=" * 72)
        print(
            f"PEEC local index = {solver_index}"
        )
        print(
            f"Gmsh node tag    = {gmsh_tag}"
        )
        print(
            "xyz              = "
            f"{xyz[0]:.12e} "
            f"{xyz[1]:.12e} "
            f"{xyz[2]:.12e}"
        )
        print("=" * 72)

        redraw_marker(
            pick_mode
        )

        update_status()
        plotter.render()

    def print_selected() -> None:
        print()
        print("=" * 72)
        print("SELECTED PEEC NODES")
        print("=" * 72)

        for key in (
            "strike",
            "return",
            "shunt_a",
            "shunt_b",
        ):
            value = selected[
                key
            ]

            if value is None:
                print(
                    f"{labels[key]:8s}: not selected"
                )
                continue

            solver_index, gmsh_tag, xyz = value

            print(
                (
                    f"{labels[key]:8s}: "
                    f"PEEC={solver_index}, "
                    f"Gmsh={gmsh_tag}, "
                    f"xyz="
                    f"{xyz[0]:.12e} "
                    f"{xyz[1]:.12e} "
                    f"{xyz[2]:.12e}"
                )
            )

        print()
        print(
            "Ready-to-copy PeecSolver INDEX CLI:"
        )
        print()

        if selected[
            "strike"
        ] is not None:
            solver_index, _, _ = selected[
                "strike"
            ]
            print(
                f"--strike-node {solver_index} \\"
            )

        if selected[
            "return"
        ] is not None:
            solver_index, _, _ = selected[
                "return"
            ]
            print(
                f"--return-node {solver_index} \\"
            )

        if selected[
            "shunt_a"
        ] is not None:
            solver_index, _, _ = selected[
                "shunt_a"
            ]
            print(
                f"--shunt-a-node {solver_index} \\"
            )

        if selected[
            "shunt_b"
        ] is not None:
            solver_index, _, _ = selected[
                "shunt_b"
            ]
            print(
                f"--shunt-b-node {solver_index}"
            )

        print("=" * 72)

    def clear_all() -> None:
        for key in selected:
            selected[
                key
            ] = None

            for name in (
                actor_names[key],
                actor_names[key] + "_label",
            ):
                try:
                    plotter.remove_actor(
                        name
                    )
                except Exception:
                    pass

        update_status()
        plotter.render()

    plotter.add_key_event(
        "s",
        lambda: set_mode(
            "strike"
        ),
    )

    plotter.add_key_event(
        "r",
        lambda: set_mode(
            "return"
        ),
    )

    plotter.add_key_event(
        "a",
        lambda: set_mode(
            "shunt_a"
        ),
    )

    plotter.add_key_event(
        "b",
        lambda: set_mode(
            "shunt_b"
        ),
    )

    plotter.add_key_event(
        "p",
        print_selected,
    )

    plotter.add_key_event(
        "c",
        clear_all,
    )

    plotter.enable_surface_point_picking(
        callback=on_pick,
        show_message=False,
        show_point=False,
        left_clicking=True,
    )

    update_status()
    plotter.show_axes()

    set_camera(
        plotter,
        surface,
    )

    print()
    print("=" * 72)
    print("PEEC EXACT C-NODE PICKER")
    print("=" * 72)
    print(
        f"Mesh              : {mesh_file}"
    )
    print(
        f"PEEC nodes        : {len(mesh.xyz)}"
    )
    print(
        f"Quad cells        : {len(mesh.quads)}"
    )
    print()
    print(
        "Node indices are reconstructed directly from $Nodes"
    )
    print(
        "using the same order as src/mesh/mesh_gmsh.c."
    )
    print()
    print("S : select STRIKE")
    print("R : select RETURN")
    print("A : select SHUNT A")
    print("B : select SHUNT B")
    print("P : print ready-to-copy PEEC indices")
    print("C : clear")
    print("Left click : select nearest PEEC node")
    print("=" * 72)

    plotter.show()


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Exact PEEC C-solver node picker. "
            "Local node indices are reconstructed directly "
            "from Gmsh 4.x ASCII $Nodes exactly like mesh_gmsh.c."
        )
    )

    parser.add_argument(
        "--mesh",
        required=True,
        help="Gmsh 4.x ASCII .msh file.",
    )

    parser.add_argument(
        "--point-size",
        type=float,
        default=7.0,
    )

    parser.add_argument(
        "--surface-opacity",
        type=float,
        default=0.28,
    )

    args = parser.parse_args()

    if args.point_size <= 0.0:
        raise ValueError(
            "--point-size must be > 0."
        )

    if not (
        0.0
        <= args.surface_opacity
        <= 1.0
    ):
        raise ValueError(
            "--surface-opacity must be between 0 and 1."
        )

    run_node_picker(
        mesh_file=Path(
            args.mesh
        ),
        point_size=args.point_size,
        surface_opacity=args.surface_opacity,
    )


if __name__ == "__main__":
    main()
