from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

try:
    import gmsh
except ImportError as exc:
    raise SystemExit(
        "Python module 'gmsh' is required.\n"
        "Install it with:\n"
        "    python3 -m pip install gmsh pyvista numpy\n"
    ) from exc


# ============================================================
# SOLVER NODE ORDER
# ============================================================

def read_solver_node_order(
    mesh_file: Path,
) -> tuple[np.ndarray, dict[int, int]]:
    """
    Reproduce exactly the node ordering used by the C solver.

    mesh_gmsh.c assigns local indices in the order in which node
    tags appear inside the $Nodes entity blocks:

        local_index = 0, 1, 2, ...

    Gmsh's Python API is free to return nodes in a different order,
    therefore np.arange(len(gmsh.model.mesh.getNodes())) must NOT be
    interpreted as a PEEC solver node index.
    """

    lines = mesh_file.read_text(
        encoding="utf-8",
    ).splitlines()

    try:
        start = lines.index("$Nodes")
    except ValueError as exc:
        raise RuntimeError(
            "No $Nodes section in Gmsh file."
        ) from exc

    header = lines[
        start + 1
    ].split()

    if len(header) != 4:
        raise RuntimeError(
            "Invalid Gmsh $Nodes header."
        )

    n_blocks = int(
        header[0]
    )

    n_nodes = int(
        header[1]
    )

    cursor = start + 2

    solver_tags: list[int] = []

    for _ in range(
        n_blocks
    ):
        block = lines[
            cursor
        ].split()

        cursor += 1

        if len(block) != 4:
            raise RuntimeError(
                "Invalid Gmsh node block header."
            )

        parametric = int(
            block[2]
        )

        n_block_nodes = int(
            block[3]
        )

        if parametric != 0:
            raise RuntimeError(
                "Parametric Gmsh nodes are not supported "
                "by the current C solver."
            )

        block_tags: list[int] = []

        for _ in range(
            n_block_nodes
        ):
            block_tags.append(
                int(
                    lines[cursor].strip()
                )
            )

            cursor += 1

        solver_tags.extend(
            block_tags
        )

        # Skip coordinates. The C solver reads them in exactly
        # the same block/tag order.
        cursor += n_block_nodes

    if len(solver_tags) != n_nodes:
        raise RuntimeError(
            f"Expected {n_nodes} solver nodes, "
            f"read {len(solver_tags)}."
        )

    tag_to_solver_index = {
        tag: index
        for index, tag in enumerate(
            solver_tags
        )
    }

    return (
        np.asarray(
            solver_tags,
            dtype=np.int64,
        ),
        tag_to_solver_index,
    )


# ============================================================
# GMSH LOADING
# ============================================================

def load_gmsh_surface(
    mesh_file: Path,
) -> tuple[
    pv.PolyData,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """
    Read a Gmsh .msh file with the Gmsh Python API.

    Returns:
        surface   : PyVista PolyData containing 2D elements
        node_tags : Gmsh node tags, shape (N,)
        points    : node coordinates, shape (N, 3)

    The point order in `surface.points`, `node_tags`, and `points`
    is identical.
    """

    if not mesh_file.exists():
        raise FileNotFoundError(mesh_file)

    _, tag_to_solver_index = read_solver_node_order(
        mesh_file
    )

    gmsh.initialize()

    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(mesh_file))

        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()

        node_tags = np.asarray(
            node_tags,
            dtype=np.int64,
        )

        points = np.asarray(
            node_coords,
            dtype=float,
        ).reshape(-1, 3)

        if len(node_tags) == 0:
            raise RuntimeError(
                "The Gmsh mesh contains no nodes."
            )

        solver_indices = np.asarray(
            [
                tag_to_solver_index[
                    int(tag)
                ]
                for tag in node_tags
            ],
            dtype=np.int64,
        )

        tag_to_index = {
            int(tag): i
            for i, tag in enumerate(node_tags)
        }

        faces: list[int] = []

        element_types, element_tags, element_node_tags = (
            gmsh.model.mesh.getElements(dim=2)
        )

        del element_tags

        for element_type, connectivity in zip(
            element_types,
            element_node_tags,
        ):
            name, dim, order, num_nodes, local_coords, num_primary = (
                gmsh.model.mesh.getElementProperties(
                    int(element_type)
                )
            )

            del name, dim, order, local_coords

            connectivity = np.asarray(
                connectivity,
                dtype=np.int64,
            ).reshape(-1, int(num_nodes))

            # For higher-order elements only the primary corner nodes
            # are used to draw the surface.
            n_corner = int(num_primary)

            if n_corner not in (3, 4):
                continue

            for elem in connectivity:
                corner_tags = elem[:n_corner]

                try:
                    corner_indices = [
                        tag_to_index[int(tag)]
                        for tag in corner_tags
                    ]
                except KeyError as exc:
                    raise RuntimeError(
                        "Element references an unknown Gmsh node tag."
                    ) from exc

                faces.append(n_corner)
                faces.extend(corner_indices)

        if not faces:
            raise RuntimeError(
                "No supported 2D triangle/quad elements were found "
                "in the Gmsh mesh."
            )

        surface = pv.PolyData(
            points,
            np.asarray(
                faces,
                dtype=np.int64,
            ),
        )

        # Useful metadata for optional inspection in PyVista.
        surface.point_data[
            "gmsh_node_tag"
        ] = node_tags

        surface.point_data[
            "peec_solver_index"
        ] = solver_indices

        return (
            surface,
            node_tags,
            solver_indices,
            points,
        )

    finally:
        gmsh.finalize()


# ============================================================
# VISUAL HELPERS
# ============================================================

def model_scale(
    dataset: pv.DataSet,
) -> float:
    bounds = np.asarray(
        dataset.bounds,
        dtype=float,
    )

    scale = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )

    return float(scale)


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
# NODE PICKER
# ============================================================

def run_node_picker(
    mesh_file: Path,
    point_size: float,
    surface_opacity: float,
) -> None:
    surface, node_tags, solver_indices, points = load_gmsh_surface(
        mesh_file
    )

    selected: dict[str, tuple[int, int, int, np.ndarray] | None] = {
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
        title=f"PEEC Gmsh node picker: {mesh_file.name}"
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
        points
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
    ) -> tuple[int, int, int, np.ndarray]:
        p = np.asarray(
            picked_point,
            dtype=float,
        )

        delta = (
            points
            -
            p[None, :]
        )

        distance2 = np.einsum(
            "ij,ij->i",
            delta,
            delta,
        )

        local_index = int(
            np.argmin(
                distance2
            )
        )

        gmsh_tag = int(
            node_tags[local_index]
        )

        solver_index = int(
            solver_indices[local_index]
        )

        xyz = points[
            local_index
        ].copy()

        return (
            local_index,
            gmsh_tag,
            solver_index,
            xyz,
        )

    def status_text() -> str:
        lines = [
            "PEEC GMSH NODE PICKER",
            "",
            f"mode: {labels[pick_mode]}",
            "",
            "S : strike",
            "R : return",
            "A : shunt A",
            "B : shunt B",
            "P : print CLI values",
            "C : clear",
            "",
            "Left click: nearest Gmsh node",
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
            else:
                api_index, gmsh_tag, solver_index, xyz = value

                lines.append(
                    (
                        f"{labels[key]:8s}: "
                        f"PEEC idx={solver_index}, "
                        f"Gmsh tag={gmsh_tag}"
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

        try:
            plotter.remove_actor(
                actor
            )
        except Exception:
            pass

        try:
            plotter.remove_actor(
                actor + "_label"
            )
        except Exception:
            pass

        value = selected[
            key
        ]

        if value is None:
            return

        api_index, gmsh_tag, solver_index, xyz = value

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
                    f"PEEC idx {solver_index}\n"
                    f"Gmsh tag {gmsh_tag}"
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

        api_index, gmsh_tag, solver_index, xyz = nearest_node(
            picked_point
        )

        selected[
            pick_mode
        ] = (
            api_index,
            gmsh_tag,
            solver_index,
            xyz,
        )

        print()
        print("=" * 72)
        print(
            f"{labels[pick_mode]} SELECTED"
        )
        print("=" * 72)

        print(
            f"Gmsh node tag = {gmsh_tag}"
        )

        print(
            f"PEEC index    = {solver_index}"
        )

        print(
            f"Gmsh API idx  = {api_index}"
        )

        print(
            (
                "xyz           = "
                f"{xyz[0]:.12e} "
                f"{xyz[1]:.12e} "
                f"{xyz[2]:.12e}"
            )
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
        print("SELECTED PEEC GMSH NODES")
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

            api_index, gmsh_tag, solver_index, xyz = value

            print(
                (
                    f"{labels[key]:8s}: "
                    f"peec_index={solver_index}, "
                    f"gmsh_tag={gmsh_tag}, "
                    f"gmsh_api_index={api_index}, "
                    f"xyz="
                    f"{xyz[0]:.12e} "
                    f"{xyz[1]:.12e} "
                    f"{xyz[2]:.12e}"
                )
            )

        print()
        print("Ready-to-copy PeecSolver coordinate CLI:")
        print()

        if selected["strike"] is not None:
            _, _, _, xyz = selected["strike"]

            print(
                (
                    "--strike "
                    f"{xyz[0]:.12g} "
                    f"{xyz[1]:.12g} "
                    f"{xyz[2]:.12g} \\"
                )
            )

        if selected["return"] is not None:
            _, _, _, xyz = selected["return"]

            print(
                (
                    "--return "
                    f"{xyz[0]:.12g} "
                    f"{xyz[1]:.12g} "
                    f"{xyz[2]:.12g} \\"
                )
            )

        if selected["shunt_a"] is not None:
            _, _, _, xyz = selected["shunt_a"]

            print(
                (
                    "--shunt-a "
                    f"{xyz[0]:.12g} "
                    f"{xyz[1]:.12g} "
                    f"{xyz[2]:.12g} \\"
                )
            )

        if selected["shunt_b"] is not None:
            _, _, _, xyz = selected["shunt_b"]

            print(
                (
                    "--shunt-b "
                    f"{xyz[0]:.12g} "
                    f"{xyz[1]:.12g} "
                    f"{xyz[2]:.12g}"
                )
            )

        print()
        print("Ready-to-copy PeecSolver INDEX CLI:")
        print()

        if selected["strike"] is not None:
            _, _, solver_index, _ = selected["strike"]
            print(f"--strike-node {solver_index} \\")

        if selected["return"] is not None:
            _, _, solver_index, _ = selected["return"]
            print(f"--return-node {solver_index} \\")

        if selected["shunt_a"] is not None:
            _, _, solver_index, _ = selected["shunt_a"]
            print(f"--shunt-a-node {solver_index} \\")

        if selected["shunt_b"] is not None:
            _, _, solver_index, _ = selected["shunt_b"]
            print(f"--shunt-b-node {solver_index}")

        print()
        print("PEEC / Gmsh node identifiers:")
        print()

        for key in (
            "strike",
            "return",
            "shunt_a",
            "shunt_b",
        ):
            value = selected[
                key
            ]

            if value is not None:
                _, gmsh_tag, solver_index, _ = value

                print(
                    (
                        f"{labels[key]:8s}: "
                        f"PEEC={solver_index}, Gmsh={gmsh_tag}"
                    )
                )

        print("=" * 72)

    def clear_all() -> None:
        for key in selected:
            selected[
                key
            ] = None

            try:
                plotter.remove_actor(
                    actor_names[key]
                )
            except Exception:
                pass

            try:
                plotter.remove_actor(
                    actor_names[key]
                    + "_label"
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
    print("PEEC GMSH NODE PICKER")
    print("=" * 72)
    print(f"Mesh: {mesh_file}")
    print(f"Nodes: {len(points)}")
    print(f"Surface cells: {surface.n_cells}")
    print()
    print("S : select STRIKE")
    print("R : select RETURN")
    print("A : select SHUNT A")
    print("B : select SHUNT B")
    print("P : print PEEC indices / coordinates / Gmsh tags")
    print("C : clear")
    print("Left click : select nearest Gmsh node")
    print("=" * 72)

    plotter.show()


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Interactive Gmsh mesh node picker for PEEC. "
            "Select lightning strike/return nodes and shunt terminals."
        )
    )

    parser.add_argument(
        "--mesh",
        required=True,
        help="Gmsh .msh file.",
    )

    parser.add_argument(
        "--point-size",
        type=float,
        default=7.0,
        help="Displayed mesh-node point size.",
    )

    parser.add_argument(
        "--surface-opacity",
        type=float,
        default=0.28,
        help="Surface opacity in [0,1].",
    )

    args = parser.parse_args()

    if args.point_size <= 0.0:
        raise ValueError(
            "--point-size must be > 0."
        )

    if not (
        0.0 <= args.surface_opacity <= 1.0
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
