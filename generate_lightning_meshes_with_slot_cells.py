from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import gmsh


# ============================================================
# DATA
# ============================================================

@dataclass(frozen=True)
class Geometry:
    length: float
    width: float
    height: float

    opening_w: float
    opening_h: float

    @property
    def x0(self) -> float:
        return 0.0

    @property
    def x1(self) -> float:
        return self.length

    @property
    def y0(self) -> float:
        return -0.5 * self.width

    @property
    def y1(self) -> float:
        return +0.5 * self.width

    @property
    def z0(self) -> float:
        return 0.0

    @property
    def z1(self) -> float:
        return self.height

    @property
    def open_x0(self) -> float:
        return 0.5 * (self.length - self.opening_w)

    @property
    def open_x1(self) -> float:
        return 0.5 * (self.length + self.opening_w)

    @property
    def open_z0(self) -> float:
        return 0.5 * (self.height - self.opening_h)

    @property
    def open_z1(self) -> float:
        return 0.5 * (self.height + self.opening_h)


@dataclass
class BuiltGeometry:
    body_surfaces: list[int]
    front_surface: int
    opening_curves: list[int]
    cover_surface: int | None


@dataclass
class ParsedMesh:
    node_xyz: dict[int, tuple[float, float, float]]

    # Gmsh node tags in exactly the order in which the C reader
    # stores them in mesh.xyz.
    #
    # IMPORTANT:
    #     internal PEEC node index != gmsh_tag - 1 in general.
    #
    # The C reader assigns:
    #
    #     internal_index = position in $Nodes file order.
    node_order: list[int]
    tag_to_internal: dict[int, int]

    # Quads in exactly the order in which they occur in $Elements.
    # Each item:
    #
    #     (surface_entity_tag, [n0,n1,n2,n3])
    #
    quads: list[tuple[int, tuple[int, int, int, int]]]


# ============================================================
# BASIC GMSH GEOMETRY
# ============================================================

def add_point(
    x: float,
    y: float,
    z: float,
    h: float,
) -> int:
    return gmsh.model.geo.addPoint(
        x,
        y,
        z,
        h,
    )


def add_quad_surface(
    p0: int,
    p1: int,
    p2: int,
    p3: int,
) -> int:
    l0 = gmsh.model.geo.addLine(p0, p1)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p0)

    loop = gmsh.model.geo.addCurveLoop(
        [l0, l1, l2, l3]
    )

    return gmsh.model.geo.addPlaneSurface(
        [loop]
    )


def build_body_geometry(
    geometry: Geometry,
    h: float,
    with_cover: bool,
) -> BuiltGeometry:
    """
    Build the SAME enclosure geometry in two variants.

    OPEN:
        front wall contains the aperture.

    CLOSED:
        same front wall + a temporary PEC surface TEMP_COVER
        filling exactly the aperture.

    The closed and open models are therefore separate meshes,
    but the aperture boundary has the same physical coordinates.
    """

    g = geometry

    # --------------------------------------------------------
    # Main box corner points
    # --------------------------------------------------------
    p000 = add_point(g.x0, g.y0, g.z0, h)
    p100 = add_point(g.x1, g.y0, g.z0, h)
    p110 = add_point(g.x1, g.y1, g.z0, h)
    p010 = add_point(g.x0, g.y1, g.z0, h)

    p001 = add_point(g.x0, g.y0, g.z1, h)
    p101 = add_point(g.x1, g.y0, g.z1, h)
    p111 = add_point(g.x1, g.y1, g.z1, h)
    p011 = add_point(g.x0, g.y1, g.z1, h)

    body_surfaces: list[int] = []

    # --------------------------------------------------------
    # z = z0
    # --------------------------------------------------------
    body_surfaces.append(
        add_quad_surface(
            p000,
            p100,
            p110,
            p010,
        )
    )

    # --------------------------------------------------------
    # z = z1
    # --------------------------------------------------------
    body_surfaces.append(
        add_quad_surface(
            p001,
            p011,
            p111,
            p101,
        )
    )

    # --------------------------------------------------------
    # y = y0
    # --------------------------------------------------------
    body_surfaces.append(
        add_quad_surface(
            p000,
            p001,
            p101,
            p100,
        )
    )

    # --------------------------------------------------------
    # x = x0
    # --------------------------------------------------------
    body_surfaces.append(
        add_quad_surface(
            p000,
            p010,
            p011,
            p001,
        )
    )

    # --------------------------------------------------------
    # x = x1
    # --------------------------------------------------------
    body_surfaces.append(
        add_quad_surface(
            p100,
            p101,
            p111,
            p110,
        )
    )

    # --------------------------------------------------------
    # Front wall y = y1
    #
    # Outer loop.
    # --------------------------------------------------------
    l_outer_bottom = gmsh.model.geo.addLine(
        p010,
        p110,
    )

    l_outer_right = gmsh.model.geo.addLine(
        p110,
        p111,
    )

    l_outer_top = gmsh.model.geo.addLine(
        p111,
        p011,
    )

    l_outer_left = gmsh.model.geo.addLine(
        p011,
        p010,
    )

    outer_loop = gmsh.model.geo.addCurveLoop(
        [
            l_outer_bottom,
            l_outer_right,
            l_outer_top,
            l_outer_left,
        ]
    )

    # --------------------------------------------------------
    # Aperture rectangle.
    # --------------------------------------------------------
    po0 = add_point(
        g.open_x0,
        g.y1,
        g.open_z0,
        h,
    )

    po1 = add_point(
        g.open_x1,
        g.y1,
        g.open_z0,
        h,
    )

    po2 = add_point(
        g.open_x1,
        g.y1,
        g.open_z1,
        h,
    )

    po3 = add_point(
        g.open_x0,
        g.y1,
        g.open_z1,
        h,
    )

    li0 = gmsh.model.geo.addLine(po0, po1)
    li1 = gmsh.model.geo.addLine(po1, po2)
    li2 = gmsh.model.geo.addLine(po2, po3)
    li3 = gmsh.model.geo.addLine(po3, po0)

    opening_curves = [
        li0,
        li1,
        li2,
        li3,
    ]

    # Hole orientation opposite to outer loop.
    inner_loop = gmsh.model.geo.addCurveLoop(
        [
            -li3,
            -li2,
            -li1,
            -li0,
        ]
    )

    front_surface = gmsh.model.geo.addPlaneSurface(
        [
            outer_loop,
            inner_loop,
        ]
    )

    body_surfaces.append(
        front_surface
    )

    # --------------------------------------------------------
    # Temporary PEC cover.
    #
    # IMPORTANT:
    # it uses exactly the SAME aperture curves. Therefore the
    # closed mesh is conforming at the wall/cover interface.
    # --------------------------------------------------------
    cover_surface: int | None = None

    if with_cover:
        cover_loop = gmsh.model.geo.addCurveLoop(
            [
                li0,
                li1,
                li2,
                li3,
            ]
        )

        cover_surface = gmsh.model.geo.addPlaneSurface(
            [cover_loop]
        )

    gmsh.model.geo.synchronize()

    # --------------------------------------------------------
    # Physical groups.
    #
    # They are useful for inspection in Gmsh even though the
    # current PEEC C reader does not need them for the solve.
    # --------------------------------------------------------
    body_group = gmsh.model.addPhysicalGroup(
        2,
        body_surfaces,
    )

    gmsh.model.setPhysicalName(
        2,
        body_group,
        "BODY",
    )

    aperture_boundary_group = gmsh.model.addPhysicalGroup(
        1,
        opening_curves,
    )

    gmsh.model.setPhysicalName(
        1,
        aperture_boundary_group,
        "APERTURE_BOUNDARY",
    )

    if cover_surface is not None:
        cover_group = gmsh.model.addPhysicalGroup(
            2,
            [cover_surface],
        )

        gmsh.model.setPhysicalName(
            2,
            cover_group,
            "TEMP_COVER",
        )

    return BuiltGeometry(
        body_surfaces=body_surfaces,
        front_surface=front_surface,
        opening_curves=opening_curves,
        cover_surface=cover_surface,
    )


# ============================================================
# MESH SETTINGS
# ============================================================

def configure_mesh(
    built: BuiltGeometry,
    h: float,
) -> None:
    gmsh.option.setNumber(
        "Mesh.MeshSizeMin",
        h,
    )

    gmsh.option.setNumber(
        "Mesh.MeshSizeMax",
        h,
    )

    gmsh.option.setNumber(
        "Mesh.Algorithm",
        8,
    )

    gmsh.option.setNumber(
        "Mesh.RecombinationAlgorithm",
        1,
    )

    gmsh.option.setNumber(
        "Mesh.RecombineAll",
        1,
    )

    surfaces = list(
        built.body_surfaces
    )

    if built.cover_surface is not None:
        surfaces.append(
            built.cover_surface
        )

    for surface in surfaces:
        gmsh.model.mesh.setRecombine(
            2,
            surface,
        )


# ============================================================
# GENERATE ONE FILE
# ============================================================

def generate_one_mesh(
    output: Path,
    geometry: Geometry,
    h: float,
    with_cover: bool,
    model_name: str,
) -> int | None:
    gmsh.clear()
    gmsh.model.add(
        model_name
    )

    built = build_body_geometry(
        geometry,
        h,
        with_cover,
    )

    configure_mesh(
        built,
        h,
    )

    gmsh.model.mesh.generate(
        2
    )

    gmsh.model.mesh.removeDuplicateNodes()

    gmsh.option.setNumber(
        "Mesh.MshFileVersion",
        4.1,
    )

    gmsh.option.setNumber(
        "Mesh.Binary",
        0,
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    gmsh.write(
        str(output)
    )

    return built.cover_surface


# ============================================================
# SMALL GMSH 4.1 ASCII PARSER
# ============================================================
#
# We parse the generated files again instead of relying on
# in-memory Gmsh ordering.
#
# This gives us exactly the node/quad ordering that the C solver
# will see in the .msh file.
# ============================================================

def _find_section(
    lines: list[str],
    begin: str,
    end: str,
) -> tuple[int, int]:
    try:
        i0 = lines.index(begin)
        i1 = lines.index(end)
    except ValueError as exc:
        raise RuntimeError(
            f"Missing Gmsh section {begin}"
        ) from exc

    return i0 + 1, i1


def parse_gmsh_41_ascii(
    path: Path,
) -> ParsedMesh:
    lines = [
        line.strip()
        for line in path.read_text(
            encoding="utf-8"
        ).splitlines()
    ]

    # --------------------------------------------------------
    # Nodes
    # --------------------------------------------------------
    n0, n1 = _find_section(
        lines,
        "$Nodes",
        "$EndNodes",
    )

    header = lines[n0].split()

    if len(header) != 4:
        raise RuntimeError(
            "Unsupported $Nodes header."
        )

    n_blocks = int(
        header[0]
    )

    cursor = n0 + 1

    node_xyz: dict[
        int,
        tuple[float, float, float],
    ] = {}

    node_order: list[int] = []

    for _ in range(n_blocks):
        block_header = lines[cursor].split()
        cursor += 1

        entity_dim = int(
            block_header[0]
        )

        _entity_tag = int(
            block_header[1]
        )

        parametric = int(
            block_header[2]
        )

        count = int(
            block_header[3]
        )

        tags: list[int] = []

        for _ in range(count):
            tags.append(
                int(lines[cursor])
            )
            cursor += 1

        # This is the exact order in which mesh_gmsh.c assigns
        # local PEEC node indices.
        node_order.extend(
            tags
        )

        for tag in tags:
            values = lines[cursor].split()
            cursor += 1

            if len(values) < 3:
                raise RuntimeError(
                    "Invalid node coordinate line."
                )

            node_xyz[tag] = (
                float(values[0]),
                float(values[1]),
                float(values[2]),
            )

            # Parametric coordinates, when present, are on the
            # same line after xyz and can simply be ignored.
            _ = entity_dim
            _ = parametric

    # --------------------------------------------------------
    # Elements
    # --------------------------------------------------------
    e0, e1 = _find_section(
        lines,
        "$Elements",
        "$EndElements",
    )

    header = lines[e0].split()

    if len(header) != 4:
        raise RuntimeError(
            "Unsupported $Elements header."
        )

    e_blocks = int(
        header[0]
    )

    cursor = e0 + 1

    quads: list[
        tuple[
            int,
            tuple[int, int, int, int],
        ]
    ] = []

    for _ in range(e_blocks):
        block_header = lines[cursor].split()
        cursor += 1

        entity_dim = int(
            block_header[0]
        )

        entity_tag = int(
            block_header[1]
        )

        element_type = int(
            block_header[2]
        )

        count = int(
            block_header[3]
        )

        for _ in range(count):
            values = [
                int(value)
                for value in lines[cursor].split()
            ]

            cursor += 1

            # Gmsh type 3 = 4-node quadrangle.
            if entity_dim == 2 and element_type == 3:
                if len(values) != 5:
                    raise RuntimeError(
                        "Unexpected quadrangle element format."
                    )

                quads.append(
                    (
                        entity_tag,
                        (
                            values[1],
                            values[2],
                            values[3],
                            values[4],
                        ),
                    )
                )

    if not node_xyz:
        raise RuntimeError(
            f"No nodes in {path}"
        )

    if not quads:
        raise RuntimeError(
            f"No quadrangles in {path}"
        )

    tag_to_internal = {
        tag: index
        for index, tag in enumerate(node_order)
    }

    return ParsedMesh(
        node_xyz=node_xyz,
        node_order=node_order,
        tag_to_internal=tag_to_internal,
        quads=quads,
    )


# ============================================================
# INTERNAL NODE INDEX
# ============================================================
#
# The C Gmsh reader does NOT assume:
#
#     internal_index = gmsh_tag - 1.
#
# It scans the $Nodes entity blocks in file order and assigns:
#
#     0, 1, 2, ...
#
# in exactly that order.
#
# Therefore all generated auxiliary files must use the mapping
# reconstructed from the actual written .msh file.
# ============================================================

def internal_node(
    mesh: ParsedMesh,
    gmsh_tag: int,
) -> int:
    try:
        return mesh.tag_to_internal[
            gmsh_tag
        ]
    except KeyError as exc:
        raise RuntimeError(
            f"Unknown Gmsh node tag {gmsh_tag}."
        ) from exc


def validate_node_index_mapping(
    mesh: ParsedMesh,
    path: Path,
) -> None:
    if len(mesh.node_order) != len(mesh.node_xyz):
        raise RuntimeError(
            f"{path}: duplicate or missing node tags in $Nodes."
        )

    if len(mesh.tag_to_internal) != len(mesh.node_xyz):
        raise RuntimeError(
            f"{path}: invalid Gmsh-tag -> PEEC-index mapping."
        )

    expected = set(
        range(len(mesh.node_order))
    )

    actual = set(
        mesh.tag_to_internal.values()
    )

    if actual != expected:
        raise RuntimeError(
            f"{path}: internal PEEC node mapping is not contiguous."
        )


# ============================================================
# TOPOLOGY RECONSTRUCTION
# ============================================================

def build_global_edges(
    mesh: ParsedMesh,
) -> tuple[
    list[tuple[int, int]],
    dict[tuple[int, int], int],
]:
    """
    Reproduce the CURRENT PEEC topology exactly.

    The C code scans quads in file order and local edges:

        q0 -> q1
        q1 -> q2
        q2 -> q3
        q3 -> q0

    The canonical pair (min,max) is ONLY the lookup key.

    The stored global edge orientation is the orientation of the
    first local occurrence, exactly as in topology.c.
    """

    edges: list[
        tuple[int, int]
    ] = []

    edge_to_id: dict[
        tuple[int, int],
        int,
    ] = {}

    for _entity_tag, quad in mesh.quads:
        q = [
            internal_node(
                mesh,
                tag,
            )
            for tag in quad
        ]

        local_pairs = [
            (q[0], q[1]),
            (q[1], q[2]),
            (q[2], q[3]),
            (q[3], q[0]),
        ]

        for a, b in local_pairs:
            key = (
                min(a, b),
                max(a, b),
            )

            if key not in edge_to_id:
                edge_to_id[key] = len(
                    edges
                )

                # IMPORTANT:
                # keep the first-appearance orientation a -> b.
                edges.append(
                    (a, b)
                )

    return edges, edge_to_id


# ============================================================
# APERTURE BOUNDARY
# ============================================================

def is_close(
    a: float,
    b: float,
    tolerance: float,
) -> bool:
    return abs(a - b) <= tolerance


def aperture_boundary_tags(
    mesh: ParsedMesh,
    geometry: Geometry,
    tolerance: float,
) -> list[int]:
    g = geometry

    result: list[int] = []

    for tag, xyz in mesh.node_xyz.items():
        x, y, z = xyz

        if not is_close(
            y,
            g.y1,
            tolerance,
        ):
            continue

        inside_x = (
            g.open_x0 - tolerance
            <= x
            <= g.open_x1 + tolerance
        )

        inside_z = (
            g.open_z0 - tolerance
            <= z
            <= g.open_z1 + tolerance
        )

        if not (
            inside_x
            and inside_z
        ):
            continue

        on_boundary = (
            is_close(
                x,
                g.open_x0,
                tolerance,
            )
            or is_close(
                x,
                g.open_x1,
                tolerance,
            )
            or is_close(
                z,
                g.open_z0,
                tolerance,
            )
            or is_close(
                z,
                g.open_z1,
                tolerance,
            )
        )

        if on_boundary:
            result.append(
                tag
            )

    result.sort()

    return result


# ============================================================
# COORDINATE MAPPING
# ============================================================

def distance_squared(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return (
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def build_node_mapping(
    closed_mesh: ParsedMesh,
    open_mesh: ParsedMesh,
    geometry: Geometry,
    tolerance: float,
) -> tuple[list[dict[str, int]], float]:
    closed_tags = aperture_boundary_tags(
        closed_mesh,
        geometry,
        tolerance,
    )

    open_tags = aperture_boundary_tags(
        open_mesh,
        geometry,
        tolerance,
    )

    if not closed_tags:
        raise RuntimeError(
            "No closed-mesh aperture boundary nodes found."
        )

    if not open_tags:
        raise RuntimeError(
            "No open-mesh aperture boundary nodes found."
        )

    remaining_open = set(
        open_tags
    )

    mapping: list[
        dict[str, int]
    ] = []

    max_error = 0.0

    for closed_tag in closed_tags:
        closed_xyz = closed_mesh.node_xyz[
            closed_tag
        ]

        best_open_tag: int | None = None
        best_d2 = math.inf

        for open_tag in remaining_open:
            d2 = distance_squared(
                closed_xyz,
                open_mesh.node_xyz[
                    open_tag
                ],
            )

            if d2 < best_d2:
                best_d2 = d2
                best_open_tag = open_tag

        if best_open_tag is None:
            raise RuntimeError(
                "Aperture node mapping failed."
            )

        error = math.sqrt(
            best_d2
        )

        if error > tolerance:
            raise RuntimeError(
                "Closed/open aperture meshes are not compatible: "
                f"nearest-node error = {error:.6e} m, "
                f"tolerance = {tolerance:.6e} m."
            )

        remaining_open.remove(
            best_open_tag
        )

        mapping.append(
            {
                "closed_node": internal_node(
                    closed_mesh,
                    closed_tag,
                ),
                "open_node": internal_node(
                    open_mesh,
                    best_open_tag,
                ),
            }
        )

        max_error = max(
            max_error,
            error,
        )

    if remaining_open:
        raise RuntimeError(
            "Open aperture boundary has extra nodes. "
            f"closed={len(closed_tags)}, open={len(open_tags)}."
        )

    return mapping, max_error


# ============================================================
# COVER EDGES
# ============================================================

def build_cover_edges(
    closed_mesh: ParsedMesh,
    cover_surface_entity: int,
) -> tuple[
    list[int],
    int,
]:
    """
    Find every global PEEC edge belonging to TEMP_COVER.

    This includes:
        perimeter edges
        cut edges
        internal cover edges.

    aperture_coupling.c later classifies them automatically.
    """

    _, edge_to_id = build_global_edges(
        closed_mesh
    )

    cover_quads = [
        quad
        for entity_tag, quad in closed_mesh.quads
        if entity_tag == cover_surface_entity
    ]

    if not cover_quads:
        raise RuntimeError(
            "TEMP_COVER has no quadrangle elements."
        )

    cover_edge_ids: set[int] = set()

    for quad in cover_quads:
        q = [
            internal_node(
                closed_mesh,
                tag,
            )
            for tag in quad
        ]

        local_pairs = [
            (q[0], q[1]),
            (q[1], q[2]),
            (q[2], q[3]),
            (q[3], q[0]),
        ]

        for a, b in local_pairs:
            key = (
                min(a, b),
                max(a, b),
            )

            edge_id = edge_to_id.get(
                key
            )

            if edge_id is None:
                raise RuntimeError(
                    "Internal topology reconstruction error."
                )

            cover_edge_ids.add(
                edge_id
            )

    return (
        sorted(
            cover_edge_ids
        ),
        len(
            cover_quads
        ),
    )


# ============================================================
# DISTRIBUTED SLOT CELLS
# ============================================================

def _slot_longitudinal_axis(
    geometry: Geometry,
    requested: str,
) -> str:
    if requested in ("x", "z"):
        return requested

    # auto: the longer aperture dimension is treated as the slot
    # longitudinal direction. For a square aperture choose x
    # deterministically.
    if geometry.opening_w >= geometry.opening_h:
        return "x"

    return "z"


def _slot_side_tags(
    mesh: ParsedMesh,
    geometry: Geometry,
    axis: str,
    side: str,
    tolerance: float,
) -> list[int]:
    """
    Return one of the two opposite aperture-boundary sides.

    axis == "x":
        longitudinal coordinate = x,
        cells bridge z=open_z0 -> z=open_z1.

    axis == "z":
        longitudinal coordinate = z,
        cells bridge x=open_x0 -> x=open_x1.
    """

    result: list[int] = []

    for tag, xyz in mesh.node_xyz.items():
        x, y, z = xyz

        if not is_close(
            y,
            geometry.y1,
            tolerance,
        ):
            continue

        if axis == "x":
            target = (
                geometry.open_z0
                if side == "a"
                else geometry.open_z1
            )

            if not is_close(
                z,
                target,
                tolerance,
            ):
                continue

            if not (
                geometry.open_x0 - tolerance
                <= x
                <= geometry.open_x1 + tolerance
            ):
                continue

        else:
            target = (
                geometry.open_x0
                if side == "a"
                else geometry.open_x1
            )

            if not is_close(
                x,
                target,
                tolerance,
            ):
                continue

            if not (
                geometry.open_z0 - tolerance
                <= z
                <= geometry.open_z1 + tolerance
            ):
                continue

        result.append(
            tag
        )

    coordinate = (
        (lambda tag: mesh.node_xyz[tag][0])
        if axis == "x"
        else
        (lambda tag: mesh.node_xyz[tag][2])
    )

    result.sort(
        key=coordinate
    )

    return result


def build_slot_cells(
    open_mesh: ParsedMesh,
    geometry: Geometry,
    tolerance: float,
    axis_request: str,
) -> tuple[
    list[tuple[int, int, float]],
    dict[str, float | int | str],
]:
    """
    Build distributed equivalent aperture cells directly from the
    OPEN mesh.

    Every cell connects matching nodes on two opposite sides of
    the aperture.

    The third value is the longitudinal dual support length dl_m.

    End nodes receive half of the adjacent interval; interior
    nodes receive half of the interval on each side. Therefore:

        sum_m dl_m = slot longitudinal length.
    """

    axis = _slot_longitudinal_axis(
        geometry,
        axis_request,
    )

    side_a = _slot_side_tags(
        open_mesh,
        geometry,
        axis,
        "a",
        tolerance,
    )

    side_b = _slot_side_tags(
        open_mesh,
        geometry,
        axis,
        "b",
        tolerance,
    )

    if len(side_a) < 2 or len(side_b) < 2:
        raise RuntimeError(
            "Not enough nodes on opposite slot sides."
        )

    coordinate = (
        (lambda tag: open_mesh.node_xyz[tag][0])
        if axis == "x"
        else
        (lambda tag: open_mesh.node_xyz[tag][2])
    )

    # Match the two sides by longitudinal coordinate.
    remaining_b = set(
        side_b
    )

    pairs: list[
        tuple[int, int, float]
    ] = []

    max_pair_error = 0.0

    for tag_a in side_a:
        s_a = coordinate(
            tag_a
        )

        best_b: int | None = None
        best_error = math.inf

        for tag_b in remaining_b:
            error = abs(
                coordinate(tag_b)
                - s_a
            )

            if error < best_error:
                best_error = error
                best_b = tag_b

        if best_b is None:
            raise RuntimeError(
                "Slot-side pairing failed."
            )

        if best_error > tolerance:
            raise RuntimeError(
                "Opposite slot sides have incompatible node positions: "
                f"error={best_error:.9e} m, "
                f"tolerance={tolerance:.9e} m."
            )

        remaining_b.remove(
            best_b
        )

        max_pair_error = max(
            max_pair_error,
            best_error,
        )

        pairs.append(
            (
                tag_a,
                best_b,
                s_a,
            )
        )

    if remaining_b:
        raise RuntimeError(
            "Opposite slot sides contain different numbers of nodes: "
            f"A={len(side_a)}, B={len(side_b)}."
        )

    pairs.sort(
        key=lambda item: item[2]
    )

    positions = [
        item[2]
        for item in pairs
    ]

    support: list[float] = []

    for i, s in enumerate(positions):
        if i == 0:
            dl = 0.5 * (
                positions[1]
                - positions[0]
            )
        elif i + 1 == len(positions):
            dl = 0.5 * (
                positions[-1]
                - positions[-2]
            )
        else:
            dl = 0.5 * (
                positions[i + 1]
                - positions[i - 1]
            )

        if dl <= 0.0:
            raise RuntimeError(
                "Non-positive slot-cell support length."
            )

        support.append(
            dl
        )

    cells: list[
        tuple[int, int, float]
    ] = []

    for (tag_a, tag_b, _), dl in zip(
        pairs,
        support,
    ):
        cells.append(
            (
                internal_node(
                    open_mesh,
                    tag_a,
                ),
                internal_node(
                    open_mesh,
                    tag_b,
                ),
                dl,
            )
        )

    longitudinal_length = (
        geometry.opening_w
        if axis == "x"
        else geometry.opening_h
    )

    gap_width = (
        geometry.opening_h
        if axis == "x"
        else geometry.opening_w
    )

    total_support = sum(
        dl
        for _, _, dl in cells
    )

    if abs(total_support - longitudinal_length) > max(
        tolerance,
        1.0e-10 * longitudinal_length,
    ):
        raise RuntimeError(
            "Slot-cell dual supports do not cover the full slot length: "
            f"sum(dl)={total_support:.9e}, "
            f"length={longitudinal_length:.9e}."
        )

    diagnostics: dict[
        str,
        float | int | str,
    ] = {
        "axis": axis,
        "n_cells": len(cells),
        "slot_length_m": longitudinal_length,
        "slot_gap_width_m": gap_width,
        "support_sum_m": total_support,
        "max_pair_error_m": max_pair_error,
    }

    return cells, diagnostics


def write_slot_cells(
    output: Path,
    open_mesh: ParsedMesh,
    geometry: Geometry,
    tolerance: float,
    axis_request: str,
) -> dict[str, float | int | str]:
    cells, diagnostics = build_slot_cells(
        open_mesh,
        geometry,
        tolerance,
        axis_request,
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    axis = diagnostics[
        "axis"
    ]

    lines = [
        "# ============================================================",
        "# PEEC distributed equivalent slot cells",
        "#",
        "# node_a node_b support_length_m",
        "#",
        "# Node indices are INTERNAL 0-based PEEC indices reconstructed",
        "# exactly from the Gmsh $Nodes file order.",
        "#",
        f"# longitudinal_axis = {axis}",
        f"# slot_length_m      = {diagnostics['slot_length_m']:.17e}",
        f"# slot_gap_width_m   = {diagnostics['slot_gap_width_m']:.17e}",
        "# ============================================================",
        "",
    ]

    for node_a, node_b, dl in cells:
        lines.append(
            f"{node_a} {node_b} {dl:.17e}"
        )

    output.write_text(
        "\n".join(lines)
        + "\n",
        encoding="utf-8",
    )

    return diagnostics


# ============================================================
# JSON
# ============================================================

def write_aperture_map(
    output: Path,
    closed_mesh_path: Path,
    open_mesh_path: Path,
    geometry: Geometry,
    h: float,
    closed_mesh: ParsedMesh,
    open_mesh: ParsedMesh,
    cover_surface_entity: int,
    tolerance: float,
) -> dict[str, float | int]:
    node_mapping, max_mapping_error = (
        build_node_mapping(
            closed_mesh,
            open_mesh,
            geometry,
            tolerance,
        )
    )

    cover_edges, cover_quads = (
        build_cover_edges(
            closed_mesh,
            cover_surface_entity,
        )
    )

    # --------------------------------------------------------
    # Compact JSON format consumed by the current C solver.
    #
    # Only the data required for aperture coupling are written:
    #
    #   node_mapping : closed boundary node -> open boundary node
    #   cover_edges  : all PEEC edges belonging to TEMP_COVER
    #
    # Geometry and diagnostics remain available in the generator
    # itself and are printed to the terminal, but are deliberately
    # not stored in aperture_map.json.
    # --------------------------------------------------------
    data = {
        "node_mapping": node_mapping,
        "cover_edges": cover_edges,
    }

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output.write_text(
        json.dumps(
            data,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "boundary_mappings": len(node_mapping),
        "cover_quads": cover_quads,
        "n_cover_edges": len(cover_edges),
        "max_mapping_error_m": max_mapping_error,
    }


# ============================================================
# DIAGNOSTICS
# ============================================================

def print_mesh_summary(
    title: str,
    path: Path,
    mesh: ParsedMesh,
) -> None:
    edges, _ = build_global_edges(
        mesh
    )

    print()
    print(title)
    print("=" * len(title))

    print(
        f"File  : {path}"
    )

    print(
        f"Nodes : {len(mesh.node_xyz)}"
    )

    print(
        f"Edges : {len(edges)}"
    )

    print(
        f"Quads : {len(mesh.quads)}"
    )

    print(
        "Euler : "
        f"{len(mesh.node_xyz) - len(edges) + len(mesh.quads)}"
    )


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a matched PEEC lightning mesh pair: "
            "closed body with temporary PEC cover, open body "
            "with aperture, aperture_map.json, and distributed "
            "slot_cells.dat using exact PEEC internal indices."
        )
    )

    parser.add_argument(
        "--closed",
        default="meshes/body_closed.msh",
        help="Closed Stage-1 mesh.",
    )

    parser.add_argument(
        "--open",
        dest="open_file",
        default="meshes/body_open.msh",
        help="Open Stage-2 mesh.",
    )

    parser.add_argument(
        "--map",
        dest="map_file",
        default="cases/aperture_map.json",
        help="Generated aperture map JSON.",
    )


    parser.add_argument(
        "--slot-cells",
        dest="slot_cells_file",
        default="cases/slot_cells.dat",
        help=(
            "Generated distributed slot-cell file for the "
            "one-stage equivalent-aperture model."
        ),
    )

    parser.add_argument(
        "--slot-axis",
        choices=[
            "auto",
            "x",
            "z",
        ],
        default="auto",
        help=(
            "Longitudinal direction of the equivalent slot cells. "
            "auto chooses the longer aperture dimension."
        ),
    )

    parser.add_argument(
        "--h",
        type=float,
        default=0.25,
        help="Target mesh size [m].",
    )

    parser.add_argument(
        "--length",
        type=float,
        default=7.0,
    )

    parser.add_argument(
        "--width",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--height",
        type=float,
        default=2.0,
    )

    parser.add_argument(
        "--opening-w",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--opening-h",
        type=float,
        default=1.0,
    )

    parser.add_argument(
        "--mapping-tol",
        type=float,
        default=None,
        help=(
            "Closed/open aperture-node matching tolerance [m]. "
            "Default: max(1e-10, 1e-8*h)."
        ),
    )

    parser.add_argument(
        "--gui",
        choices=[
            "none",
            "closed",
            "open",
        ],
        default="none",
        help="Open one generated mesh in the Gmsh GUI.",
    )

    args = parser.parse_args()

    if args.h <= 0.0:
        raise SystemExit(
            "--h must be > 0."
        )

    geometry = Geometry(
        length=args.length,
        width=args.width,
        height=args.height,
        opening_w=args.opening_w,
        opening_h=args.opening_h,
    )

    if (
        geometry.length <= 0.0
        or geometry.width <= 0.0
        or geometry.height <= 0.0
    ):
        raise SystemExit(
            "Body dimensions must be > 0."
        )

    if (
        geometry.opening_w <= 0.0
        or geometry.opening_h <= 0.0
        or geometry.opening_w >= geometry.length
        or geometry.opening_h >= geometry.height
    ):
        raise SystemExit(
            "Opening must fit strictly inside the front wall."
        )

    closed_path = Path(
        args.closed
    )

    open_path = Path(
        args.open_file
    )

    map_path = Path(
        args.map_file
    )

    slot_cells_path = Path(
        args.slot_cells_file
    )

    tolerance = (
        args.mapping_tol
        if args.mapping_tol is not None
        else max(
            1.0e-10,
            1.0e-8 * args.h,
        )
    )

    gmsh.initialize()

    try:
        print()
        print(
            "Generating OPEN Stage-2 mesh ..."
        )

        generate_one_mesh(
            open_path,
            geometry,
            args.h,
            with_cover=False,
            model_name="body_open",
        )

        print()
        print(
            "Generating CLOSED Stage-1 mesh ..."
        )

        cover_surface_entity = generate_one_mesh(
            closed_path,
            geometry,
            args.h,
            with_cover=True,
            model_name="body_closed",
        )

        if cover_surface_entity is None:
            raise RuntimeError(
                "Closed mesh has no TEMP_COVER surface."
            )

        # ----------------------------------------------------
        # Parse files exactly as written.
        # ----------------------------------------------------
        open_mesh = parse_gmsh_41_ascii(
            open_path
        )

        closed_mesh = parse_gmsh_41_ascii(
            closed_path
        )

        validate_node_index_mapping(
            open_mesh,
            open_path,
        )

        validate_node_index_mapping(
            closed_mesh,
            closed_path,
        )

        # ----------------------------------------------------
        # Build JSON automatically.
        # ----------------------------------------------------
        diagnostics = write_aperture_map(
            map_path,
            closed_path,
            open_path,
            geometry,
            args.h,
            closed_mesh,
            open_mesh,
            cover_surface_entity,
            tolerance,
        )


        slot_diagnostics = write_slot_cells(
            slot_cells_path,
            open_mesh,
            geometry,
            tolerance,
            args.slot_axis,
        )

        print_mesh_summary(
            "OPEN BODY",
            open_path,
            open_mesh,
        )

        print_mesh_summary(
            "CLOSED BODY",
            closed_path,
            closed_mesh,
        )

        print()
        print("APERTURE MAP")
        print("============")

        print(
            f"File              : {map_path}"
        )

        print(
            "Boundary mappings : "
            f"{diagnostics['boundary_mappings']}"
        )

        print(
            "Cover quads       : "
            f"{diagnostics['cover_quads']}"
        )

        print(
            "Cover edges       : "
            f"{diagnostics['n_cover_edges']}"
        )

        print(
            "Max map error     : "
            f"{diagnostics['max_mapping_error_m']:.9e} m"
        )

        print()
        print("SLOT CELLS")
        print("==========")

        print(
            f"File              : {slot_cells_path}"
        )

        print(
            "Axis              : "
            f"{slot_diagnostics['axis']}"
        )

        print(
            "Cells             : "
            f"{slot_diagnostics['n_cells']}"
        )

        print(
            "Slot length       : "
            f"{slot_diagnostics['slot_length_m']:.9e} m"
        )

        print(
            "Gap width         : "
            f"{slot_diagnostics['slot_gap_width_m']:.9e} m"
        )

        print(
            "Support sum       : "
            f"{slot_diagnostics['support_sum_m']:.9e} m"
        )

        print(
            "Max pair error    : "
            f"{slot_diagnostics['max_pair_error_m']:.9e} m"
        )

        print()
        print("Generated successfully:")
        print(
            f"  Stage 1: {closed_path}"
        )
        print(
            f"  Stage 2: {open_path}"
        )
        print(
            f"  Map    : {map_path}"
        )

        print(
            f"  Slot   : {slot_cells_path}"
        )

        # ----------------------------------------------------
        # Optional GUI.
        # ----------------------------------------------------
        if args.gui != "none":
            gui_path = (
                closed_path
                if args.gui == "closed"
                else open_path
            )

            gmsh.clear()

            gmsh.open(
                str(gui_path)
            )

            gmsh.fltk.run()

    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
