#!/usr/bin/env python3
"""End-to-end Fortran checks and comparisons against the original C executable.
Only the Python standard library is required. Run from the repository root.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def run(exe, args, log, success=True):
    proc = subprocess.run([str(exe), *map(str, args)], cwd=ROOT, text=True,
                          capture_output=True, timeout=240)
    log.write_text(shlex.join([str(exe), *map(str, args)]) + "\n" +
                   proc.stdout + "\n" + proc.stderr, encoding="utf-8")
    if (proc.returncode == 0) != success:
        raise AssertionError(f"unexpected exit {proc.returncode}: {log}\n{proc.stdout}\n{proc.stderr}")
    return proc.stdout


def fixture(path):
    xyz = [(x / 2, y / 2, 0) for y in range(3) for x in range(3)]
    quads = [(3*y+x, 3*y+x+1, 3*(y+1)+x+1, 3*(y+1)+x)
             for y in range(2) for x in range(2)]
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat",
             "$Nodes", "1 9 1 9", "2 1 0 9",
             *map(str, range(1, 10)),
             *(" ".join(map(str, p)) for p in xyz),
             "$EndNodes", "$Elements", "1 4 1 4", "2 1 3 4"]
    lines += [f"{i+1} " + " ".join(str(v+1) for v in q) for i, q in enumerate(quads)]
    path.write_text("\n".join(lines + ["$EndElements", ""]), encoding="ascii")
    edges = []
    for q in quads:
        for a, b in zip(q, q[1:] + q[:1]):
            if (a, b) not in edges and (b, a) not in edges:
                edges.append((a, b))
    return edges


def matrix(path):
    lines = path.read_text().splitlines()
    header = lines[0].split()
    data = [line.split() for line in lines[1:] if line.strip() and not line.startswith("%")]
    rows, cols = map(int, data[0][:2])
    a = [0.0] * (rows*cols)
    if header[2] == "array":
        values = [float(v) for row in data[1:] for v in row]
        assert len(values) == len(a), path
        for j in range(cols):
            for i in range(rows):
                a[i*cols+j] = values[j*rows+i]
    else:
        assert len(data)-1 == int(data[0][2]), path
        for i, j, v in data[1:]:
            a[(int(i)-1)*cols+int(j)-1] += float(v)
    assert all(map(math.isfinite, a)), path
    return rows, cols, a


def vtk(path):
    words = path.read_text().split()
    result = {}
    n = 0
    i = 0
    while i < len(words):
        if words[i] in ("POINT_DATA", "CELL_DATA"):
            n = int(words[i+1]); i += 2
        elif words[i] == "SCALARS":
            name = words[i+1]
            assert words[i+4:i+6] == ["LOOKUP_TABLE", "default"], path
            result[name] = list(map(float, words[i+6:i+6+n]))
            assert len(result[name]) == n, path
            assert all(map(math.isfinite, result[name])), path
            i += 6+n
        else:
            i += 1
    assert result, path
    return result


def close(actual, expected, tol, name):
    assert len(actual) == len(expected), name
    scale = max([abs(x) for x in expected] + [1e-280])
    error = max([abs(a-b) for a, b in zip(actual, expected)] + [0]) / scale
    assert error <= tol, f"{name}: relative max error {error:.6g} > {tol}"
    return error


def scenario(exe, name, task, common, directory, config=False):
    out = directory / name
    args = ["--task", task, *common, "--parallel", "2", "--vtk", "--vtk-dir", str(out),
            "--save-matrices", "--matrix-dir", str(out / "matrices")]
    if task == "rcs":
        args += ["--output", str(out / "rcs.csv")]
    if config:
        cfg = directory / (name + ".cfg")
        cfg.write_text("# Generated Fortran parser check\n" + shlex.join(list(map(str, args))) + "\n",
                       encoding="utf-8")
        run(exe, ["--config", cfg], directory / (name + ".log"))
    else:
        run(exe, args, directory / (name + ".log"))
    assert list(out.rglob("*.vtk")), out
    for p in out.rglob("*.vtk"):
        vtk(p)
    for p in out.rglob("*.mtx"):
        matrix(p)
    return out


def compare_outputs(actual, expected, name):
    worst_matrix = worst_field = 0.0
    for ref in expected.rglob("*.mtx"):
        got = actual / ref.relative_to(expected)
        r, c, values = matrix(ref)
        gr, gc, gvalues = matrix(got)
        assert (r, c) == (gr, gc), got
        worst_matrix = max(worst_matrix, close(gvalues, values, 2e-11, str(got)))
    for ref in expected.rglob("*.vtk"):
        got = actual / ref.relative_to(expected)
        a, b = vtk(got), vtk(ref)
        assert set(a) == set(b), (got, set(a), set(b))
        for key in b:
            # One norm for each complex pair avoids relative error on a component
            # that vanishes by symmetry. The absolute fields are compared directly.
            if "_re_" in key or "_im_" in key:
                stem, suffix = key.rsplit("_re_", 1) if "_re_" in key else key.rsplit("_im_", 1)
                magnitude = stem + "_abs_" + suffix
                scale = max(b.get(magnitude, [0.0]))
                floor = 2e-7 * max(scale, 1e-280)
                assert max(abs(x-y) for x, y in zip(a[key], b[key])) <= floor, (got, key)
            else:
                worst_field = max(worst_field, close(a[key], b[key], 2e-7, f"{got}:{key}"))
    ref = expected / "rcs.csv"
    if ref.exists():
        with ref.open() as f:
            b = list(csv.DictReader(f))
        with (actual / "rcs.csv").open() as f:
            a = list(csv.DictReader(f))
        assert len(a) == len(b)
        close([float(r["sigma_m2"]) for r in a], [float(r["sigma_m2"]) for r in b],
              2e-7, name + " RCS")
    for ref in expected.glob("*_shunt_*.dat"):
        def values(p):
            return {r.split()[0]: float(r.split()[1]) for r in p.read_text().splitlines()
                    if r.split()[0] in ("time", "shunt_current_A", "shunt_resistance_ohm")}
        a, b = values(actual/ref.name), values(ref)
        for k in b:
            close([a[k]], [b[k]], 2e-7, name+" "+k)
    print(f"PASS {name}: matrices {worst_matrix:.3g}; fields {worst_field:.3g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fortran", required=True, type=Path)
    parser.add_argument("--c-solver", type=Path)
    parser.add_argument("--smoke-only", action="store_true")
    args = parser.parse_args()
    exe = args.fortran.resolve()
    parent = ROOT / "build" / "fortran-tests"
    parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="run-", dir=parent))
    mesh = work / "tiny.msh"
    edges = fixture(mesh)
    cover = [i for i, ab in enumerate(edges) if 4 in ab]
    mapping = {"node_mapping": [{"closed_node": v, "open_node": v} for v in (1, 3, 5, 7)],
               "cover_edges": cover}
    aperture = work / "aperture.json"
    aperture.write_text(json.dumps(mapping))
    cells = work / "slot.dat"
    cells.write_text("0 8 0.5\n2 6 0.5\n")
    trans = ["--dt", "1e-7", "--t-end", "3e-7", "--time-order", "2",
             "--strike-node", "0", "--return-node", "8", "--lightning-K", "1"]
    shunt = ["--shunt-R", "1000", "--shunt-a-node", "1", "--shunt-b-node", "7"]
    tests = [
        ("scattering", "scattering", ["--mesh", mesh, "--frequency", "1e6", "--phi", "23"]),
        ("rcs", "rcs", ["--mesh", mesh, "--ka", "0.02", "--a", "1", "--rcs-order", "4",
                        "--phi-start", "0", "--phi-end", "180", "--phi-step", "45"]),
        ("closed", "lightning", ["--lightning-case", "closed", "--closed-mesh", mesh, *trans]),
        ("two-stage", "lightning", ["--lightning-case", "two-stage", "--closed-mesh", mesh,
                                  "--open-mesh", mesh, "--aperture-map", aperture, *trans, *shunt]),
        ("slot-cells", "lightning", ["--lightning-case", "slot-cells", "--open-mesh", mesh,
                                   "--slot-cells", cells, "--slot-width", ".1",
                                   "--slot-wall-span", "2", "--slot-wall-thickness", ".002",
                                   *trans, *shunt]),
        ("incident-pulse", "lightning", ["--lightning-case", "closed", "--closed-mesh", mesh,
                                        "--excitation", "incident-pulse", "--pulse-K", "1",
                                        "--dt", "1e-7", "--t-end", "3e-7", "--time-order", "1"]),
    ]
    for name, task, common in tests:
        got = scenario(exe, "f-"+name, task, common, work, config=(name=="scattering"))
        print("PASS Fortran end-to-end " + name)
        if args.c_solver and not args.smoke_only:
            ref = scenario(args.c_solver.resolve(), "c-"+name, task, common, work)
            compare_outputs(got, ref, name)
    rejects = [
        ["--task", "rcs", "--mesh", mesh, "--rcs-order", "7"],
        ["--task", "rcs", "--mesh", mesh, "--phi-step", "0"],
        ["--task", "scattering", "--mesh", mesh, "--frequency", "nan"],
        ["--task", "mesh-info", "--mesh", mesh, "--physics", "retarded"],
        ["--task", "mesh-info", "--unknown-option"],
        ["--task", "lightning", "--closed-mesh", mesh, *trans, "--strike-node", "900"],
        ["--task", "lightning", "--closed-mesh", mesh, *trans, "--t-end", "3.5e-7"],
    ]
    bad = work / "bad-aperture.json"
    badmap = json.loads(json.dumps(mapping))
    badmap["node_mapping"][0]["open_node"] = 0
    bad.write_text(json.dumps(badmap))
    rejects.append(["--task", "lightning", "--lightning-case", "two-stage",
                    "--closed-mesh", mesh, "--open-mesh", mesh, "--aperture-map", bad, *trans])
    for i, case in enumerate(rejects):
        run(exe, case, work/f"reject-{i}.log", success=False)
    print(f"PASS {len(rejects)} invalid-input cases")
    if args.c_solver and not args.smoke_only:
        common = ["--mesh", ROOT/"meshes/plate.msh", "--frequency", "1e6", "--phi", "23"]
        got = scenario(exe, "f-plate", "scattering", common, work)
        ref = scenario(args.c_solver.resolve(), "c-plate", "scattering", common, work)
        compare_outputs(got, ref, "original plate mesh")
    print("Reports: " + str(work))


if __name__ == "__main__":
    main()
