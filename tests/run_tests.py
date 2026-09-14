#!/usr/bin/env python3
"""Run PEEC checks independently; preserve every log; return nonzero on failure.

The default and sanitizer groups use only the Python standard library.
The separate spice group additionally requires NumPy and SciPy.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def project_path(value: str | None, default: str) -> Path:
    path = Path(value or default)
    return (path if path.is_absolute() else ROOT / path).resolve()


class Suite:
    def __init__(self, args: argparse.Namespace):
        self.args = args
        output = ROOT / "build" / "verification" / "runs"
        output.mkdir(parents=True, exist_ok=True)
        self.output = Path(tempfile.mkdtemp(prefix="run-", dir=output))
        self.results: list[dict] = []
        self.core = project_path(args.core, "build/verification/test_core")
        self.solver = project_path(args.solver, "build/verification/PeecSolver")
        self.config = project_path(args.config_parser, "build/verification/test_config_cli")
        self.asan = project_path(args.asan, "build/verification/asan/test_config_lifetime")

    def process(self, name: str, command: list[str]) -> subprocess.CompletedProcess:
        log = self.output / (name + ".log")
        heading = "$ " + shlex.join(command) + "\n\n"
        try:
            result = subprocess.run(command, cwd=ROOT, capture_output=True,
                                    text=True, errors="replace", timeout=120)
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout or b""
            stderr = exc.stderr or b""
            if isinstance(stdout, bytes):
                stdout = stdout.decode(errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode(errors="replace")
            log.write_text(heading + stdout + stderr + "\nTIMEOUT (120 s)\n", encoding="utf-8")
            raise RuntimeError("Process exceeded 120 seconds; see its log.") from exc
        log.write_text(heading + result.stdout + "\n[stderr]\n" + result.stderr
                       + f"\n[exit] {result.returncode}\n", encoding="utf-8")
        return result

    def check(self, name: str, action) -> None:
        if self.args.only and name not in self.args.only:
            return
        start = time.monotonic()
        try:
            passed, detail = action()
        except Exception as exc:
            passed, detail = False, f"{type(exc).__name__}: {exc}"
        log = self.output / (name + ".log")
        if not log.exists():
            log.write_text(detail + "\n", encoding="utf-8")
        entry = {"name": name, "status": "PASS" if passed else "FAIL",
                 "detail": detail, "seconds": time.monotonic() - start,
                 "log": str(log.relative_to(ROOT))}
        self.results.append(entry)
        print(f"{entry['status']:4} {name}: {detail}", flush=True)

    def core_tests(self) -> None:
        listing = self.process("core_list", [str(self.core), "--list"])
        if listing.returncode != 0 or not listing.stdout.strip():
            self.check("core_discovery", lambda: (False, "Cannot list C tests."))
            return
        for name in listing.stdout.splitlines():
            name = name.strip()
            if not re.fullmatch(r"[a-z0-9_]+", name):
                raise RuntimeError("Unexpected C test name: " + repr(name))
            def run(name=name):
                result = self.process(name, [str(self.core), name])
                lines = [line for line in result.stdout.splitlines()
                         if line and not line.startswith(("PASS ", "FAIL "))]
                detail = lines[-1] if lines else f"exit={result.returncode}"
                return result.returncode == 0 and f"PASS {name}" in result.stdout, detail
            self.check(name, run)

    def harmonic(self, name: str, frequency: str) -> tuple[bool, str]:
        result = self.process(name, [str(self.solver), "--task", "scattering",
            "--mesh", "meshes/plate.msh", "--frequency", frequency, "--parallel", "1"])
        if result.returncode != 0:
            errors = [line for line in result.stderr.splitlines() if "ERROR" in line]
            return False, errors[0] if errors else f"exit={result.returncode}"
        match = re.search(r"Backward error\s*:\s*(\S+)", result.stdout)
        if not match:
            return False, "No numerical residual reported."
        value = float(match.group(1))
        ok = math.isfinite(value) and value <= 1e-9 and "ERROR:" not in result.stderr
        return ok, f"f={frequency} Hz; backward error={value:.3e}; limit=1e-9"

    def rcs(self, name: str, order: int, expect_success: bool) -> tuple[bool, str]:
        csv_path = self.output / (name + ".csv")
        command = [str(self.solver), "--task", "rcs", "--mesh", "meshes/plate.msh",
            "--ka", "0.02", "--a", "1", "--parallel", "1", "--rcs-dual",
            "--phi-start", "0", "--phi-end", "90", "--phi-step", "90",
            "--rcs-order", str(order), "--output", str(csv_path)]
        result = self.process(name, command)
        if not expect_success:
            diagnostic = (result.stdout + result.stderr).lower()
            ok = result.returncode > 0 and ("order" in diagnostic or "quadrature" in diagnostic)
            return ok, f"unsupported order={order}; exit={result.returncode}; expected diagnostic and nonzero exit"
        if result.returncode != 0:
            return False, f"valid order={order}; exit={result.returncode}"
        with csv_path.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        values = [float(row["sigma_m2"]) for row in rows]
        db_values = [float(row["sigma_dbsm"]) for row in rows]
        ok = len(rows) == 2 and all(math.isfinite(x) and x > 0 for x in values)
        ok = ok and all(math.isfinite(x) for x in db_values)
        ok = ok and "ERROR:" not in result.stderr
        return ok, f"valid RCS: rows={len(rows)}, sigma={values}"

    def config_case(self, name: str, filename: str) -> tuple[bool, str]:
        arguments = shlex.split((ROOT / filename).read_text(encoding="utf-8"), comments=True)
        result = self.process(name, [str(self.config), *arguments])
        errors = [line for line in result.stderr.splitlines() if "ERROR" in line]
        return result.returncode == 0, errors[0] if errors else f"{filename}: exit={result.returncode}"

    def file_config(self):
        result = self.process("config_file_cli", [str(self.solver), "--config",
                              "tests/fixtures/file_mode.cfg"])
        ok = result.returncode == 0 and "Done." in result.stdout and "ERROR:" not in result.stderr
        return ok, f"PeecSolver --config: exit={result.returncode}"

    def cli_tests(self) -> None:
        self.check("config_file_cli", self.file_config)
        for name, frequency in [("harmonic_1MHz", "1e6"), ("harmonic_1Hz", "1"),
                                ("harmonic_1kHz", "1000")]:
            self.check(name, lambda name=name, frequency=frequency: self.harmonic(name, frequency))
        self.check("rcs_valid_order", lambda: self.rcs("rcs_valid_order", 4, True))
        self.check("rcs_reject_order7", lambda: self.rcs("rcs_reject_order7", 7, False))
        for name, filename in [("case_scattering", "cases/run_scattering.cfg"),
                               ("case_rcs", "cases/run_rcs.cfg")]:
            self.check(name, lambda name=name, filename=filename: self.config_case(name, filename))

    def asan_test(self) -> None:
        def run():
            result = self.process("config_file_lifetime", [str(self.asan)])
            combined = result.stdout + result.stderr
            ok = result.returncode == 0 and "PASS config_file_lifetime" in result.stdout
            ok = ok and "runtime error:" not in combined and "ERROR: AddressSanitizer" not in combined
            if "heap-use-after-free" in combined:
                detail = "AddressSanitizer detected heap-use-after-free after config parsing."
            else:
                detail = f"sanitizer probe exit={result.returncode} (must be 0)"
            return ok, detail
        self.check("config_file_lifetime", run)

    def spice_tests(self) -> None:
        # Import only for this explicitly requested group. Missing dependencies fail.
        try:
            import numpy as np
            from scipy.io import mmread
        except ImportError as exc:
            self.check("spice_dependencies", lambda: (False, f"NumPy and SciPy required: {exc}"))
            return

        matrices = self.output / "matrices"
        name = "spice_assemble"
        result = self.process(name, [str(self.solver), "--task", "scattering",
            "--mesh", "meshes/plate.msh", "--frequency", "1e6", "--parallel", "1",
            "--save-matrices", "--matrix-dir", str(matrices)])
        self.check(name, lambda: (result.returncode == 0, f"fresh matrix assembly exit={result.returncode}"))
        if result.returncode != 0:
            return

        def read_matrix(filename):
            value = mmread(matrices / filename)
            return np.asarray(value.toarray() if hasattr(value, "toarray") else value, dtype=float)

        l, p, a, r = [read_matrix(filename) for filename in ["L.mtx", "P.mtx", "A.mtx", "R.mtx"]]
        spec = importlib.util.spec_from_file_location("peec_export_under_test", ROOT / "export_peec_spice.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        netlist = self.output / "model.cir"
        try:
            module.validate_inputs(l, p, a, r)
            module.write_netlist(output=netlist, L=l, P=p, A=a, R=r,
                strike_node=0, return_node=1, shunt_a_node=None, shunt_b_node=None,
                shunt_R=1000, K_amp=1, alpha=1500, beta=1e6, delay=0,
                tstep=1e-7, tstop=1e-6, r_floor=1e-12, cap_cutoff=0, k_cutoff=0)
        except Exception as exc:
            # Explicitly rejecting an inexact export is acceptable after a fix,
            # but require a relevant diagnostic, not any unrelated exception.
            message = str(exc).lower()
            relevant = isinstance(exc, (ValueError, RuntimeError)) and any(
                word in message for word in ("negative", "capacit", "passiv", "exact"))
            self.check("spice_exact_or_reject", lambda: (relevant, f"export rejected: {exc}"))
            return

        text = netlist.read_text(encoding="utf-8")
        def exact_export():
            reference = np.linalg.solve(p, np.eye(p.shape[0]))
            reconstructed = np.zeros_like(reference)
            for line in text.splitlines():
                parts = line.split()
                if not parts or not parts[0].upper().startswith("CPEEC"):
                    continue
                _, na, nb, value = parts[:4]
                value = float(value)
                i = None if na == "0" else int(na[1:])
                j = None if nb == "0" else int(nb[1:])
                if i is not None:
                    reconstructed[i, i] += value
                if j is not None:
                    reconstructed[j, j] += value
                if i is not None and j is not None:
                    reconstructed[i, j] -= value
                    reconstructed[j, i] -= value
            error = np.linalg.norm(reconstructed-reference, np.inf)/np.linalg.norm(reference, np.inf)
            return bool(np.isfinite(error) and error <= 1e-10), f"exported C relative infinity-norm error={error:.9g}; limit=1e-10"
        self.check("spice_exact_or_reject", exact_export)

        def command_order():
            active = [line.strip().lower() for line in text.splitlines()
                      if line.strip() and not line.lstrip().startswith("*")]
            simulations = [i for i, line in enumerate(active)
                           if line == "run" or line.startswith("tran ")]
            saves = [i for i, line in enumerate(active) if line.startswith("wrdata ")]
            ok = bool(simulations and saves and min(simulations) < min(saves))
            return ok, "simulation command must precede wrdata in the emitted control script"
        self.check("spice_run_before_wrdata", command_order)

    def finish(self) -> int:
        passed = sum(result["status"] == "PASS" for result in self.results)
        failed = len(self.results) - passed
        if not self.results:
            print("FAIL no tests selected", flush=True)
            return 2
        if self.args.only:
            missing = set(self.args.only) - {r["name"] for r in self.results}
            if missing:
                print("FAIL unknown/unavailable tests: " + ", ".join(sorted(missing)))
                failed += len(missing)
        report = {"group": self.args.group, "passed": passed, "failed": failed,
                  "results": self.results}
        path = self.output / "report.json"
        path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"\nRESULT: {passed} PASS, {failed} FAIL")
        print(f"Report: {path.relative_to(ROOT)}")
        return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=["default", "core", "cli", "asan", "spice"], default="default")
    parser.add_argument("--core")
    parser.add_argument("--solver")
    parser.add_argument("--config-parser")
    parser.add_argument("--asan")
    parser.add_argument("--only", action="append", help="Run just this test; may be repeated.")
    args = parser.parse_args()
    suite = Suite(args)
    if args.group in ("default", "core"):
        suite.core_tests()
    if args.group in ("default", "cli"):
        suite.cli_tests()
    if args.group == "asan":
        suite.asan_test()
    if args.group == "spice":
        suite.spice_tests()
    return suite.finish()


if __name__ == "__main__":
    raise SystemExit(main())
