#!/usr/bin/env python3
"""Dependency-free unit and CLI smoke tests. Run from any working directory."""
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    subprocess.run(["make", "-C", str(ROOT)], check=True)
    with tempfile.TemporaryDirectory(prefix="containment-test-") as tmp:
        tmp = Path(tmp)
        unit = tmp / "unit"
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) + [
            "-O2", "-std=c++20", "-I" + str(ROOT / "include"),
            str(ROOT / "test/test_containment.cpp"), "-o", str(unit)
        ], check=True)
        subprocess.run([str(unit)], check=True)
        output = tmp / "results.tsv"

        def run(algo="alphamaxgeom", **options):
            flags = dict(t=0.5, metric="containment", seeds=200, steps=3,
                         base_n=1024, size_multiplier=40, algo=algo, alpha=0.45,
                         k=2000, scale=0.1, out=str(output))
            flags.update(options)
            args = [str(ROOT / "bin/expt_growth")]
            for key, value in flags.items():
                args.extend(["--" + key, str(value)])
            result = subprocess.run(args, capture_output=True, text=True, check=True)
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            return rows, result

        for algo in ("alphamaxgeom", "maxgeom", "bottomk", "fracminhash"):
            rows, _ = run(algo)
            assert len(rows) == 3
            for step, row in enumerate(rows):
                assert row["metric"] == "containment"
                assert int(row["|A|"]) == 1024 * 2**step
                assert int(row["|B|"]) == 40 * int(row["|A|"])
                assert float(row["true_sim"]) == 0.5
                assert int(row["valid_seeds"]) >= 190
                assert abs(float(row["mean_est"]) - 0.5) < 0.07, (algo, row)
                assert 0 <= float(row["mse"]) < 0.09, (algo, row)
            print(algo, [(r["mean_est"], r["mse"], r["valid_seeds"]) for r in rows])
            # Exact endpoints, reverse size ratio, and x10 growth.
            for target, multiplier in ((0, 40), (1, 40), (0.25, 0.5)):
                rows, _ = run(algo, t=target, size_multiplier=multiplier,
                              base_n=128, seeds=60, steps=2, growth="x10")
                for row in rows:
                    assert float(row["true_sim"]) == target
                    assert int(row["|B|"]) == multiplier * int(row["|A|"])
                    assert abs(float(row["mean_est"]) - target) < 0.08
            # Integer rounding is reflected in true_sim, not hidden by --t.
            rows, _ = run(algo, base_n=101, t=0.333, size_multiplier=1.5,
                          seeds=5, steps=1)
            assert (int(rows[0]["|A|"]), int(rows[0]["|B|"])) == (101, 152)
            assert abs(float(rows[0]["true_sim"]) - 34 / 101) < 1e-6

        rows, _ = run("fracminhash", base_n=101, size_multiplier=3,
                      seeds=5, steps=1, scale=1)
        assert float(rows[0]["mse"]) == 0
        rows, result = run("fracminhash", scale=1e-30, base_n=10,
                           seeds=5, steps=1)
        assert rows[0]["valid_seeds"] == "0"
        assert math.isnan(float(rows[0]["mean_est"])) and math.isnan(float(rows[0]["mse"]))
        assert "WARNING" in result.stderr
        rows, result = run("bottomk", k=1, seeds=100, steps=1)
        assert 0 < int(rows[0]["valid_seeds"]) < 100
        assert "WARNING" in result.stderr
        first, _ = run("bottomk", seeds=20, steps=1)
        assert first == run("bottomk", seeds=20, steps=1)[0]
        assert first != run("bottomk", seeds=20, steps=1, seed=43)[0]

        for options in (dict(size_multiplier=0), dict(size_multiplier=-1),
                        dict(size_multiplier="nan"), dict(size_multiplier="inf"),
                        dict(size_multiplier="40oops"), dict(t=1.1), dict(t="nan"),
                        dict(t=0.75, size_multiplier=0.5), dict(base_n=0),
                        dict(base_n=-1), dict(steps=0), dict(seeds=0),
                        dict(steps=64), dict(size_multiplier=1e308),
                        dict(base_n=10, size_multiplier=0.001, t=0),
                        dict(metric="jaccard"), dict(algo="minhash"),
                        dict(algo="bad"), dict(growth="x3")):
            try:
                run(**options)
            except subprocess.CalledProcessError as err:
                assert err.returncode > 0 and "ERROR:" in err.stderr, (options, err)
            else:
                raise AssertionError(f"Invalid arguments accepted: {options}")
        print("Containment CLI smoke tests passed")


if __name__ == "__main__":
    main()
