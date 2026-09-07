#!/usr/bin/env python3
"""
compare_rmsd_dirs.py — run `foldcomp rmsd` on every matching filename between
two directories and summarize backbone/all-atom RMSD statistics.

Usage:
    compare_rmsd_dirs.py <dir1> <dir2> [--csv OUT.csv] [--foldcomp PATH]
"""

import argparse
import csv
import statistics
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir1", type=Path, help="First directory of structure files")
    parser.add_argument("dir2", type=Path, help="Second directory of structure files")
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("rmsd_comparison.csv"),
        help="Output CSV path (default: rmsd_comparison.csv)",
    )
    parser.add_argument(
        "--foldcomp",
        type=Path,
        default=Path("build_gpu/foldcomp"),
        help="Path to the foldcomp binary (default: build_gpu/foldcomp)",
    )
    args = parser.parse_args()

    for d in (args.dir1, args.dir2):
        if not d.is_dir():
            sys.exit(f"error: {d} is not a directory")
    if not args.foldcomp.is_file():
        sys.exit(f"error: foldcomp binary not found at {args.foldcomp}")

    files1 = {p.name: p for p in args.dir1.iterdir() if p.is_file()}
    files2 = {p.name: p for p in args.dir2.iterdir() if p.is_file()}
    common = sorted(set(files1) & set(files2))

    if not common:
        sys.exit(f"error: no matching filenames between {args.dir1} and {args.dir2}")

    rows = []
    errors = []
    for name in common:
        f1, f2 = files1[name], files2[name]
        try:
            result = subprocess.run(
                [str(args.foldcomp), "rmsd", str(f1), str(f2)],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            errors.append((name, e.stderr.strip() or e.stdout.strip()))
            continue
        fields = result.stdout.strip().split("\t")
        if len(fields) != 6:
            errors.append((name, f"unexpected output: {result.stdout.strip()!r}"))
            continue
        _, _, n_residue, n_atom, backbone_rmsd, all_atom_rmsd = fields
        rows.append(
            {
                "file1": str(f1),
                "file2": str(f2),
                "n_residue": int(n_residue),
                "n_atom": int(n_atom),
                "backbone_rmsd": float(backbone_rmsd),
                "all_atom_rmsd": float(all_atom_rmsd),
            }
        )

    if not rows:
        sys.exit("error: no successful comparisons")

    with args.csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "file1",
                "file2",
                "n_residue",
                "n_atom",
                "backbone_rmsd",
                "all_atom_rmsd",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    def summarize(key):
        values = [r[key] for r in rows]
        return {
            "min": min(values),
            "max": max(values),
            "mean": statistics.mean(values),
            "median": statistics.median(values),
        }

    print(f"Compared {len(rows)} file pairs ({len(errors)} failed).")
    print(f"Results saved to {args.csv}")
    for key, label in (
        ("backbone_rmsd", "Backbone RMSD"),
        ("all_atom_rmsd", "All-atom RMSD"),
    ):
        s = summarize(key)
        print(f"\n{label}:")
        print(f"  min:    {s['min']:.6f}")
        print(f"  max:    {s['max']:.6f}")
        print(f"  mean:   {s['mean']:.6f}")
        print(f"  median: {s['median']:.6f}")

    if errors:
        print(f"\n{len(errors)} file(s) failed comparison:", file=sys.stderr)
        for name, msg in errors:
            print(f"  {name}: {msg}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
