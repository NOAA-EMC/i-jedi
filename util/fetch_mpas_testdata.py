#!/usr/bin/env python3
"""Fetch and stage MPAS test data for ijedi geometry tests.

This runs at CMake configure time. It sparse-clones the required path from
mpas-jedi-data and stages only required files into the test working dir.
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys


REPO_URL = "https://github.com/JCSDA-internal/mpas-jedi-data.git"
BRANCH = "develop"
SPARSE_PATH = "testinput_tier_1/480km/bg"


def run_cmd(cmd: list[str]) -> None:
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed (rc={result.returncode}): {' '.join(cmd)}")


def ensure_symlink(src: pathlib.Path, dst: pathlib.Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        return
    dst.symlink_to(src)
    print(f"Staged (symlink): {dst}")


def ensure_clone(clone_dir: pathlib.Path) -> None:
    git_dir = clone_dir / ".git"
    if git_dir.exists():
        print(f"mpas-jedi-data clone already present at {clone_dir}, skipping clone.")
        return

    if clone_dir.exists() and any(clone_dir.iterdir()):
        raise RuntimeError(f"CLONE_DIR exists but is not a git checkout: {clone_dir}")

    if clone_dir.exists():
        shutil.rmtree(clone_dir)

    print("Cloning mpas-jedi-data (sparse, no blobs)...")
    run_cmd([
        "git",
        "clone",
        "--depth",
        "1",
        "--branch",
        BRANCH,
        "--filter=blob:none",
        "--no-checkout",
        REPO_URL,
        str(clone_dir),
    ])
    run_cmd(["git", "-C", str(clone_dir), "sparse-checkout", "init", "--cone"])
    run_cmd(["git", "-C", str(clone_dir), "sparse-checkout", "set", SPARSE_PATH])
    run_cmd(["git", "-C", str(clone_dir), "checkout", BRANCH])

    print(f"Fetching LFS objects for {SPARSE_PATH}...")
    run_cmd(["git", "-C", str(clone_dir), "lfs", "pull", f"--include={SPARSE_PATH}/**"])


def write_partition_file(part_file: pathlib.Path) -> None:
    if part_file.exists():
        print(f"Already present: {part_file}")
        return

    print(f"Generating block partition: {part_file}")
    n_cells = 2562
    n_parts = 6
    per_part = 427

    part_file.parent.mkdir(parents=True, exist_ok=True)
    with part_file.open("w", encoding="ascii") as f:
        for cell in range(1, n_cells + 1):
            part = (cell - 1) // per_part
            if part >= n_parts:
                part = n_parts - 1
            f.write(f"{part}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clone-dir", required=True)
    parser.add_argument("--dest-dir", required=True)
    parser.add_argument("--data-link-source", required=False)
    args = parser.parse_args()

    clone_dir = pathlib.Path(args.clone_dir)
    dest_dir = pathlib.Path(args.dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    if args.data_link_source:
        data_link_source = pathlib.Path(args.data_link_source)
        if not data_link_source.exists():
            raise RuntimeError(f"DATA_LINK_SOURCE does not exist: {data_link_source}")
        ensure_symlink(data_link_source, dest_dir / "data")

    ensure_clone(clone_dir)

    src = clone_dir / "testinput_tier_1" / "480km" / "bg" / "restart.2018-04-15_00.00.00.nc"
    if not src.exists():
        raise RuntimeError(f"Expected file not found after clone: {src}")

    dst = dest_dir / "Data" / "480km" / "bg" / "restart.2018-04-15_00.00.00.nc"
    ensure_symlink(src, dst)

    write_partition_file(dest_dir / "x1.2562.graph.info.part.6")
    print(f"MPAS test data ready in {dest_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # pragma: no cover
        print(str(exc), file=sys.stderr)
        sys.exit(1)
