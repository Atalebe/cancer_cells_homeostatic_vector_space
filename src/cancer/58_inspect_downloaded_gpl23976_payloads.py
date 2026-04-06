#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    args = ap.parse_args()

    indir = Path(args.indir)
    print("[info] files:")
    for p in sorted(indir.glob("*")):
        if p.is_file():
            print(f"  {p}  ({p.stat().st_size} bytes)")
            try:
                head = p.read_text(encoding="utf-8", errors="replace")[:500]
                print("----- preview start -----")
                print(head.replace("\x00", " "))
                print("----- preview end -----")
            except Exception as e:
                print(f"[warn] preview failed: {e}")


if __name__ == "__main__":
    main()
