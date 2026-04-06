# src/cancer/61_peek_manifest_header_region.py
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import zipfile
from pathlib import Path


def safe_decode(b: bytes) -> str:
    for enc in ("utf-8", "utf-8-sig", "latin-1", "cp1252"):
        try:
            return b.decode(enc)
        except Exception:
            pass
    return b.decode("utf-8", errors="replace")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--member", default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--n", type=int, default=80)
    args = ap.parse_args()

    path = Path(args.input)

    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path, "r") as zf:
            member = args.member
            if member is None:
                members = [m for m in zf.namelist() if not m.endswith("/")]
                print("[zip members]")
                for m in members:
                    print(m)
                print("\n[using first member]")
                member = members[0]
            raw = zf.read(member)
            text = safe_decode(raw)
            print(f"[member] {member}")
    else:
        raw = path.read_bytes()
        text = safe_decode(raw)

    lines = text.splitlines()
    end = min(len(lines), args.start + args.n)

    print(f"[total_lines] {len(lines)}")
    print(f"[showing] lines {args.start} to {end-1}")
    for i in range(args.start, end):
        print(f"{i:05d}: {lines[i]}")


if __name__ == "__main__":
    main()
