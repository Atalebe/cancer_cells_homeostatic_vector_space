# src/cancer/60_inspect_manual_illumina_manifest_inputs.py
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import pandas as pd


def safe_decode(b: bytes) -> str:
    for enc in ("utf-8", "utf-8-sig", "latin-1", "cp1252"):
        try:
            return b.decode(enc)
        except Exception:
            pass
    return b.decode("utf-8", errors="replace")


def inspect_plain_text(path: Path, n_lines: int = 120) -> dict:
    raw = path.read_bytes()
    text = safe_decode(raw)
    lines = text.splitlines()

    head = lines[:n_lines]
    nonempty = [x for x in head if x.strip()]

    candidate_header_idx = None
    for i, line in enumerate(head):
        low = line.lower()
        if any(
            key in low
            for key in [
                "ilmnid",
                "name,",
                "name\t",
                "mapinfo",
                "chromosome",
                "chr,",
                "ucsc_refgene_name",
                "relation_to_ucsc_cpg_island",
            ]
        ):
            candidate_header_idx = i
            break

    sep_votes = {
        "comma": sum(line.count(",") for line in nonempty[:40]),
        "tab": sum(line.count("\t") for line in nonempty[:40]),
        "semicolon": sum(line.count(";") for line in nonempty[:40]),
        "pipe": sum(line.count("|") for line in nonempty[:40]),
    }

    return {
        "path": str(path),
        "kind": "plain_text",
        "size_bytes": path.stat().st_size,
        "n_preview_lines": len(head),
        "candidate_header_idx": candidate_header_idx,
        "separator_votes": sep_votes,
        "preview_lines": head,
    }


def inspect_zip(path: Path, n_lines: int = 80) -> dict:
    out = {
        "path": str(path),
        "kind": "zip",
        "size_bytes": path.stat().st_size,
        "members": [],
    }

    with zipfile.ZipFile(path, "r") as zf:
        for member in zf.namelist():
            rec = {"member": member}
            if member.endswith("/"):
                rec["is_dir"] = True
                out["members"].append(rec)
                continue

            info = zf.getinfo(member)
            rec["is_dir"] = False
            rec["file_size"] = info.file_size
            rec["compress_size"] = info.compress_size

            try:
                raw = zf.read(member)
                text = safe_decode(raw[:400000])
                lines = text.splitlines()[:n_lines]

                candidate_header_idx = None
                for i, line in enumerate(lines):
                    low = line.lower()
                    if any(
                        key in low
                        for key in [
                            "ilmnid",
                            "name,",
                            "name\t",
                            "mapinfo",
                            "chromosome",
                            "chr,",
                            "ucsc_refgene_name",
                            "relation_to_ucsc_cpg_island",
                        ]
                    ):
                        candidate_header_idx = i
                        break

                rec["candidate_header_idx"] = candidate_header_idx
                rec["preview_lines"] = lines[:30]
            except Exception as e:
                rec["read_error"] = repr(e)

            out["members"].append(rec)

    return out


def inspect_csv_guess(path: Path, sep: str, skiprows: int = 0) -> dict:
    try:
        df = pd.read_csv(path, sep=sep, skiprows=skiprows, nrows=5, engine="python")
        return {
            "ok": True,
            "sep": repr(sep),
            "skiprows": skiprows,
            "columns": [str(c) for c in df.columns],
            "n_cols": int(df.shape[1]),
            "preview_rows": df.head(3).astype(str).to_dict(orient="records"),
        }
    except Exception as e:
        return {
            "ok": False,
            "sep": repr(sep),
            "skiprows": skiprows,
            "error": repr(e),
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = []

    for p in args.inputs:
        path = Path(p)
        if not path.exists():
            summary.append({"path": str(path), "exists": False})
            continue

        rec = {"path": str(path), "exists": True}

        if zipfile.is_zipfile(path):
            rec["inspection"] = inspect_zip(path)
        else:
            rec["inspection"] = inspect_plain_text(path)

            suffix = "".join(path.suffixes).lower()
            if suffix.endswith(".csv") or suffix.endswith(".txt") or suffix.endswith(".tsv"):
                csv_trials = []
                for skip in [0, 1, 5, 10, 20, 30, 40, 50]:
                    for sep in [",", "\t", ";", "|"]:
                        csv_trials.append(inspect_csv_guess(path, sep=sep, skiprows=skip))
                rec["csv_trials"] = csv_trials

        summary.append(rec)

    with open(outdir / "manual_manifest_input_inspection.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    rows = []
    for rec in summary:
        base = {
            "path": rec["path"],
            "exists": rec["exists"],
        }
        if not rec["exists"]:
            rows.append(base)
            continue
        insp = rec["inspection"]
        base["kind"] = insp.get("kind")
        base["size_bytes"] = insp.get("size_bytes")
        if insp.get("kind") == "plain_text":
            base["candidate_header_idx"] = insp.get("candidate_header_idx")
            base["separator_votes"] = json.dumps(insp.get("separator_votes", {}))
        else:
            base["n_members"] = len(insp.get("members", []))
        rows.append(base)

    pd.DataFrame(rows).to_csv(outdir / "manual_manifest_input_inventory.csv", index=False)

    print("[ok] wrote inspection json:", outdir / "manual_manifest_input_inspection.json")
    print("[ok] wrote inventory csv:", outdir / "manual_manifest_input_inventory.csv")
    print("[info] inspect these:")
    print(" ", outdir / "manual_manifest_input_inventory.csv")
    print(" ", outdir / "manual_manifest_input_inspection.json")


if __name__ == "__main__":
    main()
