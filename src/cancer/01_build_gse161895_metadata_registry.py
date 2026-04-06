#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
import re
import tarfile
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml
from xml.etree import ElementTree as ET


REPO_ROOT = Path.home() / "cancer_cells_hrsm_sims"
CONFIG_PATH = REPO_ROOT / "configs" / "gse161895.yaml"


def read_config() -> dict:
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_soft(soft_gz: Path) -> pd.DataFrame:
    rows = []
    current = None

    with gzip.open(soft_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                gsm = line.split("=", 1)[1].strip()
                current = {"gsm": gsm}
            elif current and line.startswith("!Sample_"):
                key, val = line.split("=", 1)
                key = key.replace("!Sample_", "").strip()
                val = val.strip()
                if key in current:
                    current[key] = f"{current[key]} | {val}"
                else:
                    current[key] = val

    if current:
        rows.append(current)

    return pd.DataFrame(rows)


def parse_miniml(miniml_tgz: Path) -> pd.DataFrame:
    rows = []
    with tarfile.open(miniml_tgz, "r:gz") as tar:
        xml_members = [m for m in tar.getmembers() if m.name.endswith(".xml")]
        if not xml_members:
            return pd.DataFrame()
        xml_file = tar.extractfile(xml_members[0])
        if xml_file is None:
            return pd.DataFrame()
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for sample in root.findall(".//{*}Sample"):
            row = {}
            iid = sample.attrib.get("iid")
            if iid:
                row["gsm"] = iid
            title = sample.findtext("{*}Title")
            if title:
                row["title"] = title

            for ch in sample.findall("{*}Channel"):
                for char in ch.findall("{*}Characteristics"):
                    tag = char.attrib.get("tag", "characteristics")
                    text = (char.text or "").strip()
                    if text:
                        if tag in row:
                            row[tag] = f"{row[tag]} | {text}"
                        else:
                            row[tag] = text

            for s in sample.findall("{*}Supplementary-Data"):
                text = (s.text or "").strip()
                if text:
                    row["supplementary_data"] = f"{row.get('supplementary_data', '')} | {text}".strip(" |")

            if row:
                rows.append(row)

    return pd.DataFrame(rows)


def extract_series_matrix_headers(matrix_gz: Path) -> pd.DataFrame:
    meta = {}
    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_"):
                key, rest = line.split("\t", 1)
                key = key.replace("!Sample_", "").strip()
                fields = [x.strip().strip('"') for x in rest.rstrip("\n").split("\t")]
                meta[key] = fields

    if not meta:
        return pd.DataFrame()

    # infer width from first field
    keys = list(meta.keys())
    n = len(meta[keys[0]])
    rows = []
    for i in range(n):
        row = {}
        for k, vals in meta.items():
            row[k] = vals[i] if i < len(vals) else None
        rows.append(row)

    return pd.DataFrame(rows)


def harmonize_registry(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    out = df.copy()
    out.columns = [re.sub(r"[^a-zA-Z0-9]+", "_", c.strip()).strip("_").lower() for c in out.columns]

    if "gsm" not in out.columns:
        for c in out.columns:
            if "sample" == c or "geo_accession" == c:
                out = out.rename(columns={c: "gsm"})
                break

    # heuristic donor/patient/cell inference
    source_cols = [c for c in out.columns if c != "gsm"]

    def pick(*patterns: str, text: str) -> str | None:
        for p in patterns:
            m = re.search(p, text, flags=re.I)
            if m:
                return m.group(1)
        return None

    patient_ids = []
    donor_types = []
    treatment_states = []

    for _, row in out.iterrows():
        blob = " | ".join(str(row[c]) for c in source_cols if pd.notna(row[c]))
        patient = pick(r"(Patient\d+)", r"(NormalDonor\d+)", text=blob) or pick(r"(ETP[_ -]?\d+)", text=blob)
        patient_ids.append(patient)

        donor_type = None
        if re.search(r"normal", blob, flags=re.I):
            donor_type = "normal_donor"
        elif re.search(r"patient", blob, flags=re.I):
            donor_type = "patient"
        donor_types.append(donor_type)

        treat = None
        if re.search(r"untreated", blob, flags=re.I):
            treat = "untreated"
        elif re.search(r"treated", blob, flags=re.I):
            treat = "treated"
        treatment_states.append(treat)

    out["patient_or_donor_id"] = patient_ids
    out["source_class"] = donor_types
    out["treatment_state"] = treatment_states

    return out


def main() -> None:
    cfg = read_config()
    raw_dir = REPO_ROOT / cfg["raw_dir"]
    meta_dir = REPO_ROOT / cfg["metadata_dir"]
    mkdir(meta_dir)

    soft = raw_dir / "GSE161895_family.soft.gz"
    miniml = raw_dir / "GSE161895_family.xml.tgz"
    matrix = raw_dir / "GSE161895_series_matrix.txt.gz"
    super_matrix = raw_dir / "GSE161901_series_matrix.txt.gz"

    candidates = []
    if soft.exists():
        df = parse_soft(soft)
        if not df.empty:
            candidates.append(("soft", df))
    if miniml.exists():
        df = parse_miniml(miniml)
        if not df.empty:
            candidates.append(("miniml", df))
    if matrix.exists():
        df = extract_series_matrix_headers(matrix)
        if not df.empty:
            candidates.append(("matrix", df))
    if super_matrix.exists():
        df = extract_series_matrix_headers(super_matrix)
        if not df.empty:
            candidates.append(("super_matrix", df))

    if not candidates:
        raise RuntimeError("No metadata source could be parsed for GSE161895.")

    source_name, registry = max(candidates, key=lambda x: len(x[1].columns))
    registry = harmonize_registry(registry)

    out_csv = meta_dir / "cell_metadata_registry.csv"
    out_parquet = meta_dir / "cell_metadata_registry.parquet"
    out_json = meta_dir / "cell_metadata_registry_summary.json"

    registry.to_csv(out_csv, index=False)
    registry.to_parquet(out_parquet, index=False)

    summary = {
        "source_used": source_name,
        "n_rows": int(len(registry)),
        "n_columns": int(len(registry.columns)),
        "columns": registry.columns.tolist(),
        "source_class_counts": registry["source_class"].value_counts(dropna=False).to_dict() if "source_class" in registry.columns else {},
        "treatment_state_counts": registry["treatment_state"].value_counts(dropna=False).to_dict() if "treatment_state" in registry.columns else {},
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[ok] wrote {out_csv}")
    print(f"[ok] wrote {out_parquet}")
    print(f"[ok] wrote {out_json}")


if __name__ == "__main__":
    main()
