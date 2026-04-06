from __future__ import annotations

from pathlib import Path
import pandas as pd


def load_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


TARGETS = ["SAMPLE 24", "SAMPLE 35", "SAMPLE 341"]


def main() -> None:
    outdir = Path("results/gse240704/condition_axis_tests_curated")
    outdir.mkdir(parents=True, exist_ok=True)

    registry = pd.read_csv("data/metadata/gse240704/sample_registry.csv")
    ann_cur = load_table("data/processed/gse240704/sample_annotations_curated.parquet")
    ann_raw = load_table("data/processed/gse240704/sample_annotations.parquet")

    files = [
        "results/gse240704/metadata_inspection/selected_sample_metadata_table.csv",
        "results/gse240704/metadata_inspection/matched_sample_metadata_table.csv",
        "results/gse240704/inspection/annotations/gse240704_metadata_mapping_template.csv",
    ]

    outputs = []

    for f in files:
        path = Path(f)
        if not path.exists():
            continue
        df = pd.read_csv(path)
        candidate_cols = [c for c in df.columns if "sample" in c.lower() or "mgmt" in c.lower() or "tumor" in c.lower() or "condition" in c.lower()]
        candidate_cols = list(dict.fromkeys(candidate_cols))
        subset = pd.DataFrame()
        matched = False

        for col in [c for c in df.columns if "sample" in c.lower()]:
            s = df[col].astype(str)
            hit = df[s.isin(TARGETS)].copy()
            if len(hit):
                subset = hit[candidate_cols] if candidate_cols else hit
                matched = True
                break

        if not matched and "sample_id" in df.columns:
            hit = df[df["sample_id"].astype(str).isin(TARGETS)].copy()
            if len(hit):
                subset = hit[candidate_cols] if candidate_cols else hit

        if len(subset):
            out = outdir / f"{path.stem}_target_samples.csv"
            subset.to_csv(out, index=False)
            outputs.append(str(out))

    reg_sub = registry[registry["sample_id"].astype(str).isin(TARGETS)].copy()
    reg_sub.to_csv(outdir / "sample_registry_target_samples.csv", index=False)

    raw_sub = ann_raw[ann_raw["sample_id"].astype(str).isin(TARGETS)].copy()
    raw_sub.to_csv(outdir / "sample_annotations_raw_target_samples.csv", index=False)

    cur_sub = ann_cur[ann_cur["sample_id"].astype(str).isin(TARGETS)].copy()
    cur_sub.to_csv(outdir / "sample_annotations_curated_target_samples.csv", index=False)

    print("[ok] wrote target-sample audit tables")
    print("[info] inspect:")
    print("  results/gse240704/condition_axis_tests_curated/sample_registry_target_samples.csv")
    print("  results/gse240704/condition_axis_tests_curated/sample_annotations_raw_target_samples.csv")
    print("  results/gse240704/condition_axis_tests_curated/sample_annotations_curated_target_samples.csv")
    for o in outputs:
        print(" ", o)


if __name__ == "__main__":
    main()
