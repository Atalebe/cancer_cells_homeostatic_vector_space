import argparse
import json
from pathlib import Path

import pandas as pd


def infer_caption(filename: str) -> str:
    name = filename.replace(".png", "")

    caption_map = {
        "candidate_core_score_boxplot_by_population": "Candidate core score distribution by population.",
        "candidate_core_score_boxplot_by_sector": "Candidate core score distribution by branch sector.",
        "candidate_core_score_vs_H_by_population": "Candidate core score versus reserve, colored by population.",
        "candidate_core_score_vs_H_by_sector": "Candidate core score versus reserve, colored by branch sector.",
        "candidate_core_score_vs_M_by_population": "Candidate core score versus malignant commitment, colored by population.",
        "candidate_core_score_vs_M_by_sector": "Candidate core score versus malignant commitment, colored by branch sector.",
        "candidate_core_score_vs_R_by_population": "Candidate core score versus refined recoverability, colored by population.",
        "candidate_core_score_vs_R_by_sector": "Candidate core score versus refined recoverability, colored by branch sector.",
        "candidate_core_score_vs_S_by_population": "Candidate core score versus refined structural stability, colored by population.",
        "candidate_core_score_vs_S_by_sector": "Candidate core score versus refined structural stability, colored by branch sector.",
        "candidate_core_heatmap_by_population": "Candidate core gene heatmap across populations.",
        "candidate_core_heatmap_by_population_rowz": "Row-standardized candidate core gene heatmap across populations.",
        "candidate_core_heatmap_by_sector": "Candidate core gene heatmap across branch sectors.",
        "candidate_core_heatmap_by_sector_rowz": "Row-standardized candidate core gene heatmap across branch sectors.",
        "hrsm_H_vs_R_refined": "Refined HRSM geometry, reserve versus recoverability.",
        "hrsm_M_vs_R_refined": "Refined HRSM geometry, commitment versus recoverability.",
        "hrsm_S_vs_R_refined": "Refined HRSM geometry, structural stability versus recoverability.",
        "pca_pc1_vs_pc2_by_population": "PCA geometry colored by population.",
        "boxplot_H_by_population": "Reserve by population.",
        "boxplot_S_by_population": "Structural stability by population.",
        "boxplot_M_by_population": "Commitment by population.",
        "boxplot_R_by_population": "Recoverability by population.",
    }
    return caption_map.get(name, f"Placeholder caption for {filename}.")


def infer_section(filename: str) -> str:
    if filename.startswith("candidate_core_"):
        return "candidate_core_results"
    if filename.startswith("hrsm_") or filename.startswith("boxplot_") or filename.startswith("pca_"):
        return "refined_state_geometry"
    return "appendix_misc"


def infer_logbook_entry(filename: str) -> str:
    if filename.startswith("candidate_core_"):
        return "Entry 9"
    if filename.startswith("hrsm_") or filename.startswith("boxplot_") or filename.startswith("pca_"):
        return "Entry 5"
    return "unknown"


def make_label(dataset_id: str, filename: str) -> str:
    stem = filename.replace(".png", "")
    return f"fig:{dataset_id}_{stem}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--figures-dir", required=True)
    parser.add_argument("--output-dir", default="results/tables")
    args = parser.parse_args()

    figures_dir = Path(args.figures_dir)
    output_dir = Path(args.output_dir) / args.dataset_id
    output_dir.mkdir(parents=True, exist_ok=True)

    pngs = sorted(figures_dir.glob("*.png"))
    rows = []

    for i, path in enumerate(pngs, start=1):
        filename = path.name
        rows.append(
            {
                "figure_number": i,
                "dataset_id": args.dataset_id,
                "png_file": str(path),
                "png_filename": filename,
                "figure_key": filename.replace(".png", ""),
                "latex_label": make_label(args.dataset_id, filename),
                "caption_draft": infer_caption(filename),
                "result_section": infer_section(filename),
                "logbook_entry": infer_logbook_entry(filename),
            }
        )

    df = pd.DataFrame(rows)

    csv_path = output_dir / "appendix_figure_map.csv"
    json_path = output_dir / "appendix_figure_map.json"
    df.to_csv(csv_path, index=False)

    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)

    summary = {
        "dataset_id": args.dataset_id,
        "figures_dir": str(figures_dir),
        "n_figures": int(len(rows)),
        "output_csv": str(csv_path),
        "output_json": str(json_path),
    }

    with open(output_dir / "appendix_figure_map_summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
