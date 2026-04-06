from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    out = ROOT / "logbook" / "run_ledger.csv"
    df = pd.DataFrame(columns=[
        "timestamp_utc",
        "dataset_id",
        "branch",
        "script_name",
        "status",
        "major_result",
        "next_action",
        "notes",
    ])
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"[ok] wrote run ledger template to {out}")


if __name__ == "__main__":
    main()
