from pathlib import Path
from src.common.registry_utils import write_registry

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    out = ROOT / "data" / "registries" / "dataset_registry.csv"
    write_registry(out)
    print(f"[ok] wrote empty dataset registry to {out}")


if __name__ == "__main__":
    main()
