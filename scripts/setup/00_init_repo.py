from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    marker = ROOT / ".repo_initialized"
    marker.write_text("initialized\n", encoding="utf-8")
    print(f"[ok] repo marker written: {marker}")


if __name__ == "__main__":
    main()
