import argparse
import json

from src.cancer.downloaders import dry_run_download_plan, execute_geo_series_download


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--download-files", action="store_true")
    args = parser.parse_args()

    if args.dry_run:
        plan = dry_run_download_plan(args.config)
        print(json.dumps(plan, indent=2))
        return

    manifest = execute_geo_series_download(
        args.config,
        download_files=args.download_files,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
