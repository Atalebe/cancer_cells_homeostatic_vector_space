from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass
class RunMetadata:
    script_name: str
    dataset_id: str | None
    branch: str | None
    status: str
    parameters: dict[str, Any]
    input_files: list[str]
    output_files: list[str]
    notes: str = ""


def utc_now_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def write_text(path: str | Path, text: str) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def start_run_log(base_dir: str | Path, script_name: str) -> Path:
    run_dir = ensure_dir(Path(base_dir) / script_name / utc_now_stamp())
    return run_dir


def finalize_run(run_dir: str | Path, metadata: RunMetadata) -> None:
    run_dir = Path(run_dir)
    write_json(run_dir / "run_meta.json", asdict(metadata))
    write_json(run_dir / "outputs_manifest.json", {"output_files": metadata.output_files})
