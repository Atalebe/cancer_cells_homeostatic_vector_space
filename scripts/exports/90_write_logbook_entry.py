from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]

TEMPLATE = r"""\section*{Entry %(entry_number)s: %(title)s}
\addcontentsline{toc}{section}{Entry %(entry_number)s: %(title)s}

\subsection*{What is being tested}
%(what_tested)s

\subsection*{Dataset}
%(dataset)s

\subsection*{Methods}
%(methods)s

\subsection*{Files produced}
%(files_produced)s

\subsection*{Results}
%(results)s

\subsection*{Interpretation}
%(interpretation)s

\subsection*{Next steps}
%(next_steps)s

\subsection*{Tables}
%(tables)s

\subsection*{Figures}
%(figures)s
"""


def main() -> None:
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = ROOT / "logbook" / "entries" / f"entry_template_{ts}.tex"
    payload = {
        "entry_number": "X",
        "title": "Milestone title",
        "what_tested": "Describe the exact branch or diagnostic under test.",
        "dataset": "List dataset accession, branch definition, and key metadata scope.",
        "methods": "Summarize realization, filters, proxies, geometry, overlays, and validation used.",
        "files_produced": "List the main tables, figures, and manifests produced.",
        "results": "State the clean results, not the wishful ones.",
        "interpretation": "Explain what the results mean and what they do not mean.",
        "next_steps": "State the immediate next tests.",
        "tables": "List tables to include or update.",
        "figures": "List figures to include or update.",
    }
    out.write_text(TEMPLATE % payload, encoding="utf-8")
    print(f"[ok] wrote logbook entry template to {out}")


if __name__ == "__main__":
    main()
