# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#   "marimo",
#   "simple_parsing",
# ]
# ///

import marimo

app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    from typing import List
    import os
    import re

    @dataclass
    class Settings:
        input_mgf: str = field(
            default="scratch/consolidated_spectra.mgf",
            metadata={"help": "Path to input MGF file."},
        )
        output_tsv: str = field(
            default="scratch/batch_gnps.tsv",
            metadata={"help": "Path to output TSV."},
        )
        dry_run: bool = field(
            default=False,
            metadata={"help": "If True, do not write output file; just report counts."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args():
        try:
            import marimo as mo  # type: ignore

            if mo.running_in_notebook():
                return Settings()
        except Exception:
            pass
        return parser.parse_args().settings

    settings = parse_args()


# ---------------- Internal Utilities ---------------- #
@app.function
def read_text_fallback(path: str) -> List[str]:
    try:
        with open(path, encoding="utf-8") as f:
            return f.readlines()
    except UnicodeDecodeError:
        with open(path, encoding="latin-1") as f:
            return f.readlines()


@app.function
def write_text_utf8(path: str, lines: List[str]):
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)


@app.function
def parse_mgf_file(path: str):
    """Parse MGF into list of dictionaries with fields and peaks."""
    lines = read_text_fallback(path)
    blocks = []
    cur_block = {}
    peaks = []
    inside = False

    for line in lines:
        line = line.strip()
        if line.startswith("BEGIN IONS"):
            inside = True
            cur_block = {}
            peaks = []
        elif line.startswith("END IONS"):
            inside = False
            cur_block["PEAKS"] = peaks
            blocks.append(cur_block)
        elif inside:
            if "=" in line:
                k, v = line.split("=", 1)
                cur_block[k.strip()] = v.strip()
            else:
                peaks.append(line)
    return blocks


# ---------------- Core Processing ---------------- #
@app.function
def mgf_to_tsv(settings: Settings):
    blocks = parse_mgf_file(settings.input_mgf)
    if not blocks:
        return {"error": "No spectra found in MGF."}

    # TSV Columns (fixed)
    tsv_columns = [
        "FILENAME",
        "SEQ",
        "COMPOUND_NAME",
        "MOLECULEMASS",
        "INSTRUMENT",
        "IONSOURCE",
        "EXTRACTSCAN",
        "SMILES",
        "SELFIES",
        "INCHI",
        "INCHIAUX",
        "CHARGE",
        "IONMODE",
        "ACQUISITION",
        "EXACTMASS",
        "DATACOLLECTOR",
        "DATACURATOR",
        "ADDUCT",
        "LIBQUALITY",
        "PI",
        "SPECIES",
        "CASNUMBER",
        "PUBMED",
        "STRAIN",
        "INTEREST",
        "GENUS",
    ]

    na_defaults = {"SPECIES", "CASNUMBER", "PUBMED", "STRAIN", "INTEREST", "GENUS"}

    tsv_lines = ["\t".join(tsv_columns) + "\n"]

    for blk in blocks:
        # Normalize field names first
        normalized = {}
        for k, v in blk.items():
            key = k.upper()
            if key == "DATA_COLLECTOR":
                key = "DATACOLLECTOR"
            elif key == "DATA_CURATOR":
                key = "DATACURATOR"
            normalized[key] = v

        row = []
        for col in tsv_columns:
            if col == "LIBQUALITY":
                val = "1"
            elif col == "SEQ":
                val = "*..*"
            elif col == "INSTRUMENT":
                val = normalized.get(
                    "INSTRUMENT_TYPE", normalized.get("INSTRUMENT", "")
                )
            elif col == "CHARGE":
                raw_val = normalized.get("CHARGE", "")
                m = re.search(r"-?\d+", raw_val)
                val = m.group(0) if m else ""
            elif col == "IONMODE":
                raw_val = normalized.get("IONMODE", "")
                val = raw_val.strip().capitalize() if raw_val else ""
                if val not in {"Positive", "Negative"}:
                    val = ""
            elif col in na_defaults:
                val = normalized.get(col, "N/A")
            else:
                val = normalized.get(col, "")
            row.append(val)
        tsv_lines.append("\t".join(row) + "\n")

    if not settings.dry_run:
        os.makedirs(os.path.dirname(settings.output_tsv), exist_ok=True)
        write_text_utf8(settings.output_tsv, tsv_lines)

    return {"spectra_total": len(blocks), "output_tsv": settings.output_tsv}


# ------------- CLI entry ------------- #
if __name__ == "__main__":
    result = mgf_to_tsv(settings)
    if "error" in result:
        print("ERROR:", result["error"])
    else:
        print(result)
