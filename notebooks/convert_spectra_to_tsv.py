# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "simple_parsing",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    import marimo as mo
    from typing import List
    import os
    import re

    @dataclass
    class Settings:
        input_mgf: str = field(
            default="data/multims2_spectra.mgf",
            metadata={"help": "Path to input MGF file."},
        )
        output_tsv: str = field(
            default="data/MULTIMS2.tsv",
            metadata={"help": "Path to output TSV."},
        )
        dry_run: bool = field(
            default=False,
            metadata={"help": "If True, do not write output file."},
        )
        change_mzml_to_mzxml: bool = field(
            default=False,
            metadata={"help": "If True, change .mzml extension to .mzxml in FILENAME."},
        )
        split: int = field(
            default=0,
            metadata={
                "help": "Split the output TSV into multiple files with up to this many data lines each."
            },
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        if mo.running_in_notebook():
            return Settings()
        else:
            return parser.parse_args().settings

    settings = parse_args()


@app.function
def read_text_fallback(path: str) -> List[str]:
    """Read text file with encoding fallback."""
    try:
        with open(path, encoding="utf-8") as f:
            return f.readlines()
    except UnicodeDecodeError:
        with open(path, encoding="latin-1") as f:
            return f.readlines()


@app.function
def parse_mgf_file(path: str) -> list[dict]:
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


@app.function
def mgf_to_tsv(settings: Settings) -> dict:
    """Convert MGF to GNPS batch TSV format."""
    if not os.path.exists(settings.input_mgf):
        return {"error": f"Input file not found: {settings.input_mgf}"}

    blocks = parse_mgf_file(settings.input_mgf)
    if not blocks:
        return {"error": "No spectra found in MGF."}

    # TSV columns (GNPS format)
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
        "FRAGMENTATION_METHOD",
        "COLLISION_ENERGY",
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
    header = "\t".join(tsv_columns) + "\n"
    tsv_lines = []

    for blk in blocks:
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
            elif col == "FILENAME":
                val = normalized.get("FILENAME", "")
                if settings.change_mzml_to_mzxml and val.lower().endswith(".mzml"):
                    val = val[:-5] + ".mzXML"
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

    if settings.dry_run:
        return {"spectra_total": len(blocks), "output_tsv": None}

    os.makedirs(os.path.dirname(settings.output_tsv), exist_ok=True)

    base, ext = os.path.splitext(settings.output_tsv)
    total_lines = len(tsv_lines)
    split_size = settings.split
    output_files = []

    # Disable splitting if split <= 0
    if split_size <= 0 or total_lines <= split_size:
        with open(settings.output_tsv, "w", encoding="utf-8") as f:
            f.write(header)
            f.writelines(tsv_lines)
        output_files.append(settings.output_tsv)
    else:
        num_splits = (total_lines + split_size - 1) // split_size
        for i in range(num_splits):
            start = i * split_size
            end = min((i + 1) * split_size, total_lines)
            split_path = f"{base}-PARTITION-{i+1}{ext}"
            with open(split_path, "w", encoding="utf-8") as f:
                f.write(header)
                f.writelines(tsv_lines[start:end])
            output_files.append(split_path)

    return {
        "spectra_total": len(blocks),
        "output_tsv": output_files,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## Convert Spectra to TSV Settings

    - **Input MGF**: `{settings.input_mgf}`
    - **Output TSV**: `{settings.output_tsv}`
    - **Dry run**: {settings.dry_run}
    - **Change .mzml to .mzXML**: {settings.change_mzml_to_mzxml}

    Converts MGF spectra to GNPS batch upload TSV format.
    """)
    return


@app.cell
def run_conversion():
    result = mgf_to_tsv(settings)

    if "error" in result:
        mo.md(f"**Error:** {result['error']}")
    else:
        mo.md(f"""
        ### Conversion Complete
        - **Spectra processed**: {result['spectra_total']:,}
        - **Output file**: `{result['output_tsv'] or 'None (dry run)'}`
        """)
    return


if __name__ == "__main__":
    app.run()
