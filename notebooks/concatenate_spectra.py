# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "matchms",
#     "polars",
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
    import glob
    import os
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf
    import polars as pl

    @dataclass
    class Settings:
        mgf_dir: str = field(
            default="scratch/",
            metadata={"help": "Directory containing .mgf files."},
        )
        output_tsv: str = field(
            default="scratch/spectra_metadata.tsv",
            metadata={"help": "TSV output file for spectrum metadata."},
        )
        output_mgf: str = field(
            default="scratch/all_spectra.mgf",
            metadata={"help": "MGF output file for all concatenated spectra."},
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
def extract_numeric_fields(meta: dict, field_name: str, converter) -> any:
    """Extract and convert numeric fields from spectrum metadata."""
    val = meta.get(field_name)
    if val is not None and val != "":
        try:
            return converter(val)
        except Exception:
            return None
    return None


@app.function
def concatenate_spectra(settings: Settings) -> dict:
    """Concatenate all MGF files in directory and export metadata."""
    mgf_files = sorted(glob.glob(os.path.join(settings.mgf_dir, "*.mgf")))

    if not mgf_files:
        return {"error": f"No MGF files found in {settings.mgf_dir}"}

    all_spectra = []
    rows = []

    numeric_fields = [
        ("retention_time", float),
        ("feature_ms1_height", float),
        ("precursor_purity", float),
        ("num_peaks", int),
        ("quality_explained_intensity", float),
        ("quality_explained_signals", float),
    ]

    for mgf_file in mgf_files:
        print(f"Loading {mgf_file}...")
        filename = os.path.basename(mgf_file)

        for spectrum in load_from_mgf(mgf_file):
            meta = spectrum.metadata
            row = dict(meta)
            row["mgf_filename"] = filename

            for field_name, converter in numeric_fields:
                converted = extract_numeric_fields(meta, field_name, converter)
                if converted is not None:
                    row[field_name] = converted

            rows.append(row)
            all_spectra.append(spectrum)

    # Export metadata
    if rows:
        all_keys = set()
        for r in rows:
            all_keys.update(r.keys())

        columns = ["mgf_filename"] + sorted(k for k in all_keys if k != "mgf_filename")
        flat_rows = [
            {
                col: str(r.get(col, "")) if r.get(col) is not None else ""
                for col in columns
            }
            for r in rows
        ]

        df = pl.DataFrame(flat_rows)
        os.makedirs(os.path.dirname(settings.output_tsv), exist_ok=True)
        df.write_csv(settings.output_tsv, separator="\t")
        print(f"Wrote {len(rows)} spectra metadata rows to {settings.output_tsv}")

    # Export concatenated MGF
    os.makedirs(os.path.dirname(settings.output_mgf), exist_ok=True)
    save_as_mgf(all_spectra, settings.output_mgf)
    print(f"Exported {len(all_spectra)} spectra to {settings.output_mgf}")

    return {
        "total_spectra": len(all_spectra),
        "total_files": len(mgf_files),
        "output_tsv": settings.output_tsv,
        "output_mgf": settings.output_mgf,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## Concatenate Spectra Settings

    - **MGF directory**: `{settings.mgf_dir}`
    - **Output metadata TSV**: `{settings.output_tsv}`
    - **Output concatenated MGF**: `{settings.output_mgf}`
    """)
    return


@app.cell
def run_concatenation():
    result = concatenate_spectra(settings)
    if "error" in result:
        mo.md(f"**Error:** {result['error']}")
    else:
        mo.md(f"""
        ### Results
        - **Total spectra**: {result['total_spectra']:,}
        - **MGF files processed**: {result['total_files']}
        - **Metadata saved to**: `{result['output_tsv']}`
        - **Spectra saved to**: `{result['output_mgf']}`
        """)
    return result


if __name__ == "__main__":
    app.run()
