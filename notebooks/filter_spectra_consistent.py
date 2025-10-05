# /// script
# requires-python = "<3.13,>=3.12"
# dependencies = [
#     "marimo",
#     "simple_parsing",
#     "matchms",
#     "polars",
# ]
# ///

import marimo

__generated_with = "0.16.3"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    import marimo as mo
    import glob
    import os
    from typing import Optional, Tuple
    from simple_parsing import ArgumentParser
    import polars as pl
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf

    @dataclass
    class Settings:
        mgf_dir: str = field(
            default="scratch/",
            metadata={"help": "Directory containing .mgf files."},
        )
        output_tsv: str = field(
            default="scratch/spectra_metadata.tsv",
            metadata={"help": "TSV output file for spectrum metadata (not aggregated)."},
        )
        output_mgf_all: str = field(
            default="scratch/all_spectra.mgf",
            metadata={"help": "MGF output file for all spectra (concatenated, unfiltered)."},
        )
        output_mgf_final: str = field(
            default="scratch/filtered_spectra.mgf",
            metadata={"help": "MGF output file for spectra passing all post-metadata filters."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        if mo.running_in_notebook():
            return Settings()
        else:
            args = parser.parse_args()
            return args.settings

    settings = parse_args()

@app.function
def export_full_and_final_filtered(
    settings: Settings,
) -> str:
    """
    Write all spectra metadata (one row per spectrum, not aggregated) to TSV,
    output all spectra as concatenated MGF, and output only the final filtered spectra as MGF.
    Filtering is based on the unique combination of
    (collision_energy, description, fragmentation_method, inchi_aux, spectrum_id).
    If output files already exist, they are removed before writing.
    """
    # Remove output files if they exist
    for out_path in [settings.output_tsv, settings.output_mgf_all, settings.output_mgf_final]:
        if os.path.exists(out_path):
            print(f"Removing existing output: {out_path}")
            os.remove(out_path)

    mgf_files = glob.glob(os.path.join(settings.mgf_dir, "*.mgf"))
    all_spectra = []
    rows = []
    spectrum_keys = []  # Will hold (collision_energy, description, fragmentation_method, inchi_aux, spectrum_id)

    numeric_fields = [
        ("retention_time", float),
        ("feature_ms1_height", float),
        ("precursor_purity", float),
        ("num_peaks", int),
        ("quality_explained_intensity", float),
        ("quality_explained_signals", float),
    ]

    for mgf_file in mgf_files:
        print(f"Loading {mgf_file} ...")
        filename = os.path.basename(mgf_file)
        for spectrum in load_from_mgf(mgf_file):
            meta = spectrum.metadata
            row = dict(meta)
            row["mgf_filename"] = filename
            # Numeric conversion for filtering, but keep original names
            for orig, conv in numeric_fields:
                val = meta.get(orig)
                if val is not None and val != "":
                    try:
                        row[orig] = conv(val)
                    except Exception:
                        pass  # leave as is if conversion fails
            # Save tuple for filtering (include inchi_aux)
            spectrum_key = (
                str(meta.get("collision_energy", "")),
                str(meta.get("description", "")),
                str(meta.get("fragmentation_method", "")),
                str(meta.get("inchi_aux", "")),
                str(meta.get("spectrum_id", "")),
            )
            spectrum_keys.append(spectrum_key)
            rows.append(row)
            all_spectra.append(spectrum)

    # CSV export: Keep only original keys found in the input, plus mgf_filename
    if rows:
        all_keys = set()
        for r in rows:
            all_keys.update(r.keys())
        columns = ["mgf_filename"] + sorted([k for k in all_keys if k != "mgf_filename"])
        flat_rows = [
            {col: str(r.get(col, "")) if r.get(col, "") is not None else "" for col in columns}
            for r in rows
        ]
        df = pl.DataFrame(flat_rows)
        df.write_csv(settings.output_tsv, separator="\t")
        print(f"Wrote {len(rows)} spectra metadata rows to {settings.output_tsv}")
    else:
        df = pl.DataFrame([])

    # Write out the full concatenated spectra to output_mgf_all
    save_as_mgf(all_spectra, settings.output_mgf_all)
    print(f"Exported {len(all_spectra)} spectra to {settings.output_mgf_all}")

    # ----------- Advanced filtering -----------
    df_meta = df

    # Thresholds
    min_precursor_height = 1000.0
    min_precursor_purity = 0.9
    min_signals = 3
    min_explained_intensity = 0.4
    min_explained_signals = 0.05
    min_modalities = 3
    min_intensity_ratio = 0.8
    min_signals_ratio = 0.4

    filter_cols = [
        "adduct",
        "feature_ms1_height",
        "precursor_purity",
        "collision_energy",
        "description",
        "fragmentation_method",
        "inchi_aux",
        "num_peaks",
        "quality_explained_intensity",
        "quality_explained_signals",
        "spectrum_id",
    ]
    df_dis = df_meta.unique(subset=filter_cols)

    print("Filtering now...")
    df_clean = df_dis.filter(
        (pl.col("feature_ms1_height").cast(pl.Float64, strict=False) >= min_precursor_height) &
        (pl.col("precursor_purity").cast(pl.Float64, strict=False) >= min_precursor_purity) &
        (pl.col("num_peaks").cast(pl.Int64, strict=False) >= min_signals) &
        (pl.col("quality_explained_intensity").cast(pl.Float64, strict=False) >= min_explained_intensity) &
        (pl.col("quality_explained_signals").cast(pl.Float64, strict=False) >= min_explained_signals)
    ).unique()
    print(df_clean)

    # Only add 'connectivity' to the DataFrame for filtering; do not export to CSV
    df_clean = df_clean.with_columns(
        pl.col("inchi_aux").str.replace(r"-.*", "", literal=False).alias("connectivity")
    )
    print(df_clean)
    group_cols = ["adduct", "connectivity"]

    # --- MODALITIES COUNT LOGIC ---
    modalities = (
        df_clean
        .unique(subset=["adduct", "collision_energy", "fragmentation_method", "inchi_aux", "connectivity"])
    )
    modalities_n = (
        modalities
        .group_by(group_cols)
        .agg(n=pl.len())
    )
    modalities = modalities.join(modalities_n, on=group_cols)
    modalities = modalities.filter(pl.col("n") >= min_modalities)
    modalities = modalities.unique(subset=["adduct", "connectivity", "inchi_aux"])

    df_fin = df_clean.join(modalities, on=["adduct", "connectivity", "inchi_aux"], how="inner", suffix="_joined")
    # Remove any duplicated columns from join, keep only original names (except 'n' which is needed)
    cols_to_drop = [col for col in df_fin.columns if col.endswith("_joined") and col != "n"]
    if cols_to_drop:
        df_fin = df_fin.drop(cols_to_drop)
    print(df_fin)

    group_cols2 = ["adduct", "collision_energy", "fragmentation_method", "inchi_aux"]
    df_grouped = (
        df_fin
        .group_by(group_cols2)
        .agg([
            pl.col("quality_explained_intensity").cast(pl.Float64, strict=False).max().alias("qual_int_max"),
            pl.col("quality_explained_signals").cast(pl.Float64, strict=False).max().alias("qual_sig_max"),
            pl.len().alias("n_spe")
        ])
    )
    print(df_grouped)
    df_fin2 = df_fin.join(df_grouped, on=group_cols2, how="left", suffix="_joined2")
    # Drop any duplicate columns from this join as well
    cols_to_drop2 = [col for col in df_fin2.columns if col.endswith("_joined2")]
    if cols_to_drop2:
        df_fin2 = df_fin2.drop(cols_to_drop2)
    print(df_fin2)
    df_final = (
        df_fin2
        .filter(
            (pl.col("quality_explained_intensity").cast(pl.Float64, strict=False) >= min_intensity_ratio * pl.col("qual_int_max")) &
            (pl.col("quality_explained_signals").cast(pl.Float64, strict=False) >= min_signals_ratio * pl.col("qual_sig_max"))
        )
        .unique()
    )
    print(df_final)

    spectrum_key_cols = ["collision_energy", "description", "fragmentation_method", "inchi_aux", "spectrum_id"]
    keep_keys = set(
        tuple(df_final.get_column(col)[i] for col in spectrum_key_cols)
        for i in range(df_final.height)
    )

    # Select spectra: only keep spectra whose unique key is in keep_keys
    final_spectra = []
    for spectrum, key in zip(all_spectra, spectrum_keys):
        key_str = tuple(str(x) for x in key)
        if key_str in keep_keys:
            final_spectra.append(spectrum)

    save_as_mgf(final_spectra, settings.output_mgf_final)
    print(f"Exported {len(final_spectra)} final spectra to {settings.output_mgf_final}")
    return (
        f"Wrote {len(rows)} (one per spectrum) to {settings.output_tsv}, "
        f"all spectra to {settings.output_mgf_all}, "
        f"final spectra to {settings.output_mgf_final}"
    )

@app.cell
def run_app():
    export_full_and_final_filtered(settings)
    return

if __name__ == "__main__":
    app.run()