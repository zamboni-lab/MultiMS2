# /// script
# requires-python = ">=3.13,<4"
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
    import os
    from simple_parsing import ArgumentParser
    import polars as pl
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf

    @dataclass
    class Settings:
        input_mgf: str = field(
            default="scratch/all_spectra.mgf",
            metadata={
                "help": "MGF input file for all spectra (concatenated, unfiltered)."
            },
        )
        output_mgf_final: str = field(
            default="scratch/filtered_spectra.mgf",
            metadata={
                "help": "MGF output file for spectra passing all post-metadata filters."
            },
        )
        min_modalities: int = field(
            default=3,
            metadata={
                "help": "Minimum number of unique modalities per inchi_aux to keep spectra for that molecule. Modality is defined as (adduct, collision_energy, fragmentation_method) if include_adduct_in_modality is True, else (collision_energy, fragmentation_method)."
            },
        )
        include_adduct_in_modality: bool = field(
            default=False,
            metadata={
                "help": "If True, include adduct in modality definition for filtering; if False, only use collision_energy and fragmentation_method."
            },
        )
        min_precursor_height: float = field(
            default=1000.0,
            metadata={"help": "Minimum precursor MS1 height."},
        )
        min_precursor_purity: float = field(
            default=0.9,
            metadata={"help": "Minimum precursor purity."},
        )
        min_signals: int = field(
            default=3,
            metadata={"help": "Minimum number of fragment signals."},
        )
        min_explained_intensity: float = field(
            default=0.4,
            metadata={"help": "Minimum explained intensity."},
        )
        min_explained_signals: float = field(
            default=0.05,
            metadata={"help": "Minimum explained signals."},
        )
        min_intensity_ratio: float = field(
            default=0.8,
            metadata={"help": "Minimum ratio to group max for explained intensity."},
        )
        min_signals_ratio: float = field(
            default=0.4,
            metadata={"help": "Minimum ratio to group max for explained signals."},
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
def filter_spectra(
    settings: Settings,
) -> str:
    """
    Filter spectra based on metadata and quality criteria.
    Reads a concatenated MGF file (produced by concat_spectra.py).
    Outputs only the final filtered spectra as MGF.
    """
    import polars as pl
    from matchms.importing import load_from_mgf
    from matchms.exporting import save_as_mgf
    import os

    # Read concatenated MGF
    mgf_input_path = settings.input_mgf
    if not os.path.exists(mgf_input_path):
        raise FileNotFoundError(f"Concatenated MGF not found: {mgf_input_path}")
    all_spectra = list(load_from_mgf(mgf_input_path))

    # Build metadata DataFrame from spectra
    rows = []
    spectrum_key_cols = [
        "adduct",
        "collision_energy",
        "description",
        "fragmentation_method",
        "inchi_aux",
        "spectrum_id",
    ]
    numeric_fields = [
        ("feature_ms1_height", float),
        ("precursor_purity", float),
        ("num_peaks", int),
        ("quality_explained_intensity", float),
        ("quality_explained_signals", float),
    ]
    for spectrum in all_spectra:
        meta = spectrum.metadata
        row = dict(meta)
        for orig, conv in numeric_fields:
            val = meta.get(orig)
            if val is not None and val != "":
                try:
                    row[orig] = conv(val)
                except Exception:
                    pass
        rows.append(row)
    if rows:
        all_keys = set()
        for r in rows:
            all_keys.update(r.keys())
        columns = sorted(all_keys)
        flat_rows = [
            {
                col: str(r.get(col, "")) if r.get(col, "") is not None else ""
                for col in columns
            }
            for r in rows
        ]
        df = pl.DataFrame(flat_rows)
    else:
        df = pl.DataFrame([])

    def norm(x):
        if isinstance(x, str):
            return x.strip().replace("[", "").replace("]", "")
        return str(x).replace("[", "").replace("]", "").strip()

    spectrum_keys = [
        tuple(norm(row.get(col, "")) for col in spectrum_key_cols) for row in rows
    ]

    # ----------- Strict Filtering Steps -----------
    # 1. Remove spectra that do not have the minimal precursor height
    df_step = df.filter(
        pl.col("feature_ms1_height").cast(pl.Float64, strict=False)
        >= settings.min_precursor_height
    )
    print(
        f"After min_precursor_height: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 2. Remove spectra that do not have the minimal precursor purity
    df_step = df_step.filter(
        pl.col("precursor_purity").cast(pl.Float64, strict=False)
        >= settings.min_precursor_purity
    )
    print(
        f"After min_precursor_purity: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 3. Remove spectra that do not have the minimal signals
    df_step = df_step.filter(
        pl.col("num_peaks").cast(pl.Int64, strict=False) >= settings.min_signals
    )
    print(
        f"After min_signals: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 4. Remove spectra that do not have the minimal explained intensity
    df_step = df_step.filter(
        pl.col("quality_explained_intensity").cast(pl.Float64, strict=False)
        >= settings.min_explained_intensity
    )
    print(
        f"After min_explained_intensity: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 5. Remove spectra that do not have the minimal explained signals
    df_step = df_step.filter(
        pl.col("quality_explained_signals").cast(pl.Float64, strict=False)
        >= settings.min_explained_signals
    )
    print(
        f"After min_explained_signals: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 6. For each (inchi_aux, adduct, modality) group, calculate max explained intensity and max explained signals
    # Define modality column according to settings
    if settings.include_adduct_in_modality:
        modality_expr = (
            pl.col("adduct").cast(pl.Utf8)
            + "|"
            + pl.col("collision_energy").cast(pl.Utf8)
            + "|"
            + pl.col("fragmentation_method").cast(pl.Utf8)
        )
    else:
        modality_expr = (
            pl.col("collision_energy").cast(pl.Utf8)
            + "|"
            + pl.col("fragmentation_method").cast(pl.Utf8)
        )
    df_step = df_step.with_columns([modality_expr.alias("modality")])

    group_cols = ["inchi_aux", "adduct", "modality"]
    df_grouped = df_step.group_by(group_cols).agg(
        [
            pl.col("quality_explained_intensity")
            .cast(pl.Float64, strict=False)
            .max()
            .alias("max_intensity"),
            pl.col("quality_explained_signals")
            .cast(pl.Float64, strict=False)
            .max()
            .alias("max_signals"),
        ]
    )
    df_step = df_step.join(df_grouped, on=group_cols, how="left")

    # 7. Remove spectra below min_intensity_ratio * max_intensity for their (inchi_aux, adduct) group
    df_step = df_step.filter(
        pl.col("quality_explained_intensity").cast(pl.Float64, strict=False)
        >= settings.min_intensity_ratio * pl.col("max_intensity")
    )
    print(
        f"After min_intensity_ratio: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 8. Remove spectra below min_signals_ratio * max_signals for their (inchi_aux, adduct) group
    df_step = df_step.filter(
        pl.col("quality_explained_signals").cast(pl.Float64, strict=False)
        >= settings.min_signals_ratio * pl.col("max_signals")
    )
    print(
        f"After min_signals_ratio: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 9. Remove spectra that do not have the minimal explained intensity
    df_step = df_step.filter(
        pl.col("quality_explained_intensity").cast(pl.Float64, strict=False)
        >= settings.min_explained_intensity
    )
    print(
        f"After min_explained_intensity: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 10. Remove spectra that do not have the minimal explained signals
    df_step = df_step.filter(
        pl.col("quality_explained_signals").cast(pl.Float64, strict=False)
        >= settings.min_explained_signals
    )
    print(
        f"After min_explained_signals: {df_step.shape[0]} spectra, {df_step.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct)"
    )

    # 11. For each (inchi_aux, adduct) pair, count unique modalities (collision_energy + fragmentation_method)
    modality_expr = (
        pl.col("collision_energy").cast(pl.Utf8)
        + "|"
        + pl.col("fragmentation_method").cast(pl.Utf8)
    )
    df_modalities = df_step.with_columns(
        [
            modality_expr.alias("modality"),
            pl.col("inchi_aux").cast(pl.Utf8),
            pl.col("adduct").cast(pl.Utf8),
        ]
    )
    df_modalities_unique = df_modalities.unique(
        subset=["inchi_aux", "adduct", "modality"]
    )
    modalities_count = df_modalities_unique.group_by(["inchi_aux", "adduct"]).agg(
        pl.col("modality").count().alias("n_modalities")
    )
    keep_pairs = (
        modalities_count.filter(pl.col("n_modalities") >= settings.min_modalities)
        .select(["inchi_aux", "adduct"])
        .to_dict(as_series=False)
    )
    keep_pairs_set = set(zip(keep_pairs["inchi_aux"], keep_pairs["adduct"]))
    df_final = df_step.filter(
        pl.struct(["inchi_aux", "adduct"]).map_elements(
            lambda row: (row["inchi_aux"], row["adduct"]) in keep_pairs_set,
            return_dtype=pl.Boolean,
        )
    )
    print(
        f"After min_modalities per (inchi_aux, adduct): {df_final.select(['inchi_aux']).unique().shape[0]} inchi_aux, {df_final.select(['inchi_aux', 'adduct']).unique().shape[0]} unique (inchi_aux, adduct), {df_modalities_unique.shape[0]} unique (inchi_aux, adduct, modality)"
    )

    # Normalize key fields in df_final and drop duplicates
    def normalize_key(row):
        # Remove brackets, whitespace, and cast to string
        def norm(x):
            if isinstance(x, str):
                return x.strip().replace("[", "").replace("]", "")
            return str(x).replace("[", "").replace("]", "").strip()

        return tuple(norm(row[col]) for col in spectrum_key_cols)

    # Add a normalized key column for duplicate checking
    df_final = df_final.with_columns(
        [pl.struct(spectrum_key_cols).map_elements(normalize_key).alias("_norm_key")]
    )
    n_before = df_final.shape[0]
    df_final = df_final.unique(subset=["_norm_key"])
    n_after = df_final.shape[0]
    if n_after < n_before:
        print(
            f"Dropped {n_before - n_after} duplicate spectra based on normalized key."
        )

    # Build keep_keys from normalized keys (convert lists to tuples)
    keep_keys = set(
        tuple(x) if isinstance(x, list) else x for x in df_final["_norm_key"].to_list()
    )

    # Select spectra: only keep spectra whose normalized key is in keep_keys
    final_spectra = []
    seen_keys = set()
    for spectrum, key in zip(all_spectra, spectrum_keys):
        # Normalize the key for comparison
        def norm(x):
            if isinstance(x, str):
                return x.strip().replace("[", "").replace("]", "")
            return str(x).replace("[", "").replace("]", "").strip()

        key_str = tuple(norm(x) for x in key)
        if key_str in keep_keys:
            if key_str not in seen_keys:
                final_spectra.append(spectrum)
                seen_keys.add(key_str)
            else:
                print(
                    f"Warning: duplicate spectrum key {key_str} detected and skipped."
                )
    print(f"Final spectra selected for output: {len(final_spectra)}")

    # No need for a final robust chimeric filter, as all chimeric filtering is done in the DataFrame

    # Sort spectra: positive polarity first, then by description, then by increasing spectrum length
    def get_sort_tuple(spectrum):
        meta = spectrum.metadata
        ionmode = str(meta.get("ionmode", meta.get("IONMODE", "")).lower())
        polarity_order = 0 if ionmode == "positive" else 1
        description = str(meta.get("description", ""))
        num_peaks = (
            len(spectrum.peaks.mz)
            if hasattr(spectrum, "peaks") and hasattr(spectrum.peaks, "mz")
            else 0
        )
        return (polarity_order, description, num_peaks)

    final_spectra_sorted = sorted(final_spectra, key=get_sort_tuple)

    save_as_mgf(final_spectra_sorted, settings.output_mgf_final)
    print(
        f"Exported {len(final_spectra_sorted)} final spectra to {settings.output_mgf_final}"
    )
    return (
        f"Read {df.shape[0]} (one per spectrum) from {settings.input_mgf}, "
        f"final spectra to {settings.output_mgf_final}"
    )


@app.cell
def run_app():
    filter_spectra(settings)
    return


if __name__ == "__main__":
    from simple_parsing import ArgumentParser

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")
    args = parser.parse_args()
    filter_spectra(args.settings)
