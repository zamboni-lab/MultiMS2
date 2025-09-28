# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "marimo",
#     "polars",
#     "matchms",
#     "rdkit",
#     "altair",
#     "pandas",
#     "pyarrow",
#     "vegafusion",
#     "vl-convert-python",
#     "cmcrameri",
# ]
# ///

import marimo

__generated_with = "0.14.17"
app = marimo.App(width="full")

with app.setup:
    import glob
    import os
    import re
    from pathlib import Path
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import altair as alt
    import cmcrameri.cm as cmc
    import matplotlib as mpl
    import numpy as np
    import polars as pl
    from matchms.importing import load_from_mgf
    from multiprocessing import cpu_count
    from rdkit import Chem


@app.function
def smiles_to_formula(smiles: str) -> str | None:
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.rdMolDescriptors.CalcMolFormula(mol) if mol else None
    except Exception:
        return None


@app.function
def parse_spectrum_id(idx: int, spectrum) -> tuple[str, float | None, float | None]:
    mz = spectrum.get("precursor_mz")
    rt = spectrum.get("retention_time")

    mz = float(mz) if mz is not None else None
    rt = float(rt) if rt is not None else None

    if mz is not None and rt is not None:
        sid = f"{idx}_mz_{mz:.4f}_rt_{rt:.2f}"
    else:
        sid = None

    return sid, mz, rt


@app.function
def read_formula_results(tsv_path: str) -> pl.DataFrame | None:
    try:
        return pl.read_csv(tsv_path, separator="\t", infer_schema_length=1000)
    except Exception:
        return None


@app.function
def best_hit_row(df: pl.DataFrame, formula: str) -> dict | None:
    if df is None or formula is None or "formula" not in df.columns:
        return None

    # Use lazy evaluation and more efficient filtering
    sub = df.lazy().filter(pl.col("formula") == formula)

    if "rank" in df.columns:
        sub = sub.sort("rank")

    # Collect only the first row
    result = sub.limit(1).collect()

    if result.is_empty():
        return None

    return result.row(0, named=True)


@app.function
def batch_check_directories(base_paths: list[str]) -> dict[str, bool]:
    """Batch check directory existence to reduce I/O calls"""
    return {path: os.path.isdir(path) for path in base_paths}


@app.function
def process_mgf(mgf_path: str, msbuddy_root: str) -> pl.DataFrame:
    mgf_name = os.path.splitext(os.path.basename(mgf_path))[0]
    msbuddy_dir = os.path.join(msbuddy_root, mgf_name)

    # Load all spectra at once
    spectra = list(load_from_mgf(mgf_path))

    # Pre-build all potential paths and check existence in batch
    potential_paths = []
    spectrum_info = []

    for idx, spectrum in enumerate(spectra, start=1):
        sid, pepmass, rt = parse_spectrum_id(idx, spectrum)
        smiles = spectrum.metadata.get("smiles")
        adduct = spectrum.metadata.get("adduct")
        expected_formula = smiles_to_formula(smiles)

        # Build potential directory paths
        if sid:
            sid_suffix = re.sub(r"^\d+_", "", sid)
            glob_pattern = os.path.join(msbuddy_dir, f"*_{sid_suffix}")
            matches = glob.glob(glob_pattern)
            spec_dir = matches[0] if matches else os.path.join(msbuddy_dir, sid)
        else:
            spec_dir = None

        potential_paths.append(spec_dir)
        spectrum_info.append(
            {
                "idx": idx,
                "sid": sid,
                "pepmass": pepmass,
                "rt": rt,
                "smiles": smiles,
                "adduct": adduct,
                "expected_formula": expected_formula,
                "spec_dir": spec_dir,
            },
        )

    # Batch check directory existence
    dir_exists = batch_check_directories([p for p in potential_paths if p is not None])

    # Process records with reduced I/O
    recs = []
    for i, info in enumerate(spectrum_info):
        spec_dir = info["spec_dir"]
        formula_found = dir_exists.get(spec_dir, False) if spec_dir else False

        rec = {
            "mgf": mgf_name,
            "spectrum_index": info["idx"],
            "spectrum_id": info["sid"],
            "pepmass": info["pepmass"],
            "rt_seconds": info["rt"],
            "smiles": info["smiles"],
            "adduct": info["adduct"],
            "formula_expected": info["expected_formula"],
            "formula_found": formula_found,
            "has_results_tsv": False,
            "formula_match": False,
            "rank": None,
            "estimated_prob": None,
            "normalized_estimated_prob": None,
            "estimated_fdr": None,
            "mz_error_ppm": None,
            "ms1_isotope_similarity": None,
            "explained_ms2_peak": None,
            "total_valid_ms2_peak": None,
            "explained_peak_fraction": None,
            "explained_intensity": None,
            "ms2_explanation_idx": None,
            "ms2_explanation": None,
        }

        if formula_found and spec_dir:
            tsv_path = os.path.join(spec_dir, "formula_results.tsv")
            if os.path.exists(tsv_path):
                rec["has_results_tsv"] = True
                df = read_formula_results(tsv_path)
                row = best_hit_row(df, info["expected_formula"])
                if row is not None:
                    rec["formula_match"] = True
                    # Batch update all fields from row
                    field_mapping = {
                        "rank": "rank",
                        "estimated_prob": "estimated_prob",
                        "normalized_estimated_prob": "normalized_estimated_prob",
                        "estimated_fdr": "estimated_fdr",
                        "mz_error_ppm": "mz_error_ppm",
                        "ms1_isotope_similarity": "ms1_isotope_similarity",
                        "explained_ms2_peak": "explained_ms2_peak",
                        "total_valid_ms2_peak": "total_valid_ms2_peak",
                        "ms2_explanation_idx": "ms2_explanation_idx",
                        "ms2_explanation": "ms2_explanation",
                    }

                    for rec_key, row_key in field_mapping.items():
                        value = row.get(row_key)
                        if rec_key in [
                            "explained_ms2_peak",
                            "total_valid_ms2_peak",
                        ]:
                            rec[rec_key] = (
                                int(value) if value not in (None, "NA") else None
                            )
                        elif rec_key == "ms1_isotope_similarity":
                            rec[rec_key] = None if value in (None, "NA") else value
                        else:
                            rec[rec_key] = value

                    # Calculate explained peak fraction efficiently
                    if rec["explained_ms2_peak"] and rec["total_valid_ms2_peak"]:
                        rec["explained_peak_fraction"] = (
                            rec["explained_ms2_peak"] / rec["total_valid_ms2_peak"]
                        )

                    # Compute explained_intensity using polars
                    ms2_explanation_idx = rec["ms2_explanation_idx"]
                    if ms2_explanation_idx not in (None, "NA"):
                        mse_path = os.path.join(spec_dir, "ms2_preprocessed.tsv")
                        if os.path.exists(mse_path):
                            try:
                                mse_df = pl.read_csv(mse_path, separator="\t")
                                # indices can be comma-separated string
                                if isinstance(ms2_explanation_idx, str):
                                    indices = [
                                        int(i) for i in ms2_explanation_idx.split(",")
                                    ]
                                elif isinstance(ms2_explanation_idx, int):
                                    indices = [ms2_explanation_idx]
                                else:
                                    indices = list(ms2_explanation_idx)

                                total_intensity = mse_df["intensity"].sum()
                                explained_sum = mse_df.filter(
                                    pl.col("raw_idx").is_in(indices),
                                )["intensity"].sum()

                                if total_intensity > 0:
                                    rec["explained_intensity"] = (
                                        explained_sum / total_intensity
                                    )
                                else:
                                    rec["explained_intensity"] = None
                            except Exception:
                                rec["explained_intensity"] = None

        recs.append(rec)

    return pl.DataFrame(recs)


@app.function
def process_mgf_wrapper(args):
    """Wrapper function for multiprocessing"""
    mgf_path, msbuddy_root = args
    return process_mgf(mgf_path, msbuddy_root)


@app.function
def process_directory_parallel(root_dir: str, max_workers: int = None) -> pl.DataFrame:
    """Process directory with parallel MGF processing"""
    spectra_dir = os.path.join(root_dir, "spectra")
    msbuddy_dir = os.path.join(root_dir, "msbuddy")
    mgfs = glob.glob(os.path.join(spectra_dir, "*.mgf"))

    if not mgfs:
        return pl.DataFrame()

    # Determine optimal number of workers
    if max_workers is None:
        max_workers = min(len(mgfs), os.cpu_count() or 1)

    dfs = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_mgf = {
            executor.submit(process_mgf_wrapper, (mgf, msbuddy_dir)): mgf
            for mgf in mgfs
        }

        # Collect results as they complete
        for future in as_completed(future_to_mgf):
            try:
                df = future.result()
                if not df.is_empty():
                    dfs.append(df)
            except Exception as e:
                mgf = future_to_mgf[future]
                print(f"Error processing {mgf}: {e}")

    return pl.concat(dfs, how="vertical") if dfs else pl.DataFrame()


@app.function
def process_directory(root_dir: str) -> pl.DataFrame:
    """Fallback to original sequential processing"""
    spectra_dir = os.path.join(root_dir, "spectra")
    msbuddy_dir = os.path.join(root_dir, "msbuddy")
    mgfs = glob.glob(os.path.join(spectra_dir, "*.mgf"))
    dfs = [process_mgf(mgf, msbuddy_dir) for mgf in mgfs]
    return pl.concat(dfs, how="vertical") if dfs else pl.DataFrame()


@app.function
def make_status_plot(df: pl.DataFrame):
    if df.is_empty():
        return (
            alt.Chart(pl.DataFrame({"note": ["No data"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )

    # Use lazy evaluation and more efficient aggregation
    df_summary = (
        df.lazy()
        .select(["mgf", "formula_found", "formula_match"])
        .with_columns(
            pl.when(~pl.col("formula_found"))
            .then(pl.lit("Formula not found"))
            .when(pl.col("formula_found") & ~pl.col("formula_match"))
            .then(pl.lit("Found but no match"))
            .when(pl.col("formula_found") & pl.col("formula_match"))
            .then(pl.lit("Match"))
            .alias("status"),
        )
        .group_by(["mgf", "status"])
        .agg(
            [
                pl.len().alias("count"),
            ],
        )
        .collect()
        .to_pandas()
    )

    alt.data_transformers.enable("vegafusion")
    return (
        alt.Chart(df_summary)
        .mark_bar()
        .encode(
            x="mgf:N",
            y="count:Q",
            color=alt.Color(
                "status:N",
                scale=alt.Scale(
                    domain=[
                        "Formula not found",
                        "Found but no match",
                        "Match",
                    ],
                    range=["#004488", "#DDAA33", "#BB5566"],
                ),
            ),
            tooltip=["mgf:N", "status:N", alt.Tooltip("count:Q", title="n")],
        )
        .properties(title="MSBuddy: Formula Identification")
    )


@app.function
def cmap_to_hex_list(cmap, n_colors: int = 256) -> list[str]:
    """Convert colormap to hex list efficiently"""
    return [mpl.colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n_colors)]


@app.function
def make_match_plot(df: pl.DataFrame, metric: str = "estimated_prob"):
    if df.is_empty() or metric not in df.columns:
        return (
            alt.Chart(
                pl.DataFrame({"note": ["No matches or invalid metric"]}).to_pandas(),
            )
            .mark_text()
            .encode(text="note:N")
        )

    # Use lazy evaluation for better performance
    df_matches = (
        df.lazy()
        .filter((pl.col("formula_found") == True) & (pl.col("formula_match") == True))
        .filter(pl.col(metric).is_not_null())
        .group_by("mgf")
        .agg(
            [
                pl.col(metric).mean().alias(f"mean_{metric}"),
                pl.len().alias("count"),
            ],
        )
        .collect()
        .to_pandas()
    )

    if df_matches.empty:
        return (
            alt.Chart(pl.DataFrame({"note": ["No valid matches"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )

    min_val, max_val = (
        df_matches[f"mean_{metric}"].min(),
        df_matches[f"mean_{metric}"].max(),
    )

    return (
        alt.Chart(df_matches)
        .mark_bar()
        .encode(
            x="mgf:N",
            y="count:Q",
            color=alt.Color(
                f"mean_{metric}:Q",
                scale=alt.Scale(
                    domain=[min_val, max_val],
                    range=cmap_to_hex_list(cmc.batlow, 256),
                ),
            ),
            tooltip=[
                "mgf:N",
                alt.Tooltip(f"mean_{metric}:Q", format=".3f", title=f"Mean {metric}"),
                alt.Tooltip("count:Q", title="n"),
            ],
        )
        .properties(title=f"MSBuddy: Match Quality by {metric}")
    )


@app.cell
def _():
    root_dir = Path("/Volumes/T7/data/zeno_lib_v2")
    out_name = "summary_test.csv"
    use_parallel = True  # Toggle parallel processing
    return root_dir, use_parallel


@app.cell
def _process(root_dir, use_parallel):
    if use_parallel:
        n_cpus = cpu_count()
        df = process_directory_parallel(root_dir=root_dir, max_workers=n_cpus)
    else:
        df = process_directory(root_dir=root_dir)
    return (df,)


@app.cell
def _(df):
    df
    return


@app.cell
def _filter(df):
    # Use lazy evaluation for filtering
    df_filtered = (
        df.lazy()
        .filter(
            pl.col("smiles").is_not_null() & (~pl.col("adduct").str.contains("2M")),
            # & (~pl.col("adduct").str.contains("-2H"))
            # & (
            #     (~pl.col("adduct").str.starts_with("[M-"))
            #     | (pl.col("adduct") == "[M-H]-")
            # )
        )
        .collect()
    )

    def summary(df: pl.DataFrame):
        """Optimized summary function using polars native operations"""
        summary_data = []
        for col in df.columns:
            dtype = df[col].dtype
            null_count = df[col].null_count()

            if dtype == pl.Boolean:
                true_count = df[col].sum()
                false_count = df.shape[0] - true_count - null_count
                summary_data.append(
                    f"{col} ({dtype}): {null_count} nulls, {true_count} True, {false_count} False",
                )
            else:
                unique_count = df[col].n_unique()
                summary_data.append(
                    f"{col} ({dtype}): {null_count} nulls, {unique_count} unique",
                )

        for line in summary_data:
            print(line)

    summary(df_filtered)
    df_filtered
    return (df_filtered,)


@app.cell
def _plot(df_filtered):
    # Use lazy writing for better performance
    df_filtered.write_csv("results.tsv")

    charts = {
        "status": make_status_plot(df_filtered),
        "prob": make_match_plot(df_filtered, "estimated_prob"),
        "ms2peak": make_match_plot(df_filtered, "explained_ms2_peak"),
        "peakfrac": make_match_plot(df_filtered, "explained_peak_fraction"),
        "intensity": make_match_plot(df_filtered, "explained_intensity"),
    }

    # Parallel chart saving could be added here if needed
    for name, chart in charts.items():
        chart.save(f"{name}.svg", format="svg")
        chart.save(f"{name}.png", format="png")
    return


if __name__ == "__main__":
    app.run()
