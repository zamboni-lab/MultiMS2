# /// script
# requires-python = ">=3.13,<4"
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

__generated_with = "0.16.3"
app = marimo.App(width="full")

with app.setup:
    import os
    import re
    from pathlib import Path
    import altair as alt
    import cmcrameri.cm as cmc
    import matplotlib as mpl
    import numpy as np
    import polars as pl
    from matchms.importing import load_from_mgf
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
def process_mgf(mgf_path: str, msbuddy_root: str) -> pl.DataFrame:
    """Process single MGF mapping spectra to flat msbuddy result dirs using feature_id only."""

    def _meta_lookup_local(spectrum, *keys):
        for k in keys:
            v = spectrum.metadata.get(k)
            if v not in (None, ""):
                return v
        return None

    def _sanitize_group_value_local(v):
        if v is None:
            return "na"
        v = str(v).strip().lower()
        v = re.sub(r"[^a-z0-9.+-]", "_", v)
        v = re.sub(r"__+", "_", v).strip("_")
        return v or "na"

    mgf_file = os.path.basename(mgf_path)
    spectra = list(load_from_mgf(mgf_path))
    try:
        featureid_to_dir = {
            d.split("_")[0]: os.path.join(msbuddy_root, d)
            for d in os.listdir(msbuddy_root)
            if os.path.isdir(os.path.join(msbuddy_root, d))
        }
    except FileNotFoundError:
        featureid_to_dir = {}

    recs = []
    for idx, spectrum in enumerate(spectra, start=1):
        # Use feature_id from metadata if available, else fallback to index
        feature_id = (
            spectrum.metadata.get("feature_id")
            or str(idx)
        )
        ce = _meta_lookup_local(
            spectrum,
            "collision_energy",
            "collisionenergy",
            "collisionenergy_ev",
            "ce",
        )
        frag = _meta_lookup_local(
            spectrum,
            "fragmentation_method",
            "fragmentationmethod",
            "frag_method",
            "fragmode",
        )
        ion = _meta_lookup_local(spectrum, "ionmode", "ion_mode", "polarity")
        group_label = "_".join(
            [
                _sanitize_group_value_local(frag),
                _sanitize_group_value_local(ce),
                _sanitize_group_value_local(ion),
            ]
        )

        spec_dir = featureid_to_dir.get(str(feature_id))

        has_results_tsv = False
        formula_match = False
        formula_found = spec_dir is not None
        expected_formula = smiles_to_formula(spectrum.metadata.get("smiles"))

        metric_map = {}
        if spec_dir:
            tsv_path = os.path.join(spec_dir, "formula_results.tsv")
            if os.path.exists(tsv_path):
                has_results_tsv = True
                df = read_formula_results(tsv_path)
                row = best_hit_row(df, expected_formula)
                if row is not None:
                    formula_match = True
                    metric_map = {
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

        rec = {
            "mgf_file": mgf_file,
            "mgf": group_label,
            "spectrum_index": idx,
            "feature_id": feature_id,
            "pepmass": spectrum.get("precursor_mz") or spectrum.get("pepmass"),
            "rt_seconds": spectrum.get("retention_time"),
            "smiles": spectrum.metadata.get("smiles"),
            "adduct": spectrum.metadata.get("adduct"),
            "collision_energy": ce,
            "fragmentation_method": frag,
            "ionmode": ion,
            "formula_expected": expected_formula,
            "formula_found": formula_found,
            "has_results_tsv": has_results_tsv,
            "formula_match": formula_match,
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
            "result_dir": os.path.basename(spec_dir) if spec_dir else None,
        }

        if formula_match and spec_dir and has_results_tsv:
            tsv_path = os.path.join(spec_dir, "formula_results.tsv")
            df = read_formula_results(tsv_path)
            row = best_hit_row(df, expected_formula)
            if row:
                for rec_key, row_key in metric_map.items():
                    val = row.get(row_key)
                    if rec_key in ("explained_ms2_peak", "total_valid_ms2_peak"):
                        rec[rec_key] = int(val) if val not in (None, "NA") else None
                    elif rec_key == "ms1_isotope_similarity":
                        rec[rec_key] = None if val in (None, "NA") else val
                    else:
                        rec[rec_key] = val
                if rec["explained_ms2_peak"] and rec["total_valid_ms2_peak"]:
                    rec["explained_peak_fraction"] = (
                        rec["explained_ms2_peak"] / rec["total_valid_ms2_peak"]
                    )
                ms2_explanation_idx = rec["ms2_explanation_idx"]
                if ms2_explanation_idx not in (None, "NA"):
                    mse_path = os.path.join(spec_dir, "ms2_preprocessed.tsv")
                    if os.path.exists(mse_path):
                        try:
                            mse_df = pl.read_csv(mse_path, separator="\t")
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
                                pl.col("raw_idx").is_in(indices)
                            )["intensity"].sum()
                            if total_intensity > 0:
                                rec["explained_intensity"] = (
                                    explained_sum / total_intensity
                                )
                        except Exception:
                            pass
        recs.append(rec)

    df = pl.DataFrame(recs)
    return df


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
    sub = df.lazy().filter(pl.col("formula") == formula)
    if "rank" in df.columns:
        sub = sub.sort("rank")
    result = sub.limit(1).collect()
    if result.is_empty():
        return None
    return result.row(0, named=True)


@app.function
def cmap_to_hex_list(cmap, n_colors: int = 256) -> list[str]:
    return [mpl.colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n_colors)]


@app.function
def make_status_plot(df: pl.DataFrame):
    if df.is_empty():
        return (
            alt.Chart(pl.DataFrame({"note": ["No data"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )

    statuses = ["Formula not found", "Found but no match", "Match"]

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
        .agg([pl.len().alias("count")])
        .collect()
    )

    # Generate all mgf × status combinations
    mgfs = df_summary.select(pl.col("mgf")).unique()
    all_combos = mgfs.join(pl.DataFrame({"status": statuses}), how="cross")

    # Left join summary to fill missing combinations with zero
    df_summary = (
        all_combos.join(df_summary, on=["mgf", "status"], how="left")
        .with_columns(
            pl.col("count").fill_null(0),
            pl.col("count").sum().over("mgf").alias("total_count"),
        )
        .to_pandas()  # Altair needs pandas
    )

    alt.data_transformers.enable("vegafusion")
    return (
        alt.Chart(df_summary)
        .mark_bar()
        .encode(
            x=alt.X(
                "mgf:N",
                sort=alt.EncodingSortField(field="total_count", order="descending"),
            ),
            y="count:Q",
            color=alt.Color(
                "status:N",
                scale=alt.Scale(
                    domain=statuses,
                    range=["#004488", "#DDAA33", "#BB5566"],
                ),
            ),
            tooltip=[
                "mgf:N",
                "status:N",
                alt.Tooltip("count:Q", title="n"),
                alt.Tooltip("total_count:Q", title="total"),
            ],
        )
        .properties(title="MSBuddy: Formula Identification (ordered by total count)")
    )


@app.function
def make_match_plot(df: pl.DataFrame, metric: str = "estimated_prob"):
    if df.is_empty() or metric not in df.columns:
        return (
            alt.Chart(
                pl.DataFrame({"note": ["No matches or invalid metric"]}).to_pandas()
            )
            .mark_text()
            .encode(text="note:N")
        )
    df_matches = (
        df.lazy()
        .filter((pl.col("formula_found") == True) & (pl.col("formula_match") == True))
        .filter(pl.col(metric).is_not_null())
        .group_by("mgf")
        .agg([pl.col(metric).mean().alias(f"mean_{metric}"), pl.len().alias("count")])
        .collect()
    )

    # Ensure all mgfs appear
    all_mgfs = df.select(pl.col("mgf")).unique()
    df_matches = (
        all_mgfs.join(df_matches, on="mgf", how="left")
        .with_columns(
            pl.col(f"mean_{metric}").fill_null(0), pl.col("count").fill_null(0)
        )
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
            x=alt.X(
                "mgf:N",
                sort=alt.EncodingSortField(field="count", order="descending"),
            ),
            y="count:Q",
            color=alt.Color(
                f"mean_{metric}:Q",
                scale=alt.Scale(
                    domain=[min_val, max_val], range=cmap_to_hex_list(cmc.batlow, 256)
                ),
            ),
            tooltip=[
                "mgf:N",
                alt.Tooltip(f"mean_{metric}:Q", format=".3f", title=f"Mean {metric}"),
                alt.Tooltip("count:Q", title="n"),
            ],
        )
        .properties(title=f"MSBuddy: Match Quality by {metric} (ordered by count)")
    )


@app.cell
def _():
    # Configuration for single-MGF workflow
    mgf_path = Path("scratch/consolidated_spectra.mgf")
    msbuddy_root = Path("scratch/msbuddy")
    if not mgf_path.exists():
        print(f"Warning: MGF file missing: {mgf_path}")
    if not msbuddy_root.exists():
        print(f"Warning: MSBuddy result root missing: {msbuddy_root}")
    return mgf_path, msbuddy_root


@app.cell
def _process(mgf_path, msbuddy_root):
    df = process_mgf(str(mgf_path), str(msbuddy_root))
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
    df_filtered.write_csv("scratch/msbuddy_summary_results.tsv")

    charts = {
        "status": make_status_plot(df_filtered),
        "prob": make_match_plot(df_filtered, "estimated_prob"),
        "ms2peak": make_match_plot(df_filtered, "explained_ms2_peak"),
        "peakfrac": make_match_plot(df_filtered, "explained_peak_fraction"),
        "intensity": make_match_plot(df_filtered, "explained_intensity"),
    }

    os.makedirs("figures", exist_ok=True)
    for name, chart in charts.items():
        chart.save(f"figures/{name}.svg", format="svg")
        chart.save(f"figures/{name}.pdf", format="pdf")
        chart.save(f"figures/{name}.png", format="png")
    return


if __name__ == "__main__":
    app.run()
