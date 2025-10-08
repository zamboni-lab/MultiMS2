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
    """Output one row per unique FEATURE_ID, using only the first spectrum's metadata and the msbuddy result for that feature."""
    mgf_file = os.path.basename(mgf_path)
    spectra = list(load_from_mgf(mgf_path))
    # Build a mapping from FEATURE_ID to the first spectrum with that FEATURE_ID
    featureid_to_spectrum = {}
    for idx, spectrum in enumerate(spectra, start=1):
        feature_id = spectrum.metadata.get("feature_id") or str(idx)
        if feature_id not in featureid_to_spectrum:
            featureid_to_spectrum[feature_id] = spectrum
    # Map FEATURE_ID (first part of dir name) to msbuddy result dir
    try:
        featureid_to_dir = {
            d.split("_")[0]: os.path.join(msbuddy_root, d)
            for d in os.listdir(msbuddy_root)
            if os.path.isdir(os.path.join(msbuddy_root, d))
        }
    except FileNotFoundError:
        featureid_to_dir = {}
    recs = []
    for feature_id, spectrum in featureid_to_spectrum.items():
        spec_dir = featureid_to_dir.get(str(feature_id))
        smiles = spectrum.metadata.get("smiles")
        adduct = spectrum.metadata.get("adduct")
        pepmass = spectrum.get("precursor_mz") or spectrum.get("pepmass")
        rt = spectrum.get("retention_time")
        ce = (
            spectrum.metadata.get("collision_energy")
            or spectrum.metadata.get("collisionenergy")
            or spectrum.metadata.get("collisionenergy_ev")
            or spectrum.metadata.get("ce")
        )
        frag = (
            spectrum.metadata.get("fragmentation_method")
            or spectrum.metadata.get("fragmentationmethod")
            or spectrum.metadata.get("frag_method")
            or spectrum.metadata.get("fragmode")
        )
        ion = (
            spectrum.metadata.get("ionmode")
            or spectrum.metadata.get("ion_mode")
            or spectrum.metadata.get("polarity")
        )
        group_label = "_".join(
            [
                str(frag).strip().lower() if frag else "na",
                str(ce).strip().lower() if ce else "na",
                str(ion).strip().lower() if ion else "na",
            ]
        )
        has_results_tsv = False
        formula_match = False
        formula_found = spec_dir is not None
        expected_formula = smiles_to_formula(smiles) if smiles else None
        metric_map = {
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
        if spec_dir:
            tsv_path = os.path.join(spec_dir, "formula_results.tsv")
            if os.path.exists(tsv_path):
                has_results_tsv = True
                df = read_formula_results(tsv_path)
                row = best_hit_row(df, expected_formula)
                if row is not None:
                    formula_match = True
                    for rec_key in metric_map:
                        val = row.get(rec_key)
                        if rec_key in ("explained_ms2_peak", "total_valid_ms2_peak"):
                            metric_map[rec_key] = (
                                int(val) if val not in (None, "NA") else None
                            )
                        elif rec_key == "ms1_isotope_similarity":
                            metric_map[rec_key] = None if val in (None, "NA") else val
                        else:
                            metric_map[rec_key] = val
                    if (
                        metric_map["explained_ms2_peak"]
                        and metric_map["total_valid_ms2_peak"]
                    ):
                        metric_map["explained_peak_fraction"] = (
                            metric_map["explained_ms2_peak"]
                            / metric_map["total_valid_ms2_peak"]
                        )
                    ms2_explanation_idx = metric_map["ms2_explanation_idx"]
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
                                    metric_map["explained_intensity"] = (
                                        explained_sum / total_intensity
                                    )
                            except Exception:
                                pass
        rec = {
            "mgf_file": mgf_file,
            "mgf": group_label,
            "feature_id": feature_id,
            "pepmass": pepmass,
            "rt_seconds": rt,
            "smiles": smiles,
            "adduct": adduct,
            "collision_energy": ce,
            "fragmentation_method": frag,
            "ionmode": ion,
            "formula_expected": expected_formula,
            "formula_found": formula_found,
            "has_results_tsv": has_results_tsv,
            "formula_match": formula_match,
            "rank": metric_map["rank"],
            "estimated_prob": metric_map["estimated_prob"],
            "normalized_estimated_prob": metric_map["normalized_estimated_prob"],
            "estimated_fdr": metric_map["estimated_fdr"],
            "mz_error_ppm": metric_map["mz_error_ppm"],
            "ms1_isotope_similarity": metric_map["ms1_isotope_similarity"],
            "explained_ms2_peak": metric_map["explained_ms2_peak"],
            "total_valid_ms2_peak": metric_map["total_valid_ms2_peak"],
            "explained_peak_fraction": metric_map["explained_peak_fraction"],
            "explained_intensity": metric_map["explained_intensity"],
            "ms2_explanation_idx": metric_map["ms2_explanation_idx"],
            "ms2_explanation": metric_map["ms2_explanation"],
            "result_dir": os.path.basename(spec_dir) if spec_dir else None,
        }
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

    statuses = ["Formula found, correct", "Formula found, incorrect", "Formula not found"]

    df_summary = (
      df.lazy()
      .select(["mgf", "formula_found", "formula_match"])
      .with_columns(
          pl.when(~pl.col("formula_found"))
          .then(pl.lit("Formula not found"))
          .when(pl.col("formula_found") & ~pl.col("formula_match"))
          .then(pl.lit("Formula found, incorrect"))
          .when(pl.col("formula_found") & pl.col("formula_match"))
          .then(pl.lit("Formula found, correct"))
          .alias("status"),
          )
        .group_by(["mgf", "status"])
        .agg([pl.len().alias("count")])
        .collect()
    )

    # Generate all mgf × status combinations
    mgf_labels = df_summary.select(pl.col("mgf")).unique()
    all_statuses = pl.DataFrame({"status": statuses})
    all_combos = mgf_labels.join(all_statuses, how="cross")

    # Left join summary to fill missing combinations with zero
    df_summary = (
        all_combos.join(df_summary, on=["mgf", "status"], how="left")
        .with_columns(
            pl.col("count").fill_null(0),
            pl.col("count").sum().over("mgf").alias("total_count"),
        )
        .with_columns(
            (pl.col("count") / pl.col("total_count")).alias("percent")
        )
        .sort(["mgf", "status"])
        .with_columns(
            pl.col("count").cum_sum().over("mgf").alias("cumsum_count"),
        )
        .with_columns(
            (pl.col("cumsum_count") - pl.col("count") / 2).alias("x_center")
        )
        .to_pandas()  # Altair needs pandas
    )

    alt.data_transformers.enable("vegafusion")
    chart = alt.Chart(df_summary).mark_bar().encode(
        y=alt.Y(
            "mgf:N",
            sort=alt.EncodingSortField(field="total_count", order="descending"),
            title="Group",
        ),
        x=alt.X("count:Q", title = "Count", stack="zero"),
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
            alt.Tooltip("percent:Q", title="%", format=".0%"),
        ],
    )
    text = alt.Chart(df_summary).mark_text(align="center", color = "white", size=8).encode(
        y=alt.Y("mgf:N", sort=alt.EncodingSortField(field="total_count", order="descending")),
        x=alt.X("x_center:Q"),
        detail="status:N",
        text=alt.Text("percent:Q", format=".0%")
    ).transform_filter("datum.percent >= 0.05")
    return (
        (chart + text)
        .properties(title="MSBuddy: Formula identification")
    )


@app.function
def make_match_plot(df: pl.DataFrame, metric: str = "estimated_prob", n_bins: int = 10):
    import pandas as pd
    if df.is_empty() or metric not in df.columns:
        return (
            alt.Chart(
                pl.DataFrame({"note": ["No matches or invalid metric"]}).to_pandas()
            )
            .mark_text()
            .encode(text="note:N")
        )
    # Filter to matches only
    df_matches = (
        df.lazy()
        .filter((pl.col("formula_found") == True) & (pl.col("formula_match") == True))
        .filter(pl.col(metric).is_not_null())
        .select(["mgf", metric])
        .collect()
    )
    if df_matches.is_empty():
        return (
            alt.Chart(pl.DataFrame({"note": ["No valid matches"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )
    # Use pandas for robust binning
    df_matches_pd = df_matches.to_pandas()
    values = df_matches_pd[metric].values
    # Regular bins between min and max, rounded to 1 decimal if in [0,1]
    min_val, max_val = float(np.nanmin(values)), float(np.nanmax(values))
    if min_val >= 0 and max_val <= 1:
        bins = np.round(np.linspace(0, 1, n_bins + 1), 2)
    else:
        bins = np.linspace(min_val, max_val, n_bins + 1)
    bin_labels = [f"{bins[i]:.2f}–{bins[i+1]:.2f}" for i in range(len(bins)-1)]
    df_matches_pd["metric_bin"] = pd.cut(values, bins=bins, labels=bin_labels, include_lowest=True, right=False)
    # Count per group/bin
    df_binned = (
        df_matches_pd.groupby(["mgf", "metric_bin"]).size().reset_index(name="count")
    )
    # Ensure all group/bin combos appear
    all_mgfs = df_matches_pd["mgf"].unique()
    all_bins = bin_labels
    all_combos = pd.MultiIndex.from_product([all_mgfs, all_bins], names=["mgf", "metric_bin"]).to_frame(index=False)
    df_binned = all_combos.merge(df_binned, on=["mgf", "metric_bin"], how="left").fillna({"count": 0})
    # Add total_count and percent columns
    df_binned["total_count"] = df_binned.groupby("mgf")["count"].transform("sum")
    df_binned["percent"] = df_binned["count"] / df_binned["total_count"]
    # Compute x_center for label centering
    df_binned = df_binned.sort_values(["mgf", "metric_bin"])
    df_binned["cumsum_count"] = df_binned.groupby("mgf")["count"].cumsum()
    df_binned["x_center"] = df_binned["cumsum_count"] - df_binned["count"] / 2
    # Use cmcrameri colors for bins
    from cmcrameri import cm
    bin_palette = cmap_to_hex_list(cm.batlow, n_bins)
    # Horizontal stacked bar chart with text labels
    alt.data_transformers.enable("vegafusion")
    chart = alt.Chart(df_binned).mark_bar().encode(
        y=alt.Y("mgf:N", title="Group", sort=alt.EncodingSortField(field="count", op="sum", order="descending")),
        x=alt.X("count:Q", title="Count", stack="zero"),
        color=alt.Color("metric_bin:N", title=metric.replace("_", " "), scale=alt.Scale(range=bin_palette)),
        tooltip=["mgf:N", "metric_bin:N", alt.Tooltip("count:Q", title="n"), alt.Tooltip("total_count:Q", title="total"), alt.Tooltip("percent:Q", title="%", format=".0%")],
    )
    text = alt.Chart(df_binned).mark_text(align="center", color = "white", size=8).encode(
        y=alt.Y("mgf:N", sort=alt.EncodingSortField(field="count", op="sum", order="descending")),
        x=alt.X("x_center:Q"),
        detail="metric_bin:N",
        text=alt.Text("percent:Q", format=".0%")
    ).transform_filter("datum.percent >= 0.05")
    return (
        (chart + text)
        .properties(title=f"MSBuddy: {metric.replace('_', ' ')}")
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
            pl.col("smiles").is_not_null()
            # & (~pl.col("adduct").str.contains("2M"))
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
