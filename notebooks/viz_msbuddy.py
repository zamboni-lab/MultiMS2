# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "altair",
#     "cmcrameri",
#     "marimo",
#     "matchms",
#     "pandas",
#     "polars",
#     "pyarrow",
#     "rdkit",
#     "vegafusion",
#     "vl-convert-python",
# ]
# ///

import marimo

__generated_with = "0.18.1"
app = marimo.App(width="full")

with app.setup:
    import os
    import altair as alt
    import matplotlib as mpl
    import numpy as np
    import polars as pl
    from matchms.importing import load_from_mgf
    from pathlib import Path
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
def safe_div(numerator, denominator, default=None):
    try:
        return float(numerator) / float(denominator) if denominator else default
    except Exception:
        return default


@app.function
def pd_unique_preserve_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


@app.function
def cmap_to_hex_list(cmap, n_colors: int = 256) -> list[str]:
    return [mpl.colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, n_colors)]


@app.function
def get_mgf_order(df_pd):
    """Return ordered MGF groups with custom priority."""
    desired_order = [
        "cid_[20.0]_positive",
        "cid_[40.0]_positive",
        "cid_[60.0]_positive",
        "ead_[12.0]_positive",
        "ead_[16.0]_positive",
        "ead_[24.0]_positive",
        "cid_[20.0]_negative",
        "cid_[40.0]_negative",
        "cid_[60.0]_negative",
        "ead_[12.0]_negative",
        "ead_[16.0]_negative",
        "ead_[24.0]_negative",
    ]
    try:
        values = list(df_pd["mgf"].astype(str).values)
    except Exception:
        values = list(df_pd["mgf"].values)
    mgf_order = [m for m in desired_order if m in values]
    rest = [m for m in pd_unique_preserve_order(values) if m not in mgf_order]
    mgf_order.extend(rest)
    return mgf_order


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
    return None if result.is_empty() else result.row(0, named=True)


@app.function
def process_mgf(mgf_path: str, msbuddy_root: str) -> pl.DataFrame:
    """Output one row per unique FEATURE_ID."""
    mgf_file = os.path.basename(mgf_path)
    spectra = list(load_from_mgf(mgf_path))
    featureid_to_spectrum = {
        s.metadata.get("feature_id") or str(i): s for i, s in enumerate(spectra, 1)
    }

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
        metadata = spectrum.metadata
        pepmass = spectrum.get("precursor_mz") or spectrum.get("pepmass")
        rt = spectrum.get("retention_time")

        ce = (
            metadata.get("collision_energy")
            or metadata.get("collisionenergy")
            or metadata.get("collisionenergy_ev")
            or metadata.get("ce")
        )
        frag = (
            metadata.get("fragmentation_method")
            or metadata.get("fragmentationmethod")
            or metadata.get("frag_method")
            or metadata.get("fragmode")
        )
        ion = (
            metadata.get("ionmode")
            or metadata.get("ion_mode")
            or metadata.get("polarity")
        )
        smiles = metadata.get("smiles")
        adduct = metadata.get("adduct")

        group_label = "_".join(
            [
                str(frag).strip().lower() if frag else "na",
                str(ce).strip().lower() if ce else "na",
                str(ion).strip().lower() if ion else "na",
            ]
        )

        expected_formula = smiles_to_formula(smiles) if smiles else None
        has_results_tsv = False
        formula_match = False
        formula_found = spec_dir is not None
        metric_map = {
            key: None
            for key in [
                "rank",
                "estimated_prob",
                "normalized_estimated_prob",
                "estimated_fdr",
                "mz_error_ppm",
                "ms1_isotope_similarity",
                "explained_ms2_peak",
                "total_valid_ms2_peak",
                "explained_peak_fraction",
                "explained_intensity",
                "ms2_explanation_idx",
                "ms2_explanation",
            ]
        }

        if spec_dir:
            tsv_path = os.path.join(spec_dir, "formula_results.tsv")
            if os.path.exists(tsv_path):
                has_results_tsv = True
                df_tsv = read_formula_results(tsv_path)
                row = best_hit_row(df_tsv, expected_formula)
                if row is not None:
                    formula_match = True
                    for key in metric_map:
                        val = row.get(key)
                        if key in ("explained_ms2_peak", "total_valid_ms2_peak"):
                            try:
                                metric_map[key] = (
                                    int(val) if val not in (None, "NA", "") else None
                                )
                            except Exception:
                                metric_map[key] = None
                        else:
                            metric_map[key] = (
                                val if val not in (None, "NA", "") else None
                            )

                    total_peak = metric_map["total_valid_ms2_peak"]
                    ms2_peak = metric_map["explained_ms2_peak"]
                    metric_map["explained_peak_fraction"] = safe_div(
                        ms2_peak, total_peak, default=0.0
                    )
                    ms2_expl_idx = metric_map["ms2_explanation_idx"]
                    mse_path = os.path.join(spec_dir, "ms2_preprocessed.tsv")
                    if ms2_expl_idx not in (None, "NA", "") and os.path.exists(
                        mse_path
                    ):
                        try:
                            mse_df = pl.read_csv(mse_path, separator="\t")
                            indices = [
                                int(i)
                                for i in str(ms2_expl_idx).split(",")
                                if i.strip().isdigit()
                            ]
                            total_intensity = mse_df["intensity"].sum()
                            explained_sum = mse_df.filter(
                                pl.col("raw_idx").is_in(indices)
                            )["intensity"].sum()
                            metric_map["explained_intensity"] = safe_div(
                                explained_sum, total_intensity, default=0.0
                            )
                        except Exception:
                            metric_map["explained_intensity"] = None

        recs.append(
            {
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
                **metric_map,
                "result_dir": os.path.basename(spec_dir) if spec_dir else None,
            }
        )
    return pl.DataFrame(recs)


@app.function
def make_status_plot(df: pl.DataFrame, panel_label: str = "A", chart_width: int = 500):
    if df.is_empty():
        return (
            alt.Chart(pl.DataFrame({"note": ["No data"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )
    import pandas as pd

    statuses = [
        "Formula found, correct",
        "Formula found, incorrect",
        "Formula not found",
    ]
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
            .alias("status")
        )
        .group_by(["mgf", "status"])
        .agg([pl.len().alias("count")])
        .collect()
    )
    mgf_labels = df_summary.select(pl.col("mgf")).unique()
    all_combos = mgf_labels.join(pl.DataFrame({"status": statuses}), how="cross")
    df_summary = (
        all_combos.join(df_summary, on=["mgf", "status"], how="left")
        .with_columns(
            pl.col("count").fill_null(0),
            pl.col("count").sum().over("mgf").alias("total_count"),
        )
        .with_columns((pl.col("count") / pl.col("total_count")).alias("percent"))
        .sort(["mgf", "status"])
        .with_columns(pl.col("count").cum_sum().over("mgf").alias("cumsum_count"))
        .with_columns((pl.col("cumsum_count") - pl.col("count") / 2).alias("x_center"))
        .to_pandas()
    )

    # build mgf order and set categorical on the pandas DF so the order travels with data
    mgf_order = get_mgf_order(df_summary)
    df_summary["mgf"] = pd.Categorical(
        df_summary["mgf"], categories=mgf_order, ordered=True
    )

    alt.data_transformers.enable("vegafusion")
    status_colors = ["#4477AA", "#DDAA33", "#BB5566"]
    chart = (
        alt.Chart(df_summary, width=chart_width, height=300)
        .mark_bar(stroke="#4d4d4d", strokeWidth=1)
        .encode(
            y=alt.Y(
                "mgf:N",
                sort=mgf_order,
                title="Group",
                axis=alt.Axis(
                    labelExpr='replace(replace(replace(datum.value, "[", ""), "]", ""), /_/g, " ")',
                    labelFontSize=11 * 2,
                    titleFontSize=13 * 2,
                    labelFont="Arial",
                    titleFont="Arial",
                    ticks=False,
                    domain=False,
                    labelPadding=6,
                ),
                scale=alt.Scale(domain=mgf_order),
            ),
            x=alt.X(
                "count:Q",
                title="Count",
                stack="zero",
                axis=alt.Axis(
                    labelFontSize=11 * 2,
                    titleFontSize=13 * 2,
                    labelFont="Arial",
                    titleFont="Arial",
                    gridColor="#4d4d4d",
                    gridOpacity=0.5,
                ),
            ),
            color=alt.Color(
                "status:N",
                scale=alt.Scale(domain=statuses, range=status_colors),
                title="",
                legend=alt.Legend(
                    labelFontSize=11 * 2,
                    titleFontSize=12 * 2,
                    labelFont="Arial",
                    titleFont="Arial",
                    labelLimit=1000,
                    titleLimit=1000,
                ),
            ),
            tooltip=[
                "mgf:N",
                alt.Tooltip("status:N", title="Status"),
                alt.Tooltip("count:Q", title="n"),
                alt.Tooltip("total_count:Q", title="total"),
                alt.Tooltip("percent:Q", title="%", format=".0%"),
            ],
        )
    )

    max_count = df_summary["count"].max()
    text = (
        alt.Chart(df_summary, width=chart_width, height=300)
        .mark_text(align="center", color="#767676", fontSize=10 * 2, font="Arial")
        .encode(
            y=alt.Y("mgf:N", sort=mgf_order, scale=alt.Scale(domain=mgf_order)),
            x=alt.X("x_center:Q"),
            detail="status:N",
            text=alt.Text("percent:Q", format=".0%"),
        )
        .transform_filter(
            "datum.percent >= 0.05 && datum.count >= %f" % (0.10 * max_count)
        )
    )

    panel_label_chart = (
        alt.Chart(pd.DataFrame({"label": [panel_label]}))
        .mark_text(
            align="left",
            baseline="top",
            fontSize=18 * 2,
            fontWeight="bold",
            dx=-200,
            dy=-30,
            font="Arial",
            color="#767676",
        )
        .encode(x=alt.value(0), y=alt.value(0), text="label:N")
    )

    # Add subtitle
    title_chart = (
        alt.Chart(pd.DataFrame({"title": ["Formula identification"]}))
        .mark_text(
            align="center",
            fontSize=14 * 2,
            fontWeight="bold",
            dy=-15,
            font="Arial",
            color="#767676",
        )
        .encode(x=alt.value(chart_width / 2), y=alt.value(0), text="title:N")
    )

    return panel_label_chart + title_chart + chart + text


@app.function
def make_match_plot(
    df: pl.DataFrame,
    metric: str = "estimated_prob",
    n_bins: int = 10,
    panel_label: str = "B",
    chart_width: int = 500,
    show_legend: bool = True,
):
    import pandas as pd

    if df.is_empty() or metric not in df.columns:
        return (
            alt.Chart(pl.DataFrame({"note": ["No data"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )
    df_matches = (
        df.lazy()
        .filter((pl.col("formula_found") == True) & (pl.col("formula_match") == True))
        .filter(pl.col(metric).is_not_null())
        .select(["mgf", metric])
        .collect()
    )
    if df_matches.is_empty():
        return (
            alt.Chart(pl.DataFrame({"note": ["No data"]}).to_pandas())
            .mark_text()
            .encode(text="note:N")
        )
    df_matches_pd = df_matches.to_pandas()
    non_null_values = df_matches_pd[metric].dropna().values
    unique_non_null = np.unique(non_null_values)
    if len(unique_non_null) < n_bins:
        n_bins = max(2, len(unique_non_null))
    all_ints = np.all(np.isclose(non_null_values, np.round(non_null_values), atol=1e-8))
    if all_ints and len(non_null_values) > 0:
        min_val, max_val = (
            int(np.nanmin(non_null_values)),
            int(np.nanmax(non_null_values)),
        )
        bin_edges = np.unique(np.linspace(0, max_val + 1, n_bins + 1, dtype=int))
        if len(bin_edges) < 2:
            bin_edges = np.array([0, max_val + 1])
        pad = max(len(str(bin_edges[0])), len(str(bin_edges[-1])))
        bin_labels = [
            f"{str(bin_edges[i]).zfill(pad)}–{str(bin_edges[i+1]-1).zfill(pad)}"
            for i in range(len(bin_edges) - 1)
        ]
        binned = pd.cut(
            df_matches_pd[metric],
            bins=bin_edges,
            labels=bin_labels,
            include_lowest=True,
            right=False,
        )
    else:
        if len(non_null_values) > 0:
            min_val, max_val = (
                float(np.nanmin(non_null_values)),
                float(np.nanmax(non_null_values)),
            )
            bins = (
                np.round(np.linspace(0, 1, n_bins + 1), 2)
                if min_val >= 0 and max_val <= 1
                else np.linspace(min(0, min_val), max_val, n_bins + 1)
            )
            bin_labels = [
                f"{bins[i]:.1f}–{bins[i+1]:.1f}" for i in range(len(bins) - 1)
            ]
            binned = pd.cut(
                df_matches_pd[metric],
                bins=bins,
                labels=bin_labels,
                include_lowest=True,
                right=True,
            )
        else:
            bin_labels = ["0.0–1.0"]
            binned = pd.Series(
                ["0.0–1.0"] * len(df_matches_pd), index=df_matches_pd.index
            )
    df_matches_pd["metric_bin"] = pd.Categorical(
        binned.astype(object), categories=bin_labels, ordered=True
    )
    df_binned = (
        df_matches_pd.groupby(["mgf", "metric_bin"], observed=True)
        .size()
        .reset_index(name="count")
    )
    all_combos = pd.MultiIndex.from_product(
        [df_matches_pd["mgf"].unique(), bin_labels], names=["mgf", "metric_bin"]
    ).to_frame(index=False)
    df_binned = all_combos.merge(
        df_binned, on=["mgf", "metric_bin"], how="left"
    ).fillna({"count": 0})
    df_binned["total_count"] = df_binned.groupby("mgf")["count"].transform("sum")
    df_binned["percent"] = df_binned["count"] / df_binned["total_count"]
    df_binned["cumsum_count"] = df_binned.groupby("mgf")["count"].cumsum()
    df_binned["x_center"] = df_binned["cumsum_count"] - df_binned["count"] / 2

    # build mgf order and make sure mgf column is a pandas categorical to carry the order
    mgf_order = get_mgf_order(df_binned)
    df_binned["mgf"] = pd.Categorical(
        df_binned["mgf"].astype(str), categories=mgf_order, ordered=True
    )

    from cmcrameri import cm

    bin_palette = cmap_to_hex_list(cm.batlow, len(bin_labels))
    max_count = df_binned["count"].max()
    alt.data_transformers.enable("vegafusion")
    METRIC_DISPLAY_NAMES = {
        "estimated_prob": "Estimated Probability",
        "explained_ms2_peak": "Number of Explained MS2 Peaks",
        "explained_peak_fraction": "Explained Peak Fraction",
        "explained_intensity": "Explained Intensity Fraction",
    }
    metric_display = METRIC_DISPLAY_NAMES.get(metric, metric.replace("_", " ").title())

    legend_config = (
        alt.Legend(
            title="",
            labelFontSize=11 * 2,
            titleFontSize=12 * 2,
            labelFont="Arial",
            titleFont="Arial",
            labelLimit=1000,
            titleLimit=1000,
        )
        if show_legend
        else None
    )

    chart = (
        alt.Chart(df_binned, width=chart_width, height=300)
        .mark_bar(stroke="#4d4d4d", strokeWidth=1)
        .encode(
            y=alt.Y(
                "mgf:N",
                title="Group",
                sort=mgf_order,
                axis=alt.Axis(
                    labelExpr='replace(replace(replace(datum.value, "[", ""), "]", ""), /_/g, " ")',
                    labelFontSize=11 * 2,
                    titleFontSize=13 * 2,
                    labelFont="Arial",
                    titleFont="Arial",
                    ticks=False,
                    domain=False,
                    labelPadding=6,
                ),
                scale=alt.Scale(domain=mgf_order),
            ),
            x=alt.X(
                "count:Q",
                title="Count",
                stack="zero",
                axis=alt.Axis(
                    labelFontSize=11 * 2,
                    titleFontSize=13 * 2,
                    labelFont="Arial",
                    titleFont="Arial",
                    gridColor="#4d4d4d",
                    gridOpacity=0.5,
                ),
            ),
            color=alt.Color(
                "metric_bin:N",
                title=metric_display,
                scale=alt.Scale(range=bin_palette),
                sort="ascending",
                legend=legend_config,
            ),
            order=alt.Order("metric_bin:N", sort="ascending"),
            tooltip=[
                "mgf:N",
                alt.Tooltip("metric_bin:N", title=metric_display),
                alt.Tooltip("count:Q", title="n"),
                alt.Tooltip("total_count:Q", title="total"),
                alt.Tooltip("percent:Q", title="%", format=".0%"),
            ],
        )
    )
    text = (
        alt.Chart(df_binned, width=chart_width, height=300)
        .mark_text(align="center", color="#767676", fontSize=10 * 2, font="Arial")
        .encode(
            y=alt.Y("mgf:N", sort=mgf_order, scale=alt.Scale(domain=mgf_order)),
            x=alt.X("x_center:Q"),
            detail="metric_bin:N",
            text=alt.Text("percent:Q", format=".0%"),
        )
        .transform_filter(
            "datum.percent >= 0.05 && datum.count >= %f" % (0.10 * max_count)
        )
    )
    panel_label_chart = (
        alt.Chart(pd.DataFrame({"label": [panel_label]}))
        .mark_text(
            align="left",
            baseline="top",
            fontSize=18 * 2,
            fontWeight="bold",
            dx=-200,
            dy=-30,
            font="Arial",
            color="#767676",
        )
        .encode(x=alt.value(0), y=alt.value(0), text="label:N")
    )

    # Add subtitle
    title_chart = (
        alt.Chart(pd.DataFrame({"title": [metric_display]}))
        .mark_text(
            align="center",
            fontSize=14 * 2,
            fontWeight="bold",
            dy=-15,
            font="Arial",
            color="#767676",
        )
        .encode(x=alt.value(chart_width / 2), y=alt.value(0), text="title:N")
    )

    return panel_label_chart + title_chart + chart + text


@app.function
def build_combined_msbuddy_figure(df: pl.DataFrame):
    chart_width = 450
    chart_a = make_match_plot(
        df,
        "estimated_prob",
        panel_label="A",
        chart_width=chart_width,
    )
    chart_b = make_match_plot(
        df,
        "explained_intensity",
        panel_label="B",
        chart_width=chart_width,
    )
    chart_c = make_match_plot(
        df,
        "explained_peak_fraction",
        panel_label="C",
        chart_width=chart_width,
    )
    chart_supp_1 = make_match_plot(
        df,
        "explained_ms2_peak",
        panel_label="",
        chart_width=chart_width,
    )
    chart_supp_2 = make_status_plot(
        df,
        panel_label="",
        chart_width=chart_width,
    )

    main_figure = (
        alt.vconcat(chart_a, chart_b, chart_c)
        .properties(background="none")
        .configure_view(strokeWidth=0, fill=None)
        .configure_axis(
            labelFontSize=11,
            titleFontSize=13,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            gridColor="#4d4d4d",
            gridOpacity=0.5,
        )
        .configure_legend(
            labelFontSize=11,
            titleFontSize=12,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            orient="bottom",
            columns=4,
        )
    )

    supp_figure_1 = (
        chart_supp_1.properties(background="none")
        .configure_view(strokeWidth=0, fill=None)
        .configure_axis(
            labelFontSize=11,
            titleFontSize=13,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            gridColor="#4d4d4d",
            gridOpacity=0.5,
        )
        .configure_legend(
            labelFontSize=11,
            titleFontSize=12,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            orient="bottom",
            columns=4,
        )
    )

    supp_figure_2 = (
        chart_supp_2.properties(background="none")
        .configure_view(strokeWidth=0, fill=None)
        .configure_axis(
            labelFontSize=11,
            titleFontSize=13,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            gridColor="#4d4d4d",
            gridOpacity=0.5,
        )
        .configure_legend(
            labelFontSize=11,
            titleFontSize=12,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            orient="bottom",
            columns=1,
        )
    )
    return main_figure, supp_figure_1, supp_figure_2


@app.cell
def _():
    mgf_path = Path("data/multims2_spectra.mgf")
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
    df_filtered = df.lazy().filter(pl.col("smiles").is_not_null()).collect()

    print(f"Total rows before filtering: {len(df)}")
    print(f"Total rows after filtering (smiles not null): {len(df_filtered)}")
    print(f"Rows with null smiles: {len(df) - len(df_filtered)}")

    def summary(df: pl.DataFrame):
        for col in df.columns:
            dtype = df[col].dtype
            null_count = df[col].null_count()
            if dtype == pl.Boolean:
                true_count = df[col].sum()
                false_count = df.shape[0] - true_count - null_count
                print(
                    f"{col} ({dtype}): {null_count} nulls, {true_count} True, {false_count} False"
                )
            else:
                unique_count = df[col].n_unique()
                print(f"{col} ({dtype}): {null_count} nulls, {unique_count} unique")

    summary(df_filtered)
    df_filtered
    return (df_filtered,)


@app.cell
def _plot(df_filtered):
    df_filtered.write_csv("scratch/msbuddy_summary_results.tsv", separator="\t")
    main_figure, supp_figure_1, supp_figure_2 = build_combined_msbuddy_figure(
        df_filtered
    )
    os.makedirs("figures", exist_ok=True)
    main_figure.save("figures/msbuddy_main.svg", format="svg", background=None)
    main_figure.save("figures/msbuddy_main.pdf", format="pdf", background=None)
    main_figure.save("figures/msbuddy_main.png", format="png", background=None)
    supp_figure_1.save("figures/msbuddy_supp_1.svg", format="svg", background=None)
    supp_figure_1.save("figures/msbuddy_supp_1.pdf", format="pdf", background=None)
    supp_figure_1.save("figures/msbuddy_supp_1.png", format="png", background=None)
    supp_figure_2.save("figures/msbuddy_supp_2.svg", format="svg", background=None)
    supp_figure_2.save("figures/msbuddy_supp_2.pdf", format="pdf", background=None)
    supp_figure_2.save("figures/msbuddy_supp_2.png", format="png", background=None)
    return main_figure, supp_figure_1, supp_figure_2


@app.cell
def _show_main(main_figure):
    main_figure
    return


@app.cell
def _show_supp(supp_figure_1):
    supp_figure_1
    return


@app.cell
def _show_supp(supp_figure_2):
    supp_figure_2
    return


if __name__ == "__main__":
    app.run()
