# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "polars",
#     "pyarrow",
#     "numpy",
#     "matchms",
#     "rdkit",
#     "simple_parsing",
#     "vl-convert-python",
#     "altair",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    import marimo as mo
    import polars as pl
    import numpy as np
    import os
    import logging
    import altair as alt
    from simple_parsing import ArgumentParser
    from matchms.importing import load_from_mgf
    from rdkit import Chem

    parser = ArgumentParser()

    @dataclass
    class Settings:
        mgf_path: str = field(
            default="data/multims2_spectra.mgf",
            metadata={"help": "Path to consolidated spectra MGF file"},
        )
        top_n_sets: int = field(
            default=12,
            metadata={
                "help": "Maximum number of sets (groups) to retain before intersection ranking (unused if all)"
            },
        )
        top_n_intersections: int = field(
            default=20,
            metadata={"help": "Show only the largest N intersections"},
        )

    parser.add_arguments(Settings, dest="settings")

    def parse_args():
        if mo.running_in_notebook():
            return Settings()
        else:
            args = parser.parse_args()
            return args.settings

    settings = parse_args()

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )


@app.function
def smiles_to_inchikey_first_layer(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            inchikey = Chem.MolToInchiKey(mol)
            return inchikey[:14]
        return None
    except Exception:
        return None


@app.function
def read_consolidated_mgf(settings: Settings):
    """Read a single consolidated MGF and derive group memberships.

    Group label logic replicates viz_msbuddy: group = f"{fragmentation_method}_{collision_energy}_{ionmode}" (sanitized).
    Returns dictionaries for upset:
      group_inchikeys: group -> {InChIKey14}
      all_inchikeys: sorted list of unique InChIKey14
      group_adduct_inchikey: group -> {(adduct, InChIKey14)}
      all_adduct_inchikey: sorted list of all (adduct, InChIKey14) pairs
    """

    # Local sanitize (avoid execution-order issues if global _sanitize not yet defined)
    def _sanitize_local(v):
        if v is None:
            return "na"
        v = str(v).strip().lower()
        import re as _re

        v = _re.sub(r"[^a-z0-9.+-]", "_", v)
        v = _re.sub(r"__+", "_", v).strip("_")
        return v or "na"

    mgf_path = settings.mgf_path
    if not os.path.isfile(mgf_path):
        logging.error("MGF file not found: %s", mgf_path)
        return {}, [], {}, []

    spectra = list(load_from_mgf(mgf_path))
    group_inchikeys: dict[str, set[str]] = {}
    group_adduct_inchikey: dict[str, set[tuple[str | None, str]]] = {}

    smiles2ik: dict[str, str | None] = {}

    # metadata synonyms
    ce_keys = ["collision_energy", "collisionenergy", "collisionenergy_ev", "ce"]
    frag_keys = [
        "fragmentation_method",
        "fragmentationmethod",
        "frag_method",
        "fragmode",
    ]
    ion_keys = ["ionmode", "ion_mode", "polarity"]

    def meta_lookup(md: dict, keys):
        for k in keys:
            v = md.get(k)
            if v not in (None, ""):
                return v
        return None

    triplet_set = set()
    for spectrum in spectra:
        md = getattr(spectrum, "metadata", {}) or {}
        smiles = md.get("smiles")
        adduct = md.get("adduct")
        ce = meta_lookup(md, ce_keys)
        frag = meta_lookup(md, frag_keys)
        ion = meta_lookup(md, ion_keys)
        group = "_".join(
            [_sanitize_local(frag), _sanitize_local(ce), _sanitize_local(ion)]
        )

        group_inchikeys.setdefault(group, set())
        group_adduct_inchikey.setdefault(group, set())

        if smiles:
            if smiles not in smiles2ik:
                smiles2ik[smiles] = smiles_to_inchikey_first_layer(smiles)
            ik14 = smiles2ik[smiles]
            if ik14:
                group_inchikeys[group].add(ik14)
                # add (adduct, ik14) for plotting
                if adduct:
                    group_adduct_inchikey[group].add((adduct, ik14))
                # add (adduct, ik14, ce) for summary
                if adduct and ce:
                    triplet_set.add((adduct, ik14, str(ce)))

    all_inchikeys = sorted({ik for ik in smiles2ik.values() if ik})
    all_adduct_inchikey = sorted(
        {(a, k) for pairs in group_adduct_inchikey.values() for (a, k) in pairs},
        key=lambda x: (x[0] or "", x[1]),
    )

    return (
        group_inchikeys,
        all_inchikeys,
        group_adduct_inchikey,
        all_adduct_inchikey,
        triplet_set,
    )


@app.function
def create_upset_data(group_items, all_items, item_label):
    data_dict = {}
    group_names = list(group_items.keys())
    group_sizes = []
    for group_name in group_names:
        s = group_items[group_name]
        data_dict[group_name] = [1 if item in s else 0 for item in all_items]
        group_sizes.append((group_name, sum(data_dict[group_name])))
    # Sort by size descending, then group name
    group_sizes_sorted = sorted(group_sizes, key=lambda x: (-x[1], x[0]))
    print(f"TOTAL: {len(all_items)} unique {item_label}(s)")
    for group_name, size in group_sizes_sorted:
        print(f"{group_name}: {size} {item_label}(s)")
    # Print number of unique (adduct, connectivity, energy) triplets if possible
    if all_items and isinstance(all_items[0], tuple) and len(all_items[0]) >= 3:
        # Assume tuple: (adduct, connectivity, energy) or similar
        triplets = set()
        for item in all_items:
            # Only use first three elements
            triplets.add((item[0], item[1], item[2]))
        print(f"TOTAL: {len(triplets)} unique Adduct-Connectivity-Energy triplet(s)")
    return pl.DataFrame(data_dict), group_names


@app.function
def filter_upset_data(data: pl.DataFrame, group_names: list[str], top_n: int):
    pdf = data.to_pandas()
    set_sums = pdf.sum(axis=0)
    # Take top N sets by membership size
    top_sets = set_sums.nlargest(min(top_n, len(set_sums))).index
    pdf_sub = pdf[top_sets]
    group_names_sub = [g for g in group_names if g in top_sets]
    return pdf_sub, group_names_sub


@app.function
def membership_top_intersections(pdf, top_n: int):
    """Compute top N intersections from a pandas membership matrix (items x sets)."""
    # pdf: rows = items, columns = sets (0/1)
    counts = {}
    cols = list(pdf.columns)
    arr = pdf.to_numpy(dtype=int)
    for row in arr:
        idxs = np.nonzero(row)[0]
        if len(idxs) == 0:
            continue
        key = tuple(sorted(idxs))
        counts[key] = counts.get(key, 0) + 1
    # build sorted list
    ranked = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    if top_n > 0:
        ranked = ranked[:top_n]
    # active set indices
    active_indices = sorted({i for key, _ in ranked for i in key})
    active_sets = [cols[i] for i in active_indices]
    # Build intersection dataframe
    records = []
    for rank_pos, (key, ct) in enumerate(ranked, start=1):
        set_names = [cols[i] for i in key]
        intersection_label = " ∩ ".join(set_names)
        records.append(
            {
                "intersection": intersection_label,
                "count": ct,
                "order": rank_pos,
                "sets": set_names,
            },
        )
    return records, active_sets


@app.function
def build_single_upset_chart(
    pdf, top_n_intersections: int, panel_label: str, title_text: str
):
    """Build a single upset chart with a panel label."""
    records, active_sets = membership_top_intersections(pdf, top_n_intersections)
    if not records:
        return (
            alt.Chart(pl.DataFrame({"note": ["No intersections"]}).to_pandas())
            .mark_text()
            .encode(text="note")
        )
    import pandas as _pd

    inter_df = _pd.DataFrame(records)
    matrix_rows = []
    for r in records:
        for s in active_sets:
            matrix_rows.append(
                {
                    "intersection": r["intersection"],
                    "set": s,
                    "present": 1 if s in r["sets"] else 0,
                    "count": r["count"],
                },
            )
    matrix_df = _pd.DataFrame(matrix_rows)
    set_sizes_full = pdf.sum(axis=0).reset_index()
    set_sizes_full.columns = ["set", "size"]
    set_sizes = set_sizes_full[set_sizes_full["set"].isin(active_sets)].copy()

    # Define custom order: cid 20, 40, 60 pos, ead 12, 16, 24 pos, then same for neg
    desired_order = [
        "cid_[20.0]_negative",
        "cid_[20.0]_positive",
        "cid_[40.0]_negative",
        "cid_[40.0]_positive",
        "cid_[60.0]_negative",
        "cid_[60.0]_positive",
        "ead_[12.0]_negative",
        "ead_[12.0]_positive",
        "ead_[16.0]_negative",
        "ead_[16.0]_positive",
        "ead_[24.0]_negative",
        "ead_[24.0]_positive",
    ]

    # Filter to only include sets that are active, maintaining desired order
    set_order = [s for s in desired_order if s in active_sets]
    # Add any active sets not in desired_order at the end
    set_order.extend([s for s in active_sets if s not in set_order])

    # Don't sort set_sizes by size - keep it in the custom order
    # Create a categorical column to preserve the order
    set_sizes["set"] = _pd.Categorical(
        set_sizes["set"], categories=set_order, ordered=True
    )
    set_sizes = set_sizes.sort_values("set")

    # Dynamic width and heights
    inter_width = max(400, len(records) * 28)
    row_height = 24
    matrix_height = max(80, row_height * len(set_order))

    # Paul Tol color scheme - vibrant blue with light neutral background
    base_color = "#4477AA"  # Paul Tol vibrant blue
    grid_color = "#4d4d4d"  # soft grid lines

    max_size = int(set_sizes["size"].max()) if len(set_sizes) else 0
    set_sizes = set_sizes.assign(maxsize=max_size)

    x_inter = alt.X(
        "intersection:N",
        sort=alt.SortField(field="count", order="descending"),
        axis=alt.Axis(labels=False, ticks=False, title=None),
    )

    # Bars (top)
    bars = (
        alt.Chart(inter_df, width=inter_width, height=180)
        .mark_bar(color=base_color, stroke="white", strokeWidth=1)
        .encode(
            x=x_inter,
            y=alt.Y(
                "count:Q",
                title="Intersection Size",
                axis=alt.Axis(
                    titleFontSize=13 * 2,
                    labelFontSize=11 * 2,
                    titleFont="Arial",
                    labelFont="Arial",
                    labelColor="#4d4d4d",
                    titleColor="#4d4d4d",
                ),
            ),
            tooltip=["intersection:N", alt.Tooltip("count:Q", title="size")],
        )
    )

    bar_labels = (
        alt.Chart(inter_df, width=inter_width, height=180)
        .mark_text(
            dy=-5, color="#4d4d4d", fontSize=11 * 1.5, font="Arial", fontWeight="normal"
        )
        .encode(x=x_inter, y=alt.Y("count:Q"), text=alt.Text("count:Q"))
    )

    # Panel label for top section - positioned at very left
    panel_label_chart = (
        alt.Chart(_pd.DataFrame({"label": [panel_label], "x": [0], "y": [0]}))
        .mark_text(
            align="left",
            baseline="top",
            fontSize=18 * 2,
            fontWeight="bold",
            dx=-200,
            dy=-48,
            font="Arial",
            color="#4d4d4d",
        )
        .encode(x=alt.value(0), y=alt.value(0), text="label:N")
    )

    # Title - positioned at the top center of the bars
    title_chart = (
        alt.Chart(_pd.DataFrame({"title": [title_text]}))
        .mark_text(
            align="center",
            fontSize=14 * 2,
            fontWeight="bold",
            dy=-30,
            font="Arial",
            color="#4d4d4d",
        )
        .encode(x=alt.value(inter_width / 2), y=alt.value(0), text="title:N")
    )

    # Shared y encoding
    y_sets = alt.Y(
        "set:N",
        sort=set_order,
        axis=alt.Axis(
            labelExpr='replace(replace(replace(datum.value, "[", ""), "]", ""), /_/g, " ")',
            title="Group",
            labelPadding=6,
            ticks=False,
            domain=False,
            titleFontSize=13 * 2,
            labelFontSize=11 * 2,
            titleFont="Arial",
            labelFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
        ),
    )

    matrix = (
        alt.Chart(matrix_df, width=inter_width, height=matrix_height)
        .mark_rect(stroke=grid_color, strokeWidth=0.8)
        .encode(
            x=x_inter,
            y=y_sets,
            color=alt.condition(
                alt.datum.present == 1,
                alt.value(base_color),
                alt.value("none"),
            ),
            tooltip=["set:N", "intersection:N", "present:Q"],
        )
    )

    # Background full-width bars for consistent row grid lines
    bg_bars = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_bar(fill="none", stroke=grid_color, strokeWidth=0.8)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X(
                "maxsize:Q",
                title="Group Size",
                scale=alt.Scale(domain=[0, max_size]),
                axis=alt.Axis(
                    titleFontSize=13 * 2,
                    labelFontSize=11 * 2,
                    titleFont="Arial",
                    labelFont="Arial",
                    labelColor="#4d4d4d",
                    titleColor="#4d4d4d",
                ),
            ),
        )
    )

    size_bars = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_bar(color=base_color, stroke="white", strokeWidth=1)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X(
                "size:Q",
                title="Group Size",
                scale=alt.Scale(domain=[0, max_size]),
                axis=alt.Axis(
                    titleFontSize=13 * 2,
                    labelFontSize=11 * 2,
                    titleFont="Arial",
                    labelFont="Arial",
                    labelColor="#4d4d4d",
                    titleColor="#4d4d4d",
                ),
            ),
            tooltip=["set:N", alt.Tooltip("size:Q", title="set size")],
        )
    )

    size_labels = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_text(
            align="left",
            dx=4,
            color="#4d4d4d",
            fontSize=11 * 2,
            font="Arial",
            fontWeight="normal",
        )
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X("size:Q"),
            text=alt.Text("size:Q"),
        )
    )

    right_panel = bg_bars + size_bars + size_labels

    lower = alt.hconcat(matrix, right_panel, spacing=0).resolve_scale(y="shared")

    chart = alt.vconcat(
        panel_label_chart + title_chart + bars + bar_labels, lower
    ).resolve_scale(x="shared", color="independent")
    return chart


@app.function
def build_combined_upset_figure(
    pd_sub_inchikeys, pd_sub_pairs, top_n_intersections: int
):
    """Build combined figure with panels A and B."""
    # Define grid color globally for consistency
    grid_color = "#4d4d4d"

    chart_a = build_single_upset_chart(
        pd_sub_inchikeys,
        top_n_intersections,
        "A",
        f"Connectivities: Top {min(top_n_intersections, len(pd_sub_inchikeys))} Intersections",
    )

    chart_b = build_single_upset_chart(
        pd_sub_pairs,
        top_n_intersections,
        "B",
        f"Adduct–Connectivities: Top {min(top_n_intersections, len(pd_sub_pairs))} Intersections",
    )

    combined = (
        alt.vconcat(chart_a, chart_b, spacing=0)
        .properties(background="none")
        .configure_view(strokeWidth=0, fill=None)
        .configure_axis(
            labelFontSize=11 * 2,
            titleFontSize=13 * 2,
            labelFont="Arial",
            titleFont="Arial",
            labelColor="#4d4d4d",
            titleColor="#4d4d4d",
            gridColor=grid_color,
            gridOpacity=0.5,
        )
        .configure_title(font="Arial", fontSize=14 * 2, color="#4d4d4d")
    )

    return combined


@app.function
def log_triplet_count(triplet_set):
    print(f"TOTAL: {len(triplet_set)} unique Adduct-Connectivity-Energy triplet(s)")
    return len(triplet_set)


@app.cell
def _read():
    (
        group_inchikeys,
        all_inchikeys,
        group_adduct_inchikey,
        all_adduct_inchikey,
        triplet_set,
    ) = read_consolidated_mgf(settings)
    log_triplet_count(triplet_set)
    return (
        all_adduct_inchikey,
        all_inchikeys,
        group_adduct_inchikey,
        group_inchikeys,
    )


@app.cell
def _prepare(
    all_adduct_inchikey,
    all_inchikeys,
    group_adduct_inchikey,
    group_inchikeys,
):
    data_inchikeys, names_inchikeys = create_upset_data(
        group_items=group_inchikeys,
        all_items=all_inchikeys,
        item_label="InChIKey",
    )
    data_pairs, names_pairs = create_upset_data(
        group_items=group_adduct_inchikey,
        all_items=all_adduct_inchikey,
        item_label="Adduct–Connectivity pair",
    )
    return data_inchikeys, data_pairs, names_inchikeys, names_pairs


@app.cell
def _filter(data_inchikeys, data_pairs, names_inchikeys, names_pairs):
    pd_sub_inchikeys, names_sub_inchikeys = filter_upset_data(
        data=data_inchikeys,
        group_names=names_inchikeys,
        top_n=settings.top_n_sets,
    )
    pd_sub_pairs, names_sub_pairs = filter_upset_data(
        data=data_pairs,
        group_names=names_pairs,
        top_n=settings.top_n_sets,
    )
    return pd_sub_inchikeys, pd_sub_pairs


@app.cell
def _plot(pd_sub_inchikeys, pd_sub_pairs):
    combined_chart = build_combined_upset_figure(
        pd_sub_inchikeys, pd_sub_pairs, settings.top_n_intersections
    )

    os.makedirs("figures", exist_ok=True)
    combined_chart.save("figures/combined_upset.svg", format="svg", background=None)
    combined_chart.save("figures/combined_upset.pdf", format="pdf", background=None)
    combined_chart.save("figures/combined_upset.png", format="png", background=None)
    return (combined_chart,)


@app.cell
def _show_combined(combined_chart):
    combined_chart
    return


if __name__ == "__main__":
    app.run()
