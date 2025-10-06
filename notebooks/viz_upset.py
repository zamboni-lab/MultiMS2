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

__generated_with = "0.16.3"
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
            default="scratch/consolidated_spectra.mgf",
            metadata={"help": "Path to consolidated spectra MGF file"},
        )
        top_n_sets: int = field(
            default=12,
            metadata={
                "help": "Maximum number of sets (groups) to retain before intersection ranking (unused if all)"
            },
        )
        top_n_intersections: int = field(
            default=32,
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

    for spectrum in spectra:
        md = getattr(spectrum, "metadata", {}) or {}
        smiles = md.get("smiles")
        adduct = md.get("adduct")
        # Build group label
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
                # only add (adduct, ik14) if adduct exists
                if adduct:
                    group_adduct_inchikey[group].add((adduct, ik14))

    all_inchikeys = sorted({ik for ik in smiles2ik.values() if ik})
    all_adduct_inchikey = sorted(
        {(a, k) for pairs in group_adduct_inchikey.values() for (a, k) in pairs},
        key=lambda x: (x[0] or "", x[1]),
    )

    return group_inchikeys, all_inchikeys, group_adduct_inchikey, all_adduct_inchikey


@app.function
def create_upset_data(group_items, all_items, item_label):
    data_dict = {}
    group_names = list(group_items.keys())
    for group_name in group_names:
        s = group_items[group_name]
        data_dict[group_name] = [1 if item in s else 0 for item in all_items]
        print(f"{group_name}: {sum(data_dict[group_name])} {item_label}(s)")
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
def build_upset_charts(pdf, top_n_intersections: int, title_prefix: str):
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
    set_sizes = set_sizes_full[set_sizes_full["set"].isin(active_sets)]
    set_sizes = set_sizes.sort_values(["size", "set"], ascending=[False, True])
    set_order = list(set_sizes["set"])

    # Dynamic width and heights
    inter_width = max(400, len(records) * 26)
    row_height = 22
    matrix_height = max(80, row_height * len(set_order))

    # Unified palette
    base_color = "#2E5B88"  # consistent deep blue
    absent_color = "#F1F4F7"

    max_size = int(set_sizes["size"].max()) if len(set_sizes) else 0
    set_sizes = set_sizes.assign(maxsize=max_size)

    x_inter = alt.X(
        "intersection:N",
        sort=alt.SortField(field="count", order="descending"),
        axis=alt.Axis(labels=False, ticks=False, title=None),
    )

    # Bars (top)
    bars = (
        alt.Chart(inter_df, width=inter_width, height=200)
        .mark_bar(color=base_color, stroke="white", strokeWidth=0.6)
        .encode(
            x=x_inter,
            y=alt.Y("count:Q", title="Intersection Size"),
            tooltip=["intersection:N", alt.Tooltip("count:Q", title="size")],
        )
        .properties(title=f"{title_prefix}: Top {len(records)} Intersections")
    )

    bar_labels = (
        alt.Chart(inter_df, width=inter_width, height=200)
        .mark_text(dy=-4, color="#222", fontSize=10)
        .encode(x=x_inter, y=alt.Y("count:Q"), text=alt.Text("count:Q"))
    )

    # Shared y encoding
    y_sets = alt.Y(
        "set:N",
        sort=set_order,
        axis=alt.Axis(title="Sets", labelPadding=4, ticks=False, domain=False),
    )

    matrix = (
        alt.Chart(matrix_df, width=inter_width, height=matrix_height)
        .mark_rect(stroke="#d0d5da", strokeWidth=0.5)
        .encode(
            x=x_inter,
            y=y_sets,
            color=alt.condition(
                alt.datum.present == 1,
                alt.value(base_color),
                alt.value(absent_color),
            ),
            tooltip=["set:N", "intersection:N", "present:Q"],
        )
    )

    # Background full-width bars for consistent row grid lines
    bg_bars = (
        alt.Chart(set_sizes, width=200, height=matrix_height)
        .mark_bar(fill=absent_color, stroke="#d0d5da", strokeWidth=0.5)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X(
                "maxsize:Q",
                title="Set Size",
                scale=alt.Scale(domain=[0, max_size]),
            ),
        )
    )

    size_bars = (
        alt.Chart(set_sizes, width=200, height=matrix_height)
        .mark_bar(color=base_color, stroke="white", strokeWidth=0.6)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X(
                "size:Q",
                title="Set Size",
                scale=alt.Scale(domain=[0, max_size]),
            ),
            tooltip=["set:N", alt.Tooltip("size:Q", title="set size")],
        )
    )

    size_labels = (
        alt.Chart(set_sizes, width=200, height=matrix_height)
        .mark_text(align="left", dx=3, color="#222", fontSize=10)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X("size:Q"),
            text=alt.Text("size:Q"),
        )
    )

    right_panel = bg_bars + size_bars + size_labels

    lower = alt.hconcat(matrix, right_panel, spacing=16).resolve_scale(y="shared")

    chart = (
        alt.vconcat(bars + bar_labels, lower)
        .resolve_scale(x="shared", color="independent")
        .configure_view(strokeWidth=0)
        .configure_axis(labelFontSize=11, titleFontSize=12)
    )
    return chart


@app.cell
def _read(settings):
    group_inchikeys, all_inchikeys, group_adduct_inchikey, all_adduct_inchikey = (
        read_consolidated_mgf(settings)
    )
    return group_inchikeys, all_inchikeys, group_adduct_inchikey, all_adduct_inchikey


@app.cell
def _prepare(
    group_inchikeys, all_inchikeys, group_adduct_inchikey, all_adduct_inchikey
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
def _filter(data_inchikeys, data_pairs, names_inchikeys, names_pairs, settings):
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
    return names_sub_inchikeys, names_sub_pairs, pd_sub_inchikeys, pd_sub_pairs


@app.cell
def _plot(
    names_sub_inchikeys, names_sub_pairs, pd_sub_inchikeys, pd_sub_pairs, settings
):
    chart_inchikeys = build_upset_charts(
        pd_sub_inchikeys, settings.top_n_intersections, "Connectivities"
    )
    chart_pairs = build_upset_charts(
        pd_sub_pairs, settings.top_n_intersections, "Adduct–Connectivities"
    )
    os.makedirs("figures", exist_ok=True)
    chart_inchikeys.save("figures/connectivities_upset.svg", format="svg")
    chart_pairs.save("figures/pairs_upset.svg", format="svg")
    chart_inchikeys.save("figures/connectivities_upset.pdf", format="pdf")
    chart_pairs.save("figures/pairs_upset.pdf", format="pdf")
    chart_inchikeys.save("figures/connectivities_upset.png", format="png")
    chart_pairs.save("figures/pairs_upset.png", format="png")
    return chart_inchikeys, chart_pairs


@app.cell
def _show_inchikey(chart_inchikeys):
    chart_inchikeys
    return


@app.cell
def _show_pairs(chart_pairs):
    chart_pairs
    return


if __name__ == "__main__":
    app.run()
