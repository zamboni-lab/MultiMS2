# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "marimo",
#     "altair",
#     "matchms",
#     "numpy",
#     "pandas",
#     "polars",
#     "pyarrow",
#     "rdkit",
#     "simple_parsing",
#     "vl-convert-python",
# ]
# ///

import marimo

__generated_with = "0.18.1"
app = marimo.App(width="full")

with app.setup:
    import logging
    import marimo
    import re
    import altair as alt
    import numpy as np
    import polars as pl
    from dataclasses import dataclass, field
    from matchms.importing import load_from_mgf
    from pathlib import Path
    from rdkit import Chem
    from simple_parsing import ArgumentParser
    from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

    logger = logging.getLogger("multims_upset")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"),
        )
        logger.addHandler(handler)

    @dataclass
    class Settings:
        mgf_path: str = field(
            default="data/multims2_spectra.mgf",
            metadata={"help": "Path to consolidated spectra MGF file"},
        )
        top_n_sets: int = field(
            default=12,
            metadata={
                "help": "Maximum number of sets (groups) to retain before intersection ranking (unused if all)",
            },
        )
        top_n_intersections: int = field(
            default=20,
            metadata={"help": "Show only the largest N intersections"},
        )

    parser = ArgumentParser()

    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        """Return settings from CLI or default Settings when running in notebook.

        Returns
        -------
        Settings
            Parsed runtime settings, or defaults when running in a notebook.
        """
        if marimo.running_in_notebook():
            logger.info("Running in notebook: using default Settings()")
            return Settings()
        args = parser.parse_args()
        return args.settings

    settings = parse_args()

    SANITIZE_RE = re.compile(r"[^a-z0-9.+-]", flags=re.IGNORECASE)
    MULT_UNDER_RE = re.compile(r"_+")


@app.function
def sanitize_local(v: Any) -> str:
    """Sanitize a metadata value into a normalized group token (lowercase, safe).

    Parameters
    ----------
    v : Any
        V.

    Returns
    -------
    str
        Sanitized lowercase token safe for group labeling.
    """
    if v is None:
        return "na"
    s = str(v).strip().lower()
    s = SANITIZE_RE.sub("_", s)
    s = MULT_UNDER_RE.sub("_", s).strip("_")
    return s or "na"


@app.function
def meta_lookup(md: dict, keys: Iterable[str]) -> Optional[str]:
    """Return the first non-empty metadata value matching any key in `keys`.

    Parameters
    ----------
    md : dict
        Md.
    keys : Iterable[str]
        Keys.

    Returns
    -------
    Optional[str]
        First non-empty metadata value for the provided keys, or ``None``.
    """
    for k in keys:
        v = md.get(k)
        if v not in (None, ""):
            return v
    return None


@app.function
def smiles_to_inchikey_first_layer(smiles: str) -> Optional[str]:
    """Convert SMILES to first 14 chars of InChIKey (or None on failure).

    Parameters
    ----------
    smiles : str
        Smiles.

    Returns
    -------
    Optional[str]
        First-layer InChIKey (14 chars), or ``None`` when conversion fails.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        ik = Chem.MolToInchiKey(mol)
        return ik[:14]
    except Exception as exc:
        # Do not spam logs for many failures, keep debug level
        logger.debug("Failed to convert SMILES to InChIKey: %s -> %s", smiles, exc)
        return None


@app.function
def read_consolidated_mgf(settings: Settings):
    """Read a consolidated MGF and compute group memberships.

        Group label: "{fragmentation_method}_{collision_energy}_{ionmode}" (sanitized).

    Parameters
    ----------
    settings : Settings
        Settings.
    """
    mgf_path = Path(settings.mgf_path)
    if not mgf_path.is_file():
        logger.error("MGF file not found: %s", mgf_path)
        return {}, [], {}, [], set()

    logger.info("Loading spectra from %s", mgf_path)
    try:
        spectra = list(load_from_mgf(str(mgf_path)))
    except Exception as exc:
        logger.exception("Failed to load MGF: %s", exc)
        return {}, [], {}, [], set()

    # metadata synonyms
    ce_keys = ["collision_energy", "collisionenergy", "collisionenergy_ev", "ce"]
    frag_keys = [
        "fragmentation_method",
        "fragmentationmethod",
        "frag_method",
        "fragmode",
    ]
    ion_keys = ["ionmode", "ion_mode", "polarity"]

    group_inchikeys: Dict[str, Set[str]] = {}
    group_adduct_inchikey: Dict[str, Set[Tuple[Optional[str], str]]] = {}
    smiles2ik: Dict[str, Optional[str]] = {}
    triplet_set: Set[Tuple[str, str, str]] = set()

    for spec in spectra:
        md = getattr(spec, "metadata", {}) or {}
        smiles = md.get("smiles")
        adduct = md.get("adduct")
        ce = meta_lookup(md, ce_keys)
        frag = meta_lookup(md, frag_keys)
        ion = meta_lookup(md, ion_keys)

        group = "_".join(
            (sanitize_local(frag), sanitize_local(ce), sanitize_local(ion)),
        )

        # ensure group entries exist
        group_inchikeys.setdefault(group, set())
        group_adduct_inchikey.setdefault(group, set())

        if not smiles:
            continue

        if smiles not in smiles2ik:
            smiles2ik[smiles] = smiles_to_inchikey_first_layer(smiles)
        ik14 = smiles2ik[smiles]
        if not ik14:
            continue

        group_inchikeys[group].add(ik14)

        if adduct:
            group_adduct_inchikey[group].add((adduct, ik14))

        if adduct and ce:
            triplet_set.add((adduct, ik14, str(ce)))

    all_inchikeys = sorted({ik for ik in smiles2ik.values() if ik})
    all_adduct_inchikey = sorted(
        {(a, k) for pairs in group_adduct_inchikey.values() for (a, k) in pairs},
        key=lambda x: (x[0] or "", x[1]),
    )

    logger.info(
        "Parsed %d spectra -> %d groups, %d unique InChIKey14, %d adduct-connectivity pairs",
        len(spectra),
        len(group_inchikeys),
        len(all_inchikeys),
        len(all_adduct_inchikey),
    )

    return (
        group_inchikeys,
        all_inchikeys,
        group_adduct_inchikey,
        all_adduct_inchikey,
        triplet_set,
    )


@app.function
def create_upset_data(group_items: Dict[str, Set], all_items: List, item_label: str):
    """Create a membership DataFrame (polars) and ordered group names from mapping group -> set(items).

    Parameters
    ----------
    group_items : Dict[str, Set]
        Group items.
    all_items : List
        All items.
    item_label : str
        Item label.
    """
    data_dict: Dict[str, List[int]] = {}
    group_names = list(group_items.keys())
    group_sizes: List[Tuple[str, int]] = []

    for g in group_names:
        s = group_items[g]
        column = [1 if item in s else 0 for item in all_items]
        data_dict[g] = column
        group_sizes.append((g, sum(column)))

    # sort groups by size desc then name
    group_sizes_sorted = sorted(group_sizes, key=lambda x: (-x[1], x[0]))

    logger.info("TOTAL: %d unique %s(s)", len(all_items), item_label)
    for name, size in group_sizes_sorted:
        logger.info("%s: %d %s(s)", name, size, item_label)

    # If items are triplets (adduct, connectivity, energy) report their count
    if all_items and isinstance(all_items[0], tuple) and len(all_items[0]) >= 3:
        triplets = {(it[0], it[1], it[2]) for it in all_items}
        logger.info(
            "TOTAL: %d unique Adduct-Connectivity-Energy triplet(s)",
            len(triplets),
        )

    # Use polars DataFrame for downstream consumption
    return pl.DataFrame(data_dict), group_names


@app.function
def filter_upset_data(data: pl.DataFrame, group_names: List[str], top_n: int):
    """Select the top_n groups by membership size and return the pandas slice and filtered group names.

    Parameters
    ----------
    data : pl.DataFrame
        Data.
    group_names : List[str]
        Group names.
    top_n : int
        Top n.
    """
    pdf = data.to_pandas()
    set_sums = pdf.sum(axis=0)
    # top N indices (handles when top_n > available)
    n = min(top_n, len(set_sums)) if top_n > 0 else len(set_sums)
    top_sets = set_sums.nlargest(n).index.tolist()
    pdf_sub = pdf[top_sets]
    group_names_sub = [g for g in group_names if g in top_sets]
    return pdf_sub, group_names_sub


@app.function
def membership_top_intersections(pdf, top_n: int):
    """Compute top N intersections from membership matrix (pandas DataFrame: items x sets).

    Parameters
    ----------
    pdf : Any
        Pdf.
    top_n : int
        Top n.
    """
    # pdf: rows = items, columns = sets (0/1)
    counts: dict[tuple[int, ...], int] = {}
    cols = list(pdf.columns)
    arr = pdf.to_numpy(dtype=int)
    for row in arr:
        idxs = tuple(sorted(np.nonzero(row)[0].tolist()))
        if not idxs:
            continue
        counts[idxs] = counts.get(idxs, 0) + 1

    # sort by count desc
    ranked = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
    if top_n > 0:
        ranked = ranked[:top_n]

    # get active set indices and names
    active_indices = sorted({i for key, _ in ranked for i in key})
    active_sets = [cols[i] for i in active_indices]

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
    pdf,
    top_n_intersections: int,
    panel_label: str,
    title_text: str,
):
    """Build an Altair upset-like plot for the top intersections.

        The function returns an Altair Chart. It expects `pdf` as a pandas membership DataFrame
        (rows=items, columns=sets) and uses membership_top_intersections to select the top intersections.

    Parameters
    ----------
    pdf : Any
        Pdf.
    top_n_intersections : int
        Top n intersections.
    panel_label : str
        Panel label.
    title_text : str
        Title text.
    """
    records, active_sets = membership_top_intersections(pdf, top_n_intersections)
    if not records:
        # return a small text chart if empty
        return (
            alt.Chart(pl.DataFrame({"note": ["No intersections"]}).to_pandas())
            .mark_text()
            .encode(text="note")
        )

    import pandas as _pd

    inter_df = _pd.DataFrame(records)

    # Build matrix DataFrame: one row per intersection x set
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

    # compute set sizes
    set_sizes_full = pdf.sum(axis=0).reset_index()
    set_sizes_full.columns = ["set", "size"]
    set_sizes = set_sizes_full[set_sizes_full["set"].isin(active_sets)].copy()

    # Desired ordering (domain specific). Keep any sets not listed appended at end.
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
    set_order = [s for s in desired_order if s in active_sets]
    set_order.extend([s for s in active_sets if s not in set_order])

    # preserve categorical order
    set_sizes["set"] = _pd.Categorical(
        set_sizes["set"],
        categories=set_order,
        ordered=True,
    )
    set_sizes = set_sizes.sort_values("set")

    # sizes and layout
    inter_width = max(400, len(records) * 28)
    row_height = 24
    matrix_height = max(80, row_height * len(set_order))

    # Styling constants
    base_color = "#4477AA"
    grid_color = "#4d4d4d"

    max_size = int(set_sizes["size"].max()) if not set_sizes["size"].empty else 0
    set_sizes = set_sizes.assign(maxsize=max_size)

    x_inter = alt.X(
        "intersection:N",
        sort=alt.SortField(field="count", order="descending"),
        axis=alt.Axis(labels=False, ticks=False, title=None),
    )

    # Bars (top)
    bars = (
        alt.Chart(inter_df, width=inter_width, height=180)
        .mark_bar(color=base_color, stroke=grid_color, strokeWidth=1)
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
                    labelColor=grid_color,
                    titleColor=grid_color,
                ),
            ),
            tooltip=["intersection:N", alt.Tooltip("count:Q", title="size")],
        )
    )

    bar_labels = (
        alt.Chart(inter_df, width=inter_width, height=180)
        .mark_text(
            dy=-7,
            color=grid_color,
            fontSize=11 * 1.5,
            font="Arial",
            fontWeight="normal",
        )
        .encode(x=x_inter, y=alt.Y("count:Q"), text=alt.Text("count:Q"))
    )

    # Panel label (left) and title (center)
    panel_label_chart = (
        alt.Chart(_pd.DataFrame({"label": [panel_label]}))
        .mark_text(
            align="left",
            baseline="top",
            fontSize=18 * 2,
            fontWeight="bold",
            dx=-200,
            dy=-48,
            font="Arial",
            color=grid_color,
        )
        .encode(x=alt.value(0), y=alt.value(0), text="label:N")
    )

    title_chart = (
        alt.Chart(_pd.DataFrame({"title": [title_text]}))
        .mark_text(
            align="center",
            fontSize=14 * 2,
            fontWeight="bold",
            dy=-30,
            font="Arial",
            color=grid_color,
        )
        .encode(x=alt.value(inter_width / 2), y=alt.value(0), text="title:N")
    )

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
            labelColor=grid_color,
            titleColor=grid_color,
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

    # Right panel: group size bars + labels
    bg_bars = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_bar(fill="none", stroke=grid_color, strokeWidth=0.8)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X(
                "maxsize:Q",
                title="Group Size",
                scale=alt.Scale(domain=[0, max_size]),
            ),
        )
    )

    size_bars = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_bar(color=base_color, stroke=grid_color, strokeWidth=1)
        .encode(
            y=alt.Y("set:N", sort=set_order, axis=None),
            x=alt.X("size:Q", title="Group Size"),
            tooltip=["set:N", alt.Tooltip("size:Q", title="set size")],
        )
    )

    size_labels = (
        alt.Chart(set_sizes, width=220, height=matrix_height)
        .mark_text(
            align="left",
            dx=4,
            color=grid_color,
            fontSize=11 * 1.5,
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

    lower = alt.hconcat(matrix, right_panel, spacing=20).resolve_scale(y="shared")

    chart = alt.vconcat(
        panel_label_chart + title_chart + bars + bar_labels,
        lower,
    ).resolve_scale(x="shared", color="independent")

    return chart


@app.function
def build_combined_upset_figure(
    pd_sub_inchikeys,
    pd_sub_pairs,
    top_n_intersections: int,
):
    """Compose panels A and B into a single combined figure (Altair).

    Parameters
    ----------
    pd_sub_inchikeys : Any
        Pd sub inchikeys.
    pd_sub_pairs : Any
        Pd sub pairs.
    top_n_intersections : int
        Top n intersections.
    """
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
        alt.vconcat(chart_a, chart_b, spacing=50)
        .properties(background="none")
        .configure_view(strokeWidth=0, fill=None)
        .configure_axis(
            labelFontSize=11 * 2,
            titleFontSize=13 * 2,
            labelFont="Arial",
            titleFont="Arial",
            labelColor=grid_color,
            titleColor=grid_color,
            gridColor=grid_color,
            gridOpacity=0.5,
        )
        .configure_title(font="Arial", fontSize=14 * 2, color=grid_color)
    )

    return combined


@app.function
def log_triplet_count(triplet_set: Set[Tuple[str, str, str]]) -> int:
    """Log number of unique adduct-connectivity-energy triplets and return count.

    Parameters
    ----------
    triplet_set : Set[Tuple[str, str, str]]
        Triplet set.

    Returns
    -------
    int
        Number of unique adduct-connectivity-energy triplets.
    """
    count = len(triplet_set)
    logger.info("TOTAL: %d unique Adduct-Connectivity-Energy triplet(s)", count)
    return count


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
        pd_sub_inchikeys,
        pd_sub_pairs,
        settings.top_n_intersections,
    )

    outdir = Path("figures")
    outdir.mkdir(parents=True, exist_ok=True)
    # Save to multiple formats; catch errors separately so one failing format doesn't stop others
    save_path = outdir / "upset_main"
    try:
        combined_chart.save(
            str(save_path.with_suffix(".svg")),
            format="svg",
            background=None,
        )
        combined_chart.save(
            str(save_path.with_suffix(".pdf")),
            format="pdf",
            background=None,
        )
        combined_chart.save(
            str(save_path.with_suffix(".png")),
            format="png",
            background=None,
        )
        logger.info("Saved combined upset figures to %s.[svg|pdf|png]", save_path)
    except Exception as exc:
        logger.exception("Failed to save one or more chart outputs: %s", exc)
    return (combined_chart,)


@app.cell
def _show_combined(combined_chart):
    combined_chart
    return


if __name__ == "__main__":
    app.run()
