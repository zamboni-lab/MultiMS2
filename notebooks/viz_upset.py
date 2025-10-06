# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "altair_upset",
#     "polars",
#     "pyarrow",
#     "numpy",
#     "matchms",
#     "rdkit",
#     "simple_parsing",
#     "tmap-viz @ git+https://github.com/mvisani/tmap",
#     "scipy",
#     "faerun",
#     "molzip",
#     "tol-colors",
#     "cmcrameri",
#     "tqdm",
#     "vl-convert-python",
# ]
# ///

import marimo

__generated_with = "0.16.3"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass
    from dataclasses import field
    import marimo as mo
    import altair_upset as au
    import polars as pl
    import numpy as np
    import glob
    import os
    import logging
    from collections import defaultdict
    from simple_parsing import ArgumentParser
    from matchms.importing import load_from_mgf
    from rdkit import Chem
    from tqdm import tqdm
    import matplotlib.colors as mcolors
    import scipy.stats as ss
    from cmcrameri import cm
    from faerun import Faerun
    from molzip import ZipKNNGraph
    from tol_colors import high_contrast
    import tmap as tm

    parser = ArgumentParser()

    @dataclass
    class Settings:
        mgf_directory: str = field(
            default="/Volumes/T7/data/zeno_lib_v2/spectra",
            metadata={"help": "Directory containing MGF files to analyze"},
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
def extract_adduct_and_smiles_from_mgf(file_path):
    try:
        spectra = list(load_from_mgf(file_path))
        pairs = set()
        for spectrum in spectra:
            if hasattr(spectrum, "metadata"):
                smiles = spectrum.metadata.get("smiles")
                adduct = spectrum.metadata.get("adduct")
                if smiles:
                    pairs.add((adduct, smiles))
        return pairs
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return set()


@app.function
def smiles_to_inchikey_first_layer(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            inchikey = Chem.MolToInchiKey(mol)
            return inchikey[:14]
        return None
    except Exception as e:
        print(f"Error converting SMILES {smiles}: {e}")
        return None


@app.function
def read_mgf_files(settings):
    mgf_files = glob.glob(os.path.join(settings.mgf_directory, "*.mgf"))
    if not mgf_files:
        print(f"No MGF files found in {settings.mgf_directory}")

    # Precedence for TMAP grouping
    precedence = {"msmls": 3, "selleck": 2, "nexus": 1, "unknown": 0}

    def tmap_group_from_filename(filename: str) -> str:
        name = os.path.basename(filename).lower()
        if "msmls" in name:
            return "msmls"
        if "selleck" in name:
            return "selleck"
        if "nexus" in name:
            return "nexus"
        return "unknown"

    # Upset: original groups (filename stem)
    group_inchikeys = defaultdict(set)  # original_group -> {InChIKey14}
    group_adduct_smiles = defaultdict(set)  # original_group -> {(adduct, smiles)}

    # TMAP: unique SMILES + precedence group per SMILES
    unique_smiles = set()
    smiles_to_group = {}

    # Cache InChIKey conversions for efficiency
    smiles2ik = {}

    for mgf_file in mgf_files:
        original_group = os.path.splitext(os.path.basename(mgf_file))[
            0
        ]  # e.g. nexus_pos
        tmap_group = tmap_group_from_filename(mgf_file)

        for adduct, smiles in extract_adduct_and_smiles_from_mgf(mgf_file):
            # Upset: store unique (adduct, smiles) per ORIGINAL group
            group_adduct_smiles[original_group].add((adduct, smiles))

            # Upset: convert to InChIKey14 once per SMILES
            if smiles not in smiles2ik:
                smiles2ik[smiles] = smiles_to_inchikey_first_layer(smiles)
            ik14 = smiles2ik[smiles]
            if ik14:
                group_inchikeys[original_group].add(ik14)

            # TMAP: record unique SMILES and best precedence group
            unique_smiles.add(smiles)
            if (
                smiles not in smiles_to_group
                or precedence[tmap_group] > precedence[smiles_to_group[smiles]]
            ):
                smiles_to_group[smiles] = tmap_group

    # All-item universes for upset membership matrices
    all_inchikeys = sorted({ik for s, ik in smiles2ik.items() if ik})
    all_adduct_smiles = sorted(
        {(a, s) for pairs in group_adduct_smiles.values() for (a, s) in pairs},
    )

    return (
        dict(group_inchikeys),  # for upset (InChIKey14)
        all_inchikeys,
        dict(group_adduct_smiles),  # for upset (adduct, SMILES)
        all_adduct_smiles,
        smiles_to_group,  # for TMAP
        sorted(unique_smiles),  # for TMAP (unique SMILES)
    )


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
def filter_upset_data(data, group_names, group="nexus_pos"):
    pandas_df = data.to_pandas()
    set_sums = pandas_df.sum(axis=0)
    selected_cols = [col for col in set_sums.index if group in str(col)]
    if selected_cols:
        top_sets = set_sums[selected_cols].nlargest(min(6, len(selected_cols))).index
    else:
        top_sets = set_sums.nlargest(min(6, len(set_sums))).index
    pd_sub = pandas_df[top_sets]
    group_names_sub = [g for g in group_names if g in top_sets]
    return pd_sub, group_names_sub


@app.cell
def py_read_mgf_compounds():
    (
        group_inchikeys,
        all_inchikeys,
        group_adduct_smiles,
        all_adduct_smiles,
        smiles_to_group,
        all_smiles,
    ) = read_mgf_files(settings)
    return (
        all_adduct_smiles,
        all_inchikeys,
        all_smiles,
        group_adduct_smiles,
        group_inchikeys,
        smiles_to_group,
    )


@app.cell
def py_prepare_upset(
    all_adduct_smiles,
    all_inchikeys,
    group_adduct_smiles,
    group_inchikeys,
):
    data_inchikeys, names_inchikeys = create_upset_data(
        group_items=group_inchikeys,
        all_items=all_inchikeys,
        item_label="InChIKey",
    )
    data_pairs, names_pairs = create_upset_data(
        group_items=group_adduct_smiles,
        all_items=all_adduct_smiles,
        item_label="SMILES–Adduct pair",
    )
    return data_inchikeys, data_pairs, names_inchikeys, names_pairs


@app.cell
def py_filter_upset(data_inchikeys, data_pairs, names_inchikeys, names_pairs):
    pd_sub_inchikeys, names_sub_inchikey = filter_upset_data(
        data=data_inchikeys,
        group_names=names_inchikeys,  # default "nexus_pos"
    )
    pd_sub_pairs, names_sub_pairs = filter_upset_data(
        data=data_pairs,
        group_names=names_pairs,  # default "nexus_pos"
    )
    return names_sub_inchikey, names_sub_pairs, pd_sub_inchikeys, pd_sub_pairs


@app.cell
def py_plot_upset(
    names_sub_inchikey,
    names_sub_pairs,
    pd_sub_inchikeys,
    pd_sub_pairs,
):
    chart_inchikeys = au.UpSetAltair(
        data=pd_sub_inchikeys,
        sets=names_sub_inchikey,
        title="InChI Key (first layer) Intersections — Original Groups",
        height_ratio=0.9,
        width=3000,
    )
    chart_pairs = au.UpSetAltair(
        data=pd_sub_pairs,
        sets=names_sub_pairs,
        title="SMILES–Adduct Pairs Intersections — Original Groups",
        height_ratio=0.9,
        width=3000,
    )
    chart_inchikeys.chart.save("inchikeys_upset.svg", format="svg")
    chart_pairs.chart.save("pairs_upset.svg", format="svg")
    return chart_inchikeys, chart_pairs


@app.cell
def py_display_chart_inchikeys(chart_inchikeys):
    chart_inchikeys.chart
    return


@app.cell
def py_display_chart_pairs(chart_pairs):
    chart_pairs.chart
    return


@app.cell
def py_tmap_faerun(all_smiles, smiles_to_group):
    # compute descriptors on UNIQUE SMILES
    mols, hac, c_frac, ring_atom_frac, largest_ring_size = [], [], [], [], []
    for smiles in tqdm(all_smiles, desc="Processing molecules"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        atoms = mol.GetAtoms()
        size = mol.GetNumHeavyAtoms()
        n_c = sum(1 for atom in atoms if atom.GetSymbol().lower() == "c")
        n_ring_atoms = sum(1 for atom in atoms if atom.IsInRing())
        c_frac.append(n_c / size if size else 0)
        ring_atom_frac.append(n_ring_atoms / size if size else 0)
        sssr = Chem.GetSymmSSSR(mol)
        largest_ring_size.append(max((len(s) for s in sssr), default=0))
        hac.append(size)
        mols.append(mol)

    # TMAP groups (simplified)
    groups = ["nexus", "selleck", "msmls", "unknown"]
    group_to_idx = {g: i for i, g in enumerate(groups)}
    group_idx = [group_to_idx[smiles_to_group.get(s, "unknown")] for s in all_smiles]

    # MolZip graph
    zg = ZipKNNGraph()
    edge_list = zg.fit_predict(all_smiles, k=5)

    # TMAP layout
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1
    cfg.mmm_repeats = 2
    cfg.sl_repeats = 2
    x, y, s, t, _ = tm.layout_from_edge_list(
        len(all_smiles),
        edge_list,
        create_mst=True,
        config=cfg,
    )

    # carbon fraction ranks (safe normalization)
    c_arr = np.array(c_frac, dtype=float)
    denom = np.max(c_arr) if len(c_arr) and np.max(c_arr) > 0 else 1.0
    c_frac_ranked = ss.rankdata(c_arr / denom) / max(len(c_arr), 1)

    # Faerun plot
    f = Faerun(view="front", coords=False, clear_color="#ffffff")
    hc_colors = [
        high_contrast.blue,
        high_contrast.red,
        high_contrast.yellow,
        high_contrast.black,
    ]
    legend_labels_group = [(i, g) for i, g in enumerate(groups)]

    f.add_scatter(
        "attribute",
        {
            "x": list(x),
            "y": list(y),
            "c": [
                hac,
                c_frac_ranked,
                ring_atom_frac,
                largest_ring_size,
                group_idx,
            ],
            "labels": all_smiles,
        },
        shader="smoothCircle",
        point_scale=5,
        max_point_size=50,
        legend_labels=[[], [], [], [], legend_labels_group],
        categorical=[False, False, False, False, True],
        colormap=[
            cm.batlow,
            cm.batlow,
            cm.batlow,
            cm.batlow,
            mcolors.ListedColormap(hc_colors),
        ],
        series_title=[
            "HAC",
            "C Frac",
            "Ring Atom Frac",
            "Largest Ring Size",
            "Group",
        ],
        has_legend=True,
    )

    f.add_tree(
        "fragments_tree",
        {"from": list(s), "to": list(t)},
        point_helper="attribute",
        color="#e6e6e6",
    )

    f.plot("libraries_tmap", template="smiles")
    return


if __name__ == "__main__":
    app.run()
