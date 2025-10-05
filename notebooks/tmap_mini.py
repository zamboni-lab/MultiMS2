# /// script
# requires-python = "<3.13,>=3.12"
# dependencies = [
#     "marimo",
#     "tqdm",
#     "rdkit",
#     "selfies",
#     "pandas",
#     "numpy",
#     "scipy",
#     "cmcrameri",
#     "faerun",
#     "molzip",
#     "tol-colors",
#     "datasketch",
#     "map4",
#     "matplotlib",
# ]
# ///

import marimo

__generated_with = "0.16.3"
app = marimo.App(width="full")

with app.setup:
    import logging
    import time
    import tmap as tm
    import matplotlib.colors as mcolors
    import pandas as pd
    import numpy as np
    import scipy.stats as ss
    from cmcrameri import cm
    from faerun import Faerun
    from molzip import ZipKNNGraph
    from rdkit import Chem
    from tqdm import tqdm
    from tol_colors import high_contrast
    import selfies as sf
    from datasketch import MinHash, MinHashLSHForest
    from map4 import MAP4
    import marimo as mo
    import os

    # ---------------- Logging ----------------
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )

@app.function
def compute_descriptors(smiles_list):
    mols, hac, c_frac, ring_atom_frac, largest_ring_size = [], [], [], [], []
    for smiles in tqdm(smiles_list, desc="Processing molecules"):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        atoms = mol.GetAtoms()
        size = mol.GetNumHeavyAtoms()
        n_c = sum(1 for atom in atoms if atom.GetSymbol().lower() == "c")
        n_ring_atoms = sum(1 for atom in atoms if atom.IsInRing())
        c_frac.append(n_c / size if size else 0)
        ring_atom_frac.append(n_ring_atoms / size if size else 0)
        largest_ring_size.append(
            max((len(s) for s in Chem.GetSymmSSSR(mol)), default=0)
        )
        hac.append(size)
        mols.append(mol)
    return mols, hac, c_frac, ring_atom_frac, largest_ring_size

@app.function
def compute_molzip_edges(smiles_list, k=10):
    t0 = time.perf_counter()
    zg = ZipKNNGraph()
    edge_list = zg.fit_predict(smiles_list, k=k)
    elapsed = time.perf_counter() - t0
    return edge_list, elapsed

@app.function
def compute_map4_edges(smiles_list, k=5, n_perm=256):
    t0 = time.perf_counter()
    calc = MAP4(dimensions=2048, radius=2, include_duplicated_shingles=False)
    mols = [
        Chem.MolFromSmiles(sm)
        for sm in smiles_list
        if Chem.MolFromSmiles(sm) is not None
    ]
    fps = calc.calculate_many(mols)
    minhashes = []
    for fp in fps:
        m = MinHash(num_perm=n_perm)
        for bit in np.nonzero(fp)[0]:
            m.update(str(bit).encode("utf8"))
        minhashes.append(m)
    lshf = MinHashLSHForest(num_perm=n_perm)
    for i, m in enumerate(minhashes):
        lshf.add(i, m)
    lshf.index()
    edge_list = []
    for i, m in enumerate(minhashes):
        neighbors = lshf.query(m, k)
        for n in neighbors:
            if i != n:
                edge_list.append((i, n, 1.0))
    elapsed = time.perf_counter() - t0
    return edge_list, elapsed

@app.function
def compute_selfies_edges(selfies_list, k=5, ngram_size=2, n_perm=256):
    t0 = time.perf_counter()
    minhashes = []
    for s in selfies_list:
        tokens = list(sf.split_selfies(s))
        ngrams = [
            "".join(tokens[i : i + ngram_size])
            for i in range(len(tokens) - ngram_size + 1)
        ]
        m = MinHash(num_perm=n_perm)
        for ng in ngrams:
            m.update(ng.encode("utf8"))
        minhashes.append(m)
    lshf = MinHashLSHForest(num_perm=n_perm)
    for i, m in enumerate(minhashes):
        lshf.add(i, m)
    lshf.index()
    edge_list = []
    for i, m in enumerate(minhashes):
        neighbors = lshf.query(m, k)
        for n in neighbors:
            if i != n:
                edge_list.append((i, n, 1.0))
    elapsed = time.perf_counter() - t0
    return edge_list, elapsed

@app.function
def compute_tmap_layout(edge_list, n_nodes):
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1.5
    cfg.mmm_repeats = 4
    cfg.sl_repeats = 4
    x, y, s, t, _ = tm.layout_from_edge_list(
        n_nodes, edge_list, create_mst=True, config=cfg
    )
    return x, y, s, t

@app.function
def visualize_tmap(x, y, s, t, descriptors, group_idx, df, filename="tmap_out"):
    c_frac_ranked = ss.rankdata(
        np.array(descriptors["c_frac"]) / max(descriptors["c_frac"])
    ) / len(df)
    f = Faerun(view="front", coords=False, clear_color="#ffffff")
    hc_colors = [high_contrast.blue, high_contrast.red, high_contrast.yellow]
    legend_labels_group = [(i, g) for i, g in enumerate(sorted(set(df["group"])))]
    f.add_scatter(
        "attribute",
        {
            "x": list(x),
            "y": list(y),
            "c": [
                group_idx,
                descriptors["hac"],
                c_frac_ranked,
                descriptors["ring_atom_frac"],
                descriptors["largest_ring_size"],
            ],
            "labels": df["smiles"],
        },
        shader="smoothCircle",
        point_scale=5,
        max_point_size=50,
        legend_labels=[[], [], [], [], legend_labels_group],
        categorical=[False, False, False, False, True],
        colormap=[mcolors.ListedColormap(hc_colors)] + [cm.batlow] * 4,
        series_title=["Group","HAC", "C Frac", "Ring Atom Frac", "Largest Ring Size"],
        has_legend=True,
    )
    f.add_tree(
        "fragments_tree",
        {"from": list(s), "to": list(t)},
        point_helper="attribute",
        color="#e6e6e6",
    )
    f.plot(filename, template="smiles")
    logging.info(f"Saved TMAP as {filename}")

@app.function
def choose_k(df, method_name="MolZip", max_k=15):
    n_samples = len(df)
    k = max(3, int(np.sqrt(n_samples)))
    k = min(k, max_k)
    logging.info(f"{method_name}: choosing k={k} for {n_samples} molecules")
    return k

@app.cell
def select_data_file():
    mo.md("## Select a CSV file with a 'smiles' column")
    file = mo.ui.file()
    return file

@app.function
def run_pipeline(select_data_file):
    file = None
    if select_data_file is not None and hasattr(select_data_file, "value"):
        file = select_data_file.value

    # If no file, try to load the example at "../../ion-type-analysis/data/dark_smiles.csv"
    if not file:
        example_path = "metadata/nexus_metadata_pos.tsv"
        if not os.path.exists(example_path):
            return mo.md(f"No file selected, and example file not found at `{example_path}`.")
        df = pd.read_csv(example_path, sep = "\t").drop_duplicates(subset=["smiles"]).head(5000)
    else:
        import io
        df = pd.read_csv(io.BytesIO(file["content"])).drop_duplicates(subset=["smiles"]).head(5000)

    df = df[df["smiles"].apply(lambda x: isinstance(x, str) and x.strip() != "")]

    np.random.seed(42)
    groups = ["A", "B", "C"]
    df["group"] = np.random.choice(groups, size=len(df))
    group_to_idx = {g: i for i, g in enumerate(groups)}
    group_idx = [group_to_idx[g] for g in df["group"]]

    mols, hac, c_frac, ring_atom_frac, largest_ring_size = compute_descriptors(df["smiles"])
    descriptors = {
        "hac": hac,
        "c_frac": c_frac,
        "ring_atom_frac": ring_atom_frac,
        "largest_ring_size": largest_ring_size,
    }
    timings = {}

    # MolZip
    k_mz = choose_k(df, "MolZip")
    edge_list_mz, t_mz = compute_molzip_edges(df["smiles"].tolist(), k=k_mz)
    timings["MolZip"] = t_mz
    x_mz, y_mz, s_mz, t_mz_edges = compute_tmap_layout(edge_list_mz, len(df))
    visualize_tmap(
        x_mz, y_mz, s_mz, t_mz_edges, descriptors, group_idx, df, "tmap_molzip"
    )

    # MAP4
    k_map4 = choose_k(df, "MAP4")
    edge_list_map4, t_map4 = compute_map4_edges(df["smiles"].tolist(), k=k_map4)
    timings["MAP4"] = t_map4
    x_map4, y_map4, s_map4, t_map4_edges = compute_tmap_layout(edge_list_map4, len(df))
    visualize_tmap(
        x_map4, y_map4, s_map4, t_map4_edges, descriptors, group_idx, df, "tmap_map4"
    )

    # SELFIES
    k_sf = choose_k(df, "SELFIES")
    df["selfies"] = df["smiles"].apply(sf.encoder)
    edge_list_sf, t_sf = compute_selfies_edges(df["selfies"].tolist(), k=k_sf)
    timings["SELFIES"] = t_sf
    x_sf, y_sf, s_sf, t_sf_edges = compute_tmap_layout(edge_list_sf, len(df))
    visualize_tmap(
        x_sf, y_sf, s_sf, t_sf_edges, descriptors, group_idx, df, "tmap_selfies"
    )

    timing_md = "\n".join([f"- **{k}**: {v:.2f} s" for k, v in timings.items()])
    return mo.md(
        f"""
        ## TMAPs Computed and Saved

        - Output files: `tmap_molzip.html`, `tmap_map4.html`, `tmap_selfies.html` (open in browser for interactive visualization)

        ### Runtime Benchmark (seconds)
        {timing_md}
        """
    )

if __name__ == "__main__":
    run_pipeline(None)