# /// script
# requires-python = "<3.13,>=3.12"
# dependencies = [
#     "marimo",
#     "tqdm",
#     "rdkit",
#     "selfies",
#     "polars",
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
    import polars as pl
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
    from concurrent.futures import ProcessPoolExecutor

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )

    # --- MSP/MGF Parsing (Optimized) ---
    def parse_msp_smiles(path):
        smiles = []
        with open(path, "r", buffering=1024 * 1024) as f:
            content = f.read()
            blocks = content.split("\n\n")
            for block in blocks:
                for line in block.split("\n"):
                    if line.upper().startswith("SMILES:"):
                        smiles.append(line.split(":", 1)[1].strip())
                        break
        return smiles

    def parse_mgf_smiles(path):
        smiles = []
        with open(path, "r", buffering=1024 * 1024) as f:
            content = f.read()
            blocks = content.split("END IONS")
            for block in blocks:
                for line in block.split("\n"):
                    if line.upper().startswith("SMILES="):
                        smiles.append(line.split("=", 1)[1].strip())
                        break
        return smiles

    # --- Parallel descriptor computation ---
    def process_molecule_batch(smiles_batch):
        """Process a batch of SMILES in parallel"""
        results = []
        for smiles in smiles_batch:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                results.append((None, None, None, None, None, None))
                continue

            size = mol.GetNumHeavyAtoms()
            n_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol().lower() == "c")
            n_ring_atoms = sum(1 for atom in mol.GetAtoms() if atom.IsInRing())

            try:
                inchikey = Chem.MolToInchiKey(mol)
            except Exception:
                inchikey = None

            # Convert to SELFIES
            try:
                selfies_str = sf.encoder(smiles)
            except (sf.EncoderError, Exception):
                selfies_str = None

            results.append(
                (
                    size,  # hac
                    n_c / size if size else 0,  # c_frac
                    n_ring_atoms / size if size else 0,  # ring_atom_frac
                    max(
                        (len(s) for s in Chem.GetSymmSSSR(mol)), default=0
                    ),  # largest_ring
                    inchikey,
                    selfies_str,
                )
            )
        return results

    def compute_descriptors_and_conversions(smiles_list, n_workers=None):
        """Parallel computation of descriptors, mols, and SELFIES"""
        if n_workers is None:
            n_workers = min(os.cpu_count() or 4, 8)

        # Split into batches
        batch_size = max(100, len(smiles_list) // (n_workers * 4))
        batches = [
            smiles_list[i : i + batch_size]
            for i in range(0, len(smiles_list), batch_size)
        ]

        hac, c_frac, ring_atom_frac, largest_ring_size, inchikeys, selfies = (
            [],
            [],
            [],
            [],
            [],
            [],
        )

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [
                executor.submit(process_molecule_batch, batch) for batch in batches
            ]

            for future in tqdm(
                futures, desc="Processing molecules", total=len(batches), unit="batch"
            ):
                batch_results = future.result()
                for h, cf, raf, lrs, ik, sf_str in batch_results:
                    hac.append(h)
                    c_frac.append(cf)
                    ring_atom_frac.append(raf)
                    largest_ring_size.append(lrs)
                    inchikeys.append(ik)
                    selfies.append(sf_str)

        logging.info(
            f"Processed {len(smiles_list)} molecules in {len(batches)} batches"
        )
        return hac, c_frac, ring_atom_frac, largest_ring_size, inchikeys, selfies

    def convert_smiles_to_mols(smiles_list):
        """Convert SMILES to RDKit mol objects, tracking valid indices"""
        mols = []
        valid_indices = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
                valid_indices.append(i)
        logging.info(
            f"Successfully converted {len(mols)}/{len(smiles_list)} SMILES to mols"
        )
        return mols, valid_indices

    # --- Edge computation functions ---
    def compute_molzip_edges(string_list, k=10):
        """Compute edges using MolZip on any string representation"""
        t0 = time.perf_counter()
        zg = ZipKNNGraph()
        edge_list = zg.fit_predict(string_list, k=k)
        elapsed = time.perf_counter() - t0
        return edge_list, elapsed

    def compute_map4_edges(mols, k=5, n_perm=128):
        """Compute edges using MAP4 fingerprints on mol objects"""
        t0 = time.perf_counter()
        calc = MAP4(dimensions=2048, radius=2, include_duplicated_shingles=False)

        fps = calc.calculate_many(mols)

        # Create MinHashes
        minhashes = []
        for fp in fps:
            m = MinHash(num_perm=n_perm)
            nonzero_bits = np.nonzero(fp)[0]
            for bit in nonzero_bits:
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

    def compute_tmap_layout(edge_list, n_nodes):
        cfg = tm.LayoutConfiguration()
        cfg.node_size = 1.5
        cfg.mmm_repeats = 2
        cfg.sl_repeats = 2
        x, y, s, t, _ = tm.layout_from_edge_list(
            n_nodes, edge_list, create_mst=True, config=cfg
        )
        return x, y, s, t

    def visualize_tmap(x, y, s, t, descriptors, group_idx, df, filepath="tmap_out"):
        c_frac = np.array([v if v is not None else 0.0 for v in descriptors["c_frac"]])
        c_frac_ranked = ss.rankdata(
            c_frac / (max(c_frac) if np.any(c_frac) else 1.0)
        ) / len(df)
        f = Faerun(view="front", coords=False, clear_color="#ffffff")
        hc_colors = [
            high_contrast.blue,
            high_contrast.red,
            high_contrast.yellow,
            high_contrast.black,
        ]
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
            legend_labels=[legend_labels_group, [], [], [], []],
            categorical=[True, False, False, False, False],
            colormap=[mcolors.ListedColormap(hc_colors)] + [cm.batlow] * 4,
            series_title=[
                "Group",
                "HAC",
                "C Frac",
                "Ring Atom Frac",
                "Largest Ring Size",
            ],
            has_legend=True,
        )
        f.add_tree(
            "fragments_tree",
            {"from": list(s), "to": list(t)},
            point_helper="attribute",
            color="#e6e6e6",
        )

        # Split directory and filename for proper file paths
        directory = os.path.dirname(filepath) or "."
        filename = os.path.basename(filepath)
        f.plot(filename, template="smiles", path=directory)
        logging.info(f"Saved TMAP as {os.path.join(directory, filename)}.html")

    def choose_k(n_samples, method_name="Method", max_k=10):
        k = max(3, int(np.sqrt(n_samples)))
        k = min(k, max_k)
        logging.info(f"{method_name}: choosing k={k} for {n_samples} molecules")
        return k


@app.function
def run_pipeline():
    # Hardcoded file paths and their group names
    sources = [
        ("LIBGEN", "scratch/ead.msp"),
        ("LIBGEN", "scratch/hcd.msp"),
        ("LIBGEN", "scratch/uvpd.msp"),
        ("MASSBANK", "scratch/MassBank.msp_NIST"),
        ("MULTIMS2", "scratch/consolidated_spectra.mgf"),
    ]

    records = []
    for group, path in sources:
        if not os.path.exists(path):
            logging.warning(f"File not found: {path}")
            continue
        if path.endswith(".mgf"):
            smiles = parse_mgf_smiles(path)
        else:
            smiles = parse_msp_smiles(path)
        for smi in smiles:
            records.append(
                {
                    "smiles": smi,
                    "source_group": group,
                    "source_file": os.path.basename(path),
                }
            )

    if not records:
        return mo.md("No SMILES records found in the provided files.")
    df = pl.DataFrame(records)

    # --- Deduplicate SMILES early to minimize computation ---
    df_unique = df.unique(subset=["smiles"])
    smiles_unique = df_unique["smiles"].to_list()

    # Parallel computation: descriptors, mols, and SELFIES all at once
    hac, c_frac, ring_atom_frac, largest_ring_size, inchikeys, selfies = (
        compute_descriptors_and_conversions(smiles_unique)
    )

    df_unique = df_unique.with_columns(
        [
            pl.Series("hac", hac),
            pl.Series("c_frac", c_frac),
            pl.Series("ring_atom_frac", ring_atom_frac),
            pl.Series("largest_ring_size", largest_ring_size),
            pl.Series("inchikey", inchikeys),
            pl.Series("selfies", selfies),
        ]
    )

    # Remove molecules without valid inchikey and selfies
    df_unique = df_unique.filter(
      pl.col("inchikey").is_not_null() &
      pl.col("selfies").is_not_null()
      )

    # --- Merge back to full records ---
    df = df.join(
        df_unique.select(
            [
                "smiles",
                "hac",
                "c_frac",
                "ring_atom_frac",
                "largest_ring_size",
                "inchikey",
                "selfies",
            ]
        ),
        on="smiles",
        how="inner",
    )

    # --- Build group logic ---
    group_map = df.group_by("inchikey").agg(
        pl.col("source_group").unique().alias("groups")
    )
    group_map = group_map.with_columns(pl.col("groups").list.len().alias("group_count"))
    group_map = group_map.with_columns(
        pl.when(pl.col("group_count") > 1)
        .then(pl.lit("SHARED"))
        .otherwise(pl.col("groups").list.first())
        .alias("group")
    )

    df_unique_inchikeys = df.unique(subset=["inchikey"]).join(
        group_map.select(["inchikey", "group"]), on="inchikey", how="inner"
    )

    unique_groups = sorted(df_unique_inchikeys["group"].unique().to_list())
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    group_idx = [group_to_idx[g] for g in df_unique_inchikeys["group"].to_list()]

    descriptors = {
        "hac": df_unique_inchikeys["hac"].to_list(),
        "c_frac": df_unique_inchikeys["c_frac"].to_list(),
        "ring_atom_frac": df_unique_inchikeys["ring_atom_frac"].to_list(),
        "largest_ring_size": df_unique_inchikeys["largest_ring_size"].to_list(),
    }

    timings = {}

    # ============================================
    # TMAP 1: MolZip on SELFIES
    # ============================================
    selfies_list = df_unique_inchikeys["selfies"].to_list()
    # Filter out molecules that failed SELFIES conversion
    valid_selfies_indices = [i for i, s in enumerate(selfies_list) if s is not None]
    valid_selfies = [selfies_list[i] for i in valid_selfies_indices]

    if valid_selfies:
        logging.info(
            f"MolZip SELFIES: Using {len(valid_selfies)}/{len(selfies_list)} valid SELFIES"
        )

        k_selfies = choose_k(len(valid_selfies), "MolZip SELFIES")
        edge_list_selfies, t_selfies = compute_molzip_edges(valid_selfies, k=k_selfies)
        timings["MolZip_SELFIES"] = t_selfies

        # Map edges back to original indices
        edge_list_selfies_mapped = [
            (valid_selfies_indices[i], valid_selfies_indices[j], w)
            for i, j, w in edge_list_selfies
        ]

        # Filter dataframe and descriptors
        df_selfies = df_unique_inchikeys[valid_selfies_indices]
        group_idx_selfies = [group_idx[i] for i in valid_selfies_indices]
        descriptors_selfies = {
            "hac": [descriptors["hac"][i] for i in valid_selfies_indices],
            "c_frac": [descriptors["c_frac"][i] for i in valid_selfies_indices],
            "ring_atom_frac": [
                descriptors["ring_atom_frac"][i] for i in valid_selfies_indices
            ],
            "largest_ring_size": [
                descriptors["largest_ring_size"][i] for i in valid_selfies_indices
            ],
        }

        x_sf, y_sf, s_sf, t_sf = compute_tmap_layout(
            edge_list_selfies, len(valid_selfies)
        )
        visualize_tmap(
            x_sf,
            y_sf,
            s_sf,
            t_sf,
            descriptors_selfies,
            group_idx_selfies,
            df_selfies,
            "scratch/tmap_molzip_selfies",
        )
    else:
        logging.warning("No valid SELFIES found, skipping MolZip SELFIES TMAP")

    # ============================================
    # TMAP 2: MAP4
    # ============================================
    smiles_list = df_unique_inchikeys["smiles"].to_list()
    mols, valid_mol_indices = convert_smiles_to_mols(smiles_list)

    if mols:
        k_map4 = choose_k(len(mols), "MAP4")
        edge_list_map4, t_map4 = compute_map4_edges(mols, k=k_map4)
        timings["MAP4"] = t_map4

        # Map edges back to original indices
        edge_list_map4_mapped = [
            (valid_mol_indices[i], valid_mol_indices[j], w)
            for i, j, w in edge_list_map4
        ]

        # Filter dataframe and descriptors
        df_map4 = df_unique_inchikeys[valid_mol_indices]
        group_idx_map4 = [group_idx[i] for i in valid_mol_indices]
        descriptors_map4 = {
            "hac": [descriptors["hac"][i] for i in valid_mol_indices],
            "c_frac": [descriptors["c_frac"][i] for i in valid_mol_indices],
            "ring_atom_frac": [
                descriptors["ring_atom_frac"][i] for i in valid_mol_indices
            ],
            "largest_ring_size": [
                descriptors["largest_ring_size"][i] for i in valid_mol_indices
            ],
        }

        x_map4, y_map4, s_map4, t_map4 = compute_tmap_layout(edge_list_map4, len(mols))
        visualize_tmap(
            x_map4,
            y_map4,
            s_map4,
            t_map4,
            descriptors_map4,
            group_idx_map4,
            df_map4,
            "scratch/tmap_map4",
        )
    else:
        logging.warning("No valid mols found, skipping MAP4 TMAP")

    timing_md = "\n".join([f"- **{k}**: {v:.2f} s" for k, v in timings.items()])
    return mo.md(
        f"""
        ## TMAPs Computed and Saved

        - **MolZip on SELFIES**: `scratch/tmap_molzip_selfies.html`
        - **MAP4 on fingerprints**: `scratch/tmap_map4.html`

        ### Runtime Benchmark (seconds)
        {timing_md}
        """
    )


if __name__ == "__main__":
    run_pipeline()
