# /// script
# requires-python = "<3.13,>=3.12"
# dependencies = [
#     "marimo",
#     "pubchempy",
#     "simple_parsing",
#     "requests",
#     "rdkit",
# ]
# ///

import marimo

__generated_with = "0.14.17"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from typing import List, Dict
    from simple_parsing import ArgumentParser
    import marimo as mo
    import pubchempy as pcp
    import requests
    import time
    import csv
    import os
    import math
    from rdkit import Chem

    @dataclass
    class Settings:
        smiles_file: str = field(
            default="metadata/selleck_metadata_pos.tsv",
            metadata={"help": "Path to a TSV file with a 'smiles' column."},
        )
        output_tsv: str = field(
            default="metadata/selleck_metadata_with_names_pos.tsv",
            metadata={"help": "TSV output file: SMILES, Record_Name_or_InChIKey"},
        )
        batch_size: int = field(
            default=100,
            metadata={"help": "Number of SMILES per PubChem batch request (max 100)."},
        )
        sleep_time: float = field(
            default=0.22,
            metadata={"help": "Pause (seconds) between batch requests (PubChem limit ~5/sec)."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args():
        if mo.running_in_notebook():
            return Settings()
        else:
            args = parser.parse_args()
            return args.settings

    settings = parse_args()

@app.function
def get_batch_record_titles(cids: List[int]) -> Dict[int, str]:
    """Fetch record titles for multiple CIDs in one batch request."""
    if not cids:
        return {}
    
    cid_map = {}
    try:
        cid_string = ",".join(map(str, cids))
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_string}/description/JSON"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            if 'InformationList' in data and 'Information' in data['InformationList']:
                for info in data['InformationList']['Information']:
                    if 'CID' in info and 'Title' in info:
                        cid_map[info['CID']] = info['Title']
    except Exception as e:
        print(f"Warning: Batch title fetch failed: {e}")
    
    return cid_map

@app.function
def smiles_to_names(settings: "Settings"):
    with open(settings.smiles_file, newline='') as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        smiles_list = [row["smiles"].strip() for row in reader if row.get("smiles", "").strip()]
    total = len(smiles_list)
    out_rows = []
    num_batches = math.ceil(total / settings.batch_size)

    for batch_idx, batch_start in enumerate(range(0, total, settings.batch_size), 1):
        print(f"Processing batch {batch_idx} of {num_batches} ({batch_start + 1}-{min(batch_start + settings.batch_size, total)}/{total} SMILES)")
        batch = smiles_list[batch_start:batch_start + settings.batch_size]
        
        # Step 1: Get compounds from PubChem (by SMILES)
        try:
            compounds = pcp.get_compounds(batch, namespace="smiles")
        except Exception:
            compounds = []
            for smi in batch:
                try:
                    cmpd = pcp.get_compounds(smi, "smiles")
                    compounds.append(cmpd[0] if cmpd else None)
                    time.sleep(settings.sleep_time)
                except Exception:
                    compounds.append(None)
        
        # Step 2: Collect CIDs from compounds
        cids = [cmpd.cid for cmpd in compounds if cmpd and hasattr(cmpd, 'cid')]
        
        # Step 3: Batch fetch titles for all CIDs
        cid_to_title = {}
        if cids:
            cid_to_title = get_batch_record_titles(cids)
            time.sleep(settings.sleep_time)  # rate limit after title fetch
        
        # Step 4: Map results
        for smi, cmpd in zip(batch, compounds):
            if cmpd is None:
                # Use RDKit for InChIKey if no PubChem record found
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    try:
                        inchikey = Chem.inchi.MolToInchiKey(mol)
                        name = inchikey
                    except Exception:
                        name = "NOT_FOUND"
                else:
                    name = "NOT_FOUND"
            elif hasattr(cmpd, 'cid') and cmpd.cid in cid_to_title:
                name = cid_to_title[cmpd.cid]
            elif hasattr(cmpd, "inchikey") and cmpd.inchikey:
                name = cmpd.inchikey
            else:
                name = "NOT_FOUND"
            
            out_rows.append((smi, name))
        
        time.sleep(settings.sleep_time)  # between batches

    with open(settings.output_tsv, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(['SMILES', 'Record_Name_or_InChIKey'])
        writer.writerows(out_rows)
    return f"Processed {total} SMILES. Output: {settings.output_tsv}"

@app.cell
def show_settings(settings):
    mo.md(f"""
    ## SMILES to Record Name/InChIKey Batch Resolver Settings

    - **SMILES file**: `{settings.smiles_file}`
    - **Output TSV**: `{settings.output_tsv}`
    - **Batch size**: `{settings.batch_size}`
    - **Sleep per batch**: `{settings.sleep_time}` seconds

    Use `smiles_to_names(settings)` below to process your SMILES.

    **Note**: This now uses batch requests to fetch record titles efficiently. If no PubChem record is found, the InChIKey is generated using RDKit.
    """)

@app.cell
def run_smiles_to_names():
    results = smiles_to_names(settings)
    return results

if __name__ == "__main__":
    app.run()