# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "pubchempy",
#     "rdkit",
#     "simple_parsing",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    import marimo as mo
    import pubchempy as pcp
    import time
    import csv
    import os
    import math
    from rdkit import Chem

    @dataclass
    class Settings:
        smiles_file: str = field(
            default="metadata/selleck_metadata_pos.tsv",
            metadata={"help": "Path to TSV file with 'smiles' column."},
        )
        output_tsv: str = field(
            default="metadata/selleck_metadata_with_names_pos.tsv",
            metadata={"help": "TSV output: SMILES, Record_Name_or_InChIKey"},
        )
        batch_size: int = field(
            default=100,
            metadata={"help": "SMILES per PubChem batch (max 100)."},
        )
        sleep_time: float = field(
            default=0.22,
            metadata={"help": "Pause between requests (PubChem limit ~5/sec)."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        if mo.running_in_notebook():
            return Settings()
        else:
            return parser.parse_args().settings

    settings = parse_args()


@app.function
def get_batch_record_titles(cids: list[int]) -> dict[int, str]:
    """Fetch record titles for multiple CIDs in one batch."""
    if not cids:
        return {}

    cid_map = {}
    try:
        import requests

        cid_string = ",".join(map(str, cids))
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_string}/description/JSON"
        response = requests.get(url, timeout=30)

        if response.status_code == 200:
            data = response.json()
            if "InformationList" in data and "Information" in data["InformationList"]:
                for info in data["InformationList"]["Information"]:
                    if "CID" in info and "Title" in info:
                        cid_map[info["CID"]] = info["Title"]
    except Exception as e:
        print(f"Warning: Batch title fetch failed: {e}")

    return cid_map


@app.function
def smiles_to_names(settings: Settings) -> dict:
    """Batch resolve SMILES to PubChem names or InChIKeys."""
    with open(settings.smiles_file, newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        smiles_list = [
            row["smiles"].strip() for row in reader if row.get("smiles", "").strip()
        ]

    total = len(smiles_list)
    out_rows = []
    num_batches = math.ceil(total / settings.batch_size)

    for batch_idx, batch_start in enumerate(range(0, total, settings.batch_size), 1):
        print(
            f"Processing batch {batch_idx}/{num_batches} ({batch_start + 1}-{min(batch_start + settings.batch_size, total)}/{total})"
        )
        batch = smiles_list[batch_start : batch_start + settings.batch_size]

        # Get compounds from PubChem
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

        # Collect CIDs and batch fetch titles
        cids = [c.cid for c in compounds if c and hasattr(c, "cid")]
        cid_to_title = {}
        if cids:
            cid_to_title = get_batch_record_titles(cids)
            time.sleep(settings.sleep_time)

        # Map results
        for smi, cmpd in zip(batch, compounds):
            if cmpd is None:
                # Fallback to InChIKey via RDKit
                mol = Chem.MolFromSmiles(smi)
                name = Chem.MolToInchiKey(mol) if mol else "NOT_FOUND"
            elif hasattr(cmpd, "cid") and cmpd.cid in cid_to_title:
                name = cid_to_title[cmpd.cid]
            elif hasattr(cmpd, "inchikey") and cmpd.inchikey:
                name = cmpd.inchikey
            else:
                name = "NOT_FOUND"

            out_rows.append((smi, name))

        time.sleep(settings.sleep_time)

    # Write output
    os.makedirs(os.path.dirname(settings.output_tsv), exist_ok=True)
    with open(settings.output_tsv, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["SMILES", "Record_Name_or_InChIKey"])
        writer.writerows(out_rows)

    return {
        "total_processed": total,
        "output_file": settings.output_tsv,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## SMILES to Names Batch Resolver
    
    - **Input file**: `{settings.smiles_file}`
    - **Output file**: `{settings.output_tsv}`
    - **Batch size**: {settings.batch_size}
    - **Sleep per batch**: {settings.sleep_time}s
    
    Resolves SMILES to PubChem record names or InChIKeys (fallback).
    """)
    return


@app.cell
def run_resolution():
    result = smiles_to_names(settings)
    mo.md(f"""
    ### Resolution Complete
    - **Processed**: {result['total_processed']:,} SMILES
    - **Output**: `{result['output_file']}`
    """)
    return result


if __name__ == "__main__":
    app.run()
