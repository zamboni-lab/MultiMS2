# /// script
# requires-python = "<3.13,>=3.12"
# dependencies = [
#     "marimo",
#     "simple_parsing",
#     "rdkit",
#     "selfies",
# ]
# ///

import marimo

__generated_with = "0.14.17"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    import marimo as mo
    import math
    import csv
    from typing import List
    from simple_parsing import ArgumentParser
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import selfies

    @dataclass
    class Settings:
        smiles_file: str = field(
            default="metadata/nexus_metadata_pos.tsv",
            metadata={"help": "Path to a TSV file with a 'smiles' column."},
        )
        output_tsv: str = field(
            default="metadata/nexus_metadata_enhanced_pos.tsv",
            metadata={"help": "TSV output file with enhanced columns."},
        )
        batch_size: int = field(
            default=100,
            metadata={"help": "Number of SMILES to process per batch."},
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
def smiles_to_mol(smiles):
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None

@app.function
def mol_to_formula(mol):
    try:
        return Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""

@app.function    
def mol_to_mass(mol):
    try:
        return Descriptors.ExactMolWt(mol)
    except Exception:
        return float('nan')

@app.function
def mol_to_canonical_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return ""

@app.function
def mol_to_inchi(mol):
    try:
        return Chem.MolToInchi(mol)
    except Exception:
        return ""

@app.function
def mol_to_inchikey(mol):
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""

@app.function
def smiles_to_selfies(smiles):
    try:
        return selfies.encoder(smiles)
    except Exception:
        return ""

@app.function
def smiles_to_enhanced(settings: "Settings"):
    with open(settings.smiles_file, newline='') as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        smiles_list = [row["smiles"].strip() for row in reader if row.get("smiles", "").strip()]
    total = len(smiles_list)
    out_rows = []
    num_batches = math.ceil(total / settings.batch_size)

    for batch_idx, batch_start in enumerate(range(0, total, settings.batch_size), 1):
        print(f"Processing batch {batch_idx} of {num_batches} ({batch_start + 1}-{min(batch_start + settings.batch_size, total)}/{total} SMILES)")
        batch = smiles_list[batch_start:batch_start + settings.batch_size]
        
        # Step: Process each SMILES
        for smi in batch:
            mol = smiles_to_mol(smi)
            selfies_str = smiles_to_selfies(smi) if mol else ""
            formula = mol_to_formula(mol) if mol else ""
            neutral_mass = mol_to_mass(mol) if mol else float('nan')
            smiles_canonical = mol_to_canonical_smiles(mol) if mol else ""
            inchi = mol_to_inchi(mol) if mol else ""
            inchikey = mol_to_inchikey(mol) if mol else ""
            out_rows.append({
                "smiles_original": smi,
                "molecular_formula": formula,
                "mass": neutral_mass,
                "smiles_canonical": smiles_canonical,
                "selfies": selfies_str,
                "inchi": inchi,
                "inchikey": inchikey,
            })

    # Write enhanced TSV
    with open(settings.output_tsv, "w", newline="") as outfile:
        fieldnames = [
            "smiles_original", "molecular_formula", "mass", "smiles_canonical",
            "selfies", "inchi", "inchikey",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    return f"Processed {total} SMILES. Output: {settings.output_tsv}"

@app.cell
def show_settings(settings):
    mo.md(f"""
    ## SMILES to Enhanced Structure Metadata Batch Resolver Settings

    - **SMILES file**: `{settings.smiles_file}`
    - **Output TSV**: `{settings.output_tsv}`
    - **Batch size**: `{settings.batch_size}`

    Use `smiles_to_enhanced(settings)` below to process your SMILES.

    **Note**: This converts to SELFIES, computes neutral mass, formula, and generates canonical SMILES, InChI, and InChIKey.
    """)

@app.cell
def run_smiles_to_enhanced():
    results = smiles_to_enhanced(settings)
    return results

if __name__ == "__main__":
    app.run()