# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "rdkit",
#     "selfies",
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
    import csv
    import os
    import math
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import selfies

    @dataclass
    class Settings:
        smiles_file: str = field(
            default="metadata/nexus_metadata_pos.tsv",
            metadata={"help": "TSV file with 'smiles' column."},
        )
        output_tsv: str = field(
            default="metadata/nexus_metadata_enhanced_pos.tsv",
            metadata={"help": "TSV output with enhanced columns."},
        )
        batch_size: int = field(
            default=100,
            metadata={"help": "SMILES to process per batch."},
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
def smiles_to_mol(smiles: str):
    """Convert SMILES to RDKit mol object."""
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


@app.function
def mol_to_formula(mol) -> str:
    """Extract molecular formula from mol object."""
    try:
        return Chem.rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return ""


@app.function
def mol_to_mass(mol) -> float:
    """Calculate exact molecular weight."""
    try:
        return Descriptors.ExactMolWt(mol)
    except Exception:
        return float("nan")


@app.function
def mol_to_canonical_smiles(mol) -> str:
    """Convert mol to canonical SMILES."""
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return ""


@app.function
def mol_to_inchi(mol) -> str:
    """Convert mol to InChI."""
    try:
        return Chem.MolToInchi(mol)
    except Exception:
        return ""


@app.function
def mol_to_inchikey(mol) -> str:
    """Convert mol to InChIKey."""
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return ""


@app.function
def smiles_to_selfies(smiles: str) -> str:
    """Convert SMILES to SELFIES."""
    try:
        return selfies.encoder(smiles)
    except Exception:
        return ""


@app.function
def enhance_smiles(settings: Settings) -> dict:
    """Enhance SMILES with formula, mass, canonical form, SELFIES, InChI, InChIKey."""
    if not os.path.exists(settings.smiles_file):
        return {"error": f"Input file not found: {settings.smiles_file}"}

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

        for smi in batch:
            mol = smiles_to_mol(smi)

            out_rows.append(
                {
                    "smiles_original": smi,
                    "molecular_formula": mol_to_formula(mol) if mol else "",
                    "mass": mol_to_mass(mol) if mol else float("nan"),
                    "smiles_canonical": mol_to_canonical_smiles(mol) if mol else "",
                    "selfies": smiles_to_selfies(smi) if mol else "",
                    "inchi": mol_to_inchi(mol) if mol else "",
                    "inchikey": mol_to_inchikey(mol) if mol else "",
                }
            )

    # Write output
    os.makedirs(os.path.dirname(settings.output_tsv), exist_ok=True)
    with open(settings.output_tsv, "w", newline="") as outfile:
        fieldnames = [
            "smiles_original",
            "molecular_formula",
            "mass",
            "smiles_canonical",
            "selfies",
            "inchi",
            "inchikey",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    return {
        "total_processed": total,
        "output_file": settings.output_tsv,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## SMILES Enhancement Settings
    
    - **Input file**: `{settings.smiles_file}`
    - **Output file**: `{settings.output_tsv}`
    - **Batch size**: {settings.batch_size}
    
    Converts SMILES to enhanced metadata including:
    - Molecular formula
    - Exact mass
    - Canonical SMILES
    - SELFIES
    - InChI
    - InChIKey
    """)
    return


@app.cell
def run_enhancement():
    result = enhance_smiles(settings)

    if "error" in result:
        mo.md(f"**Error:** {result['error']}")
    else:
        mo.md(f"""
        ### Enhancement Complete
        - **Processed**: {result['total_processed']:,} SMILES
        - **Output**: `{result['output_file']}`
        """)
    return result


if __name__ == "__main__":
    app.run()
