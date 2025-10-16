# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#   "marimo",
#   "simple_parsing",
#   "selfies",
#   "rdkit",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    from typing import List, Optional
    import os
    import re
    import selfies
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    @dataclass
    class Settings:
        input_mgf: str = field(
            default="scratch/validated_spectra.mgf",
            metadata={"help": "Path to the input validated spectra MGF file."},
        )
        output_mgf: str = field(
            default="data/multims2_spectra.mgf",
            metadata={
                "help": "Path to write the consolidated MGF with reassigned FEATURE_IDs."
            },
        )
        add_selfies: bool = field(
            default=True,
            metadata={"help": "If True, add SELFIES (below SMILES) when missing."},
        )
        dry_run: bool = field(
            default=False,
            metadata={"help": "If True, do not write output file; only report counts."},
        )
        instrument_name: Optional[str] = field(
            default=None,
            metadata={"help": "Instrument name to add to all spectra (optional)."},
        )
        data_curator: Optional[str] = field(
            default=None,
            metadata={"help": "Data curator to add to all spectra (optional)."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args():
        try:
            import marimo as mo  # type: ignore

            if mo.running_in_notebook():
                return Settings()
        except Exception:
            pass
        return parser.parse_args().settings

    settings = parse_args()


@app.class_definition
class SpectrumBlock:
    """Container for one MGF spectrum block."""

    __slots__ = (
        "file_path",
        "lines",
        "peak_count",
        "feature_line_idx",
        "smiles_line_idx",
        "selfies_line_idx",
        "order_in_file",
    )

    def __init__(
        self,
        file_path: str,
        lines: List[str],
        peak_count: int,
        order_in_file: int,
        feature_line_idx: Optional[int],
        smiles_line_idx: Optional[int],
        selfies_line_idx: Optional[int],
    ):
        self.file_path = file_path
        self.lines = lines
        self.peak_count = peak_count
        self.order_in_file = order_in_file
        self.feature_line_idx = feature_line_idx
        self.smiles_line_idx = smiles_line_idx
        self.selfies_line_idx = selfies_line_idx


@app.function
def read_text_fallback(path: str) -> List[str]:
    try:
        with open(path, encoding="utf-8") as f:
            return f.readlines()
    except UnicodeDecodeError:
        with open(path, encoding="latin-1") as f:
            return f.readlines()


@app.function
def write_text_utf8(path: str, lines: List[str]):
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)


@app.function
def parse_mgf_file(path: str) -> List[SpectrumBlock]:
    NUMERIC_LINE_RE = re.compile(r"^\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)\s+\d")
    lines = read_text_fallback(path)
    blocks: List[SpectrumBlock] = []
    cur: List[str] = []
    feature_idx = None
    smiles_idx = None
    selfies_idx = None
    peak_count = 0
    inside = False
    order = 0

    def finalize():
        nonlocal cur, feature_idx, smiles_idx, selfies_idx, peak_count, order
        if not cur:
            return
        blocks.append(
            SpectrumBlock(
                file_path=path,
                lines=cur.copy(),
                peak_count=peak_count,
                order_in_file=order,
                feature_line_idx=feature_idx,
                smiles_line_idx=smiles_idx,
                selfies_line_idx=selfies_idx,
            )
        )
        order += 1
        cur = []
        feature_idx = None
        smiles_idx = None
        selfies_idx = None
        peak_count = 0

    for line in lines:
        if line.startswith("BEGIN IONS"):
            finalize()
            inside = True
            cur = [line]
            continue
        if inside:
            cur.append(line)
            low = line.lower()
            if low.startswith("feature_id=") and feature_idx is None:
                feature_idx = len(cur) - 1
            elif low.startswith("smiles=") and smiles_idx is None:
                smiles_idx = len(cur) - 1
            elif low.startswith("selfies=") and selfies_idx is None:
                selfies_idx = len(cur) - 1
            elif NUMERIC_LINE_RE.match(line):
                peak_count += 1
            if line.startswith("END IONS"):
                inside = False
                finalize()
    finalize()
    return blocks


@app.function
def assign_feature_ids(settings: Settings):
    """
    Assign FEATURE_ID fields based on grouping header values,
    then reorder, rename, and augment fields as requested.
    """
    if not os.path.isfile(settings.input_mgf):
        return {"error": f"Input file not found: {settings.input_mgf}"}

    blocks = parse_mgf_file(settings.input_mgf)
    if not blocks:
        return {"error": "No spectra found in input MGF."}

    # --- FEATURE_ID assignment logic ---
    grouping_fields = [
        "description",
        "adduct",
        "collision_energy",
        "fragmentation_method",
        "inchi_aux",
    ]
    canonical_key_map = {
        "description": "description",
        "adduct": "adduct",
        "collision_energy": "collision_energy",
        "collisionenergy": "collision_energy",
        "collisionenergy_ev": "collision_energy",
        "ce": "collision_energy",
        "fragmentation_method": "fragmentation_method",
        "fragmentationmethod": "fragmentation_method",
        "frag_method": "fragmentation_method",
        "inchi_aux": "inchi_aux",
        "inchi": "inchi_aux",
        "inchi_key_aux": "inchi_aux",
    }

    def extract_group_key(block: SpectrumBlock):
        values = {k: "" for k in grouping_fields}
        for line in block.lines:
            if line.startswith("BEGIN IONS"):
                continue
            if line.startswith("END IONS"):
                break
            if "=" not in line:
                continue
            key, val = line.split("=", 1)
            key_l = key.strip().lower()
            key_canon = canonical_key_map.get(key_l)
            if key_canon in values:
                values[key_canon] = val.strip()
        key_tuple = tuple(values[k] for k in grouping_fields)
        if all(v == "" for v in key_tuple):
            return ("__UNIQUE__", block.order_in_file)
        return key_tuple

    group_to_id = {}
    next_id = 1

    # --- Field order and renaming ---
    field_order = [
        "FILENAME",
        "COMPOUND_NAME",
        "MOLECULEMASS",
        "INSTRUMENT_TYPE",
        "INSTRUMENT_NAME",
        "IONSOURCE",
        "EXTRACTSCAN",
        "FORMULA",
        "SMILES",
        "SELFIES",
        "INCHI",
        "INCHIAUX",
        "CHARGE",
        "IONMODE",
        "ADDUCT",
        "EXACTMASS",
        "ACQUISITION",
        "DATA_COLLECTOR",
        "DATA_CURATOR",
        "PI",
        "DESCRIPTION",
        "FEATURE_ID",
        # "FEATURE_FULL_ID",  # REMOVE IT
        # "FEATURELIST_FEATURE_ID",  # REMOVE IT
        "FEATURE_MS1_HEIGHT",
        "SPECTYPE",
        "MERGED_ACROSS_N_SAMPLES",
        "COLLISION_ENERGY",
        "FRAGMENTATION_METHOD",
        "IMS_TYPE",
        "DATASET_ID",
        # "USI",  # REMOVE IT
        "SOURCE_SCAN_USI",
        "PRECURSOR_MZ",
        "PRECURSOR_PURITY",
        "QUALITY_CHIMERIC",
        "QUALITY_EXPLAINED_INTENSITY",
        "QUALITY_EXPLAINED_SIGNALS",
        "NUM_PEAKS",
        # "SPECTRUM_ID",  # REMOVE IT
        "RTINSECONDS",
    ]
    field_rename = {
        "FILE_NAME": "FILENAME",
        "ION_SOURCE": "IONSOURCE",
        "SCANS": "EXTRACTSCAN",
        "MOLECULEMASS": "MOLECULEMASS",
        "INSTRUMENT_TYPE": "INSTRUMENT_TYPE",
        "INSTRUMENT_NAME": "INSTRUMENT_NAME",
        "SELFIES": "SELFIES",
        "INCHI_AUX": "INCHIAUX",
        "PARENT_MASS": "EXACTMASS",
        "DATA_CURATOR": "DATA_CURATOR",
        "PRINCIPAL_INVESTIGATOR": "PI",
        "RETENTION_TIME": "RTINSECONDS",
    }

    def parse_fields(lines):
        fields = {}
        peaks = []
        for line in lines:
            if line.startswith("BEGIN IONS") or line.startswith("END IONS"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                v = v.strip()
                k = field_rename.get(k, k)
                fields[k] = v
            elif re.match(r"^\s*\d+\.?\d*", line):
                peaks.append(line)
        return fields, peaks

    def build_lines(fields, peaks):
        out = ["BEGIN IONS\n"]
        for key in field_order:
            if key == "INSTRUMENT_NAME":
                val = (
                    settings.instrument_name
                    if settings.instrument_name
                    else fields.get(key, "")
                )
            elif key == "DATA_CURATOR":
                val = (
                    settings.data_curator
                    if settings.data_curator
                    else fields.get(key, "")
                )
            elif key == "MOLECULEMASS":
                smi = fields.get("SMILES", "")
                mol = Chem.MolFromSmiles(smi) if smi else None
                val = f"{Descriptors.MolWt(mol):.4f}" if mol else ""
            elif key == "SELFIES" and settings.add_selfies:
                smi = fields.get("SMILES", "")
                val = selfies.encoder(smi) if smi else ""
            elif key == "RTINSECONDS":
                val = fields.get("RTINSECONDS", "")
            else:
                val = fields.get(key, "")
            if val != "":
                out.append(f"{key}={val}\n")
        out.extend(peaks)
        if not out or not out[-1].endswith("\n"):
            out.append("\n")
        out.append("END IONS\n\n")
        return out

    # --- Assign FEATURE_ID and reorder/rename fields ---
    for block in blocks:
        # Assign FEATURE_ID
        gkey = extract_group_key(block)
        if gkey not in group_to_id:
            group_to_id[gkey] = next_id
            next_id += 1
        fid = str(group_to_id[gkey])

        fields, peaks = parse_fields(block.lines)
        fields["FEATURE_ID"] = fid

        # Remove RETENTION_TIME if present
        fields.pop("RETENTION_TIME", None)
        # Ensure RTINSECONDS is present (copy from old RETENTION_TIME if needed)
        if "RTINSECONDS" not in fields:
            if "RETENTION_TIME" in fields and fields.get("RETENTION_TIME"):
                fields["RTINSECONDS"] = fields["RETENTION_TIME"]
        # Add SELFIES if missing and requested (handled in build_lines)

        block.lines = build_lines(fields, peaks)

    if not settings.dry_run:
        out_lines = []
        for blk in blocks:
            out_lines.extend(blk.lines)
        os.makedirs(os.path.dirname(settings.output_mgf), exist_ok=True)
        write_text_utf8(settings.output_mgf, out_lines)

    return {
        "spectra_total": len(blocks),
        "unique_feature_ids": len(group_to_id),
        "output_mgf": settings.output_mgf if not settings.dry_run else None,
        "dry_run": settings.dry_run,
    }


if __name__ == "__main__":
    app.run()
