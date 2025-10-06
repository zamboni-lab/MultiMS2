# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#   "marimo",
#   "simple_parsing",
#   "selfies",
# ]
# ///

import marimo

__generated_with = "0.16.3"
app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    from typing import List, Optional
    import os
    import re
    import selfies

    @dataclass
    class Settings:
        input_mgf: str = field(
            default="scratch/validated_spectra.mgf",
            metadata={"help": "Path to the input validated spectra MGF file."},
        )
        output_mgf: str = field(
            default="scratch/consolidated_spectra.mgf",
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

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args():  # noqa: D401
        """Parse command line arguments unless running inside a marimo notebook."""
        try:
            import marimo as mo  # type: ignore

            if mo.running_in_notebook():
                return Settings()
        except Exception:  # noqa: BLE001
            pass
        return parser.parse_args().settings

    settings = parse_args()


# ---------------- Internal Utilities ---------------- #
@app.function
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
def parse_mgf_file(path: str) -> List[SpectrumBlock]:  # keep encounter order only
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


# ---------------- Core Processing ---------------- #
@app.function
def assign_feature_ids(settings: "Settings"):
    """Assign FEATURE_ID fields based on grouping header values.

    The grouping key: (DESCRIPTION, ADDUCT, COLLISION_ENERGY, FRAGMENTATION_METHOD, INCHI_AUX).
    Empty headers do not collapse distinct blocks (unique sentinel used).
    """
    if not os.path.isfile(settings.input_mgf):
        return {"error": f"Input file not found: {settings.input_mgf}"}

    blocks = parse_mgf_file(settings.input_mgf)
    if not blocks:
        return {"error": "No spectra found in input MGF."}

    # Grouping fields for feature ID assignment (case-insensitive header matching)
    grouping_fields = [
        "description",
        "adduct",
        "collision_energy",
        "fragmentation_method",
        "inchi_aux",
    ]
    canonical_key_map = {
        # description
        "description": "description",
        # adduct
        "adduct": "adduct",
        # collision energy variants
        "collision_energy": "collision_energy",
        "collisionenergy": "collision_energy",
        "collisionenergy_ev": "collision_energy",
        "ce": "collision_energy",
        # fragmentation method variants
        "fragmentation_method": "fragmentation_method",
        "fragmentationmethod": "fragmentation_method",
        "frag_method": "fragmentation_method",
        # InChI auxiliary variants
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

    group_to_id: dict[tuple, int] = {}
    next_id = 1
    inchi_aux_index = grouping_fields.index("inchi_aux")
    unique_inchi_aux: set[str] = set()

    for block in blocks:
        gkey = extract_group_key(block)
        if gkey not in group_to_id:
            group_to_id[gkey] = next_id
            next_id += 1
        # Collect INCHI_AUX (skip sentinel unique grouping key)
        if gkey and gkey[0] != "__UNIQUE__":
            iaux_val = gkey[inchi_aux_index]
            if iaux_val:
                unique_inchi_aux.add(iaux_val)
        fid = str(group_to_id[gkey])
        if block.feature_line_idx is not None:
            block.lines[block.feature_line_idx] = f"FEATURE_ID={fid}\n"
        else:
            insert_pos = 1 if len(block.lines) > 1 else len(block.lines)
            block.lines.insert(insert_pos, f"FEATURE_ID={fid}\n")
            # Shift stored indices if needed
            if (
                block.smiles_line_idx is not None
                and block.smiles_line_idx >= insert_pos
            ):
                block.smiles_line_idx += 1
            if (
                block.selfies_line_idx is not None
                and block.selfies_line_idx >= insert_pos
            ):
                block.selfies_line_idx += 1
        # Optional SELFIES insertion
        if (
            settings.add_selfies
            and block.smiles_line_idx is not None
            and block.selfies_line_idx is None
        ):
            smi_line = block.lines[block.smiles_line_idx].strip()
            try:
                smi = smi_line.split("=", 1)[1].strip()
                if smi:
                    sf_str = selfies.encoder(smi)
                    insert_after = block.smiles_line_idx + 1
                    block.lines.insert(insert_after, f"SELFIES={sf_str}\n")
                    if (
                        block.feature_line_idx is not None
                        and block.feature_line_idx > block.smiles_line_idx
                    ):
                        block.feature_line_idx += 1
            except Exception:  # noqa: BLE001
                pass

    if not settings.dry_run:
        out_lines: List[str] = []
        for blk in blocks:
            if blk.lines and not blk.lines[-1].endswith("\n"):
                blk.lines[-1] += "\n"
            out_lines.extend(blk.lines)
            if not out_lines[-1].endswith("\n"):
                out_lines.append("\n")
            out_lines.append("\n")
        os.makedirs(os.path.dirname(settings.output_mgf), exist_ok=True)
        write_text_utf8(settings.output_mgf, out_lines)

    unique_inchi_aux_sorted = sorted(unique_inchi_aux)

    return {
        "spectra_total": len(blocks),
        "unique_feature_ids": len(group_to_id),
        "unique_inchi_aux_count": len(unique_inchi_aux_sorted),
        "output_mgf": settings.output_mgf if not settings.dry_run else None,
        "dry_run": settings.dry_run,
    }


# ------------- CLI entry ------------- #
if __name__ == "__main__":
    result = assign_feature_ids(settings)
    if isinstance(result, dict) and "error" in result:
        print("ERROR:", result["error"])  # noqa: T201
    else:
        print(result)  # noqa: T201
