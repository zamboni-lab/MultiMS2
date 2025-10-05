# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "rdkit",
#     "simple_parsing",
# ]
# ///

import marimo

app = marimo.App(width="full")

with app.setup:
    from dataclasses import dataclass, field
    from simple_parsing import ArgumentParser
    from typing import Dict, List, Tuple
    import re
    import sys
    from rdkit import Chem
    from rdkit.Chem import AllChem

    @dataclass
    class Settings:
        input: str = field(default="scratch/filtered_spectra.mgf", metadata={"help": "Path to input MGF file"})
        output: str | None = field(default="scratch/validated_spectra.mgf", metadata={"help": "Path to write validated spectra MGF"})
        invalid_out: str | None = field(default="scratch/invalidated_spectra.mgf", metadata={"help": "Path to write invalidated spectra MGF"})
        stats: bool = field(default=True, metadata={"help": "Print summary statistics when run as script"})

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        import marimo as mo
        if mo.running_in_notebook():
            return Settings()
        else:
            args = parser.parse_args()
            return args.settings

    settings = parse_args()

class MGFLossValidator:
    """Loss validation using reaction-based neutral loss rules.

    This version preserves original spectra (no structural or metadata mutation)
    and only annotates validity via a VALIDATION_NOTE downstream.
    """
    def __init__(self):
        """Initialize validator with reaction-based neutral losses."""
        # Core reactions to mimic mass-spec losses
        self.loss_reactions = {
            # 15Da: N-demethylation, O-demethylation, and S-demethylation
            "CH3_from_NOS": AllChem.ReactionFromSmarts(
                "[#7,#8,#16:1][CH3:2]>>[#7,#8,#16:1].[CH3:2]"
            ),
            # 17Da: Ammonia loss from primary amine
            "NH3": AllChem.ReactionFromSmarts(
                "[#6,#16:1][#6,#16:2][NH2:3]>>[#6,#16:1]=[#6,#16:2].[NH3:3]"
            ),
            # 18Da: Dehydration
            "H2O": AllChem.ReactionFromSmarts(
                "[#6,#7,#16:1][#6,#7,#16:2][OX2H:3]>>[#6,#7,#16:1]=[#6,#7,#16:2].[OH2:3]"
            ),
            # 28Da: CO loss from carbonyl
            "CO": AllChem.ReactionFromSmarts("[C:1][C:2]=O>>[C:1][C:2].[C]=O"),
            # 28Da: Ethylene loss
            "C2H4": AllChem.ReactionFromSmarts(
                "[C:1][C:2][C:3][C:4]>>[C:1][C:4].[C:2]=[C:3]"
            ),
            # 32Da: Methanol loss
            "CH3OH": AllChem.ReactionFromSmarts("[#6:1][O][#6:2]>>[#6:1].[O][#6:2]"),
            # 44Da: CO2 loss from carboxylic acid
            "CO2": AllChem.ReactionFromSmarts(
                "[#6:1][C](=O)[OH]>>[#6:1][H].[O]=[C]=[O]"
            ),
            # 162Da: Dehydrohexose loss
            "C6H10O5": AllChem.ReactionFromSmarts(
                "[O:1][C:2]1O[C][C][C][C][C]O1>>[O:1].[C:2]1O[C][C][C][C][C]O1"
            ),
        }

    # ---------- MGF parsing (unchanged) ---------- #
    def parse_mgf_spectrum(self, spectrum_text: str) -> Dict:
        lines = spectrum_text.strip().split("\n")
        spectrum: Dict = {}
        peaks = []
        for line in lines:
            line = line.strip()
            if line.startswith("BEGIN IONS") or line.startswith("END IONS"):
                continue
            if "=" in line and not line[:1].isdigit():
                key, value = line.split("=", 1)
                spectrum[key] = value
            elif line and line[0].isdigit():
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        peaks.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
        spectrum["peaks"] = peaks
        return spectrum

    def parse_mgf_file(self, mgf_content: str) -> List[Dict]:
        """Parse MGF content preserving original block text.
        Each returned spectrum dict includes a RAW_TEXT key with the exact
        original block (BEGIN/END IONS inclusive) to allow lossless filtering.
        """
        # Regex finds each full block
        pattern = re.compile(r"BEGIN IONS\n(.*?)END IONS", re.DOTALL)
        spectra: List[Dict] = []
        for match in pattern.finditer(mgf_content):
            block_body = match.group(1)
            full_block = f"BEGIN IONS\n{block_body}END IONS"  # original formatting retained
            spectrum = self.parse_mgf_spectrum(full_block)
            spectrum["RAW_TEXT"] = full_block
            spectra.append(spectrum)
        return spectra

    # ---------- Adduct parsing ---------- #
    def parse_adduct(self, adduct_str: str) -> Dict:
        clean = adduct_str.strip().strip("[]")
        # Remove optional leading 'M' (common notation) without touching later characters
        if clean.startswith('M') and (len(clean) == 1 or clean[1] in ['-','+']):
            clean = clean[1:]
        # Extract trailing charge symbols but we don't use charge for feasibility
        while clean.endswith('+') or clean.endswith('-'):
            clean = clean[:-1]
        losses: list[str] = []
        if self.loss_reactions:
            # Allow end-of-string as a terminator for a loss token
            loss_pattern = (
                r"-(\d*)(" + "|".join(re.escape(k) for k in self.loss_reactions) + r")(?=$|[+\-\]])"
            )
            for m in re.finditer(loss_pattern, clean):
                count = int(m.group(1)) if m.group(1) else 1
                loss_name = m.group(2)
                losses.extend([loss_name] * count)
        return {"losses": losses}

    # ---------- Core validation ---------- #
    def can_form_adduct(self, smiles: str, adduct: str) -> Tuple[bool, str]:
        """Validate feasibility of requested neutral losses.

        Strategy:
          1. Parse requested losses.
          2. Attempt reaction enumeration (greedy, shallow) for each unique loss.
          3. If enumeration underestimates capacity, fall back to functional-group
             heuristics specific to each loss type.
          4. Accept if all requested losses are supported by either enumeration
             or heuristics; reject otherwise.

        Heuristics deliberately over-estimate within reasonable chemical bounds
        to avoid false negatives for polyfunctional natural products / sugars.
        No structural mutation is performed.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES"

        info = self.parse_adduct(adduct)
        loss_list = info["losses"]
        if not loss_list:
            return True, "Adduct validated"

        # Precompute substructure matches once (lazy caches)
        cache: dict[str, int] = {}
        def count(pattern: str, key: str) -> int:
            if key in cache:
                return cache[key]
            patt = Chem.MolFromSmarts(pattern)
            if patt is None:
                cache[key] = 0
            else:
                cache[key] = len(mol.GetSubstructMatches(patt, uniquify=True))
            return cache[key]

        # Heuristic capacity calculators per loss (upper bounds)
        def capacity_H2O() -> int:
            # Hydroxyls (alcohol / phenol / carboxylic OH) as primary drivers
            oh = count('[OX2H]', 'oh')
            # Allow limited extra dehydrations for dense polyols (cyclic sugars) : every 3 OH could yield 1 additional internal dehydration
            bonus = oh // 3
            return oh + bonus
        def capacity_C2H4() -> int:
            # Broadened heuristic for ethylene eliminations:
            # 1. Adjacent sp3-sp3 single bonds (any, ring or not)
            sp3_pair = count('[CX4][CX4]', 'sp3_pair')
            # 2. Non-ring methylene pairs (legacy conservative signal of flexible chain)
            chain_pair = count('[CH2;!R][CH2;!R]', 'c2h4_chain_pair')
            # 3. Total sp3 carbon pool (each elimination needs roughly two sp3 sites)
            sp3_any = count('[CX4]', 'sp3_any')
            pool_based = sp3_any // 2
            # Combine – take the maximum plausible count; ensure at least 1 if there are >=2 sp3 carbons
            raw_cap = max(sp3_pair, chain_pair, pool_based)
            if raw_cap == 0 and sp3_any >= 2:
                raw_cap = 1
            return raw_cap
        def capacity_CH3OH() -> int:
            # Methoxy present: O-CH3 (ether or ester) groups
            return count('[OX2;!$(*=O)][CH3]', 'methoxy') + count('[C](=O)[OX2][CH3]', 'ester_methoxy')
        def capacity_CH3_from_NOS() -> int:
            return count('[#7,#8,#16]-[CH3]', 'nos_methyl')
        def capacity_NH3() -> int:
            # Primary amines: [NH2][CX4] or terminal ammonium-like sites
            return count('[NH2][CX4]', 'primary_amine')
        def capacity_CO() -> int:
            # Carbonyl groups not in amide (approx) that can decarbonylate
            return count('[CX3](=O)[#6,#7,#8,#16]', 'carbonyl_any')
        def capacity_CO2() -> int:
            # Carboxylic acids & activated carboxylates / esters
            acid = count('C(=O)[OX2H1]', 'carboxylic_acid')
            ester = count('C(=O)[OX2][CX4]', 'ester')
            # Imide / cyclic imide motif: two adjacent carbonyls bound to the same nitrogen
            imide = count('[CX3](=O)N[CX3](=O)', 'imide')
            # Generic amide (carbonyl next to nitrogen). Each match is one carbonyl; subtract imide carbonyls (2 per imide) to avoid double counting
            amide_any = count('[CX3](=O)N', 'amide_any')
            extra_amide = max(0, amide_any - imide * 2)
            # Treat each qualifying carbonyl environment (acid/ester/imide/extra amide) as a potential CO2 source heuristically
            return acid + ester + imide + extra_amide
        def capacity_C6H10O5() -> int:
            # Sugar-like ring (pyranose): O1[C][C][C][C][C]O1 coarse pattern
            return count('O1[C][C][C][C][C]O1', 'pyranose')

        heuristic_funcs = {
            'H2O': capacity_H2O,
            'C2H4': capacity_C2H4,
            'CH3OH': capacity_CH3OH,
            'CH3_from_NOS': capacity_CH3_from_NOS,
            'NH3': capacity_NH3,
            'CO': capacity_CO,
            'CO2': capacity_CO2,
            'C6H10O5': capacity_C6H10O5,
        }

        unapplied_losses: list[str] = []
        special_notes: list[str] = []

        for loss in sorted(set(loss_list)):
            requested = loss_list.count(loss)
            rxn = self.loss_reactions.get(loss)
            enumerated = 0
            if rxn is not None:
                frontier = [(mol,)]
                depth_guard = 0
                while frontier and enumerated < requested and depth_guard < requested * 4:
                    new_frontier = []
                    for tpl in frontier:
                        try:
                            result = rxn.RunReactants(tpl)
                        except Exception:
                            result = []
                        if result:
                            enumerated += 1
                            # Keep limited branching (first 2) to control explosion
                            new_frontier.extend(result[:2])
                            if enumerated >= requested:
                                break
                    frontier = new_frontier
                    depth_guard += 1
            # Heuristic fallback
            if enumerated < requested:
                cap_fn = heuristic_funcs.get(loss)
                heuristic_cap = cap_fn() if cap_fn else 0
                if requested <= heuristic_cap and heuristic_cap > 0:
                    special_notes.append(f"{loss}: heuristic capacity {heuristic_cap} used (enumerated {enumerated})")
                else:
                    unapplied_losses.append(f"{loss} (need {requested}, have {max(enumerated, heuristic_cap)})")
            elif enumerated > requested:
                special_notes.append(f"{loss}: more reactive sites ({enumerated}) than requested ({requested})")

        if unapplied_losses:
            return False, "Cannot apply losses: " + ", ".join(unapplied_losses)
        note = "Adduct validated"
        if special_notes:
            note += "; " + "; ".join(special_notes)
        return True, note

    def validate_spectrum(self, spectrum: Dict) -> Tuple[bool, str]:
        if "SMILES" not in spectrum:
            return False, "No SMILES found"
        if "ADDUCT" not in spectrum:
            return False, "No ADDUCT found"
        return self.can_form_adduct(spectrum["SMILES"], spectrum["ADDUCT"])

    def process_mgf_content(self, mgf_content: str):
        spectra = self.parse_mgf_file(mgf_content)
        valid, invalid = [], []
        stats = {"total": len(spectra), "valid": 0, "invalid": 0, "modified": 0, "reasons": {}}
        for spec in spectra:
            ok, reason = self.validate_spectrum(spec)
            # annotate in stats only; do NOT inject new metadata into RAW_TEXT
            if ok:
                valid.append(spec)  # keep original dict (with RAW_TEXT)
                stats["valid"] += 1
            else:
                invalid.append(spec)
                stats["invalid"] += 1
            if reason.startswith("Cannot apply losses"):
                reason_key = "Cannot apply losses"
            elif reason.startswith("Invalid SMILES"):
                reason_key = "Invalid SMILES"
            elif reason.startswith("No SMILES"):
                reason_key = "No SMILES"
            elif reason.startswith("No ADDUCT"):
                reason_key = "No ADDUCT"
            elif "more reactive sites" in reason:
                reason_key = "More sites than requested"
            else:
                reason_key = "Adduct validated"
            stats["reasons"][reason_key] = stats["reasons"].get(reason_key, 0) + 1
        return valid, invalid, stats

    def write_mgf_spectrum(self, spectrum: Dict) -> str:
        lines = ["BEGIN IONS"]
        for k, v in spectrum.items():
            if k == "peaks":
                continue
            lines.append(f"{k}={v}")
        for mz, inten in spectrum.get("peaks", []):
            lines.append(f"{mz} {inten}")
        lines.append("END IONS")
        return "\n".join(lines)

    def write_mgf_file(self, spectra: List[Dict]) -> str:
        """Write spectra using original untouched RAW_TEXT blocks if present."""
        blocks = []
        for spec in spectra:
            raw = spec.get("RAW_TEXT")
            if raw:
                blocks.append(raw)
            else:  # fallback (should not happen if parsed by this class)
                blocks.append(self.write_mgf_spectrum(spec))
        return "\n\n".join(blocks)

# ------------- Marimo exposed functions ------------- #
@app.function
def validate_file(settings: "Settings"):
    if not settings.input:
        return {"error": "No input file specified"}
    try:
        with open(settings.input, "r") as fh:
            content = fh.read()
    except FileNotFoundError:
        return {"error": f"Input file not found: {settings.input}"}
    validator = MGFLossValidator()
    valid, invalid, stats = validator.process_mgf_content(content)
    return {"stats": stats, "n_valid": len(valid), "n_invalid": len(invalid)}

# ------------- CLI entry (preserved) ------------- #
def _cli_main():
    validator = MGFLossValidator()
    try:
        with open(settings.input, "r") as f:
            mgf_content = f.read()
    except Exception as e:
        print(f"Error reading input: {e}", file=sys.stderr)
        sys.exit(1)
    valid, invalid, stats = validator.process_mgf_content(mgf_content)
    if settings.output:
        try:
            with open(settings.output, "w") as f:
                f.write(validator.write_mgf_file(valid))
            print(f"Valid spectra written to {settings.output}")
        except Exception as e:
            print(f"Error writing valid output: {e}", file=sys.stderr)
    if settings.invalid_out:
        try:
            with open(settings.invalid_out, "w") as f:
                f.write(validator.write_mgf_file(invalid))
            print(f"Invalid spectra written to {settings.invalid_out}")
        except Exception as e:
            print(f"Error writing invalid output: {e}", file=sys.stderr)
    if settings.stats:
        print("\nProcessing Results:")
        print(f"Total spectra: {stats['total']}")
        print(f"Valid spectra: {stats['valid']}")
        print(f"Invalid spectra: {stats['invalid']}")
        print("\nReasons:")
        for r, c in sorted(stats["reasons"].items()):
            print(f"  {r}: {c}")

if __name__ == "__main__":
    # Only run CLI main if launched as script
    _cli_main()
