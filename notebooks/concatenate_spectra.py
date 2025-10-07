# Concatenate MGF files and export metadata
import glob
import os
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf
import polars as pl
from dataclasses import dataclass, field
from simple_parsing import ArgumentParser


@dataclass
class ConcatSettings:
    mgf_dir: str = field(
        default="scratch/", metadata={"help": "Directory containing .mgf files."}
    )
    output_tsv: str = field(
        default="scratch/spectra_metadata.tsv",
        metadata={"help": "TSV output file for spectrum metadata (not aggregated)."},
    )
    output_mgf_all: str = field(
        default="scratch/all_spectra.mgf",
        metadata={
            "help": "MGF output file for all spectra (concatenated, unfiltered)."
        },
    )


def main(settings: ConcatSettings):
    mgf_files = glob.glob(os.path.join(settings.mgf_dir, "*.mgf"))
    all_spectra = []
    rows = []
    numeric_fields = [
        ("retention_time", float),
        ("feature_ms1_height", float),
        ("precursor_purity", float),
        ("num_peaks", int),
        ("quality_explained_intensity", float),
        ("quality_explained_signals", float),
    ]
    for mgf_file in mgf_files:
        print(f"Loading {mgf_file} ...")
        filename = os.path.basename(mgf_file)
        for spectrum in load_from_mgf(mgf_file):
            meta = spectrum.metadata
            row = dict(meta)
            row["mgf_filename"] = filename
            for orig, conv in numeric_fields:
                val = meta.get(orig)
                if val is not None and val != "":
                    try:
                        row[orig] = conv(val)
                    except Exception:
                        pass
            rows.append(row)
            all_spectra.append(spectrum)
    # Export metadata
    if rows:
        all_keys = set()
        for r in rows:
            all_keys.update(r.keys())
        columns = ["mgf_filename"] + sorted(
            [k for k in all_keys if k != "mgf_filename"]
        )
        flat_rows = [
            {
                col: str(r.get(col, "")) if r.get(col, "") is not None else ""
                for col in columns
            }
            for r in rows
        ]
        df = pl.DataFrame(flat_rows)
        df.write_csv(settings.output_tsv, separator="\t")
        print(f"Wrote {len(rows)} spectra metadata rows to {settings.output_tsv}")
    else:
        print("No spectra found.")
    # Export concatenated MGF
    save_as_mgf(all_spectra, settings.output_mgf_all)
    print(f"Exported {len(all_spectra)} spectra to {settings.output_mgf_all}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_arguments(ConcatSettings, dest="settings")
    args = parser.parse_args()
    main(args.settings)
