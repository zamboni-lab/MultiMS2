# /// script
# requires-python = "==3.12.*"
# dependencies = [
#     "chebai>=1.1.0", # required by chebifier
#     "chebai-graph>=1.0.0", # required by chebifier
#     "chebifier==1.2.1",
#     "chemlog-extra>=1.0.1", # required by chebifier
#     "fastobo>=0.14.1", # required by chebifier
#     "huggingface-hub>=1.8.0", # required by chebifier
#     "jsonargparse[signatures]>=4.27.7", # required by chebifier, weird error with lightning
#     "marimo",
#     "matchms==0.32.0",
#     "polars==1.39.3",
#     "simple-parsing==0.1.8",
#     "torch==2.11.0", # required by chebifier
#     "torch-geometric>=2.7.0", # required by chebifier
#     "torch-scatter>=2.1.2", # required by chebifier, needs `uv add torch-scatter --no-build-isolation`
#     "tqdm==4.67.3",
# ]
# ///

import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full")

with app.setup:
    import logging
    import os
    import time
    from dataclasses import dataclass, field
    from pathlib import Path

    import marimo as mo
    import polars as pl
    from matchms.importing import load_from_mgf
    from simple_parsing import ArgumentParser
    from tqdm import tqdm

    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )

    @dataclass
    class Settings:
        mgf_path: str = field(
            default="scratch/consolidated_spectra.mgf",
            metadata={"help": "Path to the input MGF file."},
        )
        output_csv: str = field(
            default="scratch/chebifier_results.csv",
            metadata={"help": "CSV output file for classification results."},
        )
        output_parquet: str = field(
            default="scratch/chebifier_results.parquet",
            metadata={"help": "Parquet output file for classification results."},
        )
        batch_size: int = field(
            default=256,
            metadata={"help": "Number of SMILES per chebifier batch."},
        )

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        if mo.running_in_notebook():
            return Settings()
        else:
            return parser.parse_args().settings

    settings = parse_args()

    # ============================================================
    # File Parsing
    # ============================================================

    def parse_mgf_smiles(path: str) -> list[str]:
        """Parse SMILES from an MGF file using matchms."""
        smiles = []
        for spectrum in load_from_mgf(path):
            smi = spectrum.metadata.get("smiles")
            if smi:
                smiles.append(str(smi))
        return smiles

    # ============================================================
    # Chebifier Classification
    # ============================================================

    def classify_smiles_in_batches(
        smiles_list: list[str],
        batch_size: int,
    ) -> list[list[str] | None]:
        """
        Run chebifier ensemble on SMILES in batches.

        Returns a list parallel to smiles_list where each element is:
        - A list of ChEBI class name strings, or
        - None if the molecule could not be classified.
        """
        from chebifier import BaseEnsemble

        logging.info("Initialising chebifier ensemble (downloads models on first run)…")
        t0 = time.perf_counter()
        ensemble = BaseEnsemble()
        logging.info(f"Ensemble ready in {time.perf_counter() - t0:.1f}s")

        all_predictions: list[list[str] | None] = []
        n_batches = (len(smiles_list) + batch_size - 1) // batch_size

        for i in tqdm(range(n_batches), desc="Classifying batches"):
            batch = smiles_list[i * batch_size : (i + 1) * batch_size]
            preds = ensemble.predict_smiles_list(batch)
            all_predictions.extend(preds)

        return all_predictions

    # ============================================================
    # Result Building
    # ============================================================

    def build_result_df(
        unique_smiles: list[str],
        predictions: list[list[str] | None],
    ) -> pl.DataFrame:
        """Build a Polars DataFrame from classification results."""
        rows = []
        for smiles, pred in zip(unique_smiles, predictions):
            if pred is None:
                rows.append(
                    {
                        "smiles": smiles,
                        "classified": False,
                        "n_classes": 0,
                        "classes": None,
                    }
                )
            else:
                rows.append(
                    {
                        "smiles": smiles,
                        "classified": True,
                        "n_classes": len(pred),
                        "classes": "; ".join(sorted(pred)),
                    }
                )
        return pl.DataFrame(rows)

    def top_classes(df: pl.DataFrame, n: int = 20) -> pl.DataFrame:
        """Return the top-n most frequent ChEBI classes across all molecules."""
        return (
            df.filter(pl.col("classes").is_not_null())
            .with_columns(pl.col("classes").str.split("; ").alias("class_list"))
            .explode("class_list")
            .group_by("class_list")
            .len()
            .sort("len", descending=True)
            .head(n)
            .rename({"class_list": "chebi_class", "len": "molecule_count"})
        )


@app.cell
def show_settings():
    mo.md(f"""
    ## ChEBI Classification of MGF Spectra

    Classify all unique SMILES in an MGF spectral library using the
    [chebifier](https://github.com/ChEB-AI/chebifier) ensemble model.

    | Parameter | Value |
    |-----------|-------|
    | **Input MGF** | `{settings.mgf_path}` |
    | **Output CSV** | `{settings.output_csv}` |
    | **Output Parquet** | `{settings.output_parquet}` |
    | **Batch size** | `{settings.batch_size}` |
    """)
    return


@app.function
def run_pipeline(settings: Settings) -> mo.Html | pl.DataFrame:
    """End-to-end pipeline: parse → deduplicate → classify → save."""

    # ── 1. Parse ──────────────────────────────────────────────
    if not os.path.exists(settings.mgf_path):
        return mo.md(f"**Error:** MGF file not found: `{settings.mgf_path}`")

    logging.info(f"Parsing SMILES from {settings.mgf_path} …")
    all_smiles = parse_mgf_smiles(settings.mgf_path)
    logging.info(f"  {len(all_smiles):,} SMILES parsed (including duplicates)")

    # ── 2. Deduplicate ────────────────────────────────────────
    unique_smiles = list(dict.fromkeys(s for s in all_smiles if s))
    logging.info(f"  {len(unique_smiles):,} unique non-empty SMILES")

    if not unique_smiles:
        return mo.md("**Error:** No valid SMILES found in the MGF file.")

    # ── 3. Classify ───────────────────────────────────────────
    t0 = time.perf_counter()
    predictions = classify_smiles_in_batches(unique_smiles, settings.batch_size)
    elapsed = time.perf_counter() - t0
    logging.info(f"Classification finished in {elapsed:.1f}s")

    # ── 4. Build DataFrame ────────────────────────────────────
    df = build_result_df(unique_smiles, predictions)
    n_classified = int(df["classified"].sum())
    n_failed = len(df) - n_classified

    # ── 5. Save ───────────────────────────────────────────────
    Path(settings.output_csv).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(settings.output_csv)
    df.write_parquet(settings.output_parquet)
    logging.info(f"Saved results to {settings.output_csv} and {settings.output_parquet}")

    return mo.md(f"""
    ### Classification Complete

    | Metric | Value |
    |--------|-------|
    | Total SMILES (with duplicates) | {len(all_smiles):,} |
    | Unique SMILES classified | {len(unique_smiles):,} |
    | Successfully classified | {n_classified:,} |
    | Unclassified (invalid molecules) | {n_failed:,} |
    | Runtime | {elapsed:.1f}s |
    | Output CSV | `{settings.output_csv}` |
    | Output Parquet | `{settings.output_parquet}` |
    """)


@app.cell
def run_and_show():
    result = run_pipeline(settings)
    result
    return (result,)


@app.function
def show_results(settings: Settings):
    """Load saved results and display summary tables."""
    if not os.path.exists(settings.output_parquet):
        return mo.md("_Run the pipeline above to generate results._")

    df = pl.read_parquet(settings.output_parquet)
    n_total = len(df)
    n_classified = int(df["classified"].sum())

    summary_md = mo.md(f"""
    ### Results Preview

    **{n_classified:,} / {n_total:,}** molecules successfully classified
    ({100 * n_classified / n_total:.1f}%).
    """)

    top_df = top_classes(df, n=20)

    return mo.vstack(
        [
            summary_md,
            mo.md("#### Sample of classified molecules"),
            mo.ui.table(
                df.filter(pl.col("classified")).head(50).to_pandas(),
                pagination=True,
            ),
            mo.md("#### Top-20 most frequent ChEBI classes"),
            mo.ui.table(top_df.to_pandas(), pagination=False),
        ]
    )


@app.cell
def display_results():
    show_results(settings)
    return


if __name__ == "__main__":
    app.run()
