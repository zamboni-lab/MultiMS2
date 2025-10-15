# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
#     "pooch",
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
    import pooch
    import os

    @dataclass
    class Settings:
        url: str = field(
            default="https://github.com/MassBank/MassBank-data/releases/download/2025.05.1/MassBank.msp_NIST",
            metadata={"help": "URL to download MassBank data"},
        )
        output_path: str = field(
            default="scratch/MassBank.msp_NIST",
            metadata={"help": "Local output path"},
        )
        known_hash: str = field(
            default="11662bdd2fe0bc979ae9aa5b56e8e28cf69e6fb03a564bf15273a554b5c41ed7",
            metadata={"help": "SHA256 hash for verification"},
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
def download_massbank(settings: Settings) -> dict:
    """Download MassBank data from GitHub with hash verification."""
    os.makedirs(os.path.dirname(settings.output_path) or ".", exist_ok=True)

    file_path = pooch.retrieve(
        url=settings.url,
        path=os.path.dirname(settings.output_path) or ".",
        fname=os.path.basename(settings.output_path),
        known_hash=settings.known_hash,
    )

    file_size_mb = os.path.getsize(file_path) / (1024 * 1024)

    return {
        "file_path": file_path,
        "size_mb": file_size_mb,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## Download MassBank Data Settings
    
    - **URL**: `{settings.url}`
    - **Output path**: `{settings.output_path}`
    - **Hash verification**: SHA256
    """)
    return


@app.cell
def run_download():
    result = download_massbank(settings)
    mo.md(f"""
    ### Download Complete
    - **File path**: `{result['file_path']}`
    - **Size**: {result['size_mb']:.2f} MB
    """)
    return result


if __name__ == "__main__":
    app.run()
