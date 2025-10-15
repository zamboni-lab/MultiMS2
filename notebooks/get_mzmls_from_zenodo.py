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
    import zipfile

    @dataclass
    class Settings:
        doi: str = field(
            default="10.5281/zenodo.17250693",
            metadata={"help": "Zenodo DOI for files."},
        )
        output_dir: str = field(
            default="scratch",
            metadata={"help": "Output directory for downloaded files."},
        )
        unzip: bool = field(
            default=False,
            metadata={"help": "Unzip downloaded files after downloading"},
        )
        file_filter: str | None = field(
            default=None,
            metadata={
                "help": "Optional filter pattern for files (e.g., 'selleck', 'pos')"
            },
        )

    # File list with MD5 hashes
    ZENODO_FILES = [
        ("selleck_mzml_centroided_neg_cid_20.zip", "a3533983efc09f96064b12445eda5d7f"),
        ("selleck_mzml_centroided_neg_cid_40.zip", "e2bd7e17e34ea189688f755e7ffd78a7"),
        ("selleck_mzml_centroided_neg_cid_60.zip", "843bb347bee0ac0996e15a82a5644dab"),
        ("selleck_mzml_centroided_neg_ead_12.zip", "64de8c7c0da3813e0da95d2a41c42a3d"),
        ("selleck_mzml_centroided_neg_ead_16.zip", "e87a5e427168b5c7a02d6df9f6828a5f"),
        ("selleck_mzml_centroided_neg_ead_24.zip", "33a0a6745ad4ab971690021fe04841ae"),
        ("selleck_mzml_centroided_pos_cid_20.zip", "7f5be90c1fe5e0440a6cbaa6d14c1e81"),
        ("selleck_mzml_centroided_pos_cid_40.zip", "6b04ff146759b6daf23913d367ff883b"),
        ("selleck_mzml_centroided_pos_cid_60.zip", "8ac1133a4eae63a277f46f9f21ff009f"),
        ("selleck_mzml_centroided_pos_ead_12.zip", "d0d8bd71086dd422e7b2c608f71f3df9"),
        ("selleck_mzml_centroided_pos_ead_16.zip", "bb6da1861536d53f096c71c7338b516d"),
        ("selleck_mzml_centroided_pos_ead_24.zip", "00dc8a48a371f2401c75d9676c8d8898"),
        ("nexus_mzml_centroided_neg_cid_20.zip", "2a1ab56b04c3ac78b5e1189c3e3f6301"),
        ("nexus_mzml_centroided_neg_cid_40.zip", "212ada4ccd58a3a0ec30199971dfb381"),
        ("nexus_mzml_centroided_neg_cid_60.zip", "2f4b9ee12ed3fd9436b517ebb76d844c"),
        ("nexus_mzml_centroided_neg_ead_12.zip", "f9ac5503d24df9f7b15ecf2b56e86c9f"),
        ("nexus_mzml_centroided_neg_ead_16.zip", "7c0d10dc9604a5be244d3cb7dcbce326"),
        ("nexus_mzml_centroided_neg_ead_24.zip", "1902926fb0084dea344d2ff929b458b4"),
        ("nexus_mzml_centroided_pos_cid_20.zip", "ba6bbc51acbb1b94dbb84a9ec9cfaa3e"),
        ("nexus_mzml_centroided_pos_cid_40.zip", "ee7f836e60b60ae00ff4554ba84f10c0"),
        ("nexus_mzml_centroided_pos_cid_60.zip", "28a8b99f45688a0ec607a668c8989cc0"),
        ("nexus_mzml_centroided_pos_ead_12.zip", "73a2f86d59738e161e22cf7e3a74d26d"),
        ("nexus_mzml_centroided_pos_ead_16.zip", "fe621bd229b51fa78a1f5fb1498af80b"),
        ("nexus_mzml_centroided_pos_ead_24.zip", "8a53543d1bc64afb2f1565f5208888c1"),
        ("msmls_mzml_centroided_neg_cid_20.zip", "b817193031cfa2b98ddd14fce7bf7d92"),
        ("msmls_mzml_centroided_neg_cid_40.zip", "94f3f884ca682ac1c0d09fb9cedd2768"),
        ("msmls_mzml_centroided_neg_cid_60.zip", "73ad820db2a51b49c3dc082ff72e9017"),
        ("msmls_mzml_centroided_neg_ead_16.zip", "77e229fcbb7f05f7d5c47052cf3754d2"),
        ("msmls_mzml_centroided_neg_ead_24.zip", "b28f2f86b668e59fdb83623efbf3f9f2"),
        ("msmls_mzml_centroided_pos_cid_40.zip", "9424b52347d95a3bc0eddf27e94606df"),
        ("msmls_mzml_centroided_pos_cid_60.zip", "6e3d5e8713921b42b541bf9d628ce56f"),
        ("msmls_mzml_centroided_pos_ead_16.zip", "60a2df8397a14e4035439f0a2b0bc2f4"),
        ("msmls_mzml_centroided_pos_ead_24.zip", "a03d3eb5d0389eb6d4045410e53cc8cf"),
    ]

    parser = ArgumentParser()
    parser.add_arguments(Settings, dest="settings")

    def parse_args() -> Settings:
        if mo.running_in_notebook():
            return Settings()
        else:
            return parser.parse_args().settings

    settings = parse_args()


@app.function
def download_zenodo_files(settings: Settings) -> dict:
    """Download mzML files from Zenodo with hash verification."""
    os.makedirs(settings.output_dir, exist_ok=True)

    # Filter files if requested
    files = ZENODO_FILES
    if settings.file_filter:
        files = [(f, h) for f, h in files if settings.file_filter.lower() in f.lower()]

    results = []
    errors = []

    for filename, md5 in files:
        url = f"doi:{settings.doi}/{filename}"
        file_path = os.path.join(settings.output_dir, filename)

        try:
            print(f"Downloading {filename}...")
            pooch.retrieve(
                url=url,
                known_hash=f"md5:{md5}",
                fname=filename,
                path=settings.output_dir,
                downloader=pooch.DOIDownloader(),
                progressbar=True,
            )

            msg = f"✓ {filename}"

            # Optionally unzip
            if settings.unzip and zipfile.is_zipfile(file_path):
                with zipfile.ZipFile(file_path, "r") as zip_ref:
                    zip_ref.extractall(settings.output_dir)
                msg += " (unzipped)"

            results.append(msg)
        except Exception as e:
            error_msg = f"✗ {filename}: {str(e)}"
            errors.append(error_msg)
            print(error_msg)

    return {
        "total_files": len(files),
        "successful": len(results),
        "failed": len(errors),
        "results": results,
        "errors": errors,
    }


@app.cell
def show_settings():
    file_count = len(ZENODO_FILES)
    if settings.file_filter:
        filtered = [
            f for f, _ in ZENODO_FILES if settings.file_filter.lower() in f.lower()
        ]
        file_count = len(filtered)

    mo.md(f"""
    ## Zenodo mzML Downloader Settings
    
    - **DOI**: `{settings.doi}`
    - **Output directory**: `{settings.output_dir}`
    - **Files to download**: {file_count}
    - **Unzip after download**: {settings.unzip}
    - **Filter**: {settings.file_filter or "None"}
    """)
    return


@app.cell
def run_download():
    result = download_zenodo_files(settings)

    status_lines = "\n".join(result["results"])
    error_lines = "\n".join(result["errors"]) if result["errors"] else "None"

    mo.md(f"""
    ### Download Results
    
    - **Total files**: {result['total_files']}
    - **Successful**: {result['successful']}
    - **Failed**: {result['failed']}
    
    #### Downloaded Files
    ```
    {status_lines}
    ```
    
    #### Errors
    ```
    {error_lines}
    ```
    """)
    return result


if __name__ == "__main__":
    app.run()
