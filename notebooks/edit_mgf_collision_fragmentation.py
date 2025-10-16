# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "marimo",
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
    import os

    @dataclass
    class Settings:
        mgf_path: str = field(
            default="scratch/input.mgf",
            metadata={"help": "Input MGF file path"},
        )
        frag_method: str = field(
            default="CID",
            metadata={"help": "Fragmentation method to use if missing or invalid"},
        )
        coll_energy: str = field(
            default="20.0",
            metadata={"help": "Collision energy value (e.g., 20.0)"},
        )
        export_path: str | None = field(
            default=None,
            metadata={"help": "Export path (optional, defaults to overwrite input)"},
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
def fix_mgf_fields(
    input_path: str, frag_method: str, coll_energy: str, output_path: str | None = None
) -> dict:
    """Fix/add FRAGMENTATION_METHOD and COLLISION_ENERGY in MGF files."""
    if output_path is None:
        output_path = input_path

    if not os.path.exists(input_path):
        return {"error": f"Input file not found: {input_path}"}

    with open(input_path, "r") as fin:
        lines = fin.readlines()

    new_lines = []
    in_ions = False
    frag_found = False
    coll_found = False
    block_start = None
    block_end = None

    i = 0
    while i < len(lines):
        line = lines[i]

        # Detect start of IONS block
        if line.strip() == "BEGIN IONS":
            in_ions = True
            frag_found = False
            coll_found = False
            block_start = len(new_lines)

        if in_ions:
            # Check for existing FRAGMENTATION_METHOD
            if line.startswith("FRAGMENTATION_METHOD="):
                frag_found = True
                frag_val = line.strip().split("=", 1)[1]
                if frag_val == "N.A." or not frag_val:
                    line = f"FRAGMENTATION_METHOD={frag_method}\n"

            # Check for existing COLLISION_ENERGY
            if line.startswith("COLLISION_ENERGY="):
                coll_found = True

            # Detect end of IONS block
            if line.strip() == "END IONS":
                block_end = len(new_lines)

                # Add missing FRAGMENTATION_METHOD
                if not frag_found:
                    new_lines.insert(block_end, f"FRAGMENTATION_METHOD={frag_method}\n")
                    block_end += 1

                # Add missing COLLISION_ENERGY (before FRAGMENTATION_METHOD)
                if not coll_found:
                    frag_idx = None
                    for j in range(block_start, block_end):
                        if new_lines[j].startswith("FRAGMENTATION_METHOD="):
                            frag_idx = j
                            break
                    if frag_idx is not None:
                        new_lines.insert(
                            frag_idx, f"COLLISION_ENERGY=[{coll_energy}]\n"
                        )

                in_ions = False

        new_lines.append(line)
        i += 1

    # Write output
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as fout:
        fout.writelines(new_lines)

    return {
        "input": input_path,
        "output": output_path,
        "frag_method": frag_method,
        "coll_energy": coll_energy,
    }


@app.cell
def show_settings():
    mo.md(f"""
    ## Edit MGF Collision/Fragmentation Settings

    - **Input MGF**: `{settings.mgf_path}`
    - **Fragmentation method**: `{settings.frag_method}`
    - **Collision energy**: `{settings.coll_energy}`
    - **Export path**: `{settings.export_path or settings.mgf_path + ' (overwrite)'}`
    """)
    return


@app.cell
def run_fix():
    result = fix_mgf_fields(
        settings.mgf_path,
        settings.frag_method,
        settings.coll_energy,
        settings.export_path,
    )

    if "error" in result:
        mo.md(f"**Error:** {result['error']}")
    else:
        mo.md(f"""
        ### MGF Fixed Successfully
        - **Input**: `{result['input']}`
        - **Output**: `{result['output']}`
        - **Fragmentation method**: {result['frag_method']}
        - **Collision energy**: {result['coll_energy']}
        """)
    return


if __name__ == "__main__":
    app.run()
