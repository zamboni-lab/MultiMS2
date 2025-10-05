import sys
import os
import argparse
import re


def process_mgf(input_path, frag_method, coll_energy, output_path=None):
    if output_path is None:
        output_path = input_path

    with open(input_path, "r") as fin:
        lines = fin.readlines()

    new_lines = []
    in_ions = False
    frag_found = False
    coll_found = False
    frag_idx = None
    coll_idx = None

    # We'll store indices for each BEGIN/END IONS block
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
            frag_idx = None
            coll_idx = None
            block_start = len(new_lines)
            block_end = None

        if in_ions:
            # Find FRAGMENTATION_METHOD
            if line.startswith("FRAGMENTATION_METHOD="):
                frag_found = True
                frag_idx = len(new_lines)
                frag_val = line.strip().split("=", 1)[1]
                if frag_val == "N.A." or not frag_val:
                    # Replace with ARG1
                    line = f"FRAGMENTATION_METHOD={frag_method}\n"
            # Find COLLISION_ENERGY
            if line.startswith("COLLISION_ENERGY="):
                coll_found = True
                coll_idx = len(new_lines)
            # Detect end of IONS block
            if line.strip() == "END IONS":
                block_end = len(new_lines)
                # If FRAGMENTATION_METHOD is missing, add before END IONS
                if not frag_found:
                    new_lines.insert(block_end, f"FRAGMENTATION_METHOD={frag_method}\n")
                    block_end += 1  # adjust for next insertion in this block
                # If COLLISION_ENERGY is missing, add above FRAGMENTATION_METHOD
                if not coll_found:
                    # Find the new (or old) frag line
                    # Since we may have just inserted it:
                    for j in range(block_start, block_end):
                        if new_lines[j].startswith("FRAGMENTATION_METHOD="):
                            frag_idx = j
                            break
                    new_lines.insert(frag_idx, f"COLLISION_ENERGY=[{coll_energy}]\n")
                # Reset block flags
                in_ions = False

        new_lines.append(line)
        i += 1

    with open(output_path, "w") as fout:
        fout.writelines(new_lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fix/add FRAGMENTATION_METHOD and COLLISION_ENERGY in MGF files."
    )
    parser.add_argument("mgf_path", help="Input MGF file path")
    parser.add_argument(
        "frag_method", help="Fragmentation method to use if missing or invalid"
    )
    parser.add_argument(
        "coll_energy",
        help="Collision energy (float or list for COLLISION_ENERGY, e.g. 20.0)",
    )
    parser.add_argument(
        "export_path",
        nargs="?",
        default=None,
        help="Export path (optional, default: overwrite input)",
    )

    args = parser.parse_args()
    process_mgf(args.mgf_path, args.frag_method, args.coll_energy, args.export_path)
