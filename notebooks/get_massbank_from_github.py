# /// script
# requires-python = ">=3.13,<4"
# dependencies = [
#     "pooch",
# ]
# ///

import pooch

file_path = pooch.retrieve(
    url="https://github.com/MassBank/MassBank-data/releases/download/2025.05.1/MassBank.msp_NIST",
    path=".",
    fname="scratch/MassBank.msp_NIST",
    known_hash="11662bdd2fe0bc979ae9aa5b56e8e28cf69e6fb03a564bf15273a554b5c41ed7",
)

print(f"Downloaded to: {file_path}")
