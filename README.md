# MultiMS2: A Multi-Energy, Multi-Fragmentation Spectral Library

## Overview

MultiMS2 is a comprehensive, rigorously curated mass spectrometry spectral library designed to address critical gaps in metabolomics research. This library provides high-quality MS/MS spectra across:

- **Multiple fragmentation methods**: Collision-Induced Dissociation (CID) and Electron-Activated Dissociation (EAD)
- **Multiple collision energies**: 20, 40, and 60 eV for CID; 12, 16, and 24 eV for EAD
- **Both ionization modes**: Positive and negative
- **Diverse chemical space**: Three curated compound collections (NEXUS, Selleck, MSMLS)

### Why MultiMS2?

Confident metabolite identification relies on high-quality reference spectral libraries, yet most existing resources suffer from significant limitations:

- **Limited fragmentation diversity**: Predominantly CID spectra, with EAD spectra remaining scarce despite their structural value
- **Restricted energy ranges**: Single or limited collision energy coverage
- **Narrow acquisition conditions**: Limited applicability across varied analytical workflows
- **Machine learning constraints**: Insufficient diversity for training robust, generalizable models

MultiMS2 addresses these challenges by providing a curated resource that:

- Enables reliable metabolite annotation across diverse experimental conditions
- Supports development of generalizable machine learning models for spectrum prediction and structure elucidation
- Facilitates comparative fragmentation studies between CID and EAD
- Accelerates innovation in computational metabolomics

## Data Availability

The complete dataset is publicly available through:

- **Zenodo**: [https://doi.org/10.5281/zenodo.17250693](https://doi.org/10.5281/zenodo.17250693)
- **MassIVE**: [https://doi.org/10.25345/C5GQ6RF85](https://doi.org/10.25345/C5GQ6RF85)

## Installation & Setup

### Prerequisites

- [mzmine 3+](https://mzmine.github.io/)
- [UV package manager](https://github.com/astral-sh/uv)
- [Docker](https://www.docker.com/) (for initial conversion only)

### Quick Start

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/MultiMS2.git
cd MultiMS2
```

2. **Download spectral data from Zenodo**
```bash
uv run python notebooks/get_mzmls_from_zenodo.py
```

TODO add unzip

3. **Configure mzmine batch files**

Update the metadata file path in `.mzmine/batch/*.mzbatch`:
```xml
<parameter name="Database file">
    <current_file>/path/to/your/local/msmls_metadata_neg.tsv</current_file>
</parameter>
```

## Usage

### Spectral Extraction with mzmine

The library generation uses mzmine batch processing for consistent, reproducible spectral extraction. Below are the commands for all library combinations:

#### NEXUS Library

**Negative Mode:**

```bash
# CID at different energies
mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_cid_20/*.mzML" \
  -o "scratch/nexus_neg_cid_20"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_cid_40/*.mzML" \
  -o "scratch/nexus_neg_cid_40"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_cid_60/*.mzML" \
  -o "scratch/nexus_neg_cid_60"

# EAD at different energies
mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_ead_12/*.mzML" \
  -o "scratch/nexus_neg_ead_12"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_ead_16/*.mzML" \
  -o "scratch/nexus_neg_ead_16"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" \
  -i "scratch/nexus_mzml_centroided_neg_ead_24/*.mzML" \
  -o "scratch/nexus_neg_ead_24"
```

**Positive Mode:**

```bash
# CID at different energies
mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_cid_20/*.mzML" \
  -o "scratch/nexus_pos_cid_20"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_cid_40/*.mzML" \
  -o "scratch/nexus_pos_cid_40"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_cid_60/*.mzML" \
  -o "scratch/nexus_pos_cid_60"

# EAD at different energies
mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_ead_12/*.mzML" \
  -o "scratch/nexus_pos_ead_12"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_ead_16/*.mzML" \
  -o "scratch/nexus_pos_ead_16"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" \
  -i "scratch/nexus_mzml_centroided_pos_ead_24/*.mzML" \
  -o "scratch/nexus_pos_ead_24"
```

#### Selleck Library

**Negative Mode:**

```bash
# CID energies
mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_cid_20/*.mzML" \
  -o "scratch/selleck_neg_cid_20"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_cid_40/*.mzML" \
  -o "scratch/selleck_neg_cid_40"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_cid_60/*.mzML" \
  -o "scratch/selleck_neg_cid_60"

# EAD energies
mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_ead_12/*.mzML" \
  -o "scratch/selleck_neg_ead_12"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_ead_16/*.mzML" \
  -o "scratch/selleck_neg_ead_16"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" \
  -i "scratch/selleck_mzml_centroided_neg_ead_24/*.mzML" \
  -o "scratch/selleck_neg_ead_24"
```

**Positive Mode:**

```bash
# CID energies
mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_cid_20/*.mzML" \
  -o "scratch/selleck_pos_cid_20"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_cid_40/*.mzML" \
  -o "scratch/selleck_pos_cid_40"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_cid_60/*.mzML" \
  -o "scratch/selleck_pos_cid_60"

# EAD energies
mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_ead_12/*.mzML" \
  -o "scratch/selleck_pos_ead_12"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_ead_16/*.mzML" \
  -o "scratch/selleck_pos_ead_16"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" \
  -i "scratch/selleck_mzml_centroided_pos_ead_24/*.mzML" \
  -o "scratch/selleck_pos_ead_24"
```

#### MSMLS Library

**Negative Mode:**

```bash
# CID energies
mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" \
  -i "scratch/msmls_mzml_centroided_neg_cid_20/*.mzML" \
  -o "scratch/msmls_neg_cid_20"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" \
  -i "scratch/msmls_mzml_centroided_neg_cid_40/*.mzML" \
  -o "scratch/msmls_neg_cid_40"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" \
  -i "scratch/msmls_mzml_centroided_neg_cid_60/*.mzML" \
  -o "scratch/msmls_neg_cid_60"

# EAD energies (note: 12 eV not available)
mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" \
  -i "scratch/msmls_mzml_centroided_neg_ead_16/*.mzML" \
  -o "scratch/msmls_neg_ead_16"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" \
  -i "scratch/msmls_mzml_centroided_neg_ead_24/*.mzML" \
  -o "scratch/msmls_neg_ead_24"
```

**Positive Mode:**

```bash
# CID energies (note: 20 eV not available)
mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" \
  -i "scratch/msmls_mzml_centroided_pos_cid_40/*.mzML" \
  -o "scratch/msmls_pos_cid_40"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" \
  -i "scratch/msmls_mzml_centroided_pos_cid_60/*.mzML" \
  -o "scratch/msmls_pos_cid_60"

# EAD energies (note: 12 eV not available)
mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" \
  -i "scratch/msmls_mzml_centroided_pos_ead_16/*.mzML" \
  -o "scratch/msmls_pos_ead_16"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" \
  -i "scratch/msmls_mzml_centroided_pos_ead_24/*.mzML" \
  -o "scratch/msmls_pos_ead_24"
```

**Performance Note:** For optimal performance, consider copying files to a fast local disk before processing to avoid slow network I/O.

### Quality Control and Validation

#### Metadata

Because of an issue during `.mzML` file conversion, `COLLISION_ENERGY` and `FRAGMENTATION_METHOD` are missing from negative CID files.
To fix it, run:

```bash
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/nexus_neg_cid_20_batch_library.mgf CID 20.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/nexus_neg_cid_40_batch_library.mgf CID 40.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/nexus_neg_cid_60_batch_library.mgf CID 60.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/selleck_neg_cid_20_batch_library.mgf CID 20.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/selleck_neg_cid_40_batch_library.mgf CID 40.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/selleck_neg_cid_60_batch_library.mgf CID 60.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/msmls_neg_cid_20_batch_library.mgf CID 20.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/msmls_neg_cid_40_batch_library.mgf CID 40.0
uv run python notebooks/edit_mgf_collision_fragmentation.py /Users/adrutz/Git/MultiMS2/scratch/msmls_neg_cid_60_batch_library.mgf CID 60.0
```

#### Modalities and quality filtering

After this, spectra from all sub-libraries and modalities are concatenated by running:

```bash
uv run python notebooks/filter_spectra_consistent.py
```

Only the spectra complying to the following rules are kept:

```python
# Thresholds
min_precursor_height = 1000.0
min_precursor_purity = 0.9
min_signals = 3
min_explained_intensity = 0.4
min_explained_signals = 0.05
min_modalities = 3
min_intensity_ratio = 0.8
min_signals_ratio = 0.4
```

This led to 54,996 spectra over the 164,371 initially extracted.
(Both `all` and `filtered` MGF are exported)

#### Adduct Assignment check

TODO

Spectral quality is assessed using MGF_validator to ensure accurate adduct assignments and identify chimeric spectra (spectra containing fragments from multiple co-eluting compounds).

*Documentation in progress*

#### MS-BUDDY Molecular Formula Annotation

MS-BUDDY provides molecular formula annotation and structural validation:

```bash
for f in /Volumes/T7/data/zeno_lib_v2/spectra/*.mgf; do
    base=$(basename "$f" .mgf)
    outdir="/Volumes/T7/data/zeno_lib_v2/msbuddy/${base}"
    mkdir -p "$outdir"
    uv run msbuddy -mgf "$f" -ms qtof -p -n_cpu 12 -d -hal -o "$outdir"
done
```

**Parameters:**

- `-mgf`: Input MGF spectral file
- `-ms qtof`: Specifies QTOF instrument type
- `-p`: Positive mode
- `-n_cpu 12`: Use 12 CPU cores for parallel processing
- `-d`: Enable detailed output
- `-hal`: Consider halogen atoms in formula generation
- `-o`: Output directory

## Data Processing Pipeline (here for reference)

### Initial Data Conversion

The following steps document the complete data processing workflow from raw instrument files to curated spectral libraries.

#### Step 1: Convert `.wiff` to `.mzML`

Raw AB SCIEX `.wiff` files are converted to open-format `.mzML` using ProteoWizard:

```bash
docker run -it --rm \
  -v .:/data \
  proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses \
  wine msconvert "*.wiff" \
  --ignoreUnknownInstrumentError
```

#### Step 2: Profile to Centroided Spectra

Profile mode spectra are converted to centroided format using CentroidR:

```r
files <- "/Volumes/T7/data/7600/ms2_libraries" |>
  list.files(pattern = ".mzML", recursive = TRUE, full.names = TRUE)

files |>
  purrr::walk(
    .f = CentroidR::centroid_one_file,
    pattern = "/profile/",
    replacement = "/centroided/"
  )
```

**Note:** These preprocessing steps have already been completed for the publicly available datasets on Zenodo and MassIVE.

## Library Contents

### Compound Collections

- **NEXUS**: Diverse natural product and drug-like compounds
- **Selleck**: Bioactive compound library focused on drug discovery
- **MSMLS**: Metabolomics Standards Library compounds

### Spectral Coverage

| Collection | Ionization | Fragmentation | Energies Available |
|------------|------------|---------------|-------------------|
| NEXUS | Positive | CID | 20, 40, 60 eV |
| NEXUS | Positive | EAD | 12, 16, 24 eV |
| NEXUS | Negative | CID | 20, 40, 60 eV |
| NEXUS | Negative | EAD | 12, 16, 24 eV |
| Selleck | Positive | CID | 20, 40, 60 eV |
| Selleck | Positive | EAD | 12, 16, 24 eV |
| Selleck | Negative | CID | 20, 40, 60 eV |
| Selleck | Negative | EAD | 12, 16, 24 eV |
| MSMLS | Positive | CID | 40, 60 eV* |
| MSMLS | Positive | EAD | 16, 24 eV* |
| MSMLS | Negative | CID | 20, 40, 60 eV |
| MSMLS | Negative | EAD | 16, 24 eV* |

Some energy levels not available for certain MSMLS conditions

## Applications

MultiMS2 enables:

1. **Metabolite Annotation**: High-confidence identification through multi-energy spectral matching
2. **Machine Learning Development**: Training data for spectrum prediction and structure elucidation models
3. **Fragmentation Studies**: Comparative analysis of CID vs. EAD fragmentation patterns
4. **Method Development**: Reference spectra for optimizing MS/MS acquisition parameters
5. **Quality Assessment**: Benchmarking datasets for evaluating annotation algorithms

## Citation

If you use MultiMS2 in your research, please cite:

TODO

## License

*TODO License information to be added*

## Acknowledgments

TODO
