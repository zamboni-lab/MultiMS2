# MultiMS2

## Introduction

TODO

### Initial conversion

This step is described here for transparency.

#### `.wiff` to `.mzML`

```shell
docker run -it --rm \
-v .:/data \
proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses \
wine msconvert "*.wiff" \
--ignoreUnknownInstrumentError
```

#### Profile to centroided spectra

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

The results have been archived on Zenodo (<https://doi.org/10.5281/zenodo.17250693>).

TODO ALSO MASSIVE LATER ON

## Use

### Download the files from Zenodo

```shell
uv run python notebooks/get_mzmls_from_zenodo.py
```

### Extraction

TODO small intro

There are very few things that need to be manually adjusted when not working locally:

* In `.mzmine/batch`, change

```
<parameter name="Database file">
            <current_file>/Users/adrutz/Documents/lab/Zeno/msmls_metadata_neg.tsv</current_file>
```

to a local path where the metadata file will have been copied (mzmine does not allow to pass a path for that object for now)

* Eventually, to avoid slow disk reading issues, copy the files to a fast disk.

```shell
# NEXUS

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/cid/20/*.mzML" -o "spectra/nexus_neg_cid_20"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/cid/40/*.mzML" -o "spectra/nexus_neg_cid_40"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/cid/60/*.mzML" -o "spectra/nexus_neg_cid_60"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/ead/12/*.mzML" -o "spectra/nexus_neg_ead_12"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/ead/16/*.mzML" -o "spectra/nexus_neg_ead_16"

mzmine -b ".mzmine/batch/nexus_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/neg/ead/24/*.mzML" -o "spectra/nexus_neg_ead_24"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/cid/20/*.mzML" -o "spectra/nexus_pos_cid_20"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/cid/40/*.mzML" -o "spectra/nexus_pos_cid_40"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/cid/60/*.mzML" -o "spectra/nexus_pos_cid_60"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/ead/12/*.mzML" -o "spectra/nexus_pos_ead_12"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/ead/16/*.mzML" -o "spectra/nexus_pos_ead_16"

mzmine -b ".mzmine/batch/nexus_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/nexus/mzml/centroided/pos/ead/24/*.mzML" -o "spectra/nexus_pos_ead_24"

# SELLECK

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/cid/20/*.mzML" -o "spectra/selleck_neg_cid_20"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/cid/40/*.mzML" -o "spectra/selleck_neg_cid_40"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/cid/60/*.mzML" -o "spectra/selleck_neg_cid_60"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/ead/12/*.mzML" -o "spectra/selleck_neg_ead_12"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/ead/16/*.mzML" -o "spectra/selleck_neg_ead_16"

mzmine -b ".mzmine/batch/selleck_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/neg/ead/24/*.mzML" -o "spectra/selleck_neg_ead_24"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/cid/20/*.mzML" -o "spectra/selleck_pos_cid_20"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/cid/40/*.mzML" -o "spectra/selleck_pos_cid_40"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/cid/60/*.mzML" -o "spectra/selleck_pos_cid_60"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/ead/12/*.mzML" -o "spectra/selleck_pos_ead_12"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/ead/16/*.mzML" -o "spectra/selleck_pos_ead_16"

mzmine -b ".mzmine/batch/selleck_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/selleck/mzml/centroided/pos/ead/24/*.mzML" -o "spectra/selleck_pos_ead_24"

# MSMLS

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/cid/20/*.mzML" -o "spectra/msmls_neg_cid_20"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/cid/40/*.mzML" -o "spectra/msmls_neg_cid_40"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/cid/60/*.mzML" -o "spectra/msmls_neg_cid_60"

# Not found
# mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/ead/12/*.mzML" -o "spectra/msmls_neg_ead_12"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/ead/16/*.mzML" -o "spectra/msmls_neg_ead_16"

mzmine -b ".mzmine/batch/msmls_library_generation_neg.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/neg/ead/24/*.mzML" -o "spectra/msmls_neg_ead_24"

# Not found
# mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/cid/20/*.mzML" -o "spectra/msmls_pos_cid_20"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/cid/40/*.mzML" -o "spectra/msmls_pos_cid_40"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/cid/60/*.mzML" -o "spectra/msmls_pos_cid_60"

# Not found
# mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/ead/12/*.mzML" -o "spectra/msmls_pos_ead_12"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/ead/16/*.mzML" -o "spectra/msmls_pos_ead_16"

mzmine -b ".mzmine/batch/msmls_library_generation_pos.mzbatch" -i "../7600/ms2_libraries/msmls/mzml/centroided/pos/ead/24/*.mzML" -o "spectra/msmls_pos_ead_24"
```

### Validation

#### Adducts and chimerism curation

TODO MGF_validator

#### BUDDY

TODO small intro

```shell
for f in /Volumes/T7/data/zeno_lib_v2/spectra/*.mgf; do
    base=$(basename "$f" .mgf)
    outdir="/Volumes/T7/data/zeno_lib_v2/msbuddy/${base}"
    mkdir -p "$outdir"
    uv run msbuddy -mgf "$f" -ms qtof -p -n_cpu 12 -d -hal -o "$outdir"
done
```

This should be it!

TODO