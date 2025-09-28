# ModalX-SL <img align="right" width="100" height="100" src="figures/modalx-sl.svg">

TODO

There are very few things that need to be manually adjusted when not working locally:

* In `./.mzmine/batch`, change 

```
<parameter name="Database file">
            <current_file>/Users/adrutz/Documents/lab/Zeno/msmls_metadata_neg.tsv</current_file>
```

to a local path where the metadata file will have been copied (mzmine does not allow to pass a path for that object for now)

* Eventually, to avoid slow NAS reading issues, in `./commands/02_commands_library.txt`, change

```
../7600/ms2_libraries/nexus/mzml/centroided/neg/cid/20/*.mzML"
```

to where the corresponding files have been copied

TODO MGF_validator

* Also adapt in `./commands/03_commands_buddy.txt`

This should be it!

TODO
