# Data files are not included

MIMIC-IV is released under a data use agreement that prohibits redistribution.
The analytic dataset used in this study therefore cannot be deposited here, and
no patient-level values appear anywhere in this repository.

## What is here

`cohort_all_TEMPLATE.csv` — five rows with the complete 77-column header and
every value set to `NA`. It documents the exact input structure that the
scripts expect, without disclosing any data.

## How to obtain the data

1. Complete a recognised human subjects research training course — PhysioNet
   recommends CITI "Data or Specimens Only Research".
2. Register at https://physionet.org/ and apply for credentialed access.
3. Sign the PhysioNet data use agreement for MIMIC-IV.
4. Download MIMIC-IV v3.1 from https://physionet.org/content/mimiciv/
5. Rebuild the extraction using the derived concepts (`sepsis3`, Charlson
   comorbidity, first-day SOFA and vital signs) plus `labevents`, in the layout
   described in the header of `scripts/build_dataset.py`.
6. Run `scripts/build_dataset.py` to produce `data/cohort_all.csv`.

## Why not just share the CSV

The agreement signed to obtain MIMIC-IV forbids passing the data — in whole or
in part — to anyone who has not independently completed the same training and
signed the same agreement. This applies to the derived dataset as much as to
the raw tables.
