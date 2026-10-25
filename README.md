# Prognostic indices in COPD with sepsis — analysis code

Analysis code for the manuscript:

> **Nutritional versus inflammatory indices in COPD with sepsis: the prognostic signal resides in lymphocytes rather than albumin. A head-to-head comparison of ten routine indices in the MIMIC-IV database**

Submitted to *BMC Pulmonary Medicine*.

---

## What this repository contains

The code that turns the MIMIC-IV derived extraction into every number, table and figure in the manuscript. **No patient-level data are included** — MIMIC-IV is distributed under a data use agreement that forbids redistribution. See [Data access](#data-access).

```
scripts/
  build_dataset.py     Build the analytic cohort from the MIMIC-IV extraction
  tables_all.R         Table 1 (baseline), Table 2 (AUC), Table 3 (Cox),
                       Table 4 (competing models), Table 5 (time windows)
  make_figs_R.R        Figure 1 (flow), Figure 2 (discrimination),
                       Figure 3 (forest), Figure 4 (coefficient plot),
                       Figure S1 (Kaplan-Meier by NLR tertile)
  lrt_test.R           Likelihood ratio tests on nested models
  diagnostics.R        Proportional hazards assumption (cox.zph), EPV
  vif.R                Variance inflation factors
  verify.py            Cross-check headline numbers against the source CSVs
data/
  cohort_all.csv       Column structure only (header + 5 masked rows) — see below
output/                Destination for generated tables and figures
```

---

## Data access

MIMIC-IV v3.1 is available to credentialed users from PhysioNet:

- Repository: https://physionet.org/content/mimiciv/
- Access requires completion of a recognised human subjects research training
  course (e.g. CITI "Data or Specimens Only Research") and signing of the
  PhysioNet data use agreement.

The extraction used here is built from the MIMIC-IV derived concepts
(`sepsis3`, Charlson comorbidity, first-day SOFA and vital-sign summaries)
plus `labevents` for the first albumin and lymphocyte count after ICU
admission. `scripts/build_dataset.py` expects the extracted CSVs in the
layout described in its header comment.

**The analytic dataset cannot be shared** under the terms of the data use
agreement. `data/cohort_all.csv` in this repository keeps the column
structure so that the expected input is clear, with all values masked.

---

## Reproducing the analysis

Tested with Python 3.13 and R 4.5.2 on Windows 11.

### Python

```bash
pip install pandas numpy pyyaml
python scripts/build_dataset.py        # writes cohort_all.csv to data/
python scripts/verify.py               # prints headline numbers for checking
```

### R

```r
install.packages(c("survival", "rms", "ggplot2", "survminer",
                   "dplyr", "grid", "patchwork"))

source("scripts/tables_all.R")   # Tables 1-5 -> output/
source("scripts/lrt_test.R")     # nested-model likelihood ratio tests
source("scripts/diagnostics.R")  # cox.zph, events per variable
source("scripts/vif.R")          # variance inflation factors
source("scripts/make_figs_R.R")  # Figures 1-4 and S1 -> output/figures/
```

Run order matters: `tables_all.R` must complete before `make_figs_R.R`,
which reads Tables 2–4 for the forest and coefficient plots.

---

## Cohort definition, in one paragraph

Adults (≥18 years) meeting Sepsis-3 criteria at their first ICU admission,
with an ICU stay longer than 24 hours, and COPD recorded among hospital
discharge diagnoses. COPD is identified with the Charlson comorbidity
algorithm as implemented in the MIMIC-IV derived tables, which maps ICD-9-CM
490–505, 4168, 4169, 5064, 5081, 5088 and ICD-10 J40–J47, J60–J67, I278,
I279, J684, J701, J703. Because serum albumin is measured on clinical
indication rather than routinely, patients were also required to have both a
first albumin value and a first absolute lymphocyte count obtained after ICU
admission and within 72 hours. The primary outcome is all-cause mortality
within 90 days of ICU admission, taken from the registry `dod` field.

Resulting cohort: **531 patients, 223 deaths (42.0%)**.

---

## Ten indices

| Index | Definition |
|---|---|
| PNI | prognostic nutritional index = 10 × albumin (g/dL) + 5 × lymphocytes (10⁹/L) |
| RAR | RDW-to-albumin ratio |
| Albumin | serum albumin |
| RDW | red cell distribution width |
| Lymphocyte | absolute lymphocyte count |
| NLR | neutrophil-to-lymphocyte ratio |
| SII | systemic immune-inflammation index = platelets × neutrophils / lymphocytes |
| SIRI | systemic inflammation response index = neutrophils × monocytes / lymphocytes |
| PLR | platelet-to-lymphocyte ratio |
| ACAG | albumin-corrected anion gap = anion gap + 2.5 × (4.0 − albumin) |

Each index is expressed per 1 standard deviation in the Cox models, with the
sign oriented so that higher values indicate higher risk.

---

## Statistical approach

- Discrimination: AUC with 95% CI from 2,000 bootstrap resamples.
- Association: Cox proportional hazards, Model 1 unadjusted, Model 2 adjusted
  for age, sex, race, comorbidity burden, SAPS II and renal replacement
  therapy.
- Incremental value: Harrell's C-index, plus models containing both a
  nutritional and an inflammatory term, with likelihood ratio tests on the
  nested versions.
- Assumptions: proportional hazards by scaled Schoenfeld residuals; variance
  inflation factors; events per variable.
- Sensitivity: exposure windows of 72, 48 and 24 hours.
- Shape of association: restricted cubic splines, four knots, median as
  reference.

---

## Licence

Code released under the MIT Licence. The MIMIC-IV data themselves remain
subject to the PhysioNet data use agreement and are not covered by this
licence.

## Contact

Weidong Li, MD — corresponding author — lwd10000@126.com
