# NHANES HbA1c, urinary incontinence, and depressive symptoms

Reproducible code and frozen derived results for **Association of Glycated Hemoglobin With Depressive Symptoms Among U.S. Women With Urinary Incontinence: A Complex-Survey Analysis of NHANES 2021–2023**.

`Analysis_Code_NHANES_HbA1c.R` is the analysis script; `results/` contains frozen derived tables, sample flow, diagnostics, provenance, and session information. Raw NHANES files are not redistributed. Download the August 2021–August 2023 public-use files from the [NCHS NHANES page](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2021-2023).

The analysis uses `SDMVPSU`, `SDMVSTRA`, and the HbA1c-specific phlebotomy weight `WTPH2YR`. A singleton-PSU sensitivity analysis using `adjust`, `certainty`, and `remove` is documented in the submission package.

Software: R 4.6.0 with `survey`, `haven`, and `broom`; see `results/session_info.txt`.
