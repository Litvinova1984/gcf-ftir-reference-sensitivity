# Data dictionary (proposed publication-level fields)

| File | Columns | Definition / confidentiality note |
|---|---|---|
| `data/adult_SNV_18x793.csv` | `adult_id`, `wn_*` | 18 archived normalized spectra, one individual analytical row each. A18 minor, primary excludes it. |
| `data/pediatric_SNV_25x793.csv` | `pediatric_id`, `wn_*` | Same 25 external participants for every adult model, row-wise SNV **after prior interpolation** to adult grid. |
| `data/wavenumbers_exact.csv` | `variable_index_1based`, `wavenumber_cm_1`, `spectral_column` | 793 full-precision adult target coordinates and column ordering. CSV column names round to 5 decimal places; *not* the authoritative grid. |
| `data/metadata/adult_demographics_candidate.csv` | `adult_id`,`age_years`,`sex`,`primary_adult_member`,`historical_age15_record` | 18 coded adult-archive rows (one adolescent). Exact age and sex under privacy/ethics review. |
| `data/metadata/pediatric_demographics_candidate.csv` | `pediatric_id`,`age_years`,`sex` | 25 coded children; exact age is private pending institutional approval. |

`sex` values F/M reflect the source clinical recording, without inference about gender identity. Age is age at sampling. `primary_adult_member` yes for A01–A17 only. All `wn_*` values are **dimensionless normalized intensity**; they are not raw absorbance.

No source/sample filename, original clinical code, medical record number, direct identifier, exact date of collection or clinic location is part of the proposed release files.

**Optional clinical fields (not cleared):** original periodontal/caries/anamnesis groups in the adult cohort, plus pediatric pubertal/endocrine/UCTD/periodontal categories, are provisionally parked under `RESTRICTED_REVIEW_DO_NOT_PUBLISH/`. They are *not necessary* to reproduce the paper's primary numerical PCA or age/sex summary outputs. Any decision to release them must be signed off after assessing rare category combinations and prior publication linkage.
