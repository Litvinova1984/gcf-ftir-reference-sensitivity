# Exact source and computation boundaries

1. Archived adult analytical SNV RDS: `MDPI_Data_GCF_Integrated_R_Quarto_PRIVATE_2026-10-15.zip / input_internal/adult_reference_X_snv.rds` (18×793). Converted to CSV with full binary-double round-trip decimal serialization; no smoothing, filtering, extra SNV, imputation or refit during packaging. The normalized adult rows stay in the original order; A18 is the historical 15-year-old.
2. Adult grid: same project's `adult_reference_wavenumbers.rds`, 793 full-precision values, descending order, both fixed analytical windows.
3. Pediatric 25×793 SNV: **original R-verified** output `runs.zip / runs/R_20261017_211624_53719/pediatric_on_adult_grid_SNV.csv`, not from alternative archived tables; column names matched to adult grid in source order.
4. Age/sex metadata: from existing PRIVATE R-verification inputs of 18 October, whose keyed matching against original source metadata and verified individual projections passed on Mac.
5. Frozen numerical checkpoints: from the Mac R execution of 17 October and two R extensions of 18 October. No checkpoints are training or projection inputs.

**Publication boundary:** processed SNV spectra → adult-only PCA → projection → reference-composition and K sensitivity → squared-residual localization and Q decomposition → age/sex descriptive summaries. The earliest adult DPT, OriginPro and pediatric original 1467-position absorbance plus interpolation/SNV are not included, so the repository does not reproduce those upstream operations.


## RC2 packaging lineage

RC2 is derived from `MDPI_Data_GCF_PUBLIC_CANDIDATE_REVIEW_2026-10-18.zip`, verified on Mac
in run `PUBLIC_CANDIDATE_R_20261018_220355_56691`. No clinical records or new source inputs
were added. Source data, existing R modules, numerical checkpoints, and reference figures retain
their bytes. `run_all.R` now delegates final log/pass/manifest ordering to `R/07_release_io.R`.
The prior run and its original hash manifest remain unaltered in the author's archive; RC2 does not
pretend that the old run used the newly added helper.


## RC3 notes
RC3 synchronizes the author-approved title and publication documentation dated 19 October 2026. No spectral, checkpoint, figure, or computation parameters were changed. Historical internal verification files not required for R reproducibility are retained privately in RC2.
