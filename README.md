# Reference Sensitivity in the Projection of Pediatric Gingival Crevicular Fluid Fourier-Transform Infrared Spectra

**Reproducibility package prepared for the accompanying manuscript (RC3 baseline, 19 October 2026; packaging validation completed 20 October 2026).**

The CSV-only analysis was successfully executed by the author on macOS in R 4.5.1:
`PUBLIC_CANDIDATE_R_20261018_220355_56691`. The independent audit confirmed agreement with
the verified 17–18 October results. RC3 preserves all five analytical input files, all frozen
numerical checkpoints, every existing file in `R/`, and all reference figures byte for byte.
Only run-finalization I/O, Git packaging, and documentation are revised; one new I/O helper is added.

The corresponding scientific manuscript is **Reference Sensitivity in the Projection of Pediatric Gingival Crevicular Fluid Fourier-Transform Infrared Spectra**,
author-approved scientific and language revision dated 19 October 2026. This repository does not replace or edit the manuscript.

## Minimal packaging check — no scientific recalculation

Unpack into a new folder. From its root, run:

```bash
Rscript tools/check_release.R
```

This checks package hashes, original input hashes, R syntax, and the exact run-finalization helper
on disposable non-scientific test files. It does not refit PCA or recompute projections. A successful
check produces `runs/PACKAGING_CHECK_.../PACKAGING_CHECK_PASS.txt`.

The packaging-only R self-test was completed on the author's Mac on 20 October 2026: `PACKAGING CHECK PASS`. The test validated the 60-file RC3 baseline manifest, the five analytical input hashes, the R parser, and both positive and negative run-finalization tests. It did not execute PCA or any new scientific analysis. This verification is distinct from the previously completed numerical R analyses.

## Full scientific reproduction

For a reader who wishes to reproduce the established analyses from the packaged inputs:

```bash
Rscript run_all.R
```

Base R plus the standard `shasum` or `sha256sum` utility is sufficient. No new analysis is introduced.
The final `RUN_SHA256.csv` now includes the completed log and pass flag. The latest-run pointer
is written outside the sealed run folder after verification. Details: `RUN_ON_MAC.md`.

## Inputs and analytical scope

- Historical 18 × 793 adult-series SNV matrix; A18 (age 15) is used only in the historical mixed-age sensitivity reference.
- Primary adult reference: A01–A17, n = 17; ages 21–70; nine aged 21–25 and eight aged 50–70; 13 women and four men. A04 and A07 are retained.
- Same 25 external pediatric SNV spectra on the exact 793-position adult grid in every model.
- Five reference compositions, K = 2–5, and complete leave-one-adult-out analysis of the primary reference.
- K = 3 spectral residual localization, center–subspace decomposition, and descriptive age/sex summaries.

`data/` contains the normalized spectra, exact grid, and the minimal age/sex and membership fields.
`verification/` contains output comparators, never fitting inputs. `figures/reference_verified/`
contains the unchanged original figures. `docs/DATA_DICTIONARY.md` describes the published-input proposal.
The checkpoint whose basename starts with `PRIVATE_` contains analytical IDs and computed Q values
only, **not** medical histories. Its exact Git exception is documented in `docs/RC2_PACKAGING_CHANGES.md`.

The package does not recreate raw DPT/OriginPro processing or pediatric interpolation/SNV from
pre-normalized absorbance. These upstream steps are outside its reproduction boundary. The analysis
makes no diagnostic, causal, population-reference, or universal FTIR-transfer claim.

## Optional Quarto

`analysis/report.qmd` is an optional consumer of a completed scientific run. Quarto and `knitr`
are not required by `run_all.R` or the packaging check. No rendered Quarto HTML has been verified
or promised as a required deliverable.

## Licenses, data governance, and versioned release

Original R source code is distributed under the **MIT License** (`LICENSE`). Processed coded spectra, the minimal demographic metadata, analytical output CSVs, verified figures, and authors' documentation are released under **CC BY 4.0** (`DATA_LICENSE.md`), subject to the author/rightsholder permissions and privacy safeguards described there. The license does not override ethical or consent restrictions. Raw clinical records, source `.dpt` acquisition files, personal-ID crosswalks, and additional clinical categories are not included.

The private GitHub repository and its commits are not a citable versioned release. Before changing visibility to public, the corresponding author must confirm the exact intended publication scope, including the historical 15-year-old record (A18), and approve the final public dataset. Cite the tagged release once created. A successful code check is not publication authorization. See `docs/RELEASE_GATES.md`.
