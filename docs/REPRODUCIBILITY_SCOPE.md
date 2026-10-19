# Numerical acceptance and reproduction boundary

The original RDS-based workflow and the two extensions were verified on 17–18 October 2026.
The normalized-CSV entry point was independently executed on the author's Mac in R 4.5.1,
run `PUBLIC_CANDIDATE_R_20261018_220355_56691`. It passed the 12 primary table comparisons,
the spectral diagnostics, and the demographic checks. No checkpoint is a fitting input.

RC2 changed packaging I/O only; RC3 updates public documentation and manuscript title. All eight existing `R/` modules, five input files, frozen numerical
checkpoints, and reference figures remain byte-identical to the verified candidate. The orchestrator
now calls a new helper that writes the final log and pass flag before generating the run manifest.
`tools/check_release.R` tests that helper without running any scientific analysis. Its execution on
the author's Mac remains the final packaging check; no execution is invented in the build report.

Unchanged numerical tolerances in the original code:
- Ordinary fields in the 12 primary comparisons: absolute difference <= 1e-6.
- Maximum principal angles: absolute difference <= 1e-4 degrees.
- Complete spectral diagnostic tables and individual decompositions: <= 1e-6.
- Pediatric demographic and adult age/influence summaries: <= 1e-10.

These are numerical gates, not confidence intervals, significance tests, or acceptance limits for a
clinical reference. The archived R and Python agreement is not an independent biological cohort.
The main input boundary is the normalized spectral matrices plus required demographic fields;
raw DPT/OriginPro processing and upstream pediatric absorbance interpolation/SNV are not reproduced.

The 18 October scientific result archive has been retained unchanged outside RC2. Its original
CODE_SHA256 identifies the original code, not the new packaging helper; it must not be relabelled
as a fresh RC2 run. The known old run-log hash issue is a finalization defect, not a numerical discrepancy.
