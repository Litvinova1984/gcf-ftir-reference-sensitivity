# RC2 packaging changes

## 1. Finalization order

Previously, RUN_SHA256 was written before the last log message and the global PASS flag.
`run_all.R` now calls `finalize_verified_run()` only after every unchanged numerical check.
The helper writes the pass flag, completes the log, writes RUN_SHA256, rechecks full coverage and
hashes, then updates the latest-run pointer outside the sealed directory. It removes an incomplete
pass flag if finalization fails. There are no subsequent writes inside a successfully sealed run.

`tools/check_release.R` exercises that exact helper on temporary non-scientific files. It checks
positive finalization, rejection of late log edits, unlisted files, missing pass flags, failed runs,
and cleanup following a simulated I/O failure. It also checks package/input hashes and real R syntax.
This is a packaging test, not a surrogate numerical replication. It still requires an R execution on Mac.

## 2. Git omission

An exact exception was added for:
`verification/frozen_R_20261018/PRIVATE_pediatric_Q_decomposition_K3.csv`.
The general `PRIVATE_*` rule remains active. The checkpoint contains model/analytical IDs and
computed Q/decomposition metrics only; it has no age, sex, name, medical history, or clinical-category
field. No file contents or R references were renamed. Standard Git staging is tested during assembly.
`.gitattributes` disables automatic line-ending conversion to preserve hash-pinned bytes in Git.

## 3. Documentation

README, execution instructions, provenance, reproducibility status, and release decisions now
recognize the successful 18 October CSV-only Mac run. The new manuscript title is used.
The editor note no longer incorrectly calls the n = 18 historical mixed-age model adult-only.
Quarto remains optional; its HTML has not been represented as verified. Original build checks are
retained with their historical timestamp, not rewritten as a new scientific validation.

## Preserved

All five data/grid/demographic files; all eight pre-existing R-folder modules; all 17 frozen CSV
checkpoints; all ten PNG/PDF reference-figure files. The updated orchestrator does not alter
analytical calls, numerical tolerances, model membership, K choices, normalization, or outputs.
The original manuscript, Supplement, reviewer report, source ZIP, and Mac result ZIP are not changed.

## Boundaries

The new helper has not been executed in this assembly environment because R is unavailable.
The included R packaging test is the one remaining runtime check. There is no new scientific result,
publication approval, license, DOI, or online release in this candidate.
