# RC3 author and manuscript alignment — 19 October 2026

The author group supplied authorship order, affiliations, corresponding author, CRediT roles, grant number, and declarations separately to the manuscript. RC3 synchronizes the repository title with the author-approved manuscript:

**Reference Sensitivity in the Projection of Pediatric Gingival Crevicular Fluid Fourier-Transform Infrared Spectra**

Changes relative to RC2:
- README and reproducibility documentation brought into alignment with the manuscript revised on 19 October.
- Internal assembly reports not needed for R reproducibility were removed from the intended public package (preserved in private RC2).
- Original RC2 R scripts, the five input data/grid/metadata files, all frozen numerical checkpoints, and all reference PNG/PDF figures retain their **exact bytes**.
- No PCA, leave-one-out, spectral diagnostic, or demographic computation was rerun in RC3.
- Rebuilt package SHA manifest after documentation changes.

Still required before public release: (1) run `Rscript tools/check_release.R` on the author's Mac and review its PASS report, (2) approve exact age/sex/spectra public sharing under the source permissions, (3) select data and code licenses, and (4) create and verify the real tagged GitHub release URL.

This is a technical release candidate, not a published data release.
