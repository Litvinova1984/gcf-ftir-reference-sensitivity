# Publication gates and verified package status (20 October 2026)

This checklist records the steps completed **before** a public GitHub Release. The initial GitHub push to `main` does not constitute a published or citable release.

| Item | Status | Action before public release |
|---|---|---|
| Scientific manuscript | Author-approved language and scientific version of 19 October 2026 | Verify final author declarations and repository URL before submission |
| Core numerical analysis | Original R analyses and extensions verified on Mac on 17–18 October 2026 | Do not change models, spectral inputs, or numeric outputs |
| RC3 packaging check | **PASS on Mac, 20 October 2026**; 60 baseline file hashes and five input hashes passed; R syntax and positive/negative I/O checks passed | Preserve packaging log in the author's private archive |
| Required checkpoint and Git omission | Correct exact `.gitignore` exception included and verified | Do not remove verified numerical checkpoint |
| GitHub | Private `main` branch pushed; licenses added in commit `03faf37` | Synchronize release documentation, then create a versioned GitHub Release |
| R code license | **MIT**, author selected | Maintain `LICENSE` |
| Processed spectra and research-output license | **CC BY 4.0**, author selected | Maintain `DATA_LICENSE.md` |
| Pediatric coded spectra and metadata | Scope corresponds to the previously reviewed Applied Sciences disclosure; age/sex only in the new project's metadata inputs | Corresponding author ensures existing approvals cover the actual public release |
| Adult coded spectra and metadata | Author has authorized release of adult materials | Confirm the historical adolescent record A18 (age 15) is covered by the same authorization before public exposure |
| Additional clinical records, raw `.dpt`, personal-ID crosswalks | **Not included** | Keep excluded |
| External study relationship | Adult-source IJMS 2024 and related pediatric Applied Sciences manuscript disclosed | Update related-manuscript status in Cover Letter if necessary |
| Submission link | Versioned release URL not yet available | Create tag, review release contents, insert real URL in Main, Supplement and Cover Letter as applicable |

No numerical rerun is needed merely to revise license and publication documents. The base `PACKAGE_SHA256.csv` lists exact file bytes and must be regenerated for any modified documented files and newly licensed files before repeating the packaging check. This record is a publication checklist, not an independent ethics approval.
