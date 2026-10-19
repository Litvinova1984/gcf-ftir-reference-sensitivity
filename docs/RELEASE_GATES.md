# Release decisions — RC3 private pre-submission candidate

| Item | Status | Remaining action |
|---|---|---|
| Scientific manuscript | Author-approved language revision completed 19 October 2026 | Use updated Main and Supplement dated 19 October; no new analytical cycle |
| Spectral/membership inputs | Preserved byte for byte | No redefinition of cohorts |
| CSV-only numerical reproduction | PASS on the author's Mac, 18 October 2026 | Keep original run as evidence |
| Final run manifest order | Fixed in RC2 | Execute packaging-only R self-test |
| Required checkpoint excluded by Git | Fixed; local Git staging test recorded | Keep exact exception; do not remove the general PRIVATE_* protection |
| RC2 files and Git round trip | See packaging report | R self-test is distinct from the earlier numerical PASS |
| Minimum metadata scope | Coded age/sex linked to spectra included | Confirm explicit institutional consent and disclosure approval for exact ages, including minors, before publication |
| Extra clinical metadata | Outside RC2 | Not needed to reproduce this manuscript; no automatic addition |
| Author list/roles/funding/COI | Transferred from the coauthor-provided 18 October manuscript | Corresponding author performs final sign-off |
| Data and code licenses | Not selected by this build | Authors/rightsholders approve the two licenses |
| Public release and version link | Not created | Publish only after documented approval, passing R packaging test, and selecting licenses; insert actual release URL |
| Related publications | Disclosure required | Confirm status at submission and identify reused datasets in cover letter |
| Quarto HTML | Optional and not verified | Not a prerequisite and not a promised deliverable |

This table separates technical completion from author/institutional declarations. It does not claim
that permissions are absent; approved final wording has not yet been supplied here. No medical
records or personal-ID crosswalks are included. The pre-existing demographic fields remain unchanged.


RC3 retains the same numerical verification boundary as RC2 and is NOT a claim of permissions or a published GitHub release.
