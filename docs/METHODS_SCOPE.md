# Methods-to-code map and interpretive limits

| Stage | Source inputs | R implementation | Published analysis |
|---|---|---|---|
| Adult row screening | 18-row historical SNV matrix + A18 mask | `R/00_core.R` | n=17 primary, n=18 sensitivity |
| PCA | adult-only matrix | `fit_adult` | `prcomp(center=TRUE, scale.=FALSE, rank.=5)` |
| External projection | fixed pediatric 25x793 SNV | `project_fixed` | PC scores, Q, Q_fraction, T² |
| LOO | n17 adults | `adult_loo_stats` | 17 adult models, same 25 external spectra |
| Other references | frozen adult subsets | `R/01_execute.R` | A04, A07, A04+A07, historical 18 |
| K sensitivity | K=2,3,4,5 | `R/01_execute.R` | all five reference compositions |
| Spectral localization | model residuals on 793 positions | `R/04_diagnostics.R` | mean squared residual per variable, two fixed windows |
| Center–subspace identity | baseline and reduced mean/projector | `R/04_diagnostics.R` | center only + subspace only + interaction = change in mean Q |
| Age/sex descriptions | separately coded metadata with analytical IDs | `R/05_metadata.R` | n=25 repeated descriptions, no p-values |
| Figures | numerical exports only | `R/03_figures.R`, `R/04_diagnostics.R` | three original figures + two diagnostics |

**Caveats**: The 90% adult explained-variance rationale was described **after** K=3 selection; not a preregistration. This project estimates no biological age effects, inter-instrument calibration, disease classification, oncology or FTIR domain D1–D5. Refitting adult PCA is not independent cross-cohort validation. The two groups have different collection/preparation/instrument histories (fully confounded with age group). Child observations are reused across all models; adult LOO fits overlap. Descriptive residuals should not be treated as independent inferential samples.

The wavenumber-resolved Q consists of **unweighted discrete squared SNV residuals**, not per-cm⁻¹ spectral density, biochemical abundance, or a peak with 0.01 cm⁻¹ measurement resolution. The center–subspace identity decomposes **absolute Q** (not median Q_fraction or T²), and its interaction is not causally allocated to either factor.
