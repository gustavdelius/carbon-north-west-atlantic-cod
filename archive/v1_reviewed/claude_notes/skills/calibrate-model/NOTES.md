# Project Notes — `calibrate-model`

- **Reference calibration state (`R/01_build_and_calibrate_model.R`):** The 6-species Northwest Atlantic shelf model (`Capelin`, `Herring`, `Mackerel`, `Sculpin`, `Plaice`, `Cod`) is calibrated at moderate exploitation (`F_cod = 0.22 /yr`, `F_pelagic = 0.15 /yr`, `F_mesopredator = 0.05 /yr`, `F_flatfish = 0.15 /yr`) by alternating `matchGrowth()`, `matchBiomasses()`, `steadySingleSpecies(..., keep = "biomass")`, `setBevertonHolt()`, and `tuneSteadyState()`.
- **Physiological `erepro` verification:** All calibrated species maintain `erepro << 1`, spanning `2.12e-6` (`Cod`) to `1.70e-3` (`Capelin`).
- **Beverton–Holt reproduction levels:** Prey guilds use `reproduction_level = 0.88` (allowing realistic 5.4-fold predation-release surges when adult cod collapses) and `Cod` uses `reproduction_level = 0.28` paired with the inshore spawning refuge `r_refuge` (`RefugeBevertonHoltRDD`, set to `1.1%` of reference cod recruitment).
- **Saved steady states (`models/`):**
  - `models/nw_atlantic_cod_params.rds` — Sustainable pre-collapse fishery (`F_cod = 0.15 /yr`, Cod SSB `776.3 kt`, residual `5.90e-7`).
  - `models/nw_atlantic_cod_high_f0_params.rds` — Unexploited cultivated equilibrium (`F_cod = 0.00 /yr`, Cod SSB `1018.2 kt`, residual `1.66e-7`).
  - `models/nw_atlantic_cod_moratorium_params.rds` — Collapsed moratorium trap equilibrium (`F_cod = 0.00 /yr`, Cod SSB `15.8 kt`, residual `3.08e-6`), using identical biological, interaction, resource, and rate-function parameters as `nw_atlantic_cod_high_f0_params.rds`.
