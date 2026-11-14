# Breaking the Moratorium Trap: Northwest Atlantic Cod Bistability & Ecosystem-Based Recovery (`mizer`)

This repository contains a calibrated six-species physiologically structured size-spectrum model built in [`mizer`](https://sizespectrum.org/mizer/) (v3.4), a complete scientific manuscript formatted for the *Canadian Journal of Fisheries and Aquatic Sciences* ([`paper/manuscript.pdf`](paper/manuscript.pdf)), and a self-contained interactive policy simulator ([`docs/index.html`](docs/index.html)).

Together, they address two central questions in marine ecology and fisheries management:
1. **Why did Northwest Atlantic cod (*Gadus morhua*) fail to recover despite the July 1992 commercial fishing moratorium?**
2. **How can targeted, time-limited ecosystem harvest on overabundant forage and mesopredator fishes release the juvenile cod bottleneck and restore a self-sustaining cod fishery?**

---

## Key Scientific Findings

- **Rigorous Cultivation–Depensation Bistability at $F_{\text{cod}} = 0$:** Under *identical* biological, bioenergetic, resource, and habitat-overlap parameters, continuous-time Jacobian eigenvalue analysis (`getStability()`) proves that the shelf community supports two coexisting locally asymptotically stable equilibria when directed cod fishing is zero ($F_{\text{cod}} = 0\text{ yr}^{-1}$):
  1. **Unexploited Cultivated Equilibrium:** $\text{SSB}_{\text{cod}} = 1018.2\text{ kt}$, $\max\operatorname{Re}(\lambda) = -0.1809\text{ yr}^{-1}$. Adult cod ($>1\text{ kg}$) crops down pelagic planktivores and benthic sculpins, keeping nursery predation mortality on $0.05\text{--}5\text{ g}$ cod low ($0.25\text{ yr}^{-1}$).
  2. **Collapsed Moratorium Trap Equilibrium:** $\text{SSB}_{\text{cod}} = 15.8\text{ kt}$ ($2\%$ of historical fishery SSB), $\max\operatorname{Re}(\lambda) = -0.0084\text{ yr}^{-1}$, sustained only by sheltered inshore coastal refuges (`RefugeBevertonHoltRDD`).
- **The Dual Trophic Squeeze:** Escalating otter-trawl exploitation past the saddle-node collapse fold ($F_{\text{cod}}^* = 0.250\text{ yr}^{-1}$) triggers a **5.4-fold biomass surge** in capelin, herring, mackerel, and sculpins. This guild reversal inflates instantaneous predation mortality on $0.05\text{--}5\text{ g}$ cod to **$7.33\text{ yr}^{-1}$** ($29.7\times$ the cultivated level; $10.6\times$ the historical sustainable fishery level of $0.69\text{ yr}^{-1}$) while grazing down shared zooplankton ($0.01\text{--}10\text{ g}$) and slowing juvenile cod growth through the predator pit.
- **Adaptive Ecosystem-Based Recovery:** Status-quo moratoria ($F_{\text{cod}} = 0$) remain trapped indefinitely below the $29.5\text{ kt}$ separatrix (`12.0 kt` SSB in 2066). By contrast, a **20-year selective harvest pulse (2026–2046)** on overabundant pelagic planktivores ($F_{\text{pelagic}} = 0.65\text{ yr}^{-1}$) and benthic sculpins ($F_{\text{mesopred}} = 0.45\text{ yr}^{-1}$) breaks the early-life predation bottleneck, crosses the $249.0\text{ kt}$ fishery separatrix within two decades, and allows forage/sculpin effort to **taper back to baseline** ($0.15$ and $0.05\text{ yr}^{-1}$) in 2046 while reopening a sustainable, high-value cod fishery ($F_{\text{cod}} = 0.15\text{ yr}^{-1}$).

---

## Repository Structure

```text
.
├── R/
│   ├── 00_model_functions.R                 # Smooth C^infinity inshore refuge recruitment (RefugeBevertonHoltRDD)
│   ├── 01_build_and_calibrate_model.R       # Builds & calibrates 6-species model and saves 3 steady-state MizerParams
│   ├── 02_analyse_bistability_and_scenarios.R # Jacobian eigenvalues, separatrix bisection, scanModel folds, 2026-2066 scenarios
│   └── 03_generate_figures_and_web_data.R   # Generates 6 multi-panel figures (PDF+PNG) and exports web simulation grids
├── models/
│   ├── nw_atlantic_cod_params.rds           # Sustainable pre-collapse fishery steady state (F_cod = 0.15 /yr)
│   ├── nw_atlantic_cod_high_f0_params.rds   # Unexploited cultivated high-cod steady state (F_cod = 0.00 /yr)
│   ├── nw_atlantic_cod_moratorium_params.rds# Collapsed moratorium low-cod steady state (F_cod = 0.00 /yr)
│   └── analysis_results.rds                 # Cached eigenvalues, scanModel objects, and 6 policy MizerSim trajectories
├── paper/
│   ├── manuscript.tex                       # LaTeX source (Canadian Journal of Fisheries and Aquatic Sciences format)
│   ├── references.bib                       # BibTeX bibliography
│   ├── manuscript.pdf                       # Compiled 16-page scientific paper with embedded vector figures
│   └── figures/                             # Publication figures fig1-fig6 (both .pdf and 300-dpi .png)
└── docs/
    ├── index.html                           # Self-contained interactive policy brief & mizer scenario explorer
    ├── styles.css                           # Scientific dashboard styling (responsive desktop/tablet/mobile)
    ├── app.js                               # Trilinear mizer grid interpolator & HiDPI Canvas visualizations
    ├── model_data.js                        # Precomputed 7x6x6x2 mizer simulation grid (works offline via file://)
    └── data/mizer_scenarios.json            # Raw JSON export of all model trajectories, spectra, and bifurcations
```

---

## Reproducing the Full Study

Run the three pipeline scripts sequentially from the repository root in R (`mizer >= 3.4.0`, `ggplot2`, `patchwork`, `jsonlite`):

```bash
Rscript --vanilla R/01_build_and_calibrate_model.R
Rscript --vanilla R/02_analyse_bistability_and_scenarios.R
Rscript --vanilla R/03_generate_figures_and_web_data.R
```

To recompile the manuscript PDF (from `paper/`):

```bash
pdflatex -interaction=nonstopmode manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex
```

To explore the interactive policy website locally, open [`docs/index.html`](docs/index.html) directly in any modern web browser (no local server or R backend required), or publish the `docs/` folder via GitHub Pages.
