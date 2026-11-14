# Cultivation and Depensation in a Six-Species Size-Spectrum Model with External Cod Recruitment (`mizer`)

This repository contains a theoretical six-species physiologically structured size-spectrum model built in [`mizer`](https://sizespectrum.org/mizer/) (`3.4.0.9003`), a scientific manuscript ([`paper/manuscript.pdf`](paper/manuscript.pdf)), provenance and audit ledgers ([`paper/provenance_ledger.md`](paper/provenance_ledger.md), [`models/methods_audit.json`](models/methods_audit.json), [`models/reproducibility_manifest.json`](models/reproducibility_manifest.json)), and an interactive educational scenario explorer ([`docs/index.html`](docs/index.html)).

**Scope note:** All calibration targets (`biomass_observed`) are synthetic illustrative reference values inspired by boreal continental-shelf functional guilds rather than empirical survey series from a single geographic stock boundary, and all simulations use relative simulation years ($t \in [0, 61]$ for the collapse ramp and 34-year directed-cod closure; $\tau \in [0, 40]$ for post-closure policy experiments, extended to $\tau = 340\text{ yr}$ to verify asymptotic outcomes).

---

## Key Findings

- **Coexisting Zero-Cod-Fishing Steady States ($E_{\text{cod}} = 0$):** Under one illustrative parameterization with identical biological, interaction, resource, and non-cod fishing parameters, numerical local stability analysis of the 727-dimensional consumer–resource ODE Jacobian (`getStability()`, `n_active = 727`: 505 active consumer bins + 222 resource bins) identifies two locally asymptotically stable fixed points (`residual < 1e-8`, `attractor = "fixed_point"`), separated by an unstable saddle equilibrium (`solver = "newton"`):
  1. **Zero-Cod-Fishing Cultivated State:** $\text{SSB}_{\text{cod}} = 1018.2\text{ kt}$, $\max\operatorname{Re}(\lambda) = -0.1809\text{ yr}^{-1}$ (baseline pelagic, sculpin, and flatfish fishing active).
  2. **Zero-Cod-Fishing Depensatory State:** $\text{SSB}_{\text{cod}} = 15.8\text{ kt}$, $\max\operatorname{Re}(\lambda) = -0.0084\text{ yr}^{-1}$, where the imposed exogenous cod recruitment input (`RefugeBevertonHoltRDD`, $1.1\%$ of reference recruitment) accounts for $13.5\%$ of total cod recruitment. Setting exogenous input to zero preserves bistability while reducing the lower attractor to cod extinction ($0\text{ kt}$, closed-population invasion $R_0 = 0.82 < 1$).
- **Path-Specific Basin Geometry & Predator-Resolved Bottleneck:** Equal- and inverted-SSB states with contrasting prey spectra converge to opposite attractors, demonstrating that scalar Cod SSB does not define a universal separatrix. Across the $[0.001, 15)\text{ g}$ early-life window, cumulative mortality-to-growth hazard $H = \int (\mu/g)\,\mathrm{d}w$ rises from $0.19$ to $2.82$ ($88.7\%$ direct predation increase led by mackerel, herring, and sculpin; $0.9\%$ growth slowdown alone; $10.4\%$ mortality $\times$ slower-growth synergy), cutting stationary cohort survival 14.0-fold ($83.1\%$ to $6.0\%$). Cod feeding level $f(w)$ is $1.4\%\text{--}11.7\%$ lower over $0.05\text{--}5\text{ g}$ (up to $15.0\%$ lower below $15\text{ g}$) yet higher for surviving $> 200\text{ g}$ cod; net reproductive ratio $R_0$ falls from $1.21$ to $0.88$.
- **Matched 6-Gear Policy Experiments, Bycatch Fragility, & Structural Limits:** Using a shared `method = "tr_bdf2"`, `dt = 0.05` solver across all scripts:
  - Separating `Capelin_Gear` from `Clupeid_Mackerel` shows that sparing capelin ($E_{\text{cap}} = 0.15$) while harvesting herring/mackerel ($0.65$) and sculpin ($0.45$) accelerates cod rebuilding ($745.8\text{ kt}$ vs. $433.6\text{ kt}$ at $\tau = 40\text{ yr}$; both reaching $776.3\text{ kt}$ asymptotically) while maintaining Phase-I capelin cutoff biomass $\ge 1327\text{ kt}$.
  - Rebuilding cod entails an explicit aggregate-tonnage trade-off relative to harvesting hyperabundant forage fishes, and the fixed 20-year joint pulse fails asymptotically under Phase-I cod bycatch $\ge 0.010\text{ yr}^{-1}$, whereas a state-dependent feedback harvest control rule ($\ge 280\text{ kt}$ observed SSB for 2 consecutive years) restores robustness under $25\%$ observation error.
  - Bistability is conditional: it disappears into monostable spontaneous recovery when early-life predation overlap is scaled by $s_\theta \le 0.75$ (for $r_{\text{prey}} \ge 0.70$, and for all tested $r_{\text{prey}}$ when $s_\theta \le 0.50$), under symmetric fish-larval overlap, wider feeding kernels ($1.15\times \sigma_i$), stronger cod recruitment compensation ($r_{\text{Cod}} = 0.55$), fast resource regeneration ($r_{pp} = 20$), or when exogenous cod recruitment reaches $\ge 1.5\%$ of calibrated recruitment.

---

## Repository Structure

```text
.
├── archive/v1_reviewed/                     # Preserved read-only v1 snapshot + SHA256SUMS.txt
├── R/
│   ├── 00_model_functions.R                 # RefugeBevertonHoltRDD, SaturatingShelterRDD, canonical 4/6-gear & hazard helpers
│   ├── 01_build_and_calibrate_model.R       # Builds, calibrates, and refines 4 steady states (residual < 1e-8) + t=61 state
│   ├── 02_analyse_bistability_and_scenarios.R # 727-dim stability, Newton saddles, counterexamples, 6-gear policies, bycatch & HCR
│   ├── 03_generate_figures_and_web_data.R   # Generates figures 1-6 (PDF+PNG) & exact grid-snapped web explorer data
│   └── 04_sensitivity_and_audit.R           # Solver/grid/scheme convergence, structural sweeps, fig 7, audit & manifest JSONs
├── models/
│   ├── cod_f15_params.rds                   # Baseline sustainable fishery steady state (E_cod = 0.15 /yr)
│   ├── cod_high_f0_params.rds               # Zero-cod-fishing cultivated steady state (E_cod = 0.00 /yr)
│   ├── cod_low_f0_params.rds                # Zero-cod-fishing depensatory steady state (E_cod = 0.00 /yr)
│   ├── cod_low_f15_params.rds               # Reopened-fishery low steady state (E_cod = 0.15 /yr)
│   ├── cod_saddle_f0_params.rds             # Unstable saddle fixed point at E_cod = 0.00 /yr (solver = "newton")
│   ├── cod_saddle_f15_params.rds            # Unstable saddle fixed point at E_cod = 0.15 /yr (solver = "newton")
│   ├── cod_t61_params.rds                   # Post-closure transient state at relative year t = 61
│   ├── sim_collapse_closure.rds             # Collapse + closure MizerSim trajectory (t = 0..61, saveSim)
│   ├── sims/                                # Named policy MizerSim trajectories (saveSim)
│   ├── methods_audit.json                   # Programmatic parameter, exponent, grid, & Jacobian audit
│   └── reproducibility_manifest.json        # R/package versions, mizer RemoteSha, custom rate source, & SHA-256 hashes
├── paper/
│   ├── manuscript.tex                       # Revised LaTeX manuscript (with line numbers)
│   ├── references.bib                       # Verified BibTeX bibliography
│   ├── provenance_ledger.md                 # Parameter provenance & empirical claim-to-source ledger
│   ├── referee_report.md                    # Adversarial referee report on v1
│   ├── referee_response.md                  # Point-by-point response matrix keyed to Comments 1-15
│   ├── manuscript.pdf                       # Compiled manuscript PDF
│   └── figures/                             # Figures fig1-fig7 (.pdf and 300-dpi .png)
└── docs/
    ├── index.html                           # Educational theoretical size-spectrum explorer
    ├── styles.css                           # Dashboard layout & styling
    ├── app.js                               # Exact discrete grid-node viewer (zero cross-basin interpolation)
    ├── model_data.js                        # Precomputed canonical tr_bdf2 (dt=0.05) scenario grid
    └── data/mizer_scenarios.json            # Raw JSON export of trajectories, static reference spectra, & scans
```

---

## Reproducing the Full Pipeline

Run the four scripts sequentially from the repository root in R (`mizer 3.4.0.9003`, `ggplot2`, `patchwork`, `jsonlite`, `digest`):

```bash
Rscript R/01_build_and_calibrate_model.R
Rscript R/02_analyse_bistability_and_scenarios.R
Rscript R/03_generate_figures_and_web_data.R
Rscript R/04_sensitivity_and_audit.R
```

To recompile the manuscript PDF (from `paper/`):

```bash
pdflatex -interaction=nonstopmode manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex
```
