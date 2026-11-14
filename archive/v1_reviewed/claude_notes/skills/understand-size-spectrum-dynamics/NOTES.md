# Project Notes — `understand-size-spectrum-dynamics`

- **Cultivation–depensation bistability mechanism in this model:**
  - Adult `Cod` (`> 1 kg`, `beta = 60`, `sigma = 1.6`) exerts strong top-down predation on `Capelin` (`theta = 0.45`), `Herring` (`0.43`), `Mackerel` (`0.40`), and `Sculpin` (`0.38`).
  - Conversely, adult planktivores (`Herring` and `Mackerel` `theta = 1.00`, `Capelin` `theta = 0.55`, `beta = 350–2500`) prey on pelagic `Cod` eggs and larvae (`1 mg–1 g`), while benthic `Sculpin` (`theta = 0.45`, `beta = 65`) preys on settled 0-group `Cod` (`0.5–15 g`).
- **Dual trophic squeeze in the moratorium trap (`F_cod = 0`):**
  - Mean instantaneous predation mortality on `0.05–5 g` nursery `Cod` rises from `0.25 /yr` (unexploited cultivated state) and `0.69 /yr` (`F_cod = 0.15 /yr` sustainable fishery) to `7.33 /yr` in the collapsed moratorium state (`29.7x` vs. cultivated).
  - Simultaneously, grazing by `> 4600 kt` of planktivores depresses the upper zooplankton resource (`0.01–10 g`), reducing juvenile `Cod` feeding level `f(w)` by `15–25%` while elevating `f(w)` of rare surviving `> 1 kg` adult cod.
- **Equilibrium fold vs. management horizon:** Although the collapsed branch ends at an equilibrium pelagic fishing fold near `F_pelagic = 0.175 /yr` (`F_cod = 0`), critical slowing down near the saddle-node separatrix (`29.5 kt` Cod SSB at `F_cod = 0`; `249.0 kt` at `F_cod = 0.15 /yr`) requires a stronger 20-year pulse (`F_pelagic = 0.65 /yr`, `F_mesopredator = 0.45 /yr`) to achieve full rebuilding within two decades before tapering back to baseline.

## Referee audit correction — 13 November 2026

The earlier entries summarize the manuscript rather than independently establishing the mechanism. See `paper/referee_report.md`, particularly comments 6–12. Saved-state feeding-level ratios (low/high) are about 0.883–0.987 over 0.05–5 g and 0.850–0.999 over 0.002–15 g; the claimed 15–25% reduction is not supported across those windows. Predator-specific mortality contributions have not been demonstrated by the total-mortality plot. Refuge recruitment is a constant external input at egg size, not an explicit sheltered population. The named pulse reaches about 259.1 kt SSB in 2046 and 475.3 kt in 2066; the exported grid gives 422.5 kt for the same policy with a different time step. These do not establish full rebuilding within two decades. Fold, general SSB threshold, and operational policy claims require the additional analyses specified in the report.
