# Project Notes — `analyse-stability`

- **Jacobian verification (`R/02_analyse_bistability_and_scenarios.R`):** All three equilibria of the Northwest Atlantic model are strictly locally asymptotically stable continuous-time steady states (`826 x 826` Jacobian from `getStability()`):
  - Pre-collapse fishery (`F_cod = 0.15 /yr`): `max Re(lambda) = -0.3036 /yr`.
  - Unexploited cultivated high-cod state (`F_cod = 0.00 /yr`): `max Re(lambda) = -0.1809 /yr`.
  - Moratorium low-cod trap (`F_cod = 0.00 /yr`): `max Re(lambda) = -0.0084 /yr` (slow recovery mode reflecting proximity of baseline `F_pelagic = 0.15 /yr` to the `0.175 /yr` escape fold).
- **Continuation scans (`scanModel`):** Forward/backward scans with `continuation = TRUE` locate the directed cod collapse cliff at `F_cod = 0.250 /yr` and the `F_cod = 0` pelagic escape fold at `F_pelagic = 0.175 /yr`. Note that `as.data.frame(scanModel(...))` names its columns dynamically after the scanned parameter and `value_func` (e.g. `"Effort of Cod_Trawl"`, `"Spawning stock biomass"`), so scripts index columns positionally (`d[[1]]`, `d[[2]]`).
- **Long relaxation near the collapsed attractor:** Because `max Re(lambda) = -0.0084 /yr` corresponds to a ~120-year asymptotic e-folding time, `findSteadyState()` on `params_low_f0` uses `t_max = 1200` to reach machine-tight convergence (`residual = 3.08e-6`).

## Referee audit correction — 13 November 2026

The earlier notes above repeat manuscript claims and must not be treated as independent verification. See `paper/referee_report.md`, comments 3–7, for the read-only audit of the supplied artifacts. Cached stability results contain 727 eigenvalues and `n_active = 727`, not 826. The saved low-cod object's convergence metadata is `time_limit`, `converged = FALSE`, `attractor = NA`, `years = 800`, residual `3.083853e-6`; this is below the manuscript's `1e-5` criterion but above the current script's `1e-6` criterion. It is not machine accuracy. The reported 0.250/0.175 values are coarse scan transition points, not independently located folds. The 29.5/249.0 kt numbers are finite-horizon classifications along state-space chords, not general SSB thresholds. The audit did not recompute eigenvalues or disprove bistability; source/cache reconciliation and numerical refinement remain necessary.
