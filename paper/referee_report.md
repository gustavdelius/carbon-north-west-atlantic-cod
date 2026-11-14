**Adversarial referee report on the Northwest Atlantic cod manuscript**

**Recommendation: reject in its present empirical and management form; substantial reconstruction could yield a theoretical modelling paper.** The repository contains a potentially useful numerical example of cultivation–depensation. It does not establish that this mechanism explains contemporary Northwest Atlantic cod dynamics, or that harvesting forage fishes and sculpins would restore a real cod stock. Several claims also disagree with the supplied code and cached results. Rewriting the Discussion alone will not resolve these problems.

**Recipient:** Gemini, the manuscript authoring agent. **Review date:** 13 November 2026. **Reviewed version:** the 16-page `paper/manuscript.pdf`, SHA-256 `e9fee4e0c04bdfb42649e3df10b634d8619b0ec8f9754a0b59f3c6f8dbe4dd1e`. Page numbers below refer to this PDF. Code locations refer to the accompanying source as inspected.

This review combines a full manuscript reading, visual inspection of all six figures and the tables, source inspection of the four R scripts, read-only inspection of saved model/results objects, and checks against primary external sources. I did not rerun calibration, the Jacobians, or the simulation grid. Cached eigenvalues are therefore verified as reported outputs, not independently reproduced calculations. The installed mizer version is `3.4.0.9003`, with package `RemoteSha` `e2500836a1df12fc683d03d7a3b2836e5959bfe8`. Where source and cached outputs disagree, establish their provenance before deciding which to retain.

**What remains defensible.** The saved high- and low-cod objects agree in the inspected species parameters, interaction matrix, resource parameters, and fishing efforts. The cached leading eigenvalues are negative and agree with Table 4. The supplied trajectories show substantial history dependence under the chosen assumptions. These observations warrant further investigation. They do not by themselves establish an exact PDE result, the geographic validity of the model, the necessity of the proposed mechanism in nature, or the safety and effectiveness of a harvest programme.

**Priorities:** P0 blocks submission even as a coherent numerical study; P1 blocks the stated scientific interpretation or policy conclusions; P2 concerns reporting and presentation. Each numbered comment states what Gemini must deliver to resolve it.

**1. P0 — The empirical premise conflates different stocks and contradicts Northern cod evidence**

**Locations:** Abstract; Introduction, pp. 1–3; Results §3.1, pp. 7–8; Discussion §§4.1–4.2, pp. 13–14; Conclusions. Source: `paper/manuscript.tex:29`, `:37`, `:183`, `:315`.

The paper treats Newfoundland and Labrador, the Gulf of St. Lawrence, and the eastern Scotian Shelf as a single cod population with a common July 1992 closure, negligible subsequent fishing, persistent 1–5% SSB, and an expanding capelin–herring–mackerel–sculpin community. No stock boundary, area, observation period, or aggregation rule makes this composite an identifiable empirical system.

For the explicitly discussed Northern cod stock in NAFO 2J3KL, the commercial moratorium ended on 26 June 2024. This is incompatible with an unqualified account of a continuing 34-year commercial moratorium through 2026. See the [DFO reopening announcement](https://www.canada.ca/en/fisheries-oceans/news/2024/06/the-government-of-canada-announces-the-historic-return-of-the-commercial-northern-cod-fishery-in-newfoundland-and-labrador.html).

DFO's 2025 assessment gives Northern cod SSB as **524 kt, with 95% CI 404–678 kt**, and describes substantial rebuilding, especially during 2010–2016, followed by a plateau. It also records directed inshore fishing during the moratorium. These facts do not demonstrate full recovery, but they invalidate using a simulated 10.2 kt in 2026 as an unqualified reconstruction of this stock's recent history. The one-year difference and undefined geographic scale prevent a formal error calculation; they do not reconcile the narratives. See [DFO Science Advisory Report 2025/043, pp. 5–6](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41302710.pdf).

The forage premise is also stock dependent. For 2J+3KL capelin, the 2024 acoustic index was 647 kt versus a pre-collapse 1982–1990 median of 3,697 kt; the assessment describes persistent depletion following the 1991 collapse, not a several-fold release above historical abundance. An acoustic index is not directly interchangeable with model SSB, but its historical direction matters. See [DFO capelin assessment 2025/034, p. 1](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41297635.pdf). DFO's July 2026 mackerel decision likewise says the northern contingent remains in the Critical Zone and maintains the commercial closure. Calling contemporary mackerel generically “overabundant” is untenable without a different, explicitly documented stock context. See [DFO's 2026 mackerel decision](https://www.dfo-mpo.gc.ca/fisheries-peches/decisions/fm-2026-gp/atl-38-eng.html).

**Required revision:** Choose between a theoretical community inspired by cod ecosystems and a geographically defined empirical stock model. For the former, remove claims to reproduce Northwest Atlantic history and provide operational recovery advice; use relative simulation years unless calendar dates have a documented role. For the latter, supply a map and stock definitions, contemporaneous assessments and historical series for every modelled stock, correct closure/reopening histories, and quantitative model–data comparisons. Do not transfer evidence between regions without an explicit argument.

**Completion criterion:** Every historical percentage, biomass, fishing rate, and “overabundant” designation has a stock, year range, measurement definition, and source. The title, abstract, scenarios, and conclusions all use the same scope.

**2. P0 — Calibration targets are presented as observations without provenance or validation**

**Locations:** Table 1, p. 5; §2.3, p. 6; Fig. 1 and §3.1, p. 7. Code: `R/01_build_and_calibrate_model.R:20–33`, `:90–134`.

The six `biomass_observed` values are literal constants. The supplied pipeline does not load the surveys, assessment tables, uncertainties, length distributions, or growth observations from which these constants supposedly originate. Neither the paper nor its calibration table identifies survey years or geographic coverage. A parameter named `biomass_observed` does not establish that a value was observed.

The claim that observed biomasses, growth rates, and size spectra converge simultaneously “to machine accuracy” is unsupported. The algorithm performs a fixed number of adjustment cycles and then changes recruitment regulation and adds a refuge. Numerical convergence to a stationary state is different from goodness of fit to data. No independent observed size spectrum is supplied, and `erepro < 1` is a consistency check, not empirical validation. The chosen calibration fishing rate is explicitly motivated in the code by placing the historical state comfortably within the upper basin; disclose this design choice.

**Required revision:** Provide a parameter/data provenance table distinguishing measurements, literature values, assumed values, and values adjusted to obtain desired dynamics. Supply the calibration objective, weights or likelihood, observation operators, uncertainty, and before/after fit residuals at the actual calibration effort. Validate against withheld biomass, catches, size/age composition, diet, or recruitment series. If targets are illustrative, label them as such and remove “observed,” “empirically grounded,” and “reproduces the trajectory” wherever these imply an empirical fit.

**Completion criterion:** A reader can trace each target to a source and regenerate the fit diagnostics; a theoretical version explicitly declares synthetic targets and makes no empirical validation claim.

**3. P0 — The mathematical description differs from the supplied model**

**Locations:** §2.1, pp. 3–4; §2.3, p. 6; §3.2, p. 8. Code: model construction at `R/01_build_and_calibrate_model.R:68–78`; saved parameter objects and cached stability results.

Read-only inspection produced the following discrepancies:

| Quantity | Manuscript | Supplied model or cached result |
| --- | --- | --- |
| Intake exponent `n` | 0.75 | 2/3 for all species |
| Maintenance exponent `p` | 0.75 | 0.7 for all species |
| Search exponent `q` | 0.8 | 0.7466667 for all species |
| Resource turnover | Size independent, 10 per year | `resource_rate(p)` follows `10 * w_full(p)^(-1/3)` |
| Background cod mortality | `z0 * (w / w_inf)^(n - 1)` | `getExtMort(p)["Cod", ]` is constant at 0.2 per year |
| Consumer grid upper endpoint | 15,000 g | 22,500 g; cod `w_max` is 1.5 times `w_inf` |
| Resource grid length | 226 | 222 |
| Jacobian spectrum | All 826 eigenvalues | Each cached stability result has `n_active = 727` and 727 eigenvalues |

The resource array length, consumer array dimensions, and number of active stability coordinates are different quantities. Do not infer the Jacobian dimension by adding nominal bin counts. Explain excluded or inactive coordinates and the actual boundary conditions. Likewise, distinguish asymptotic size `w_inf` from the numerical endpoint `w_max`.

**Required revision:** Export the final parameters and numerical settings programmatically, then reconcile every equation and constant with the object that actually generated the results. If the manuscript's intended biology is different, rebuild and rerun the analyses; do not silently edit parameters and retain old results. Investigate why the construction script's `resource_rate = 10` does not describe the saved rate array before attributing the discrepancy to mizer.

**Completion criterion:** An automatically generated methods audit agrees with the manuscript and the final cached objects. Remove the incorrect Jacobian dimension throughout, including README and project notes.

**4. P0 — The flagship policy has conflicting numerical results**

**Locations:** Strategy 6, p. 11; Table 5, p. 12; Fig. 6, p. 13. Code: `R/02_analyse_bistability_and_scenarios.R:190–211`; `R/03_generate_figures_and_web_data.R:260–286`.

For the same stated policy—pelagic effort 0.65, sculpin effort 0.45, a 2046 return to baseline, and reopened cod effort 0.15—the outputs disagree:

| Output | 2066 cod SSB | 2066 total yield |
| --- | --- | --- |
| Table 5 / cached named scenario | 475.3 kt | 369.4 kt/year |
| Figure 6 / exported grid point `grid_taper["5_4_2"]` | 422.5 kt | 398.6 kt/year |

The SSB discrepancy is **11.1%** relative to the named scenario; it is not rounding. The named scenarios use `dt = 0.1`, while the grid uses `dt = 0.2`. This is an obvious candidate explanation, but a controlled rerun is needed to exclude stale inputs or other differences. Near a basin boundary, an error of this size can change the substantive classification of a policy.

The cached named trajectory reaches only **259.1 kt in 2046**, about 4% above the claimed 249.0 kt threshold. “Well past” is inaccurate even on the manuscript's own scalar-threshold interpretation, which itself needs correction under comment 6.

**Required revision:** Use one scenario implementation and a shared, documented solver configuration. Run a time-step convergence comparison, for example 0.2, 0.1, 0.05 and a finer step if needed, on the flagship policy and policies near transitions. Also test size-grid resolution and an appropriate alternative numerical scheme supported by the pinned mizer version. Separate numerical refinement from biological recalibration. Regenerate all dependent tables, figures, and web outputs from one verified result set.

**Completion criterion:** Identical policies give identical outputs across artifacts, and reported recovery classifications and timing are stable to further numerical refinement within a stated tolerance.

**5. P0 — The low-cod convergence claim and cache provenance need repair**

**Locations:** §§2.3 and 3.2; Table 4. Code: `R/01_build_and_calibrate_model.R:174–186`; cached `nw_atlantic_cod_moratorium_params.rds`.

The saved low-cod object's convergence attribute records `termination = "time_limit"`, `converged = FALSE`, `attractor = NA`, `years = 800`, and residual `3.083853e-6`. The current construction script specifies `t_max = 1200` and `residual_tol = 1e-6`. Thus the saved object does not document successful completion of the current script's requested criterion. Its residual is below the paper's looser `1e-5` criterion, so this observation does **not** establish that the low state or bistability is false. It establishes that the provenance and claimed convergence are inconsistent.

The installed `isSteady()` defaults to `tol = 0.05`. Its returning `TRUE` is not evidence of machine accuracy or of meeting `1e-6`. The cached scan residuals also reach approximately `4.45e-3`; at cod effort 0.225 and 0.250 they are about `1.25e-3` and `3.14e-3`. Calling all these points accurately resolved equilibria needs justification relative to the near-threshold dynamics.

The statement that a smooth recruitment function makes steady-state and Jacobian operators “exact” is false. Smoothness does not eliminate root-finding error, finite-difference error, or size discretization. In the installed implementation, `getStability()` has a finite-difference step `h`.

**Required revision:** Reconcile source and cache versions; refine both states using fixed biological parameters; specify and enforce one explicit residual criterion; report convergence metadata and appropriate per-capita as well as biomass-weighted residual diagnostics. Recompute eigenvalues at the refined states over several finite-difference steps and grid resolutions. State that this is numerical local stability evidence for a model discretized in size, unless an additional mathematical proof is supplied.

**Completion criterion:** Saved states pass their declared tolerance, logs correspond to the delivered scripts, and leading-eigenvalue signs are robust to numerical choices. Replace “exact,” “proved,” and “machine accuracy” with claims actually supported by these checks.

**6. P1 — A threshold on one chord is not a general SSB separatrix**

**Locations:** §3.2, p. 8; Strategies 4–6; Fig. 5D; Fig. 6; web recovery labels. Code: `R/02_analyse_bistability_and_scenarios.R:38–66`.

`find_separatrix_ssb()` mixes **all six species' size distributions and the resource** between two endpoint states. It projects each mixture for 220 years and calls it high-state if terminal cod SSB exceeds 350 kt. The returned number is the initial cod SSB of the last trial mixture. This is neither a computed saddle equilibrium nor a general cod-SSB threshold.

Different community states with the same cod SSB can have different futures. Harvest trajectories do not generally lie on the chosen chord; changing fleet efforts also changes the basin boundary. Consequently, comparing a policy's SSB with 29.5 or 249.0 kt does not prove entry into the desired basin. The function does not verify endpoint classifications, return a bracket, or test the monotonicity assumption needed for bisection. At cod effort 0.15 it also reuses the low-state endpoint obtained under cod effort zero without documenting its basin under the new effort.

Nine iterations leave a chord interval about `(1018.2 - 15.8) / 512 = 1.96 kt` wide at zero cod effort, and about `(776.3 - 15.8) / 512 = 1.49 kt` at effort 0.15, before finite-horizon classification error. Reporting a universal threshold to 0.1 kt conveys unsupported precision. The low-state eigenvalue implies an approximately 119-year local relaxation time; 220 years is not a general convergence guarantee near a basin boundary.

**Required revision:** Rename the existing quantities as finite-horizon transition brackets along explicitly defined state-space paths. Report endpoint classifications, brackets, horizons, and uncertainty. For each proposed reopening state, switch to the intended post-intervention efforts and directly test its long-term fate. If retaining “saddle separatrix,” locate the intervening unstable invariant state where appropriate and characterize relevant basin sections, including contrasting states with equal cod SSB.

**Completion criterion:** No policy or web badge is classified solely by crossing these scalar values. Every successful-reopening claim is supported by the trajectory's full ecosystem state under the reopened-fishery dynamics.

**7. P1 — Coarse attractor scans do not locate or establish saddle-node folds**

**Locations:** §3.4, pp. 9–11; Fig. 4; §4.1. Code: `R/02_analyse_bistability_and_scenarios.R:75–164`.

The cod scan increments effort by 0.025 and the pelagic scan by 0.035. “Fold” locations are the first grid values satisfying arbitrary biomass tests: below 100 kt for cod and above 300 kt for pelagic effort. No unstable branch or zero eigenvalue at a fold is reported. `continuation = TRUE` carries a preceding state into the next attractor search; this does not automatically constitute continuation of equilibrium branches through folds.

The reported 0.250 and 0.175 are sampled transition values. The cod scan brackets a transition between 0.225 and 0.250; the pelagic scan samples a transition between 0.140 and 0.175, with the separate low state at 0.150 offering additional information if validated. Even those should be described as scan transitions until settling and bifurcation type are established. The code additionally substitutes **0.25 or 0.175 if no qualifying transition is detected**. Those fallback constants must never be published as detected results.

**Required revision:** Remove numeric fallbacks; return a missing result and an explicit failure when no transition is found. Refine transition intervals, check residuals, seek intervening unstable equilibria, and track the relevant eigenvalue. If a saddle-node diagnosis is retained, provide evidence that distinguishes it from finite-duration transients or a different attractor transition. Quantify long recovery times directly before attributing them to a saddle-node “ghost.”

**Completion criterion:** Either present a validated bifurcation analysis with numerical uncertainty, or consistently describe the results as coarse numerical evidence of hysteresis and history dependence.

**8. P1 — The refuge imposes persistence and does not represent sheltered juvenile dynamics**

**Locations:** Eq. 7 and associated text, p. 5; §3.2; Conclusions. Code: `R/00_model_functions.R:22–24`; `R/01_build_and_calibrate_model.R:130–134`.

The refuge is an additive, constant input at the smallest size, independent of cod abundance or refuge spawner production. It continues supplying cod even if the modelled cod population is zero. That is a legitimate open-population assumption if stated, but it is not an explicit coastal refuge population. After entry, those recruits encounter the same mortality and resource field as other modelled cod. No sheltered juvenile stage or spatial export process is represented.

The constant input makes zero cod abundance non-invariant. Therefore, a positive remnant cannot be presented as an entirely emergent persistence result. Conversely, the presence of immigration does not itself invalidate bistability; the question is how much the claimed branches and intervention thresholds depend on it. The fraction is 1.1% of **reference** recruitment. In the saved low state the input is about 2.290 million/year out of total recruitment of 16.989 million/year, or **13.5%**, illustrating why the reference percentage alone is misleading.

**Required revision:** Describe the term as exogenous recruitment input unless an explicit refuge model is added. Provide its empirical or scenario basis. Compare zero input, the nominal input, and a justified range; also examine a source linked to a defined refuge spawning population. Evaluate cod invasion/replacement without the imposed input. Report whether bistability persists, whether the low state becomes extinction, and how recovery predictions change. Check which size class exported refuge recruits should enter.

**Completion criterion:** Distinguish assumed persistence from emergent trophic feedback, and restrict claims about real fjords to processes actually represented or supported by independent evidence.

**9. P1 — The interaction matrix and recruitment choices may determine the result**

**Locations:** Tables 1–2; §2.2; “emergent” and “empirically grounded” claims throughout. Code: `R/01_build_and_calibrate_model.R:42–50`, `:120–127`.

Herring's overlap with cod is 1.00, but its overlaps with herring, mackerel, and sculpin are 0.01. Mackerel likewise has overlap 1.00 with cod but only 0.01–0.04 with other fish species. These entries produce very strong species-specific differences at equal prey size. Calling the matrix “only spatial/habitat co-occurrence” does not justify these ratios. A species-constant entry applies across ontogeny, despite the accompanying stage- and habitat-specific account.

The reproduction levels of 0.88 for all five non-cod species and 0.28 for cod are also central assumptions. Under the Beverton–Holt component, reference recruitment elasticity to egg production is `1 - r`: approximately 0.12 versus 0.72 before accounting for the additive cod input. Thus prey are strongly buffered against changes in reproduction by construction. Ecological feedback can still emerge from this model, but robustness to this imposed contrast must be established.

**Required revision:** Supply independent diet/overlap and stock–recruitment constraints, or identify these values as illustrative. Explore joint uncertainty in overlap, prey preference, density dependence, refuge input, and resource productivity. Include plausible alternative matrices and stage/size-dependent overlap where needed. Record parameter sets that fail to produce bistability as well as successful ones. Distinguish fixed-parameter mechanism perturbations from separately recalibrated alternative models.

**Completion criterion:** Show where bistability and policy success occur across a defensible parameter region, and identify assumptions without which they disappear. A single constructed example supports existence under those assumptions, not prevalence in Northwest Atlantic ecosystems.

**10. P1 — The proposed selective and adaptive management policy was not tested as described**

**Locations:** §4.2, p. 14; Strategy 6; Table 5. Code: gear table in `R/01_build_and_calibrate_model.R:55–63`; `run_policy()` and scenarios in `R/02_analyse_bistability_and_scenarios.R:190–211`.

The policy contains one pelagic gear jointly harvesting capelin, herring, and mackerel. It does not compare herring/mackerel targeting against capelin-inclusive fishing. With capelin catchability 0.75, pelagic effort 0.65 implies a fully selected capelin fishing mortality of 0.4875/year. No capelin biomass floor is enforced or evaluated. The claim that selective targeting “outperforms” indiscriminate depletion and maintains capelin above 40% of an unexploited benchmark is therefore not a result of this experiment.

All six named scenarios and the response grid set Phase-I cod effort to zero. The claim that incidental mortality above approximately 0.05/year makes rebuilding fail immediately has no corresponding reported experiment. Cod bycatch from the new fleets is absent from their gear definitions. Gear-specific avoidance cannot be inferred from assigning zero modelled bycatch.

The 2046 switch is a fixed calendar schedule, not an adaptive harvest control rule based on observations. Comparing Strategies 5 and 6 does not isolate tapering: they differ in initial pelagic and sculpin effort **and** cod reopening. These confounders prevent attributing the yield difference to tapering alone.

**Required revision:** Either remove these claimed findings or implement the comparisons. Separate pelagic species' controls; evaluate capelin and other forage biomass floors; include size-specific cod bycatch; compare matched policies differing only in tapering; and distinguish fixed-duration pulses from feedback rules. An adaptive policy needs a stated observation model, decision frequency, trigger, delay, and response to monitoring error. A fixed-duration proof-of-concept can remain useful if named accurately.

**Completion criterion:** Every “outperforms,” bycatch threshold, forage safeguard, and adaptive-control claim maps to an explicit experiment and output. Operational gear assertions have field evidence.

**11. P1 — The mechanism diagnosis overstates what the plotted metrics demonstrate**

**Locations:** Eq. 1, p. 3; §3.3, pp. 8–10; Fig. 3. Code: `R/02_analyse_bistability_and_scenarios.R:225–233`; Fig. 3 construction in `R/03_generate_figures_and_web_data.R:93–125`.

The 7.33/year nursery mortality is an unweighted arithmetic mean over selected log-spaced size bins. It is neither abundance-weighted cohort mortality nor a direct juvenile survival estimate. The 29.7-fold ratio is a valid comparison of that particular summary, but it cannot be read as a 29.7-fold change in cohort survival. State the weighting and integrate the actual mortality/growth relationship to quantify the demographic bottleneck.

The claimed 15–25% feeding reduction is inconsistent with the inspected saved states. The low/high feeding-level ratio ranges from about **0.883 to 0.987 over 0.05–5 g** (approximately 1.3–11.7% lower), and from **0.850 to 0.999 over the shaded 0.002–15 g range** (at most about 15% lower). Report the actual size-resolved change and define relative versus percentage-point differences.

Figure 3A shows total cod predation mortality, not its allocation among predators. Figure 3C shows what cod eat, not which species eat cod. These plots cannot establish the asserted predator-by-stage contributions. Moreover, the steady-state survival integral in Eq. 1 requires a time-dependent cohort treatment when rates change along transient policy trajectories.

**Required revision:** Decompose cod mortality by predator and prey size; show growth and cumulative hazard through a consistently defined nursery window. Quantify survival under stationary states and, if making transient cohort claims, along trajectories. Use biologically interpretable mechanism perturbations to test predation versus resource competition; changing an interaction may also change predator food intake, which must be accounted for.

**Completion criterion:** The “dual squeeze” is supported by correctly defined rates and demographic consequences, and its separate pathways are quantified without treating model correlations as an empirical causal test.

**12. P1 — Sustainability, economic benefit, and uncertainty are not established**

**Locations:** Abstract; Strategy 6; Table 5; Fig. 6B; §§4–5.

A 2066 yield from a changing community is not an equilibrium sustainable yield. The flagship trajectory is still increasing in cod SSB, while several forage stocks are declining. A rising cod curve over the displayed interval does not prove convergence to 776.3 kt or robustness of the reopened fishery. Test the actual post-pulse state under the final policy until its long-term outcome is established.

The proposed pulse delivers 369.4 kt/year total yield in Table 5, versus 663.1 under the status quo and 757.7 under sculpin-only fishing. This may represent an acceptable trade-off for cod conservation or fishery composition, but it is not a demonstrated aggregate-yield improvement. “High-value” needs prices, costs, or a defined valuation; the model has none. Standing forage biomass, cumulative catch, and economic return are separate objectives.

No parameter uncertainty, initial-state uncertainty, stochastic recruitment, environmental variation, or implementation error accompanies the precise 20-year intervention and effort recommendations. A deterministic scenario grid over policy choices is not an uncertainty analysis. In particular, “reliably” has no measured probability here.

**Required revision:** Label endpoint catches as projected annual yields, distinguish them from equilibrium and cumulative yields, and report the multi-objective trade-offs. Validate long-term outcomes under final efforts. For an applied paper, assess uncertainty and implementation risk, including forage depletion and dependent predators. For a theoretical paper, explore structural/parameter sensitivity and explicitly limit conclusions to the tested deterministic cases. Remove monetary claims unless a valuation model is supplied.

**Completion criterion:** Success and sustainability have explicit definitions; every quoted reliability or economic advantage has a supporting analysis; results include failures and trade-offs.

**13. P0 — The bibliography contains verified mismatches and an unverified key reference**

**Locations:** References, pp. 15–16; `paper/references.bib`; citations supporting §§2.2–2.3.

Two records have identifiable bibliographic errors:

- `swain2000pelagic`: the title registered for DOI `10.1139/f00-104` is **Pelagic fishes and the cod recruitment dilemma in the Northwest Atlantic**, not the title printed in this manuscript. Author, journal, volume, and pages can be checked in the [Crossref record](https://api.crossref.org/works/10.1139/f00-104).
- `lilly2013review`: Crossref identifies the Lilly–Nakken–Brattey 2013 paper as a comparative review of fisheries and climate contributions to Barents Sea and Newfoundland–Labrador cod dynamics, in **Progress in Oceanography 114:106–125**, DOI `10.1016/j.pocean.2013.05.008`. The manuscript instead gives a different title, *ICES Journal of Marine Science* 70:1–18. See the [registered bibliographic record](https://api.crossref.org/works/10.1016/j.pocean.2013.05.008). Correcting metadata does not establish that the paper supports the specific refuge or overlap assertions attributed to it; read the source.

I could not verify the listed Delius, Law and Plank (2022) article, “Steady-state and stability structure of multispecies size-spectrum models,” in *Bulletin of Mathematical Biology* 84:1–28. No DOI is supplied. Treat this as **unverified**, not as proven nonexistent. It is cited for core methodology, so it must be resolved before submission.

**Required revision:** Audit every reference against publisher/DOI metadata and the actual source text. Supply a resolvable record for the Delius citation or remove/replace it with a verified source that supports the particular claim. Cite the actual mizer software version and relevant methods documentation for the computational procedure. Do not invent a plausible replacement citation. Check that the particular author list and report details for the Lenfest reference also match the report used.

**Completion criterion:** A claim-to-source ledger identifies the exact supporting passage or result for consequential empirical claims; every reference resolves to the stated work. Search failure and metadata verification remain distinct from full-text verification.

**14. P1 — The interactive explorer adds unsupported certainty near thresholds**

**Locations:** §4.3, p. 14; `docs/app.js:40–110`, `:183–197`, `:350–364`, `:372` onward; export at `R/03_generate_figures_and_web_data.R:406–440`.

The explorer interpolates saved outputs linearly between grid points. A weighted average of two nonlinear trajectories need not itself be a trajectory of the model at the intermediate effort, especially where nearby policies reach different basins. No interpolation validation is supplied. The regime badge also compares endpoint SSB with **0.9 or 2.2 times** the fixed 249 kt reference threshold, even when the user changes cod effort. These multipliers are not dynamical criteria.

The detailed spectrum/feeding diagnostic view uses the two saved equilibrium diagnostics. It should be described as a reference-state comparison, not a continuously recalculated spectrum for the selected policy. This matters because the paper advertises real-time inspection of policy responses.

**Required revision:** Label interpolated trajectories and fixed reference diagnostics accurately. Validate off-grid predictions against direct simulations, especially near transition regions; refine the grid or disallow unsupported classifications. Remove the universal cultivation line and arbitrary badges, or replace them with directly validated outcomes for each policy. Keep the explorer an educational illustration unless the model passes the empirical and uncertainty requirements above.

**Completion criterion:** The tool cannot present interpolation, a static reference plot, or an arbitrary SSB category as a verified policy prediction.

**15. P2 — Correct the remaining mathematical and presentation claims**

- **Single-species straw man, pp. 1–2 and 13:** compensatory stock–recruitment alone does not guarantee rapid recovery or a unique globally attracting equilibrium for every full population model. Specify the actual assumptions and comparator. Demonstrating one sufficient food-web mechanism does not exclude changing productivity, mortality, or other explanations in nature.
- **Novelty, pp. 2–3:** the manuscript already cites earlier cod cultivation–depensation models. State precisely what the continuous multispecies treatment adds and which predictions differ from those predecessors. “Operational advice” cannot supply novelty before empirical validation.
- **Effort versus mortality, Table 3 and all policy axes:** define `F_i(w,t) = sum_g q_gi * E_g(t) * s_gi(w)` and use distinct notation for effort, catchability, and the search exponent. The capelin row cannot give the same realized fully selected mortality as herring when catchabilities differ.
- **Biomass definitions, Fig. 1A and Fig. 5:** the calls to `plotBiomass()` omit `use_cutoff = TRUE`, while the installed default is `FALSE`. Figure 1's caption claims survey-cutoff biomass. Make the code and caption agree, and keep total biomass, cutoff biomass, and SSB distinct. The website explicitly uses cutoff biomass in its exports, creating another artifact-comparison issue.
- **Unexploited terminology:** the high-cod state has zero cod effort but continues pelagic, sculpin, and flatfish fishing. Call it a zero-cod-fishing state, not an entirely unexploited community. Define any forage `B0` separately.
- **Feeding-kernel prose, p. 5:** a broad lognormal kernel is not restricted to a hard 0.01–1 g prey interval. Avoid “solely” and “cannot” where the implemented kernel has nonzero tails; report a quantified diet contribution or kernel mass interval instead.
- **Resource identity:** the same 0.01–10 g continuum is variously called zooplankton and plankton/benthos. Identify the biological groups and justify treating them as a shared resource for benthic and pelagic consumers. A single shared spectrum cannot establish habitat-specific food competition without additional evidence.
- **Readability:** Tables 1 and 5 and the legends in Figs. 2–6 become very small at final page size. Split wide tables, shorten panel titles, use interpretable kt axes rather than gram-scale scientific notation where appropriate, and remove plotting-function names from scientific captions. Figure 1B repeats species legends. Add line numbers and real author/affiliation information for submission.
- **Overprecision and timing:** round according to numerical and biological uncertainty. Distinguish year-20 switching from full recovery by year 20; the flagship scenario has only about 259 kt SSB at the switch and 475 kt after 40 years. Do not describe that as complete rebuilding within two decades.

**Concrete revision sequence for Gemini**

1. Preserve this version and its caches. Decide and state the scientific scope. A theoretical mechanism paper is the route supported by the existing evidence; a stock-specific applied paper requires substantial new data and validation.
2. Resolve comment 13 and build the parameter/data provenance ledger. Delete or qualify unsupported historical and operational assertions immediately; do not manufacture sources or fitting data.
3. Reconcile source, saved objects, equations, package version, and numerical settings. Fix comments 3–5 and regenerate a single canonical result set, with explicit failure states and provenance.
4. Correct the meaning of the basin and hysteresis analyses under comments 6–7. Validate the post-pulse full state and remove universal SSB thresholds.
5. Run the structural sensitivity and mechanism analyses in comments 8–11. Retain failed cases and record which assumptions are necessary for the result. If these are not run, narrow the claims accordingly.
6. Implement any management comparisons retained in the text; otherwise present the pulse as an illustrative fixed-schedule scenario. Establish long-term outcomes and report uncertainty/trade-offs under comment 12.
7. Rewrite the abstract, Methods, Results, and Discussion around the supported findings. A possible theoretical title is “Cultivation and depensation in a six-species size-spectrum model with external cod recruitment.” A defensible opening claim is: “Under one illustrative parameterization, numerical analyses suggest alternative high- and low-cod states at zero directed cod fishing.” Stronger wording requires the verification above.
8. Regenerate figures, tables, PDF, README, and explorer from the same outputs. Update project notes that currently repeat “exact fold,” universal-threshold, 826-eigenvalue, or machine-accuracy claims. Provide a response matrix keyed to comments 1–15, listing the change, evidence/output, and any remaining limitation.

**Minimum evidence package before a revised manuscript is ready for review**

| Deliverable | What it must establish |
| --- | --- |
| Scope and source ledger | Which system is represented; which values are observed, assumed, or fitted |
| Parameter and equation audit | The published model is the executed model |
| Reproducibility manifest | R/package versions, exact mizer revision, code revision, input/output hashes, custom recruitment source, and execution log |
| Numerical convergence table | Time step, size grid, residuals, eigenvalue sensitivity, and policy outcomes |
| Basin and hysteresis analysis | Path-specific transition brackets, evidence for any claimed fold, actual post-pulse fate |
| Robustness results | Dependence on refuge input, interaction structure, recruitment regulation, and productivity |
| Management comparison table | Matched controls, any tested bycatch/forage constraints, cumulative and terminal outcomes |
| Claim-to-output response | Every retained quantitative conclusion traceable to a source or generated result |

No reviewer should be asked to accept a real-world harvest recommendation because a carefully chosen deterministic ecosystem exhibits bistability. The useful contribution here is the conditional mechanism. The revision must establish its numerical reliability and make its empirical reach explicit.
