# ==============================================================================
# R/01_build_and_calibrate_model.R
#
# Builds, calibrates, and refines a theoretical six-species boreal-shelf
# size-spectrum model in mizer illustrating cultivation-depensation bistability
# under external cod recruitment input.
#
# Scope note:
#   All calibration targets (`biomass_observed` above `biomass_cutoff`, required
#   by `mizer::matchBiomasses()`) are synthetic, illustrative reference values
#   inspired by boreal shelf functional guilds rather than empirical survey
#   series from a single stock boundary. Simulation times are relative years
#   (t in [0, 61] for collapse + closure; tau in [0, 40] for post-closure
#   interventions).
# ==============================================================================

library(mizer)
source("R/00_model_functions.R")

dir.create("models", showWarnings = FALSE, recursive = TRUE)
dir.create("models/logs", showWarnings = FALSE, recursive = TRUE)
sink("models/logs/01_build_and_calibrate_model.log", split = TRUE)

# Remove legacy v1 cache filenames from models/ (preserved in archive/v1_reviewed/)
legacy_files <- c(
  "models/nw_atlantic_cod_params.rds",
  "models/nw_atlantic_cod_high_f0_params.rds",
  "models/nw_atlantic_cod_moratorium_params.rds"
)
for (lf in legacy_files) {
  if (file.exists(lf)) unlink(lf)
}

# ------------------------------------------------------------------------------
# 1. Species parameters: theoretical six-species boreal-shelf fish community
#    `biomass_observed` stores synthetic reference cutoff-biomass targets (g;
#    1e9 g = 1 kilotonne, kt) at calibration effort EFFORT_CAL.
# ------------------------------------------------------------------------------
sp <- data.frame(
  species          = c("Capelin", "Herring", "Mackerel", "Sculpin", "Plaice", "Cod"),
  w_inf            = c(50,        380,       850,        1100,      2600,     15000),
  w_mat            = c(7,         55,        115,        95,        180,      850),
  k_vb             = c(0.75,      0.58,      0.52,       0.42,      0.28,     0.26),
  z0               = c(0.38,      0.35,      0.34,       0.54,      0.28,     0.20),
  beta             = c(2000,      2500,      350,        65,        220,      60),
  sigma            = c(1.8,       2.2,       1.9,        1.6,       1.7,      1.6),
  a                = c(0.0060,    0.0070,    0.0080,     0.0120,    0.0090,   0.0085),
  b                = c(3.00,      3.00,      3.00,       3.00,      3.00,     3.00),
  biomass_observed = c(5.2e11,    3.8e11,    1.9e11,     1.4e11,    2.6e11,   6.8e11),
  biomass_cutoff   = c(2,         10,        15,         15,        20,       50),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. Illustrative asymmetric trophic interaction matrix (rows = predators, cols = prey)
#    Constructed to represent strong stage-reversed predation between planktivores /
#    benthic sculpins and early-life cod vs. low intra-forage larval predation.
#    Sensitivity to relaxing this asymmetry is quantified in R/04_sensitivity_and_audit.R.
# ------------------------------------------------------------------------------
inter <- matrix(0.02, nrow = 6, ncol = 6,
                dimnames = list(predator = sp$species, prey = sp$species))
# Prey order:          Capelin Herring Mackerel Sculpin Plaice  Cod
inter["Capelin", ]  <- c(0.01,   0.01,   0.01,    0.01,   0.01, 0.55)
inter["Herring", ]  <- c(0.02,   0.01,   0.01,    0.01,   0.02, 1.00)
inter["Mackerel", ] <- c(0.04,   0.02,   0.01,    0.01,   0.02, 1.00)
inter["Sculpin", ]  <- c(0.05,   0.02,   0.01,    0.02,   0.03, 0.45)
inter["Plaice", ]   <- c(0.04,   0.02,   0.01,    0.01,   0.01, 0.02)
inter["Cod", ]      <- c(0.45,   0.43,   0.40,    0.38,   0.12, 0.025)

# ------------------------------------------------------------------------------
# 3. Baseline 4-fleet gear selectivity & catchability parameters
#    Note: F_{g,i}(w) = S_{g,i}(w) * Q_{g,i} * E_g(t). Capelin catchability is
#    Q = 0.75, whereas Herring/Mackerel/Sculpin/Plaice/Cod have Q = 1.00.
# ------------------------------------------------------------------------------
gp <- data.frame(
  gear         = c("Pelagic", "Pelagic", "Pelagic", "Mesopredator", "Flatfish", "Cod_Trawl"),
  species      = c("Capelin", "Herring", "Mackerel", "Sculpin",     "Plaice",   "Cod"),
  sel_func     = "sigmoid_length",
  l50          = c(11,        20,        25,         20,            28,         52),
  l25          = c(9,         16,        20,         16,            23,         45),
  catchability = c(0.75,      1.00,      1.00,       1.00,          1.00,       1.00),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 4. Construct and calibrate the multi-species size-spectrum model
#    Note: scalar `resource_rate = 10` sets r_pp = 10 in r_R(w) = r_pp * w^(n_R - 1)
#    with n_R = 2/3, giving r_R(w) = 10 * w^(-1/3) yr^(-1).
# ------------------------------------------------------------------------------
params <- newMultispeciesParams(
  sp,
  interaction   = inter,
  gear_params   = gp,
  no_w          = 100,
  kappa         = 2.8e11,
  lambda        = 2.08,
  w_pp_cutoff   = 10,
  resource_rate = 10,
  info_level    = 0
)

params <- setColours(params, c(
  "Capelin"  = "#e69f00",
  "Herring"  = "#d95f02",
  "Mackerel" = "#b25806",
  "Sculpin"  = "#7c3aed",
  "Plaice"   = "#15805c",
  "Cod"      = "#1d5fa7",
  "Resource" = "#64748b"
))

# Calibration design choice: calibrate at E_cod = 0.22 /yr (EFFORT_CAL) so that
# the sustainable reference fishery (E_cod = 0.15 /yr, EFFORT_HIST) lies inside
# the upper cultivated basin away from the upper-branch termination (~0.235 /yr).
initial_effort(params) <- EFFORT_CAL

bio_init_uncal <- unclass(getBiomass(params, use_cutoff = TRUE)) / 1e9

params <- matchBiomasses(params, info_level = 0)
for (k in 1:10) {
  params <- steadySingleSpecies(params, keep = "biomass")
}
params <- setBevertonHolt(params, reproduction_level = 0)
params <- tuneSteadyState(params, t_max = 150, preserve = "reproduction_level",
                          progress_bar = FALSE, info_level = 0)

for (iter in 1:5) {
  params <- matchGrowth(params, info_level = 0)
  params <- matchBiomasses(params, info_level = 0)
  for (k in 1:8) {
    params <- steadySingleSpecies(params, keep = "biomass")
  }
  params <- setBevertonHolt(params, reproduction_level = 0)
  params <- tuneSteadyState(params, t_max = 150, preserve = "reproduction_level",
                            progress_bar = FALSE, info_level = 0)
}

# Impose Beverton-Holt recruitment density dependence (reproduction_level)
reproduction_level(params) <- c(
  Capelin  = 0.88,
  Herring  = 0.88,
  Mackerel = 0.88,
  Sculpin  = 0.88,
  Plaice   = 0.88,
  Cod      = 0.28
)

# Record pre-subsidy fit diagnostics at EFFORT_CAL
bio_cal_pre_sub <- unclass(getBiomass(params, use_cutoff = TRUE)) / 1e9
res_cal_pre_sub <- max(abs(getSteadyResidual(params, effort = EFFORT_CAL, measure = "biomass")))

# Register exogenous cod recruitment input (1.1% of calibrated Cod recruitment
# at EFFORT_CAL; sensitivity across 0% to 5% tested in R/04_sensitivity_and_audit.R)
rdd_cal <- getRDD(params)
species_params(params)$r_refuge <- c(0, 0, 0, 0, 0, 0.011 * unname(rdd_cal["Cod"]))
params <- setRateFunction(params, "RDD", "RefugeBevertonHoltRDD")

# Verify physiological consistency of reproductive efficiencies (all erepro << 1)
stopifnot(all(species_params(params)$erepro < 1))

# Also settle at EFFORT_CAL after adding r_refuge to document exact post-subsidy
# target deviation at calibration effort
params_cal_settled <- findSteadyState(
  params, effort = EFFORT_CAL, solver = "project", t_max = 3000,
  residual_tol = 1e-8, distance_tol = 1e-10,
  progress_bar = FALSE, info_level = 0
)
bio_cal_post_sub <- unclass(getBiomass(params_cal_settled, use_cutoff = TRUE)) / 1e9
res_cal_post_sub <- max(abs(getSteadyResidual(params_cal_settled, effort = EFFORT_CAL, measure = "biomass")))

# ------------------------------------------------------------------------------
# 5. Locate and refine the four reference steady states WITH IDENTICAL BIOLOGICAL
#    PARAMETERS to strict residual_tol = 1e-8 (fixed_point verified):
#    (1) params_f15     : Baseline sustainable cod fishery (E_cod = 0.15 /yr, high branch)
#    (2) params_high_f0 : Zero-cod-fishing cultivated equilibrium (E_cod = 0.00 /yr, high branch)
#    (3) params_low_f0  : Zero-cod-fishing depensatory equilibrium (E_cod = 0.00 /yr, low branch)
#    (4) params_low_f15 : Reopened cod fishery low equilibrium (E_cod = 0.15 /yr, low branch)
# ------------------------------------------------------------------------------
params_f15 <- findSteadyState(
  params, effort = EFFORT_HIST, solver = "project", t_max = 3500,
  residual_tol = 1e-8, distance_tol = 1e-10,
  progress_bar = FALSE, info_level = 0
)

params_high_f0 <- findSteadyState(
  params_f15, effort = EFFORT_MOR, solver = "project", t_max = 3500,
  residual_tol = 1e-8, distance_tol = 1e-10,
  progress_bar = FALSE, info_level = 0
)

# Canonical collapse (t = 0..26) + 34-year closure (t = 27..61) trajectory
# using CANONICAL_METHOD ("tr_bdf2") and CANONICAL_DT (0.05)
eff_collapse_closure <- build_collapse_closure_effort(names(EFFORT_HIST))
sim_collapse_closure <- project(
  params_f15,
  effort       = eff_collapse_closure,
  dt           = CANONICAL_DT,
  t_save       = 1,
  method       = CANONICAL_METHOD,
  progress_bar = FALSE
)
params_t61 <- finalParams(sim_collapse_closure)

params_low_f0 <- findSteadyState(
  params_t61, effort = EFFORT_MOR, solver = "project", t_max = 4500,
  residual_tol = 1e-8, distance_tol = 1e-10,
  progress_bar = FALSE, info_level = 0
)

params_low_f15 <- findSteadyState(
  params_low_f0, effort = EFFORT_HIST, solver = "project", t_max = 4500,
  residual_tol = 1e-8, distance_tol = 1e-10,
  progress_bar = FALSE, info_level = 0
)

# Enforce strict convergence & parameter-identity assertions
verify_fixed_point <- function(p_obj, eff_vec, label) {
  cv <- attr(p_obj, "convergence")
  rb <- max(abs(getSteadyResidual(p_obj, effort = eff_vec, measure = "biomass")))
  rpc_mat <- getSteadyResidual(p_obj, effort = eff_vec, measure = "per_capita")
  bio_bin <- sweep(unclass(initialN(p_obj)), 2, w(p_obj) * dw(p_obj), "*")
  rpc_active <- max(abs(rpc_mat[bio_bin > 1e-3]), na.rm = TRUE)
  rpc_all    <- max(abs(rpc_mat), na.rm = TRUE)

  if (!isTRUE(cv$converged) || !identical(cv$attractor, "fixed_point") || rb >= 1e-8) {
    stop(sprintf("Equilibrium %s failed strict convergence: conv=%s, attr=%s, res_bio=%.3e",
                 label, cv$converged, as.character(cv$attractor), rb))
  }
  data.frame(
    state            = label,
    converged        = cv$converged,
    termination      = cv$termination,
    attractor        = as.character(cv$attractor),
    years_run        = cv$years,
    residual_conv    = cv$residual,
    residual_biomass = rb,
    residual_pc_main = rpc_active,
    residual_pc_all  = rpc_all,
    cod_ssb_kt       = unname(getSSB(p_obj)["Cod"] / 1e9),
    cod_cutoff_kt    = unname(getBiomass(p_obj, use_cutoff = TRUE)["Cod"] / 1e9),
    cod_total_kt     = unname(getBiomass(p_obj, use_cutoff = FALSE)["Cod"] / 1e9),
    stringsAsFactors = FALSE
  )
}

conv_table <- rbind(
  verify_fixed_point(params_f15,     EFFORT_HIST, "High_Branch_Ecod_0.15"),
  verify_fixed_point(params_high_f0, EFFORT_MOR,  "High_Branch_Ecod_0.00"),
  verify_fixed_point(params_low_f0,  EFFORT_MOR,  "Low_Branch_Ecod_0.00"),
  verify_fixed_point(params_low_f15, EFFORT_HIST, "Low_Branch_Ecod_0.15")
)

stopifnot(isTRUE(all.equal(species_params(params_high_f0), species_params(params_low_f0))))
stopifnot(isTRUE(all.equal(gear_params(params_high_f0), gear_params(params_low_f0))))
stopifnot(isTRUE(all.equal(interaction_matrix(params_high_f0), interaction_matrix(params_low_f0))))
stopifnot(isTRUE(all.equal(resource_params(params_high_f0), resource_params(params_low_f0))))
stopifnot(identical(getRateFunction(params_high_f0), getRateFunction(params_low_f0)))
stopifnot(getSSB(params_low_f0)["Cod"] < 0.05 * getSSB(params_high_f0)["Cod"])

# ------------------------------------------------------------------------------
# 6. Attach metadata and save all MizerParams / MizerSim objects properly
# ------------------------------------------------------------------------------
meta_authors <- list(list(
  name        = "Gustav Delius",
  affiliation = "Department of Mathematics, University of York, York, UK"
))

params_f15 <- setMetadata(
  params_f15,
  title       = "Six-Species Theoretical Cod Ecosystem Model (Baseline Sustainable Fishery, E_cod = 0.15)",
  description = paste(
    "Theoretical six-species mizer size-spectrum model (Capelin, Herring,",
    "Mackerel, Sculpin, Plaice, Cod) with exogenous cod recruitment input,",
    "refined at baseline sustainable fishing effort (E_cod = 0.15 /yr, residual < 1e-8)."
  ),
  authors     = meta_authors
)

params_high_f0 <- setMetadata(
  params_high_f0,
  title       = "Six-Species Theoretical Cod Ecosystem Model (Zero-Cod-Fishing Cultivated State, E_cod = 0)",
  description = paste(
    "High-cod cultivated steady state at E_cod = 0 /yr (baseline pelagic,",
    "mesopredator, and flatfish fishing active; residual < 1e-8)."
  ),
  authors     = meta_authors
)

params_low_f0 <- setMetadata(
  params_low_f0,
  title       = "Six-Species Theoretical Cod Ecosystem Model (Zero-Cod-Fishing Depensatory State, E_cod = 0)",
  description = paste(
    "Alternative low-cod, forage-dominated steady state at E_cod = 0 /yr with",
    "identical parameters to cod_high_f0_params.rds (residual < 1e-8)."
  ),
  authors     = meta_authors
)

params_low_f15 <- setMetadata(
  params_low_f15,
  title       = "Six-Species Theoretical Cod Ecosystem Model (Low-Cod State at E_cod = 0.15)",
  description = "Low-cod steady state under baseline reopened cod fishing effort E_cod = 0.15 /yr (residual < 1e-8).",
  authors     = meta_authors
)

params_t61 <- setMetadata(
  params_t61,
  title       = "Six-Species Theoretical Cod Ecosystem Model (Transient Post-Closure State at t = 61 yr)",
  description = paste(
    "Ecosystem state after 27 relative years of intense cod overfishing (t = 0..26)",
    "followed by 34 relative years of zero directed cod fishing (t = 27..61),",
    "integrated with method = 'tr_bdf2' and dt = 0.05 yr."
  ),
  authors     = meta_authors
)

saveParams(params_f15,     "models/cod_f15_params.rds")
saveParams(params_high_f0, "models/cod_high_f0_params.rds")
saveParams(params_low_f0,  "models/cod_low_f0_params.rds")
saveParams(params_low_f15, "models/cod_low_f15_params.rds")
saveParams(params_t61,     "models/cod_t61_params.rds")
saveSim(sim_collapse_closure, "models/sim_collapse_closure.rds")

# Save plain tabular calibration & steady-state summary diagnostics
sp_final <- species_params(params_f15)
cal_target_kt <- sp$biomass_observed / 1e9
calib_summary_df <- data.frame(
  species              = sp$species,
  target_cutoff_kt     = cal_target_kt,
  cutoff_g             = sp$biomass_cutoff,
  uncal_cutoff_kt      = as.numeric(bio_init_uncal),
  cal_pre_sub_kt       = as.numeric(bio_cal_pre_sub),
  rel_err_pre_sub_pct  = 100 * (as.numeric(bio_cal_pre_sub) - cal_target_kt) / cal_target_kt,
  cal_post_sub_kt      = as.numeric(bio_cal_post_sub),
  rel_err_post_sub_pct = 100 * (as.numeric(bio_cal_post_sub) - cal_target_kt) / cal_target_kt,
  h_cal                = sp_final$h,
  gamma_cal            = sp_final$gamma,
  ks_cal               = sp_final$ks,
  erepro_cal           = sp_final$erepro,
  R_max_cal            = sp_final$R_max,
  r_refuge_cal         = sp_final$r_refuge,
  repro_level_cal      = unname(reproduction_level(params)),
  stringsAsFactors     = FALSE
)

rdd_low  <- unclass(getRDD(params_low_f0))
rdd_high <- unclass(getRDD(params_high_f0))
rdd_f15  <- unclass(getRDD(params_f15))

calibration_diagnostics <- list(
  calib_summary_df       = calib_summary_df,
  res_cal_pre_sub        = res_cal_pre_sub,
  res_cal_post_sub       = res_cal_post_sub,
  conv_table             = conv_table,
  cod_r_refuge           = sp_final$r_refuge[6],
  cod_rdd_cal            = unname(rdd_cal["Cod"]),
  cod_rdd_f15            = unname(rdd_f15["Cod"]),
  cod_rdd_high_f0        = unname(rdd_high["Cod"]),
  cod_rdd_low_f0         = unname(rdd_low["Cod"]),
  cod_sub_share_cal_pct  = 100 * sp_final$r_refuge[6] / unname(rdd_cal["Cod"]),
  cod_sub_share_high_pct = 100 * sp_final$r_refuge[6] / unname(rdd_high["Cod"]),
  cod_sub_share_low_pct  = 100 * sp_final$r_refuge[6] / unname(rdd_low["Cod"]),
  t27_cod_ssb_kt         = unname(unclass(getSSB(sim_collapse_closure))["27", "Cod"] / 1e9),
  t61_cod_ssb_kt         = unname(getSSB(params_t61)["Cod"] / 1e9)
)
saveRDS(calibration_diagnostics, "models/calibration_diagnostics.rds")

cat("=== Theoretical Six-Species Model Calibration & Equilibrium Summary ===\n")
print(calib_summary_df[, c("species", "target_cutoff_kt", "cal_pre_sub_kt",
                           "rel_err_pre_sub_pct", "cal_post_sub_kt", "erepro_cal")],
      row.names = FALSE)
cat("\n=== Refined Equilibrium Convergence Summary ===\n")
print(conv_table, row.names = FALSE)
cat(sprintf("\nCod exogenous recruitment input: %.3e /yr (%.2f%% of E_cod=0.22 cal RDD; %.2f%% of low-cod E_cod=0 RDD)\n",
            calibration_diagnostics$cod_r_refuge,
            calibration_diagnostics$cod_sub_share_cal_pct,
            calibration_diagnostics$cod_sub_share_low_pct))
cat(sprintf("Canonical collapse+closure (tr_bdf2, dt=0.05): t=27 Cod SSB = %.2f kt | t=61 Cod SSB = %.2f kt\n",
            calibration_diagnostics$t27_cod_ssb_kt,
            calibration_diagnostics$t61_cod_ssb_kt))
if (sink.number() > 0) sink()
