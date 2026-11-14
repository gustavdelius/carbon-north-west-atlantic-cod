# ==============================================================================
# R/01_build_and_calibrate_model.R
#
# Builds, calibrates, and validates a 6-species Northwest Atlantic shelf
# size-spectrum model in mizer explaining the post-1992 Atlantic cod
# non-recovery ("moratorium trap") via emergent cultivation-depensation
# bistability.
# ==============================================================================

library(mizer)
source("R/00_model_functions.R")

dir.create("models", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Species parameters: Northwest Atlantic shelf fish community
#    Biomass observations represent shelf-scale calibrated reference stocks (g;
#    1e9 g = 1 kilotonne, kt) above survey selectivity cutoffs.
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
# 2. Trophic interaction matrix (rows = predators, columns = prey)
#    Captures ontogenetic role reversal on the Northwest Atlantic shelf:
#    - Adult Cod prey strongly on Capelin, Herring, Mackerel, and Sculpin.
#    - Herring, Mackerel, and Capelin overlap with pelagic Cod eggs/larvae.
#    - Benthic Sculpin ambush newly settled 0-group / age-1 demersal Cod.
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
# 3. Fleet & gear selectivity parameters
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

# Calibrate at a moderate groundfish exploitation reference (F_cod = 0.22 /yr)
# so that historical sustainable fishing (F_cod = 0.15 /yr) sits comfortably
# inside the basin of attraction of the upper cultivated steady state.
eff_cal  <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.22)
eff_hist <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.15)
eff_mor  <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.00)

initial_effort(params) <- eff_cal

# Synchronise single-species size distributions and boundary recruitment before
# running tuneSteadyState()
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

# Register inshore coastal spawning refuge (1.1% of calibrated Cod recruitment,
# representing sheltered inshore fjords such as Smith Sound / Gilbert Bay)
rdd_cal <- getRDD(params)
species_params(params)$r_refuge <- c(0, 0, 0, 0, 0, 0.011 * unname(rdd_cal["Cod"]))
params <- setRateFunction(params, "RDD", "RefugeBevertonHoltRDD")

# Verify physiological realism of reproductive efficiencies (all erepro << 1)
stopifnot(all(species_params(params)$erepro < 1))

# ------------------------------------------------------------------------------
# 5. Locate the three key steady states WITH IDENTICAL BIOLOGICAL PARAMETERS
#    (1) Sustainable pre-collapse fishery steady state (F_cod = 0.15 /yr)
#    (2) Unexploited high-cod cultivated steady state (F_cod = 0.00 /yr)
#    (3) Post-moratorium forage-locked collapsed steady state (F_cod = 0.00 /yr)
# ------------------------------------------------------------------------------
params_hist <- findSteadyState(
  params, effort = eff_hist, t_max = 300,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)

params_high_f0 <- findSteadyState(
  params_hist, effort = eff_mor, t_max = 300,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)

# Reach the collapsed basin via historical groundfish overfishing (1965-1992)
years_hist <- 1965:1992
eff_collapse <- matrix(
  0, nrow = length(years_hist), ncol = 4,
  dimnames = list(time = years_hist, gear = names(eff_hist))
)
eff_collapse[, "Pelagic"]      <- 0.15
eff_collapse[, "Mesopredator"] <- 0.05
eff_collapse[, "Flatfish"]     <- 0.15
eff_collapse[, "Cod_Trawl"]    <- ifelse(years_hist < 1968, 0.15,
                                  ifelse(years_hist < 1981, 0.46, 0.56))

sim_collapse <- project(params_hist, effort = eff_collapse, dt = 0.1, progress_bar = FALSE)
params_1992  <- finalParams(sim_collapse)
sim_mor34    <- project(params_1992, effort = eff_mor, t_max = 34, dt = 0.1, progress_bar = FALSE)
params_2026  <- finalParams(sim_mor34)

params_low_f0 <- findSteadyState(
  params_2026, effort = eff_mor, t_max = 1200,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)

# Confirm identical species, gear, resource, interaction, and RDD parameters across
# both F_cod = 0 steady states, and confirm distinct alternative attractors!
stopifnot(isTRUE(all.equal(species_params(params_high_f0), species_params(params_low_f0))))
stopifnot(isTRUE(all.equal(interaction_matrix(params_high_f0), interaction_matrix(params_low_f0))))
stopifnot(isTRUE(all.equal(resource_params(params_high_f0), resource_params(params_low_f0))))
stopifnot(identical(getRateFunction(params_high_f0), getRateFunction(params_low_f0)))
stopifnot(getSSB(params_low_f0)["Cod"] < 0.05 * getSSB(params_high_f0)["Cod"])

# ------------------------------------------------------------------------------
# 6. Attach metadata and save models using saveParams()
# ------------------------------------------------------------------------------
meta_authors <- list(list(name = "Northwest Atlantic Size-Spectrum Study"))

params_hist <- setMetadata(
  params_hist,
  title       = "Northwest Atlantic Cod Ecosystem Model (Pre-Collapse Sustainable Regime)",
  description = paste(
    "Six-species mizer size-spectrum model of the Northwest Atlantic shelf",
    "(Capelin, Herring, Mackerel, Sculpin, American Plaice, Atlantic Cod)",
    "at the historical sustainable fishery steady state (F_cod = 0.15 /yr)."
  ),
  authors     = meta_authors
)

params_high_f0 <- setMetadata(
  params_high_f0,
  title       = "Northwest Atlantic Cod Ecosystem Model (Unexploited Cultivated Equilibrium, F_cod = 0)",
  description = paste(
    "High-cod cultivated steady state at F_cod = 0 /yr. Adult cod top-down",
    "predation suppresses forage and mesopredator abundances, maintaining low",
    "larval/juvenile cod mortality."
  ),
  authors     = meta_authors
)

params_low_f0 <- setMetadata(
  params_low_f0,
  title       = "Northwest Atlantic Cod Ecosystem Model (Collapsed Moratorium Trap Equilibrium, F_cod = 0)",
  description = paste(
    "Alternative stable low-cod, forage-dominated steady state at F_cod = 0 /yr",
    "with identical biological parameters to params_high_f0. Elevated pelagic",
    "forage and benthic sculpin abundances suppress juvenile cod survival."
  ),
  authors     = meta_authors
)

saveParams(params_hist,    "models/nw_atlantic_cod_params.rds")
saveParams(params_high_f0, "models/nw_atlantic_cod_high_f0_params.rds")
saveParams(params_low_f0,  "models/nw_atlantic_cod_moratorium_params.rds")

# Print concise verification report
cat("=== Model Calibration & Bistability Verification ===\n")
cat("erepro values:\n")
print(setNames(signif(species_params(params_hist)$erepro, 3), sp$species))
cat("\n1. Historical Sustainable Regime (F_cod = 0.15 /yr):\n")
cat("   attractor:", attr(params_hist, "convergence")$attractor,
    "| residual:", signif(attr(params_hist, "convergence")$residual, 3),
    "| isSteady:", isSteady(params_hist), "\n")
cat("   SSB (kt): ", paste(sp$species, round(getSSB(params_hist) / 1e9, 1), sep = "=", collapse = ", "), "\n")

cat("\n2. High-Cod Cultivated Equilibrium (F_cod = 0.00 /yr):\n")
cat("   attractor:", attr(params_high_f0, "convergence")$attractor,
    "| residual:", signif(attr(params_high_f0, "convergence")$residual, 3),
    "| isSteady:", isSteady(params_high_f0), "\n")
cat("   SSB (kt): ", paste(sp$species, round(getSSB(params_high_f0) / 1e9, 1), sep = "=", collapse = ", "), "\n")

cat("\n3. Collapsed Moratorium Equilibrium (F_cod = 0.00 /yr, identical parameters):\n")
cat("   attractor:", attr(params_low_f0, "convergence")$attractor,
    "| residual:", signif(attr(params_low_f0, "convergence")$residual, 3),
    "| isSteady:", isSteady(params_low_f0), "\n")
cat("   SSB (kt): ", paste(sp$species, round(getSSB(params_low_f0) / 1e9, 1), sep = "=", collapse = ", "), "\n")
cat("   1992 Cod SSB:", round(getSSB(params_1992)["Cod"] / 1e9, 1),
    "kt | 2026 Cod SSB:", round(getSSB(params_2026)["Cod"] / 1e9, 1), "kt\n")
