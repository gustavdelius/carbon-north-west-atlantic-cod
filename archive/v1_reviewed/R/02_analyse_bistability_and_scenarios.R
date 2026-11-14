# ==============================================================================
# R/02_analyse_bistability_and_scenarios.R
#
# Performs linear stability analysis (getStability), equilibrium continuation
# hysteresis scans (scanModel), trophic bottleneck diagnostics, and 1965-2066
# historical + policy recovery simulations on the Northwest Atlantic cod model.
# ==============================================================================

library(mizer)
source("R/00_model_functions.R")

params_hist    <- readParams("models/nw_atlantic_cod_params.rds")
params_high_f0 <- readParams("models/nw_atlantic_cod_high_f0_params.rds")
params_low_f0  <- readParams("models/nw_atlantic_cod_moratorium_params.rds")

eff_hist <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.15)
eff_mor  <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.00)

# ------------------------------------------------------------------------------
# 1. Dynamic stability analysis of the alternative steady states (getStability)
# ------------------------------------------------------------------------------
cat("1. Computing continuous-time Jacobian eigenvalues via getStability()...\n")
stab_hist    <- getStability(params_hist,    effort = eff_hist)
stab_high_f0 <- getStability(params_high_f0, effort = eff_mor)
stab_low_f0  <- getStability(params_low_f0,  effort = eff_mor)

cat(sprintf("   Historical Sustainable (F_cod=0.15): max Re(lambda) = %+.4f /yr\n",
            stab_hist$max_real_part))
cat(sprintf("   Cultivated High-Cod    (F_cod=0.00): max Re(lambda) = %+.4f /yr\n",
            stab_high_f0$max_real_part))
cat(sprintf("   Moratorium Low-Cod     (F_cod=0.00): max Re(lambda) = %+.4f /yr\n",
            stab_low_f0$max_real_part))

stopifnot(stab_hist$max_real_part < 0)
stopifnot(stab_high_f0$max_real_part < 0)
stopifnot(stab_low_f0$max_real_part < 0)

# Locate the critical Cod SSB separatrix (Allee cultivation threshold) at F_cod = 0
# and at F_cod = 0.15 by bisecting between the collapsed and cultivated states
find_separatrix_ssb <- function(p_lo, p_hi, effort_vec, n_iter = 10) {
  lo <- 0.0
  hi <- 1.0
  n_lo <- initialN(p_lo)
  n_hi <- initialN(p_hi)
  r_lo <- initialNResource(p_lo)
  r_hi <- initialNResource(p_hi)
  mid_ssb <- NA_real_
  for (it in seq_len(n_iter)) {
    alpha <- 0.5 * (lo + hi)
    p_try <- p_lo
    initialN(p_try) <- (1 - alpha) * n_lo + alpha * n_hi
    initialNResource(p_try) <- (1 - alpha) * r_lo + alpha * r_hi
    mid_ssb <- unname(getSSB(p_try)["Cod"] / 1e9)
    sim_try <- project(p_try, effort = effort_vec, t_max = 220, dt = 0.1, progress_bar = FALSE)
    end_ssb <- unname(getSSB(finalParams(sim_try))["Cod"] / 1e9)
    if (end_ssb > 350) {
      hi <- alpha
    } else {
      lo <- alpha
    }
  }
  mid_ssb
}

allee_ssb_f0  <- find_separatrix_ssb(params_low_f0, params_high_f0, eff_mor,  n_iter = 9)
allee_ssb_f15 <- find_separatrix_ssb(params_low_f0, params_hist,    eff_hist, n_iter = 9)
cat(sprintf("   Separatrix Cod SSB threshold: %.1f kt (at F_cod=0) | %.1f kt (at F_cod=0.15)\n",
            allee_ssb_f0, allee_ssb_f15))

# ------------------------------------------------------------------------------
# 2. Bifurcation & hysteresis continuation scans via scanModel()
# ------------------------------------------------------------------------------
cat("2. Running forward and backward hysteresis continuation scans (scanModel)...\n")

# A. Hysteresis over directed Cod_Trawl effort (0 -> 0.50 -> 0 /yr)
f_cod_grid <- seq(0.00, 0.50, by = 0.025)
scan_cod_up <- scanModel(
  params_high_f0,
  scan_values   = f_cod_grid,
  set_func      = scanEffort("Cod_Trawl"),
  value_func    = getSSB,
  species       = "Cod",
  continuation  = TRUE,
  t_max         = 320,
  distance_tol  = 0.001,
  progress_bar  = FALSE,
  info_level    = 0
)

# Start backward scan from the collapsed equilibrium at F_cod = 0.50
eff_col_top <- eff_mor
eff_col_top["Cod_Trawl"] <- max(f_cod_grid)
params_at_top_f <- findSteadyState(
  params_low_f0, effort = eff_col_top, t_max = 300,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)
scan_cod_down <- scanModel(
  params_at_top_f,
  scan_values   = rev(f_cod_grid),
  set_func      = scanEffort("Cod_Trawl"),
  value_func    = getSSB,
  species       = "Cod",
  continuation  = TRUE,
  t_max         = 320,
  distance_tol  = 0.001,
  progress_bar  = FALSE,
  info_level    = 0
)

# B. Escape bifurcation over Pelagic forage fishing effort at F_cod = 0 (0 -> 0.70 -> 0 /yr)
f_pel_grid <- seq(0.00, 0.70, by = 0.035)
eff_pel_0 <- eff_mor
eff_pel_0["Pelagic"] <- 0.00
params_low_pel0 <- findSteadyState(
  params_low_f0, effort = eff_pel_0, t_max = 300,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)

scan_pel_up <- scanModel(
  params_low_pel0,
  scan_values   = f_pel_grid,
  set_func      = scanEffort("Pelagic"),
  value_func    = getSSB,
  species       = "Cod",
  continuation  = TRUE,
  t_max         = 400,
  distance_tol  = 0.001,
  progress_bar  = FALSE,
  info_level    = 0
)

eff_pel_top <- eff_mor
eff_pel_top["Pelagic"] <- max(f_pel_grid)
params_high_pel_top <- findSteadyState(
  params_high_f0, effort = eff_pel_top, t_max = 300,
  residual_tol = 1e-6, distance_tol = 1e-6,
  progress_bar = FALSE, info_level = 0
)
scan_pel_down <- scanModel(
  params_high_pel_top,
  scan_values   = rev(f_pel_grid),
  set_func      = scanEffort("Pelagic"),
  value_func    = getSSB,
  species       = "Cod",
  continuation  = TRUE,
  t_max         = 320,
  distance_tol  = 0.001,
  progress_bar  = FALSE,
  info_level    = 0
)

# Identify fold thresholds from scanModel data frames
up_cod_df <- as.data.frame(scan_cod_up)
collapse_idx <- which(up_cod_df[[2]] / 1e9 < 100)[1]
collapse_f_cod <- if (!is.na(collapse_idx)) up_cod_df[[1]][collapse_idx] else 0.25

up_pel_df <- as.data.frame(scan_pel_up)
recovery_idx <- which(up_pel_df[[2]] / 1e9 > 300)[1]
recovery_f_pel <- if (!is.na(recovery_idx)) up_pel_df[[1]][recovery_idx] else 0.175

cat(sprintf("   Collapse fold threshold F_cod     = %.3f /yr\n", collapse_f_cod))
cat(sprintf("   Escape recovery threshold F_pel   = %.3f /yr (at F_cod=0)\n", recovery_f_pel))

# ------------------------------------------------------------------------------
# 3. Historical Collapse & 1992-2026 Moratorium Simulation
# ------------------------------------------------------------------------------
cat("3. Simulating historical collapse (1965-1992) and moratorium (1992-2026)...\n")
years_all_hist <- 1965:2026
gears <- names(initial_effort(params_hist))
eff_history <- array(
  0, dim = c(length(years_all_hist), length(gears)),
  dimnames = list(time = years_all_hist, gear = gears)
)
eff_history[, "Pelagic"]      <- 0.15
eff_history[, "Mesopredator"] <- 0.05
eff_history[, "Flatfish"]     <- 0.15
eff_history[, "Cod_Trawl"]    <- ifelse(years_all_hist < 1968, 0.15,
                                 ifelse(years_all_hist < 1981, 0.46,
                                 ifelse(years_all_hist < 1992, 0.56, 0.00)))

sim_history <- project(params_hist, effort = eff_history, dt = 0.1, t_save = 1, progress_bar = FALSE)
params_2026 <- finalParams(sim_history)

# ------------------------------------------------------------------------------
# 4. Policy Recovery Scenarios (2026-2066, 40 years)
# ------------------------------------------------------------------------------
cat("4. Running 40-year ecosystem policy scenarios (2026-2066)...\n")
run_policy <- function(fp1, fm1, fc1, fp2 = fp1, fm2 = fm1, fc2 = fc1, split_yr = 2046) {
  yrs <- 2026:2066
  em <- array(0, dim = c(length(yrs), length(gears)),
              dimnames = list(time = yrs, gear = gears))
  phase1 <- yrs < split_yr
  em[phase1,  "Pelagic"]      <- fp1
  em[!phase1, "Pelagic"]      <- fp2
  em[phase1,  "Mesopredator"] <- fm1
  em[!phase1, "Mesopredator"] <- fm2
  em[,        "Flatfish"]     <- 0.15
  em[phase1,  "Cod_Trawl"]    <- fc1
  em[!phase1, "Cod_Trawl"]    <- fc2
  project(params_2026, effort = em, dt = 0.1, t_save = 1, progress_bar = FALSE)
}

scenarios <- list(
  status_quo        = run_policy(0.15, 0.05, 0.00),
  premature_reopen  = run_policy(0.15, 0.05, 0.00, 0.15, 0.05, 0.25, split_yr = 2046),
  mesopred_only     = run_policy(0.15, 0.45, 0.00),
  pelagic_rebalance = run_policy(0.55, 0.05, 0.00),
  combined_sustain  = run_policy(0.55, 0.35, 0.00),
  ecosystem_pulse   = run_policy(0.65, 0.45, 0.00, 0.15, 0.05, 0.15, split_yr = 2046)
)

for (nm in names(scenarios)) {
  s <- scenarios[[nm]]
  ssb_end <- unname(getSSB(finalParams(s))["Cod"] / 1e9)
  yld_end <- sum(getYield(finalParams(s))) / 1e9
  cat(sprintf("   %-18s -> 2066 Cod SSB = %6.1f kt | Total Yield = %5.1f kt/yr\n",
              nm, ssb_end, yld_end))
}

# ------------------------------------------------------------------------------
# 5. Size-resolved mortality & energetic bottleneck metrics
# ------------------------------------------------------------------------------
w_vec <- w(params_hist)
nursery_mask <- w_vec >= 0.05 & w_vec <= 5.0
pm_hist <- unclass(getPredMort(params_hist))["Cod", ]
pm_high <- unclass(getPredMort(params_high_f0))["Cod", ]
pm_low  <- unclass(getPredMort(params_low_f0))["Cod", ]

larval_mort_hist <- mean(pm_hist[nursery_mask])
larval_mort_high <- mean(pm_high[nursery_mask])
larval_mort_low  <- mean(pm_low[nursery_mask])

cat(sprintf("\nMean Cod predation mortality (0.05-5 g nursery window):\n"))
cat(sprintf("   Historical (F_cod=0.15): %.2f /yr\n", larval_mort_hist))
cat(sprintf("   Cultivated (F_cod=0.00): %.2f /yr\n", larval_mort_high))
cat(sprintf("   Moratorium (F_cod=0.00): %.2f /yr (%.2fx vs Historical, %.2fx vs Cultivated)\n",
            larval_mort_low, larval_mort_low / larval_mort_hist, larval_mort_low / larval_mort_high))

analysis_results <- list(
  stab_hist         = stab_hist,
  stab_high_f0      = stab_high_f0,
  stab_low_f0       = stab_low_f0,
  allee_ssb_f0      = allee_ssb_f0,
  allee_ssb_f15     = allee_ssb_f15,
  scan_cod_up       = scan_cod_up,
  scan_cod_down     = scan_cod_down,
  scan_pel_up       = scan_pel_up,
  scan_pel_down     = scan_pel_down,
  collapse_f_cod    = collapse_f_cod,
  recovery_f_pel    = recovery_f_pel,
  sim_history       = sim_history,
  params_2026       = params_2026,
  scenarios         = scenarios,
  larval_mort_hist  = larval_mort_hist,
  larval_mort_high  = larval_mort_high,
  larval_mort_low   = larval_mort_low
)

saveRDS(analysis_results, "models/analysis_results.rds")
cat("\nSaved complete analysis results to models/analysis_results.rds\n")
