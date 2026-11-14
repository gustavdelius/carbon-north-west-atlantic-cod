# ==============================================================================
# R/02_analyse_bistability_and_scenarios.R
#
# Performs local continuous- and discrete-time stability analysis, Newton saddle
# location, path-specific basin brackets and equal/inverted-SSB counterexamples,
# hysteresis scans without fallback constants, size-resolved predator mortality
# and cumulative hazard decompositions, 4-way R0 factorials, counterfactual
# mechanism isolation, and matched 6-gear management / bycatch / feedback-HCR
# experiments using the shared canonical solver (method = "tr_bdf2", dt = 0.05).
# ==============================================================================

library(mizer)
library(parallel)
source("R/00_model_functions.R")

dir.create("models/sims", showWarnings = FALSE, recursive = TRUE)
dir.create("models/logs", showWarnings = FALSE, recursive = TRUE)
sink("models/logs/02_analyse_bistability_and_scenarios.log", split = TRUE)
NC <- min(8L, detectCores())

params_f15     <- readParams("models/cod_f15_params.rds")
params_high_f0 <- readParams("models/cod_high_f0_params.rds")
params_low_f0  <- readParams("models/cod_low_f0_params.rds")
params_low_f15 <- readParams("models/cod_low_f15_params.rds")
params_t61     <- readParams("models/cod_t61_params.rds")
sim_collapse   <- readSim("models/sim_collapse_closure.rds")

meta_authors <- list(list(
  name        = "Gustav Delius",
  affiliation = "Department of Mathematics, University of York, York, UK"
))

# ------------------------------------------------------------------------------
# 1. Local Stability Analysis: Continuous Jacobian (h sensitivity) & Discrete Step
# ------------------------------------------------------------------------------
cat("1. Computing continuous-time (getStability) & discrete-time (getDiscreteStability) spectra...\n")
flush.console()

h_grid <- c(1e-3, 1e-4, 1e-5, 1e-6)
stab_h_rows <- list()
stab_ref <- list()

for (h_val in h_grid) {
  sh <- getStability(params_f15,     effort = EFFORT_HIST, h = h_val)
  su <- getStability(params_high_f0, effort = EFFORT_MOR,  h = h_val)
  sl <- getStability(params_low_f0,  effort = EFFORT_MOR,  h = h_val)
  s15<- getStability(params_low_f15, effort = EFFORT_HIST, h = h_val)
  if (abs(h_val - 1e-5) < 1e-12) {
    stab_ref <- list(f15 = sh, high_f0 = su, low_f0 = sl, low_f15 = s15)
  }
  stab_h_rows[[length(stab_h_rows) + 1L]] <- data.frame(
    h                 = h_val,
    n_active_f15      = sh$n_active,
    n_active_high_f0  = su$n_active,
    n_active_low_f0   = sl$n_active,
    max_Re_f15        = sh$max_real_part,
    max_Re_high_f0    = su$max_real_part,
    max_Re_low_f0     = sl$max_real_part,
    max_Re_low_f15    = s15$max_real_part,
    dom_period_high   = su$dominant_period,
    osc_period_high   = su$oscillation_period,
    dom_period_low    = sl$dominant_period,
    osc_period_low    = sl$oscillation_period,
    stringsAsFactors  = FALSE
  )
}
stab_h_df <- do.call(rbind, stab_h_rows)
stopifnot(all(stab_h_df$max_Re_f15 < 0),
          all(stab_h_df$max_Re_high_f0 < 0),
          all(stab_h_df$max_Re_low_f0 < 0),
          all(stab_h_df$max_Re_low_f15 < 0))

dt_disc_grid <- c(0.20, 0.10, 0.05)
disc_rows <- lapply(dt_disc_grid, function(dt_v) {
  du <- getDiscreteStability(params_high_f0, effort = EFFORT_MOR,  dt = dt_v)
  dl <- getDiscreteStability(params_low_f0,  effort = EFFORT_MOR,  dt = dt_v)
  dh <- getDiscreteStability(params_f15,     effort = EFFORT_HIST, dt = dt_v)
  data.frame(
    dt                  = dt_v,
    spec_rad_f15        = dh$spectral_radius,
    stable_f15          = dh$stable,
    spec_rad_high_f0    = du$spectral_radius,
    stable_high_f0      = du$stable,
    spec_rad_low_f0     = dl$spectral_radius,
    stable_low_f0       = dl$stable,
    stringsAsFactors    = FALSE
  )
})
disc_stab_df <- do.call(rbind, disc_rows)

# ------------------------------------------------------------------------------
# 2. Path-Specific Chord Brackets, Unstable Saddle Equilibria (Newton), and
#    Equal / Inverted Cod-SSB Counterexamples (Comment 6)
# ------------------------------------------------------------------------------
cat("2. Bracketing state-space chords, polishing Newton saddles, and running equal-SSB counterexamples...\n")
flush.console()

bracket_state_chord <- function(p_lo, p_hi, eff_vec, n_iter = 14L, t_hor = 1000) {
  # Explicitly verify endpoint classifications
  s_lo_end <- unname(getSSB(finalParams(project(
    p_lo, effort = eff_vec, t_max = t_hor, dt = 0.10, t_save = 200,
    method = CANONICAL_METHOD, progress_bar = FALSE
  )))["Cod"] / 1e9)
  s_hi_end <- unname(getSSB(finalParams(project(
    p_hi, effort = eff_vec, t_max = t_hor, dt = 0.10, t_save = 200,
    method = CANONICAL_METHOD, progress_bar = FALSE
  )))["Cod"] / 1e9)
  if (!(s_lo_end < 50 && s_hi_end > 400)) {
    stop("Chord endpoints do not lie in distinct basins: ", s_lo_end, " vs ", s_hi_end)
  }

  # Coarse monotonicity verification across 11 evenly spaced chord fractions
  a_grid <- seq(0, 1, length.out = 11L)
  coarse_fates <- unlist(mclapply(a_grid, function(a) {
    pm <- p_lo
    initialN(pm) <- (1 - a) * initialN(p_lo) + a * initialN(p_hi)
    initialNResource(pm) <- (1 - a) * initialNResource(p_lo) + a * initialNResource(p_hi)
    unname(getSSB(finalParams(project(
      pm, effort = eff_vec, t_max = t_hor, dt = 0.10, t_save = 200,
      method = CANONICAL_METHOD, progress_bar = FALSE
    )))["Cod"] / 1e9) > 300
  }, mc.cores = NC))
  is_monotone <- all(diff(as.integer(coarse_fates)) >= 0)
  if (!is_monotone) stop("Chord fate classification is non-monotone on coarse check!")

  lo <- 0.0; hi <- 1.0
  for (it in seq_len(n_iter)) {
    am <- 0.5 * (lo + hi)
    pm <- p_lo
    initialN(pm) <- (1 - am) * initialN(p_lo) + am * initialN(p_hi)
    initialNResource(pm) <- (1 - am) * initialNResource(p_lo) + am * initialNResource(p_hi)
    s_end <- unname(getSSB(finalParams(project(
      pm, effort = eff_vec, t_max = t_hor, dt = 0.10, t_save = 200,
      method = CANONICAL_METHOD, progress_bar = FALSE
    )))["Cod"] / 1e9)
    if (s_end > 300) hi <- am else lo <- am
  }

  make_mix <- function(a) {
    pm <- p_lo
    initialN(pm) <- (1 - a) * initialN(p_lo) + a * initialN(p_hi)
    initialNResource(pm) <- (1 - a) * initialNResource(p_lo) + a * initialNResource(p_hi)
    pm
  }
  p_lo_b <- make_mix(lo)
  p_hi_b <- make_mix(hi)
  p_mid  <- make_mix(0.5 * (lo + hi))

  list(
    alpha_lo       = lo,
    alpha_hi       = hi,
    cod_ssb_lo_kt  = unname(getSSB(p_lo_b)["Cod"] / 1e9),
    cod_ssb_hi_kt  = unname(getSSB(p_hi_b)["Cod"] / 1e9),
    end_lo_ssb_kt  = s_lo_end,
    end_hi_ssb_kt  = s_hi_end,
    horizon_yr     = t_hor,
    monotone_check = is_monotone,
    p_mid          = p_mid
  )
}

chord_f0  <- bracket_state_chord(params_low_f0,  params_high_f0, EFFORT_MOR,  n_iter = 14L, t_hor = 1000)
chord_f15 <- bracket_state_chord(params_low_f15, params_f15,     EFFORT_HIST, n_iter = 14L, t_hor = 1000)

# Locate true unstable saddle fixed points via separatrix-manifold relaxation + Newton
saddle_list <- mclapply(list(
  f0  = list(chord_f0$p_mid,  EFFORT_MOR,  80, "Zero-Cod-Fishing Unstable Saddle (E_cod = 0)"),
  f15 = list(chord_f15$p_mid, EFFORT_HIST, 65, "Baseline-Fishery Unstable Saddle (E_cod = 0.15)")
), function(cfg) {
  p_guess <- finalParams(project(cfg[[1]], effort = cfg[[2]], t_max = cfg[[3]],
                                 dt = CANONICAL_DT, t_save = cfg[[3]],
                                 method = CANONICAL_METHOD, progress_bar = FALSE))
  p_newt  <- findSteadyState(p_guess, effort = cfg[[2]], solver = "newton",
                             solver_tol = 1e-10, residual_tol = 1e-8, info_level = 0)
  st      <- getStability(p_newt, effort = cfg[[2]])
  cv      <- attr(p_newt, "convergence")
  p_newt  <- setMetadata(p_newt, title = cfg[[4]],
                         description = "Unstable saddle equilibrium separating low- and high-cod basins, located via solver = 'newton'.",
                         authors = meta_authors)
  list(params = p_newt, stab = st, conv = cv,
       cod_ssb_kt = unname(getSSB(p_newt)["Cod"] / 1e9),
       n_pos_eigs = sum(Re(st$eigenvalues) > 0),
       pos_eigs   = st$eigenvalues[Re(st$eigenvalues) > 0])
}, mc.cores = 2L)

saveParams(saddle_list$f0$params,  "models/cod_saddle_f0_params.rds")
saveParams(saddle_list$f15$params, "models/cod_saddle_f15_params.rds")

saddle_summary_df <- data.frame(
  regime           = c("E_cod = 0.00", "E_cod = 0.15"),
  chord_ssb_lo_kt  = c(chord_f0$cod_ssb_lo_kt,  chord_f15$cod_ssb_lo_kt),
  chord_ssb_hi_kt  = c(chord_f0$cod_ssb_hi_kt,  chord_f15$cod_ssb_hi_kt),
  saddle_ssb_kt    = c(saddle_list$f0$cod_ssb_kt, saddle_list$f15$cod_ssb_kt),
  saddle_converged = c(saddle_list$f0$conv$converged, saddle_list$f15$conv$converged),
  saddle_residual  = c(saddle_list$f0$conv$residual,  saddle_list$f15$conv$residual),
  saddle_max_Re    = c(saddle_list$f0$stab$max_real_part, saddle_list$f15$stab$max_real_part),
  saddle_n_pos     = c(saddle_list$f0$n_pos_eigs, saddle_list$f15$n_pos_eigs),
  stringsAsFactors = FALSE
)

# Equal-SSB and Inverted-SSB counterexamples proving Cod SSB is not a scalar separatrix
make_hybrid_state <- function(p_comm, p_cod_shape, target_cod_ssb_kt) {
  p_out <- p_comm
  scale_fac <- (target_cod_ssb_kt * 1e9) / unname(getSSB(p_cod_shape)["Cod"])
  n_mat <- initialN(p_comm)
  n_mat["Cod", ] <- initialN(p_cod_shape)["Cod", ] * scale_fac
  initialN(p_out) <- n_mat
  p_out
}

ce_configs <- list(
  list(label = "Ecod0_CultPrey_22kt",  p_comm = params_high_f0, p_cod = params_low_f0,  ssb = 22.0,  eff = EFFORT_MOR),
  list(label = "Ecod0_TrapPrey_22kt",  p_comm = params_low_f0,  p_cod = params_low_f0,  ssb = 22.0,  eff = EFFORT_MOR),
  list(label = "Ecod0_TrapPrey_35kt",  p_comm = params_low_f0,  p_cod = params_low_f0,  ssb = 35.0,  eff = EFFORT_MOR),
  list(label = "Ecod15_CultPrey_180kt",p_comm = params_f15,     p_cod = params_f15,     ssb = 180.0, eff = EFFORT_HIST),
  list(label = "Ecod15_TrapPrey_180kt",p_comm = params_low_f15, p_cod = params_f15,     ssb = 180.0, eff = EFFORT_HIST)
)
ce_rows <- mclapply(ce_configs, function(cfg) {
  p0 <- make_hybrid_state(cfg$p_comm, cfg$p_cod, cfg$ssb)
  sim_ce <- project(p0, effort = cfg$eff, t_max = 800, dt = CANONICAL_DT,
                    t_save = 20, method = CANONICAL_METHOD, progress_bar = FALSE)
  p_end <- finalParams(sim_ce)
  data.frame(
    case             = cfg$label,
    E_cod            = unname(cfg$eff["Cod_Trawl"]),
    init_cod_ssb_kt  = unname(getSSB(p0)["Cod"] / 1e9),
    init_forage_kt   = sum(getBiomass(p0, use_cutoff = TRUE)[c("Capelin", "Herring", "Mackerel")]) / 1e9,
    init_sculpin_kt  = unname(getBiomass(p0, use_cutoff = TRUE)["Sculpin"] / 1e9),
    asymp_cod_ssb_kt = unname(getSSB(p_end)["Cod"] / 1e9),
    fate             = ifelse(getSSB(p_end)["Cod"] / 1e9 > 300, "High-Cod Cultivated", "Low-Cod Depensatory"),
    stringsAsFactors = FALSE
  )
}, mc.cores = min(NC, length(ce_configs)))
counterexample_df <- do.call(rbind, ce_rows)

# ------------------------------------------------------------------------------
# 3. Hysteresis Continuation Scans (scanModel) & Fine Transition Diagnostics
#    (Zero numeric fallbacks; explicit residual, eigenvalue, and escape-time tracking)
# ------------------------------------------------------------------------------
cat("3. Running hysteresis scans and fine-grid transition/escape-time diagnostics...\n")
flush.console()

f_cod_grid <- seq(0.00, 0.50, by = 0.025)
scan_cod_up <- scanModel(
  params_high_f0,
  scan_values  = f_cod_grid,
  set_func     = scanEffort("Cod_Trawl"),
  value_func   = getSSB,
  species      = "Cod",
  continuation = TRUE,
  t_max        = 1500,
  residual_tol = 1e-5,
  distance_tol = 1e-6,
  progress_bar = FALSE,
  info_level   = 0
)

eff_col_top <- EFFORT_MOR
eff_col_top["Cod_Trawl"] <- max(f_cod_grid)
params_at_top_f <- findSteadyState(
  params_low_f0, effort = eff_col_top, solver = "project", t_max = 3000,
  residual_tol = 1e-8, distance_tol = 1e-10, progress_bar = FALSE, info_level = 0
)
scan_cod_down <- scanModel(
  params_at_top_f,
  scan_values  = rev(f_cod_grid),
  set_func     = scanEffort("Cod_Trawl"),
  value_func   = getSSB,
  species      = "Cod",
  continuation = TRUE,
  t_max        = 1500,
  residual_tol = 1e-5,
  distance_tol = 1e-6,
  progress_bar = FALSE,
  info_level   = 0
)

f_pel_grid <- seq(0.00, 0.70, by = 0.035)
eff_pel_0 <- EFFORT_MOR
eff_pel_0["Pelagic"] <- 0.00
params_low_pel0 <- findSteadyState(
  params_low_f0, effort = eff_pel_0, solver = "project", t_max = 3000,
  residual_tol = 1e-8, distance_tol = 1e-10, progress_bar = FALSE, info_level = 0
)
scan_pel_up <- scanModel(
  params_low_pel0,
  scan_values  = f_pel_grid,
  set_func     = scanEffort("Pelagic"),
  value_func   = getSSB,
  species      = "Cod",
  continuation = TRUE,
  t_max        = 2000,
  residual_tol = 1e-5,
  distance_tol = 1e-6,
  progress_bar = FALSE,
  info_level   = 0
)

eff_pel_top <- EFFORT_MOR
eff_pel_top["Pelagic"] <- max(f_pel_grid)
params_high_pel_top <- findSteadyState(
  params_high_f0, effort = eff_pel_top, solver = "project", t_max = 3000,
  residual_tol = 1e-8, distance_tol = 1e-10, progress_bar = FALSE, info_level = 0
)
scan_pel_down <- scanModel(
  params_high_pel_top,
  scan_values  = rev(f_pel_grid),
  set_func     = scanEffort("Pelagic"),
  value_func   = getSSB,
  species      = "Cod",
  continuation = TRUE,
  t_max        = 1500,
  residual_tol = 1e-5,
  distance_tol = 1e-6,
  progress_bar = FALSE,
  info_level   = 0
)

# Detect coarse scan transition intervals WITHOUT fallback constants
detect_transition_bracket <- function(scan_df, comparator) {
  vals <- as.numeric(scan_df[[1]])
  ssb  <- as.numeric(scan_df[[2]]) / 1e9
  idx  <- which(comparator(ssb))[1]
  if (is.na(idx) || idx <= 1L) {
    stop("No scan transition detected in supplied scanModel output; refusing numeric fallback.")
  }
  c(last_before = vals[idx - 1L], first_after = vals[idx])
}
coarse_cod_bracket <- detect_transition_bracket(as.data.frame(scan_cod_up), function(s) s < 100)
coarse_pel_bracket <- detect_transition_bracket(as.data.frame(scan_pel_up), function(s) s > 300)

# Fine 0.005 /yr equilibrium & leading-eigenvalue scans across both transition regions
fc_fine_vals <- seq(0.215, 0.250, by = 0.005)
fine_cod_df <- do.call(rbind, mclapply(fc_fine_vals, function(fc) {
  ef <- EFFORT_MOR; ef["Cod_Trawl"] <- fc
  ps0 <- findSteadyState(params_high_f0, effort = ef, solver = "project", t_max = 2500,
                         residual_tol = 1e-6, distance_tol = 1e-8, progress_bar = FALSE, info_level = 0)
  ps  <- findSteadyState(ps0, effort = ef, solver = "newton",
                         solver_tol = 1e-10, residual_tol = 1e-8, info_level = 0)
  st <- getStability(ps, effort = ef)
  cv <- attr(ps, "convergence")
  data.frame(E_cod = fc, cod_ssb_kt = unname(getSSB(ps)["Cod"] / 1e9),
             converged = cv$converged, attractor = as.character(cv$attractor),
             residual = cv$residual, years = attr(ps0, "convergence")$years,
             max_Re = st$max_real_part, stringsAsFactors = FALSE)
}, mc.cores = NC))

fp_fine_vals <- seq(0.145, 0.185, by = 0.005)
fine_pel_df <- do.call(rbind, mclapply(fp_fine_vals, function(fp) {
  ef <- EFFORT_MOR; ef["Pelagic"] <- fp
  ps0 <- findSteadyState(params_low_f0, effort = ef, solver = "project", t_max = 3000,
                         residual_tol = 1e-6, distance_tol = 1e-8, progress_bar = FALSE, info_level = 0)
  ps  <- findSteadyState(ps0, effort = ef, solver = "newton",
                         solver_tol = 1e-10, residual_tol = 1e-8, info_level = 0)
  st <- getStability(ps, effort = ef)
  cv <- attr(ps, "convergence")
  data.frame(E_pel = fp, cod_ssb_kt = unname(getSSB(ps)["Cod"] / 1e9),
             converged = cv$converged, attractor = as.character(cv$attractor),
             residual = cv$residual, years = attr(ps0, "convergence")$years,
             max_Re = st$max_real_part, stringsAsFactors = FALSE)
}, mc.cores = NC))

idx_fine_cod <- which(fine_cod_df$cod_ssb_kt < 100)[1]
idx_fine_pel <- which(fine_pel_df$cod_ssb_kt > 300)[1]
if (is.na(idx_fine_cod) || is.na(idx_fine_pel)) stop("Fine scan failed to bracket transition!")
fine_cod_bracket <- c(last_upper = fine_cod_df$E_cod[idx_fine_cod - 1L],
                      first_lower = fine_cod_df$E_cod[idx_fine_cod])
fine_pel_bracket <- c(last_lower = fine_pel_df$E_pel[idx_fine_pel - 1L],
                      first_upper = fine_pel_df$E_pel[idx_fine_pel])

# Direct escape-time curve T_500kt(E_pel) starting from params_low_f0 at E_cod = 0
escape_pel_vals <- c(0.150, 0.160, 0.165, 0.170, 0.175, 0.185, 0.200, 0.300, 0.450, 0.550, 0.650)
escape_time_df <- do.call(rbind, mclapply(escape_pel_vals, function(fp) {
  ef <- EFFORT_MOR; ef["Pelagic"] <- fp
  sim_e <- project(params_low_f0, effort = ef, t_max = 1500, dt = 0.10, t_save = 1,
                   method = CANONICAL_METHOD, progress_bar = FALSE)
  ssb_t <- unclass(getSSB(sim_e))[, "Cod"] / 1e9
  w100  <- which(ssb_t >= 100)[1]
  w500  <- which(ssb_t >= 500)[1]
  data.frame(
    E_pel            = fp,
    t_to_100kt_yr    = if (is.na(w100)) NA_real_ else as.numeric(names(ssb_t)[w100]),
    t_to_500kt_yr    = if (is.na(w500)) NA_real_ else as.numeric(names(ssb_t)[w500]),
    ssb_at_40yr_kt   = unname(ssb_t["40"]),
    ssb_at_1500yr_kt = unname(tail(ssb_t, 1)),
    stringsAsFactors = FALSE
  )
}, mc.cores = NC))

# ------------------------------------------------------------------------------
# 4. Mechanism Diagnostics: Predator-by-Size Mortality, Feeding-Level Ratios,
#    Cumulative Hazard Decomposition, 4-Way R0 Factorial, and Counterfactuals
# ------------------------------------------------------------------------------
cat("4. Computing predator-resolved mortality, cumulative hazard, R0 factorial, and counterfactuals...\n")
flush.console()

w_vec  <- w(params_high_f0)
dw_vec <- dw(params_high_f0)

pm_pred_high <- decompose_cod_pred_mort(params_high_f0)
pm_pred_low  <- decompose_cod_pred_mort(params_low_f0)
pm_pred_f15  <- decompose_cod_pred_mort(params_f15)

haz_high_15  <- compute_cod_hazard(params_high_f0, 0.001, 15.0, EFFORT_MOR)
haz_low_15   <- compute_cod_hazard(params_low_f0,  0.001, 15.0, EFFORT_MOR)
haz_f15_15   <- compute_cod_hazard(params_f15,     0.001, 15.0, EFFORT_HIST)

haz_high_mat <- compute_cod_hazard(params_high_f0, 0.001, 850.0, EFFORT_MOR)
haz_low_mat  <- compute_cod_hazard(params_low_f0,  0.001, 850.0, EFFORT_MOR)
haz_f15_mat  <- compute_cod_hazard(params_f15,     0.001, 850.0, EFFORT_HIST)

# Exact algebraic decomposition of Delta H = H_low - H_high on [0.001, 15) g:
# m_low/g_low - m_high/g_high = (dm)/g_high + m_high*(1/g_low - 1/g_high) + dm*(1/g_low - 1/g_high)
mask_15 <- w_vec >= 0.001 & w_vec < 15.0
gh <- unclass(getEGrowth(params_high_f0))["Cod", ]
gl <- unclass(getEGrowth(params_low_f0))["Cod", ]
mh <- unclass(getMort(params_high_f0))["Cod", ]
ml <- unclass(getMort(params_low_f0))["Cod", ]
erh <- unclass(getERepro(params_high_f0))["Cod", ]
erl <- unclass(getERepro(params_low_f0))["Cod", ]

dH_mort_only <- sum(((ml[mask_15] - mh[mask_15]) / gh[mask_15]) * dw_vec[mask_15])
dH_grow_only <- sum((mh[mask_15] * (1 / gl[mask_15] - 1 / gh[mask_15])) * dw_vec[mask_15])
dH_synergy   <- sum(((ml[mask_15] - mh[mask_15]) * (1 / gl[mask_15] - 1 / gh[mask_15])) * dw_vec[mask_15])
dH_total     <- haz_low_15$H_total - haz_high_15$H_total
stopifnot(abs((dH_mort_only + dH_grow_only + dH_synergy) - dH_total) < 1e-10)

# Exact feeding-level differences across explicit windows
fl_high <- unclass(getFeedingLevel(params_high_f0))["Cod", ]
fl_low  <- unclass(getFeedingLevel(params_low_f0))["Cod", ]
fl_rel_pct <- 100 * (fl_low - fl_high) / fl_high
fl_pp_diff <- 100 * (fl_low - fl_high)

mask_005_5   <- w_vec >= 0.05  & w_vec <= 5.0
mask_0002_15 <- w_vec >= 0.002 & w_vec <= 15.0
mask_adult   <- w_vec >= 1000  & w_vec <= 10000

feeding_summary_df <- data.frame(
  window          = c("0.05-5 g", "0.002-15 g", "1-10 kg (adults)"),
  ratio_min       = c(min(fl_low[mask_005_5] / fl_high[mask_005_5]),
                      min(fl_low[mask_0002_15] / fl_high[mask_0002_15]),
                      min(fl_low[mask_adult] / fl_high[mask_adult])),
  ratio_max       = c(max(fl_low[mask_005_5] / fl_high[mask_005_5]),
                      max(fl_low[mask_0002_15] / fl_high[mask_0002_15]),
                      max(fl_low[mask_adult] / fl_high[mask_adult])),
  rel_pct_min     = c(min(fl_rel_pct[mask_005_5]), min(fl_rel_pct[mask_0002_15]), min(fl_rel_pct[mask_adult])),
  rel_pct_max     = c(max(fl_rel_pct[mask_005_5]), max(fl_rel_pct[mask_0002_15]), max(fl_rel_pct[mask_adult])),
  pp_diff_min     = c(min(fl_pp_diff[mask_005_5]), min(fl_pp_diff[mask_0002_15]), min(fl_pp_diff[mask_adult])),
  pp_diff_max     = c(max(fl_pp_diff[mask_005_5]), max(fl_pp_diff[mask_0002_15]), max(fl_pp_diff[mask_adult])),
  stringsAsFactors = FALSE
)

# Nursery unweighted vs abundance-weighted predation mortality on 0.05-5 g
pm_high_tot <- unclass(getPredMort(params_high_f0))["Cod", ]
pm_low_tot  <- unclass(getPredMort(params_low_f0))["Cod", ]
pm_f15_tot  <- unclass(getPredMort(params_f15))["Cod", ]
nh_cod <- unclass(initialN(params_high_f0))["Cod", ]
nl_cod <- unclass(initialN(params_low_f0))["Cod", ]
nf_cod <- unclass(initialN(params_f15))["Cod", ]

nursery_mort_summary <- list(
  unweighted_high = mean(pm_high_tot[mask_005_5]),
  unweighted_f15  = mean(pm_f15_tot[mask_005_5]),
  unweighted_low  = mean(pm_low_tot[mask_005_5]),
  abund_wt_high   = sum(pm_high_tot[mask_005_5] * nh_cod[mask_005_5] * dw_vec[mask_005_5]) /
                    sum(nh_cod[mask_005_5] * dw_vec[mask_005_5]),
  abund_wt_f15    = sum(pm_f15_tot[mask_005_5] * nf_cod[mask_005_5] * dw_vec[mask_005_5]) /
                    sum(nf_cod[mask_005_5] * dw_vec[mask_005_5]),
  abund_wt_low    = sum(pm_low_tot[mask_005_5] * nl_cod[mask_005_5] * dw_vec[mask_005_5]) /
                    sum(nl_cod[mask_005_5] * dw_vec[mask_005_5])
)

# 6-way demographic R0(mu, g, er) factorial
R0_factorial_df <- data.frame(
  combination = c(
    "Cultivated baseline: R0(mu_high, g_high, er_high)",
    "Juvenile growth slowdown only: R0(mu_high, g_low, er_high)",
    "Resource + forage shift on growth & adult fecundity: R0(mu_high, g_low, er_low)",
    "Early-life predation surge only: R0(mu_low, g_high, er_high)",
    "Predation surge + juvenile growth slowdown: R0(mu_low, g_low, er_high)",
    "Full depensatory state (with elevated adult feeding): R0(mu_low, g_low, er_low)"
  ),
  R0_value = c(
    compute_cod_R0(gh, mh, erh, params_high_f0),
    compute_cod_R0(gl, mh, erh, params_high_f0),
    compute_cod_R0(gl, mh, erl, params_high_f0),
    compute_cod_R0(gh, ml, erh, params_high_f0),
    compute_cod_R0(gl, ml, erh, params_high_f0),
    compute_cod_R0(gl, ml, erl, params_high_f0)
  ),
  stringsAsFactors = FALSE
)

# Dynamic mechanism counterfactuals starting from params_low_f0 at E_cod = 0
p_cf_res <- params_low_f0
initialNResource(p_cf_res) <- initialNResource(params_high_f0)
p_cf_res <- setResource(p_cf_res, resource_dynamics = "resource_constant")

p_cf_nopred <- params_low_f0
im_np <- interaction_matrix(p_cf_nopred)
im_np[c("Capelin", "Herring", "Mackerel", "Sculpin", "Plaice"), "Cod"] <- 0
interaction_matrix(p_cf_nopred) <- im_np

cf_res_out <- mclapply(list(
  res_fixed_high = p_cf_res,
  no_pred_on_cod = p_cf_nopred
), function(p_cf) {
  s_cf <- project(p_cf, effort = EFFORT_MOR, t_max = 800, dt = CANONICAL_DT,
                  t_save = 20, method = CANONICAL_METHOD, progress_bar = FALSE)
  unname(getSSB(finalParams(s_cf))["Cod"] / 1e9)
}, mc.cores = 2L)

# ------------------------------------------------------------------------------
# 5. Matched Management Scenarios (6-Gear), Capelin Benchmark, Bycatch Sweeps,
#    and State-Dependent Feedback HCR with Observation Error (Comments 4, 10, 12)
# ------------------------------------------------------------------------------
cat("5. Running matched 6-gear policies, capelin benchmark, bycatch sweeps, and feedback HCR...\n")
flush.console()

# Define zero-capelin-fishing benchmark B0_capelin in both low-cod and high-cod contexts
p6_high_cap0 <- make_6gear_params(params_high_f0)
eff6_cap0_mor <- c(Capelin_Gear = 0.00, Clupeid_Mackerel = 0.15, Mesopredator = 0.05,
                   Flatfish = 0.15, Cod_Trawl = 0.00, Cod_Bycatch = 0.00)
eff6_cap0_f15 <- c(Capelin_Gear = 0.00, Clupeid_Mackerel = 0.15, Mesopredator = 0.05,
                   Flatfish = 0.15, Cod_Trawl = 0.15, Cod_Bycatch = 0.00)
p6_high_cap0_s <- findSteadyState(p6_high_cap0, effort = eff6_cap0_mor, solver = "project",
                                  t_max = 3000, residual_tol = 1e-8, distance_tol = 1e-10,
                                  progress_bar = FALSE, info_level = 0)
p6_f15_cap0_s  <- findSteadyState(make_6gear_params(params_f15), effort = eff6_cap0_f15,
                                  solver = "project", t_max = 3000, residual_tol = 1e-8,
                                  distance_tol = 1e-10, progress_bar = FALSE, info_level = 0)
B0_capelin_high_f0_kt  <- unname(getBiomass(p6_high_cap0_s, use_cutoff = TRUE)["Capelin"] / 1e9)
B0_capelin_f15_kt      <- unname(getBiomass(p6_f15_cap0_s,  use_cutoff = TRUE)["Capelin"] / 1e9)
B_capelin_hist_f15_kt  <- unname(getBiomass(params_f15,     use_cutoff = TRUE)["Capelin"] / 1e9)
B_capelin_low_f0_kt    <- unname(getBiomass(params_low_f0,  use_cutoff = TRUE)["Capelin"] / 1e9)

p6_t61 <- make_6gear_params(params_t61)

policy_specs <- list(
  status_quo         = list(label = "1. Status-quo closure",                     args = c(0.15, 0.15, 0.05, 0.00, 0, 0.15, 0.15, 0.05, 0.00, 0)),
  premature_reopen   = list(label = "2. Unassisted cod reopening (tau=20)",      args = c(0.15, 0.15, 0.05, 0.00, 0, 0.15, 0.15, 0.05, 0.25, 0)),
  mesopred_only      = list(label = "3. Sculpin harvest only (sustained)",       args = c(0.15, 0.15, 0.45, 0.00, 0, 0.15, 0.15, 0.45, 0.00, 0)),
  pelagic_rebalance  = list(label = "4. Joint pelagic harvest (sustained)",      args = c(0.55, 0.55, 0.05, 0.00, 0, 0.55, 0.55, 0.05, 0.00, 0)),
  combined_sustain55 = list(label = "5. Combined harvest 0.55/0.35 (sustained)", args = c(0.55, 0.55, 0.35, 0.00, 0, 0.55, 0.55, 0.35, 0.00, 0)),
  joint_sustain_f0   = list(label = "6. Joint 0.65/0.45 sustained (E_cod=0)",    args = c(0.65, 0.65, 0.45, 0.00, 0, 0.65, 0.65, 0.45, 0.00, 0)),
  joint_taper_f0     = list(label = "7. Joint 0.65/0.45 tapered (E_cod=0)",      args = c(0.65, 0.65, 0.45, 0.00, 0, 0.15, 0.15, 0.05, 0.00, 0)),
  joint_sustain_f15  = list(label = "8. Joint 0.65/0.45 sustained + E_cod=0.15", args = c(0.65, 0.65, 0.45, 0.00, 0, 0.65, 0.65, 0.45, 0.15, 0)),
  joint_taper_f15    = list(label = "9. Joint 0.65/0.45 tapered + E_cod=0.15",   args = c(0.65, 0.65, 0.45, 0.00, 0, 0.15, 0.15, 0.05, 0.15, 0)),
  spare_cap_base     = list(label = "10. Capelin-sparing (E_cap=0.15) tapered",  args = c(0.15, 0.65, 0.45, 0.00, 0, 0.15, 0.15, 0.05, 0.15, 0)),
  spare_cap_zero     = list(label = "11. Capelin-closed (E_cap=0.00) tapered",   args = c(0.00, 0.65, 0.45, 0.00, 0, 0.15, 0.15, 0.05, 0.15, 0))
)

pol_runs <- mclapply(names(policy_specs), function(nm) {
  a <- policy_specs[[nm]]$args
  out <- run_policy_6gear(
    p6_t61,
    ecap1 = a[1], ehm1 = a[2], em1 = a[3], ec1 = a[4], ebyc1 = a[5],
    ecap2 = a[6], ehm2 = a[7], em2 = a[8], ec2 = a[9], ebyc2 = a[10],
    dt = CANONICAL_DT, asymp_extra_yrs = 300, asymp_dt = CANONICAL_DT
  )
  sim40 <- out$sim40
  p_end40 <- finalParams(sim40)
  p_asymp <- finalParams(out$sim_asymp)

  ssb <- unclass(getSSB(sim40)) / 1e9
  bio <- unclass(getBiomass(sim40, use_cutoff = TRUE)) / 1e9
  yld <- unclass(getYield(sim40)) / 1e9
  nr  <- nrow(ssb)

  # Trapezoidal integration of annual yield rate over tau in [0, 40] (t_save = 1 yr)
  cum_tot_yld <- sum(0.5 * (rowSums(yld[-nr, , drop = FALSE]) + rowSums(yld[-1, , drop = FALSE])))
  cum_cod_yld <- sum(0.5 * (yld[-nr, "Cod"] + yld[-1, "Cod"]))

  asymp_yld_vec <- unclass(getYield(p_asymp)) / 1e9
  asymp_cod_ssb <- unname(getSSB(p_asymp)["Cod"] / 1e9)

  list(
    sim40 = sim40,
    row = data.frame(
      id                = nm,
      strategy          = policy_specs[[nm]]$label,
      ecap_ph1          = a[1], ehm_ph1 = a[2], emeso_ph1 = a[3], ecod_ph2 = a[9],
      tapered           = !(a[1] == a[6] && a[2] == a[7] && a[3] == a[8]),
      t20_cod_ssb_kt    = unname(ssb["20", "Cod"]),
      t40_cod_ssb_kt    = unname(ssb[nr, "Cod"]),
      asymp_cod_ssb_kt  = asymp_cod_ssb,
      asymp_recovered   = asymp_cod_ssb > 400,
      min_cap_ph1_kt    = min(bio[1:21, "Capelin"]),
      min_cap_40yr_kt   = min(bio[, "Capelin"]),
      min_cap_pct_B0f15 = 100 * min(bio[, "Capelin"]) / B0_capelin_f15_kt,
      cum40_tot_yld_kt  = cum_tot_yld,
      cum40_cod_yld_kt  = cum_cod_yld,
      t40_tot_yld_kt    = sum(yld[nr, ]),
      t40_cod_yld_kt    = unname(yld[nr, "Cod"]),
      asymp_tot_yld_kt  = sum(asymp_yld_vec),
      asymp_cod_yld_kt  = unname(asymp_yld_vec["Cod"]),
      stringsAsFactors  = FALSE
    )
  )
}, mc.cores = NC)
names(pol_runs) <- names(policy_specs)

policy_comparison_df <- do.call(rbind, lapply(pol_runs, `[[`, "row"))

# Save key MizerSim objects properly using saveSim()
for (nm in c("status_quo", "premature_reopen", "mesopred_only", "pelagic_rebalance",
             "combined_sustain55", "joint_sustain_f15", "joint_taper_f15", "spare_cap_base")) {
  saveSim(pol_runs[[nm]]$sim40, sprintf("models/sims/sim_%s.rds", nm))
}

# Phase-I Cod bycatch sensitivity sweeps (juvenile l50=22cm and adult l50=52cm)
byc_vals <- c(0.000, 0.005, 0.010, 0.020, 0.030, 0.050, 0.100)
run_byc_row <- function(mode_name, ecap1, fb_juv, fb_ad) {
  out <- run_policy_6gear(
    p6_t61,
    ecap1 = ecap1, ehm1 = 0.65, em1 = 0.45, ec1 = fb_ad, ebyc1 = fb_juv,
    ecap2 = 0.15,  ehm2 = 0.15, em2 = 0.05, ec2 = 0.15,  ebyc2 = 0.00,
    dt = CANONICAL_DT, asymp_extra_yrs = 300, asymp_dt = 0.10
  )
  ssb <- unclass(getSSB(out$sim40)) / 1e9
  asymp_ssb <- unname(getSSB(finalParams(out$sim_asymp))["Cod"] / 1e9)
  data.frame(
    policy_mode      = mode_name,
    bycatch_effort   = max(fb_juv, fb_ad),
    t20_cod_ssb_kt   = unname(ssb["20", "Cod"]),
    t40_cod_ssb_kt   = unname(ssb[nrow(ssb), "Cod"]),
    asymp_cod_ssb_kt = asymp_ssb,
    asymp_recovered  = asymp_ssb > 400,
    stringsAsFactors = FALSE
  )
}

byc_jobs <- c(
  lapply(byc_vals, function(v) list("Joint_Juvenile_Bycatch",    0.65, v, 0)),
  lapply(byc_vals, function(v) list("Joint_Adult_Bycatch",       0.65, 0, v)),
  lapply(byc_vals, function(v) list("CapelinSpare_Juv_Bycatch",  0.15, v, 0)),
  lapply(byc_vals, function(v) list("CapelinSpare_Adult_Bycatch",0.15, 0, v))
)
bycatch_sweep_df <- do.call(rbind, mclapply(byc_jobs, function(j) {
  run_byc_row(j[[1]], j[[2]], j[[3]], j[[4]])
}, mc.cores = NC))

# State-dependent feedback HCR (annual monitoring, 1-yr lag, 2-yr confirmation >= 280 kt)
run_feedback_hcr <- function(trigger_kt = 280, obs_cv = 0.0, seed = 42L, spare_cap = FALSE) {
  set.seed(seed)
  p_cur     <- if (spare_cap) p6_t61 else params_t61
  switched  <- FALSE
  switch_yr <- NA_integer_
  prev_obs  <- unname(getSSB(p_cur)["Cod"] / 1e9)
  consec    <- 0L
  tot_yld   <- 0.0
  cod_yld   <- 0.0
  ssb_trace <- numeric(41L)
  ssb_trace[1L] <- prev_obs

  for (yr in 1:40) {
    obs_ssb <- if (obs_cv > 0) prev_obs * exp(rnorm(1L, -0.5 * obs_cv^2, obs_cv)) else prev_obs
    if (!switched) {
      if (obs_ssb >= trigger_kt) consec <- consec + 1L else consec <- 0L
      if (consec >= 2L) {
        switched  <- TRUE
        switch_yr <- yr - 1L
      }
    }
    eff_yr <- if (!spare_cap) {
      if (!switched) c(Pelagic = 0.65, Mesopredator = 0.45, Flatfish = 0.15, Cod_Trawl = 0.00) else EFFORT_HIST
    } else {
      if (!switched) c(Capelin_Gear = 0.15, Clupeid_Mackerel = 0.65, Mesopredator = 0.45,
                       Flatfish = 0.15, Cod_Trawl = 0.00, Cod_Bycatch = 0.00)
      else c(Capelin_Gear = 0.15, Clupeid_Mackerel = 0.15, Mesopredator = 0.05,
             Flatfish = 0.15, Cod_Trawl = 0.15, Cod_Bycatch = 0.00)
    }
    prev_obs <- unname(getSSB(p_cur)["Cod"] / 1e9)
    sim1 <- project(p_cur, effort = eff_yr, t_max = 1, dt = CANONICAL_DT, t_save = 1,
                    method = CANONICAL_METHOD, progress_bar = FALSE)
    yld_mat <- unclass(getYield(sim1)) / 1e9
    tot_yld <- tot_yld + 0.5 * (sum(yld_mat[1, ]) + sum(yld_mat[2, ]))
    cod_yld <- cod_yld + 0.5 * (yld_mat[1, "Cod"] + yld_mat[2, "Cod"])
    p_cur   <- finalParams(sim1)
    ssb_trace[yr + 1L] <- unname(getSSB(p_cur)["Cod"] / 1e9)
  }

  eff_end <- if (!spare_cap) EFFORT_HIST else
    c(Capelin_Gear = 0.15, Clupeid_Mackerel = 0.15, Mesopredator = 0.05,
      Flatfish = 0.15, Cod_Trawl = 0.15, Cod_Bycatch = 0.00)
  p_asymp <- finalParams(project(p_cur, effort = eff_end, t_max = 300, dt = 0.10,
                                 t_save = 100, method = CANONICAL_METHOD, progress_bar = FALSE))
  list(
    switch_yr        = switch_yr,
    t40_cod_ssb_kt   = ssb_trace[41L],
    asymp_cod_ssb_kt = unname(getSSB(p_asymp)["Cod"] / 1e9),
    cum40_tot_yld_kt = tot_yld,
    cum40_cod_yld_kt = cod_yld,
    ssb_trace        = ssb_trace
  )
}

hcr_det_joint <- run_feedback_hcr(280, 0.0, 42L, spare_cap = FALSE)
hcr_det_spare <- run_feedback_hcr(280, 0.0, 42L, spare_cap = TRUE)
hcr_mc_joint  <- mclapply(1:16, function(s) run_feedback_hcr(280, 0.25, 200L + s, spare_cap = FALSE), mc.cores = NC)
hcr_mc_spare  <- mclapply(1:16, function(s) run_feedback_hcr(280, 0.25, 400L + s, spare_cap = TRUE),  mc.cores = NC)

summarize_hcr_mc <- function(mc_list, label) {
  sw <- sapply(mc_list, `[[`, "switch_yr")
  t4 <- sapply(mc_list, `[[`, "t40_cod_ssb_kt")
  as <- sapply(mc_list, `[[`, "asymp_cod_ssb_kt")
  cy <- sapply(mc_list, `[[`, "cum40_tot_yld_kt")
  data.frame(
    rule               = label,
    n_reps             = length(mc_list),
    n_switched         = sum(!is.na(sw)),
    switch_yr_median   = median(sw, na.rm = TRUE),
    switch_yr_min      = min(sw, na.rm = TRUE),
    switch_yr_max      = max(sw, na.rm = TRUE),
    t40_ssb_mean_kt    = mean(t4),
    t40_ssb_min_kt     = min(t4),
    t40_ssb_max_kt     = max(t4),
    asymp_recovery_pct = 100 * mean(as > 400),
    cum40_yld_mean_kt  = mean(cy),
    stringsAsFactors   = FALSE
  )
}
hcr_summary_df <- rbind(
  summarize_hcr_mc(list(hcr_det_joint), "Feedback_Joint_Det_CV0"),
  summarize_hcr_mc(hcr_mc_joint,        "Feedback_Joint_Noisy_CV0.25"),
  summarize_hcr_mc(list(hcr_det_spare), "Feedback_SpareCap_Det_CV0"),
  summarize_hcr_mc(hcr_mc_spare,        "Feedback_SpareCap_Noisy_CV0.25")
)

# ------------------------------------------------------------------------------
# 6. Assemble plain numeric/tabular results (NO MizerParams or MizerSim inside!)
# ------------------------------------------------------------------------------
extract_scan_clean <- function(sc) {
  d <- as.data.frame(sc)
  data.frame(
    effort_val = as.numeric(d[[1]]),
    cod_ssb_kt = as.numeric(d[[2]]) / 1e9,
    ymin_kt    = if ("ymin" %in% names(d)) as.numeric(d$ymin) / 1e9 else as.numeric(d[[2]]) / 1e9,
    ymax_kt    = if ("ymax" %in% names(d)) as.numeric(d$ymax) / 1e9 else as.numeric(d[[2]]) / 1e9,
    attractor  = if ("attractor" %in% names(d)) as.character(d$attractor) else NA_character_,
    residual   = if ("residual" %in% names(d)) as.numeric(d$residual) else NA_real_,
    stringsAsFactors = FALSE
  )
}

analysis_results <- list(
  stab_h_df            = stab_h_df,
  disc_stab_df         = disc_stab_df,
  stab_eigenvalues     = list(
    f15     = stab_ref$f15$eigenvalues,
    high_f0 = stab_ref$high_f0$eigenvalues,
    low_f0  = stab_ref$low_f0$eigenvalues,
    low_f15 = stab_ref$low_f15$eigenvalues,
    sad_f0  = saddle_list$f0$stab$eigenvalues,
    sad_f15 = saddle_list$f15$stab$eigenvalues
  ),
  chord_f0             = chord_f0[setdiff(names(chord_f0), "p_mid")],
  chord_f15            = chord_f15[setdiff(names(chord_f15), "p_mid")],
  saddle_summary_df    = saddle_summary_df,
  counterexample_df    = counterexample_df,
  scan_cod_up_df       = extract_scan_clean(scan_cod_up),
  scan_cod_down_df     = extract_scan_clean(scan_cod_down),
  scan_pel_up_df       = extract_scan_clean(scan_pel_up),
  scan_pel_down_df     = extract_scan_clean(scan_pel_down),
  coarse_cod_bracket   = coarse_cod_bracket,
  coarse_pel_bracket   = coarse_pel_bracket,
  fine_cod_df          = fine_cod_df,
  fine_pel_df          = fine_pel_df,
  fine_cod_bracket     = fine_cod_bracket,
  fine_pel_bracket     = fine_pel_bracket,
  escape_time_df       = escape_time_df,
  w_vec                = w_vec,
  dw_vec               = dw_vec,
  pm_pred_high         = pm_pred_high,
  pm_pred_low          = pm_pred_low,
  pm_pred_f15          = pm_pred_f15,
  haz_high_15          = haz_high_15,
  haz_low_15           = haz_low_15,
  haz_f15_15           = haz_f15_15,
  haz_high_mat         = haz_high_mat,
  haz_low_mat          = haz_low_mat,
  haz_f15_mat          = haz_f15_mat,
  hazard_decomp_15g    = list(
    H_high       = haz_high_15$H_total,
    H_low        = haz_low_15$H_total,
    dH_total     = dH_total,
    dH_mort_only = dH_mort_only,
    dH_grow_only = dH_grow_only,
    dH_synergy   = dH_synergy
  ),
  feeding_summary_df   = feeding_summary_df,
  nursery_mort_summary = nursery_mort_summary,
  R0_factorial_df      = R0_factorial_df,
  counterfactual_ssb   = cf_res_out,
  capelin_benchmarks   = list(
    B0_capelin_high_f0_kt = B0_capelin_high_f0_kt,
    B0_capelin_f15_kt     = B0_capelin_f15_kt,
    B_capelin_hist_f15_kt = B_capelin_hist_f15_kt,
    B_capelin_low_f0_kt   = B_capelin_low_f0_kt
  ),
  policy_comparison_df = policy_comparison_df,
  bycatch_sweep_df     = bycatch_sweep_df,
  hcr_summary_df       = hcr_summary_df,
  hcr_det_traces       = list(
    joint_ssb = hcr_det_joint$ssb_trace,
    spare_ssb = hcr_det_spare$ssb_trace,
    joint_sw  = hcr_det_joint$switch_yr,
    spare_sw  = hcr_det_spare$switch_yr
  )
)

saveRDS(analysis_results, "models/analysis_results.rds")

cat("\n=== Summary of Key Bistability & Policy Analysis Results ===\n")
cat("Saddle & Chord Summary:\n")
print(saddle_summary_df, row.names = FALSE)
cat("\nEqual / Inverted Cod-SSB Counterexamples:\n")
print(counterexample_df, row.names = FALSE)
cat(sprintf("\nFine transition brackets: E_cod in [%.3f, %.3f] /yr | E_pel in [%.3f, %.3f] /yr\n",
            fine_cod_bracket[1], fine_cod_bracket[2], fine_pel_bracket[1], fine_pel_bracket[2]))
cat("\nMatched 6-Gear Policy Comparison:\n")
print(policy_comparison_df[, c("id", "t20_cod_ssb_kt", "t40_cod_ssb_kt", "asymp_cod_ssb_kt",
                               "min_cap_ph1_kt", "cum40_tot_yld_kt", "t40_tot_yld_kt", "asymp_tot_yld_kt")],
      row.names = FALSE)
cat("\nFeedback HCR Summary:\n")
print(hcr_summary_df, row.names = FALSE)
if (sink.number() > 0) sink()
