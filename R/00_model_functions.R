# ==============================================================================
# R/00_model_functions.R
#
# Shared recruitment rate functions, canonical solver constants, policy scenario
# runners, and size-structured demographic diagnostics for the theoretical
# six-species cultivation-depensation model.
# ==============================================================================

# Canonical time-stepping configuration shared across all pipeline scripts
CANONICAL_DT     <- 0.05
CANONICAL_METHOD <- "tr_bdf2"
COLLAPSE_YEARS   <- 0:61   # t = 0..26 overfishing ramp; t = 27..61 34-yr cod closure
POLICY_YEARS     <- 0:40   # tau = t - 61 in [0, 40], switch at tau = 20

EFFORT_CAL  <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.22)
EFFORT_HIST <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.15)
EFFORT_MOR  <- c(Pelagic = 0.15, Mesopredator = 0.05, Flatfish = 0.15, Cod_Trawl = 0.00)

#' Beverton-Holt recruitment with exogenous recruitment input
#'
#' Adds a constant exogenous recruitment flux `species_params$r_refuge`
#' (individuals / year entering at `w_min`) to standard `mizer::BevertonHoltRDD()`.
#' This represents an open-population recruitment subsidy (e.g. external larval
#' import or unmodelled boundary input) rather than an explicit dynamically
#' sheltered juvenile subpopulation. When `r_refuge == 0`, zero abundance is
#' strictly invariant and the function reduces to `BevertonHoltRDD()`.
#'
#' @param rdi Numeric vector of density-independent egg production rates by species.
#' @param species_params Species parameter data frame of the MizerParams object.
#' @param params MizerParams object.
#' @param t Current simulation time (years).
#' @param ... Additional arguments passed by mizer's rate pipeline.
#' @return Plain named numeric vector of recruitment fluxes (numbers / year).
RefugeBevertonHoltRDD <- function(rdi, species_params, params, t, ...) {
  BevertonHoltRDD(rdi, species_params, params, t, ...) + species_params$r_refuge
}

#' Beverton-Holt recruitment with spawner-linked saturating refuge component
#'
#' Sensitivity-analysis alternative in which the supplemental recruitment term
#' depends on modelled egg production `rdi` via a saturating Michaelis-Menten /
#' Beverton-Holt relation `r_ref_max * rdi / (rdi + k_ref)`. Unlike a constant
#' additive input, `SaturatingShelterRDD(0, ...) == 0`, so zero cod abundance is
#' an exact invariant equilibrium.
SaturatingShelterRDD <- function(rdi, species_params, params, t, ...) {
  base_bh <- BevertonHoltRDD(rdi, species_params, params, t, ...)
  r_max_s <- species_params$r_ref_max
  k_ref_s <- species_params$k_ref
  extra   <- ifelse(r_max_s > 0 & rdi > 0, r_max_s * rdi / (rdi + k_ref_s), 0)
  base_bh + extra
}

#' Build canonical collapse + closure effort array (relative years t = 0..61)
build_collapse_closure_effort <- function(gears = names(EFFORT_HIST)) {
  t_vec <- COLLAPSE_YEARS
  em <- array(0, dim = c(length(t_vec), length(gears)),
              dimnames = list(time = t_vec, gear = gears))
  em[, "Pelagic"]      <- 0.15
  em[, "Mesopredator"] <- 0.05
  em[, "Flatfish"]     <- 0.15
  em[, "Cod_Trawl"]    <- ifelse(t_vec < 3, 0.15,
                          ifelse(t_vec < 16, 0.46,
                          ifelse(t_vec < 27, 0.56, 0.00)))
  em
}

#' Run a 4-gear two-phase policy scenario from a starting state (tau = 0..40)
#' using the shared canonical solver configuration
run_policy_4gear <- function(params_start,
                             fp1, fm1, fc1 = 0.00,
                             fp2 = fp1, fm2 = fm1, fc2 = fc1,
                             split_yr = 20L,
                             yrs = POLICY_YEARS,
                             dt = CANONICAL_DT,
                             t_save = 1,
                             method = CANONICAL_METHOD) {
  gears <- names(initial_effort(params_start))
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
  project(params_start, effort = em, dt = dt, t_save = t_save,
          method = method, progress_bar = FALSE)
}

#' Construct 6-gear MizerParams splitting Capelin_Gear from Clupeid_Mackerel
#' and adding a juvenile/sub-adult Cod_Bycatch gear (l50 = 22 cm, l25 = 16 cm)
make_6gear_params <- function(params_4gear) {
  gp6 <- data.frame(
    gear         = c("Capelin_Gear", "Clupeid_Mackerel", "Clupeid_Mackerel",
                     "Mesopredator", "Flatfish",         "Cod_Trawl", "Cod_Bycatch"),
    species      = c("Capelin",      "Herring",          "Mackerel",
                     "Sculpin",      "Plaice",           "Cod",       "Cod"),
    sel_func     = "sigmoid_length",
    l50          = c(11,             20,                 25,
                     20,             28,                 52,          22),
    l25          = c(9,              16,                 20,
                     16,             23,                 45,          16),
    catchability = c(0.75,           1.00,               1.00,
                     1.00,           1.00,               1.00,        1.00),
    stringsAsFactors = FALSE
  )
  p6 <- params_4gear
  gear_params(p6) <- gp6
  initial_effort(p6) <- c(
    Capelin_Gear     = 0.15,
    Clupeid_Mackerel = 0.15,
    Mesopredator     = 0.05,
    Flatfish         = 0.15,
    Cod_Trawl        = 0.00,
    Cod_Bycatch      = 0.00
  )
  p6
}

#' Run a 6-gear two-phase policy scenario (tau = 0..40) with optional asymptotic
#' extension under constant post-switch fishing effort
run_policy_6gear <- function(params6_start,
                             ecap1, ehm1, em1, ec1 = 0.00, ebyc1 = 0.00,
                             ecap2 = ecap1, ehm2 = ehm1, em2 = em1, ec2 = ec1, ebyc2 = 0.00,
                             split_yr = 20L,
                             yrs = POLICY_YEARS,
                             dt = CANONICAL_DT,
                             t_save = 1,
                             method = CANONICAL_METHOD,
                             asymp_extra_yrs = 300,
                             asymp_dt = CANONICAL_DT) {
  gn <- names(initial_effort(params6_start))
  em <- array(0, dim = c(length(yrs), length(gn)),
              dimnames = list(time = yrs, gear = gn))
  ph1 <- yrs < split_yr
  em[ph1,  "Capelin_Gear"]     <- ecap1; em[!ph1, "Capelin_Gear"]     <- ecap2
  em[ph1,  "Clupeid_Mackerel"] <- ehm1;  em[!ph1, "Clupeid_Mackerel"] <- ehm2
  em[ph1,  "Mesopredator"]     <- em1;   em[!ph1, "Mesopredator"]     <- em2
  em[,     "Flatfish"]         <- 0.15
  em[ph1,  "Cod_Trawl"]        <- ec1;   em[!ph1, "Cod_Trawl"]        <- ec2
  em[ph1,  "Cod_Bycatch"]      <- ebyc1; em[!ph1, "Cod_Bycatch"]      <- ebyc2

  sim40 <- project(params6_start, effort = em, dt = dt, t_save = t_save,
                   method = method, progress_bar = FALSE)
  eff_final <- em[nrow(em), ]
  sim_asymp <- if (asymp_extra_yrs > 0) {
    project(finalParams(sim40), effort = eff_final, t_max = asymp_extra_yrs,
            dt = asymp_dt, t_save = max(1, asymp_extra_yrs / 6),
            method = method, progress_bar = FALSE)
  } else {
    NULL
  }
  list(sim40 = sim40, sim_asymp = sim_asymp, effort_final = eff_final)
}

#' Decompose Cod predation mortality by predator species across consumer sizes
decompose_cod_pred_mort <- function(params) {
  pr_full <- unclass(getPredRate(params))
  c_idx   <- match(w(params), w_full(params))
  stopifnot(!anyNA(c_idx))
  theta_cod <- interaction_matrix(params)[, "Cod"]
  pm_by_pred <- sweep(pr_full[, c_idx, drop = FALSE], 1, theta_cod, "*")
  dimnames(pm_by_pred) <- list(predator = rownames(interaction_matrix(params)),
                               w = dimnames(initialN(params))[[2]])
  pm_tot <- unclass(getPredMort(params))["Cod", ]
  err <- max(abs(colSums(pm_by_pred) - pm_tot))
  if (err > 1e-9) {
    stop("Predator mortality decomposition mismatch: ", err)
  }
  pm_by_pred
}

#' Compute cumulative mortality-to-growth hazard and cohort survival for Cod
compute_cod_hazard <- function(params, w_lo = 0.001, w_hi = 15.0, effort = NULL) {
  p_eval <- params
  if (!is.null(effort)) initial_effort(p_eval) <- effort
  w_vec  <- w(p_eval)
  dw_vec <- dw(p_eval)
  mask   <- w_vec >= w_lo & w_vec < w_hi

  g_cod   <- unclass(getEGrowth(p_eval))["Cod", ]
  m_cod   <- unclass(getMort(p_eval))["Cod", ]
  ext_cod <- unclass(ext_mort(p_eval))["Cod", ]
  f_cod   <- unclass(getFMort(p_eval))["Cod", ]
  pm_pred <- decompose_cod_pred_mort(p_eval)

  H_total    <- sum((m_cod[mask] / g_cod[mask]) * dw_vec[mask])
  H_by_pred  <- rowSums(sweep(pm_pred[, mask, drop = FALSE], 2, dw_vec[mask] / g_cod[mask], "*"))
  H_ext      <- sum((ext_cod[mask] / g_cod[mask]) * dw_vec[mask])
  H_fish     <- sum((f_cod[mask] / g_cod[mask]) * dw_vec[mask])
  surv_cont  <- exp(-H_total)
  surv_disc  <- prod(g_cod[mask] / (g_cod[mask] + m_cod[mask] * dw_vec[mask]))

  list(
    w_lo       = w_lo,
    w_hi       = w_hi,
    H_total    = H_total,
    H_by_pred  = H_by_pred,
    H_ext      = H_ext,
    H_fish     = H_fish,
    surv_cont  = surv_cont,
    surv_disc  = surv_disc
  )
}

#' Exact discrete steady-state net reproductive ratio R0 for Cod (without r_refuge
#' or Beverton-Holt saturation) given size-specific growth, mortality, and ERepro
compute_cod_R0 <- function(g_vec, m_vec, er_vec, params) {
  w_vec      <- w(params)
  dw_vec     <- dw(params)
  sp         <- species_params(params)
  cod_row    <- match("Cod", sp$species)
  erepro_cod <- sp$erepro[cod_row]
  w_min_cod  <- sp$w_min[cod_row]

  n_unit <- numeric(length(w_vec))
  flux   <- 1.0
  for (j in seq_along(w_vec)) {
    den <- g_vec[j] + m_vec[j] * dw_vec[j]
    if (den <= 0) break
    n_unit[j] <- flux / den
    flux      <- g_vec[j] * n_unit[j]
    if (flux <= 0) break
  }
  0.5 * erepro_cod * sum(n_unit * er_vec * dw_vec) / w_min_cod
}
