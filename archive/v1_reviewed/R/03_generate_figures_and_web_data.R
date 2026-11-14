# ==============================================================================
# R/03_generate_figures_and_web_data.R
#
# Generates 6 publication-quality multi-panel figures (PDF + PNG) using mizer's
# built-in plotting functions and exports exact mizer simulation grids and
# steady-state diagnostics for the interactive policy website (docs/).
# ==============================================================================

library(mizer)
library(ggplot2)
library(patchwork)
library(jsonlite)
source("R/00_model_functions.R")

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("docs/data",     showWarnings = FALSE, recursive = TRUE)

params_hist    <- readParams("models/nw_atlantic_cod_params.rds")
params_high_f0 <- readParams("models/nw_atlantic_cod_high_f0_params.rds")
params_low_f0  <- readParams("models/nw_atlantic_cod_moratorium_params.rds")
res            <- readRDS("models/analysis_results.rds")

sim_history <- res$sim_history
params_2026 <- res$params_2026
scenarios   <- res$scenarios

theme_pub <- theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 9.2, colour = "#334e68"),
    legend.key.size  = unit(0.42, "cm")
  )

save_pub_fig <- function(plot_obj, stem, width = 10, height = 5.5) {
  ggsave(file.path("paper/figures", paste0(stem, ".pdf")), plot_obj,
         width = width, height = height, device = "pdf")
  ggsave(file.path("paper/figures", paste0(stem, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
}

# ------------------------------------------------------------------------------
# Figure 1: Historical Collapse (1965-1992) and Moratorium Trap (1992-2026)
#           Built with mizer's plotBiomass() and plotYieldGear()
# ------------------------------------------------------------------------------
cat("Generating Figure 1 (Historical Collapse & Moratorium)...\n")
p1_bio <- plotBiomass(sim_history, highlight = "Cod") +
  geom_vline(xintercept = 1992, linetype = "dashed", colour = "#b91c1c", linewidth = 0.6) +
  annotate("text", x = 1993.2, y = 2.5e10, label = "1992 Moratorium\n(F_cod = 0)",
           hjust = 0, vjust = 0, size = 3.0, colour = "#b91c1c", fontface = "bold") +
  labs(title = "A. Community Biomass Trajectories (1965-2026)",
       subtitle = "Overfishing collapses Cod; forage & sculpin surge locks Cod in post-1992 trap") +
  theme_pub

p1_yld <- plotYieldGear(sim_history) +
  geom_vline(xintercept = 1992, linetype = "dashed", colour = "#b91c1c", linewidth = 0.5) +
  labs(title = "B. Annual Fishery Landings by Fleet (1965-2026)",
       subtitle = "Directed Cod_Trawl halted in 1992 while baseline pelagic & flatfish fleets continue") +
  theme_pub

fig1 <- (p1_bio / p1_yld) + plot_layout(heights = c(1.05, 1.0))
save_pub_fig(fig1, "fig1_collapse_moratorium", width = 10.2, height = 7.8)

# ------------------------------------------------------------------------------
# Figure 2: Alternative Stable Size Spectra at F_cod = 0
#           Built with mizer's plotSpectra2() and plotSpectraRelative()
# ------------------------------------------------------------------------------
cat("Generating Figure 2 (Alternative Stable Size Spectra)...\n")
p2_spec <- plotSpectra2(
  params_high_f0, params_low_f0,
  name1 = "Cultivated (High-Cod)", name2 = "Moratorium (Forage-Locked)",
  biomass = TRUE, per_log_size = TRUE, resource = TRUE
) +
  labs(title = "A. Alternative Stable Biomass Spectra at Identical Parameters (F_cod = 0)",
       subtitle = "Comparison of the cultivated high-cod state and the collapsed moratorium state") +
  theme_pub

p2_rel <- plotSpectraRelative(params_high_f0, params_low_f0) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  labs(title = "B. Relative Spectral Shift: Moratorium vs. Cultivated State",
       subtitle = "Symmetric relative difference between moratorium and cultivated states") +
  theme_pub

fig2 <- p2_spec + p2_rel + plot_layout(ncol = 2)
save_pub_fig(fig2, "fig2_size_spectra_cascade", width = 11.8, height = 5.2)

# ------------------------------------------------------------------------------
# Figure 3: Trophic Mechanisms of Cultivation-Depensation
#           Built with mizer's plot2(), getPredMort(), getFeedingLevel(),
#           plotDiet(), and plotGrowthCurves()
# ------------------------------------------------------------------------------
cat("Generating Figure 3 (Trophic Mechanisms & Energetics)...\n")
p3_mort <- plot2(
  getPredMort(params_low_f0), getPredMort(params_high_f0),
  "Moratorium Trap", "Cultivated Regime",
  species = "Cod"
) +
  annotate("rect", xmin = 0.002, xmax = 15, ymin = -Inf, ymax = Inf,
           fill = "#f59e0b", alpha = 0.12) +
  labs(title = "A. Early-Life Predation Bottleneck on Cod",
       subtitle = "Predation mortality mu_p(w) across egg, larval, and 0-group sizes") +
  theme_pub

p3_feed <- plot2(
  getFeedingLevel(params_low_f0), getFeedingLevel(params_high_f0),
  "Moratorium Trap", "Cultivated Regime",
  species = "Cod"
) +
  annotate("rect", xmin = 0.002, xmax = 15, ymin = -Inf, ymax = Inf,
           fill = "#f59e0b", alpha = 0.12) +
  labs(title = "B. Zooplankton Competition & Feeding Level",
       subtitle = "Planktivore grazing depresses larval/juvenile Cod feeding level f(w)") +
  theme_pub

p3_diet <- plotDiet(params_hist, species = "Cod") +
  labs(title = "C. Ontogenetic Diet Shift of Atlantic Cod",
       subtitle = "Transition from plankton (<10 g) to forage fishes & sculpin (>200 g)") +
  theme_pub

p3_grow <- plotGrowthCurves(params_hist) +
  labs(title = "D. Calibrated Somatic Growth Curves",
       subtitle = "Modelled size-at-age trajectories across the 6 shelf species") +
  theme_pub

fig3 <- (p3_mort + p3_feed) / (p3_diet + p3_grow)
save_pub_fig(fig3, "fig3_trophic_mechanisms", width = 11.4, height = 8.4)

# ------------------------------------------------------------------------------
# Figure 4: Hysteresis & Fold Bifurcations (from mizer's scanModel)
# ------------------------------------------------------------------------------
cat("Generating Figure 4 (Hysteresis & Escape Bifurcations)...\n")
extract_scan_df <- function(scan_obj, branch_label) {
  d <- as.data.frame(scan_obj)
  data.frame(
    scan_value = as.numeric(d[[1]]),
    ssb_kt     = as.numeric(d[[2]]) / 1e9,
    branch     = branch_label,
    stringsAsFactors = FALSE
  )
}
df_cod_up   <- extract_scan_df(res$scan_cod_up,   "Forward: Increasing F_cod (Collapse)")
df_cod_down <- extract_scan_df(res$scan_cod_down, "Backward: Decreasing F_cod (Moratorium)")
df_cod_both <- rbind(df_cod_up, df_cod_down)

p4_cod <- ggplot(df_cod_both, aes(x = scan_value, y = ssb_kt, colour = branch, linetype = branch)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = res$collapse_f_cod, linetype = "dotted", colour = "#b91c1c") +
  annotate("text", x = res$collapse_f_cod + 0.012, y = 880,
           label = sprintf("Collapse Fold\n(F_cod = %.2f /yr)", res$collapse_f_cod),
           hjust = 0, size = 3.0, colour = "#b91c1c", fontface = "bold") +
  scale_colour_manual(values = c("Forward: Increasing F_cod (Collapse)" = "#1d5fa7",
                                 "Backward: Decreasing F_cod (Moratorium)" = "#c53929")) +
  scale_linetype_manual(values = c("Forward: Increasing F_cod (Collapse)" = "solid",
                                   "Backward: Decreasing F_cod (Moratorium)" = "dashed")) +
  labs(title = "A. Collapse & Moratorium Hysteresis vs. Directed Cod Effort",
       subtitle = "Reducing F_cod to zero moves along the lower collapsed branch without recovery",
       x = "Directed Cod Fishing Effort F_cod (1/year)",
       y = "Equilibrium Cod SSB (kt)", colour = NULL, linetype = NULL) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical")

df_pel_up   <- extract_scan_df(res$scan_pel_up,   "Forward: Increasing F_pelagic (From Trap)")
df_pel_down <- extract_scan_df(res$scan_pel_down, "Backward: Tapering F_pelagic (Rebuilt)")
df_pel_both <- rbind(df_pel_up, df_pel_down)

p4_pel <- ggplot(df_pel_both, aes(x = scan_value, y = ssb_kt, colour = branch, linetype = branch)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = res$recovery_f_pel, linetype = "dotted", colour = "#15805c") +
  annotate("text", x = res$recovery_f_pel + 0.015, y = 280,
           label = sprintf("Equilibrium Escape Fold\n(F_pel = %.3f /yr)", res$recovery_f_pel),
           hjust = 0, size = 3.0, colour = "#15805c", fontface = "bold") +
  scale_colour_manual(values = c("Forward: Increasing F_pelagic (From Trap)" = "#d97706",
                                 "Backward: Tapering F_pelagic (Rebuilt)" = "#15805c")) +
  scale_linetype_manual(values = c("Forward: Increasing F_pelagic (From Trap)" = "solid",
                                   "Backward: Tapering F_pelagic (Rebuilt)" = "dashed")) +
  labs(title = "B. Escape Bifurcation vs. Pelagic Forage Effort (at F_cod = 0)",
       subtitle = "Forage harvest triggers jump to upper branch; recovery persists when tapered",
       x = "Pelagic Forage Fishing Effort F_pelagic (1/year)",
       y = "Equilibrium Cod SSB (kt)", colour = NULL, linetype = NULL) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical")

fig4 <- p4_cod + p4_pel + plot_layout(ncol = 2)
save_pub_fig(fig4, "fig4_hysteresis_bifurcation", width = 11.4, height = 5.5)

# ------------------------------------------------------------------------------
# Figure 5: Ecosystem-Based Recovery Strategies (2026-2066)
#           Combining mizer's plotBiomass() panels and cross-scenario Cod SSB
# ------------------------------------------------------------------------------
cat("Generating Figure 5 (Ecosystem Recovery Strategies)...\n")
p5_sq <- plotBiomass(scenarios$status_quo, highlight = "Cod") +
  labs(title = "A. Status Quo Moratorium (F_pel=0.15, F_meso=0.05)",
       subtitle = "Cod remains trapped at ~15 kt SSB despite zero directed fishing") +
  theme_pub + theme(legend.position = "none")

p5_pel <- plotBiomass(scenarios$pelagic_rebalance, highlight = "Cod") +
  labs(title = "B. Pelagic Forage Harvest (F_pel=0.55, F_cod=0)",
       subtitle = "Thinning herring, mackerel & capelin releases larval Cod bottleneck") +
  theme_pub + theme(legend.position = "none")

p5_pul <- plotBiomass(scenarios$ecosystem_pulse, highlight = "Cod") +
  geom_vline(xintercept = 2046, linetype = "dashed", colour = "#1d5fa7", linewidth = 0.55) +
  annotate("text", x = 2046.8, y = 4.5e10, label = "2046 Taper +\nReopen F_cod=0.15",
           hjust = 0, vjust = 0, size = 2.8, colour = "#1d5fa7", fontface = "bold") +
  labs(title = "C. Adaptive Pulse (2026-2046) + Taper & Reopened Fishery",
       subtitle = "Cod rebuilds and sustains F_cod=0.15 /yr after forage effort returns to baseline") +
  theme_pub

scen_labels <- c(
  status_quo        = "1. Status Quo Moratorium",
  premature_reopen  = "2. Premature Cod Reopen (2046)",
  mesopred_only     = "3. Sculpin Harvest Only",
  pelagic_rebalance = "4. Pelagic Harvest Only",
  combined_sustain  = "5. Combined Harvest (Sustained)",
  ecosystem_pulse   = "6. Adaptive Pulse + Cod F=0.15"
)
cod_comp_df <- do.call(rbind, lapply(names(scen_labels), function(k) {
  ssb_mat <- unclass(getSSB(scenarios[[k]]))
  data.frame(
    year     = as.numeric(rownames(ssb_mat)),
    cod_ssb  = as.numeric(ssb_mat[, "Cod"]) / 1e9,
    strategy = scen_labels[[k]],
    stringsAsFactors = FALSE
  )
}))

p5_comp <- ggplot(cod_comp_df, aes(x = year, y = cod_ssb, colour = strategy)) +
  geom_line(linewidth = 1.05) +
  geom_hline(yintercept = res$allee_ssb_f0, linetype = "dotted", colour = "#b91c1c") +
  scale_colour_manual(values = c(
    "1. Status Quo Moratorium"       = "#64748b",
    "2. Premature Cod Reopen (2046)" = "#c53929",
    "3. Sculpin Harvest Only"        = "#7c3aed",
    "4. Pelagic Harvest Only"        = "#d97706",
    "5. Combined Harvest (Sustained)"= "#0d9488",
    "6. Adaptive Pulse + Cod F=0.15" = "#1d5fa7"
  )) +
  labs(title = "D. Atlantic Cod Spawning Stock Biomass Across Policy Strategies",
       subtitle = "Synergistic pelagic + mesopredator harvest accelerates crossing of the cultivation threshold",
       x = "Year", y = "Cod SSB (kt)", colour = "Strategy") +
  theme_pub

fig5 <- (p5_sq + p5_pel) / (p5_pul + p5_comp)
save_pub_fig(fig5, "fig5_recovery_strategies", width = 11.8, height = 8.4)

# ------------------------------------------------------------------------------
# 6. Precompute 3D Policy Scenario Grid for Figure 6 & Interactive Website
# ------------------------------------------------------------------------------
cat("Computing 3D mizer scenario grid (sustained & adaptive taper) for Figure 6 & Web Simulator...\n")
fp_vals   <- c(0.00, 0.15, 0.35, 0.50, 0.55, 0.65, 0.80)
fm_vals   <- c(0.00, 0.05, 0.20, 0.35, 0.45, 0.60)
fcod_vals <- c(0.00, 0.10, 0.15, 0.25, 0.35, 0.50)

proj_years_full <- 2026:2066
keep_yrs        <- seq(2026, 2066, by = 2)
gears           <- names(initial_effort(params_hist))
w_vec           <- w(params_hist)
nursery_mask    <- w_vec >= 0.05 & w_vec <= 5.0

run_grid_point <- function(fp, fm, fcod, taper) {
  em <- array(0, dim = c(length(proj_years_full), length(gears)),
              dimnames = list(time = proj_years_full, gear = gears))
  phase1 <- proj_years_full < 2046
  em[phase1,  "Pelagic"]      <- fp
  em[!phase1, "Pelagic"]      <- if (taper) 0.15 else fp
  em[phase1,  "Mesopredator"] <- fm
  em[!phase1, "Mesopredator"] <- if (taper) 0.05 else fm
  em[,        "Flatfish"]     <- 0.15
  em[phase1,  "Cod_Trawl"]    <- 0.00
  em[!phase1, "Cod_Trawl"]    <- fcod

  sim <- project(params_2026, effort = em, dt = 0.2, t_save = 2, progress_bar = FALSE)
  ssb <- unclass(getSSB(sim)) / 1e9
  bio <- unclass(getBiomass(sim, use_cutoff = TRUE)) / 1e9
  yld <- unclass(getYield(sim)) / 1e9
  p_end <- finalParams(sim)
  lm_end <- mean(unclass(getPredMort(p_end))["Cod", nursery_mask])

  list(
    cod_ssb         = round(as.numeric(ssb[, "Cod"]), 1),
    forage_bio      = round(as.numeric(rowSums(bio[, c("Capelin", "Herring", "Mackerel")])), 1),
    sculpin_bio     = round(as.numeric(bio[, "Sculpin"]), 1),
    cod_yield       = round(as.numeric(yld[, "Cod"]), 1),
    pelagic_yield   = round(as.numeric(rowSums(yld[, c("Capelin", "Herring", "Mackerel")])), 1),
    other_yield     = round(as.numeric(rowSums(yld[, c("Sculpin", "Plaice")])), 1),
    larval_mort_end = round(lm_end, 3)
  )
}

grid_sustained <- list()
grid_taper     <- list()
surface_rows   <- list()

for (i in seq_along(fp_vals)) {
  for (j in seq_along(fm_vals)) {
    for (k in seq_along(fcod_vals)) {
      key <- sprintf("%d_%d_%d", i - 1L, j - 1L, k - 1L)
      g_sus <- run_grid_point(fp_vals[i], fm_vals[j], fcod_vals[k], taper = FALSE)
      g_tap <- run_grid_point(fp_vals[i], fm_vals[j], fcod_vals[k], taper = TRUE)
      grid_sustained[[key]] <- g_sus
      grid_taper[[key]]     <- g_tap

      # Record 2D surface slice for Figure 6 (Adaptive Taper + Reopened F_cod = 0.15)
      if (abs(fcod_vals[k] - 0.15) < 1e-6) {
        n_last <- length(g_tap$cod_ssb)
        surface_rows[[length(surface_rows) + 1L]] <- data.frame(
          fp         = fp_vals[i],
          fm         = fm_vals[j],
          cod_ssb_kt = g_tap$cod_ssb[n_last],
          tot_yld_kt = g_tap$cod_yield[n_last] + g_tap$pelagic_yield[n_last] + g_tap$other_yield[n_last]
        )
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Figure 6: 2D Policy Response Surface (Cod Recovery & Multispecies Yield)
# ------------------------------------------------------------------------------
cat("Generating Figure 6 (2D Policy Trade-off Surface)...\n")
surf_df <- do.call(rbind, surface_rows)
surf_df$fp_fac <- factor(sprintf("%.2f", surf_df$fp), levels = sprintf("%.2f", fp_vals))
surf_df$fm_fac <- factor(sprintf("%.2f", surf_df$fm), levels = sprintf("%.2f", fm_vals))

p6_ssb <- ggplot(surf_df, aes(x = fp_fac, y = fm_fac, fill = cod_ssb_kt)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(cod_ssb_kt)), size = 3.1, fontface = "bold",
            colour = ifelse(surf_df$cod_ssb_kt > 360, "white", "#10243e")) +
  scale_fill_gradient2(low = "#fcebe9", mid = "#93c5fd", high = "#1d5fa7", midpoint = 360,
                       name = "Cod SSB\n(kt)") +
  labs(title = "A. Year-2066 Cod Spawning Biomass (kt)",
       subtitle = "20-yr forage/mesopredator pulse (2026-2046), then tapered with F_cod = 0.15 /yr",
       x = "Phase-I Pelagic Forage Effort F_pelagic (1/year)",
       y = "Phase-I Sculpin Effort F_mesopredator (1/year)") +
  theme_pub

p6_yld <- ggplot(surf_df, aes(x = fp_fac, y = fm_fac, fill = tot_yld_kt)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(tot_yld_kt)), size = 3.1, fontface = "bold",
            colour = ifelse(surf_df$tot_yld_kt > 520, "white", "#10243e")) +
  scale_fill_gradient(low = "#fef3c7", high = "#15805c", name = "Total Yield\n(kt/yr)") +
  labs(title = "B. Year-2066 Sustainable Multispecies Yield (kt/yr)",
       subtitle = "Rebuilt Cod adds groundfish catch alongside balanced pelagic yield",
       x = "Phase-I Pelagic Forage Effort F_pelagic (1/year)",
       y = "Phase-I Sculpin Effort F_mesopredator (1/year)") +
  theme_pub

fig6 <- p6_ssb + p6_yld + plot_layout(ncol = 2)
save_pub_fig(fig6, "fig6_policy_tradeoffs", width = 11.4, height = 5.1)

# ------------------------------------------------------------------------------
# 7. Assemble and Write Interactive Website Payload (docs/model_data.js & JSON)
# ------------------------------------------------------------------------------
cat("Exporting docs/model_data.js and docs/data/mizer_scenarios.json...\n")
hist_years_keep <- seq(1965, 2025, by = 2)
if (tail(hist_years_keep, 1) != 2026) hist_years_keep <- c(hist_years_keep, 2026)
h_idx <- match(as.character(hist_years_keep), rownames(unclass(getSSB(sim_history))))

h_ssb <- unclass(getSSB(sim_history))[h_idx, , drop = FALSE] / 1e9
h_bio <- unclass(getBiomass(sim_history, use_cutoff = TRUE))[h_idx, , drop = FALSE] / 1e9
h_yld <- unclass(getYield(sim_history))[h_idx, , drop = FALSE] / 1e9

extract_spectra_diag <- function(p) {
  n_mat <- unclass(initialN(p))
  w2    <- w_vec^2
  to_log10 <- function(vec) {
    val <- as.numeric(vec) * w2
    ifelse(val >= 1e7, round(log10(val), 3), NA_real_)
  }
  list(
    spectra_log10 = list(
      Cod     = to_log10(n_mat["Cod", ]),
      Forage  = to_log10(colSums(n_mat[c("Capelin", "Herring", "Mackerel"), ])),
      Sculpin = to_log10(n_mat["Sculpin", ])
    ),
    cod_pred_mort = round(as.numeric(unclass(getPredMort(p))["Cod", ]), 4),
    cod_feeding   = round(as.numeric(unclass(getFeedingLevel(p))["Cod", ]), 4)
  )
}

precollapse_ssb  <- round(unname(getSSB(params_hist)["Cod"] / 1e9), 1)
collapse_1992    <- round(unname(unclass(getSSB(sim_history))["1992", "Cod"] / 1e9), 1)
collapse_pct_val <- round(100 * (1 - collapse_1992 / precollapse_ssb), 0)
forage_hist_bio  <- sum(getBiomass(params_hist, use_cutoff = TRUE)[1:4]) / 1e9
forage_low_bio   <- sum(getBiomass(params_low_f0, use_cutoff = TRUE)[1:4]) / 1e9
surge_ratio_val  <- round(forage_low_bio / forage_hist_bio, 1)
mort_ratio_val   <- round(res$larval_mort_low / res$larval_mort_high, 1)

make_policy_row <- function(name, fcod_txt, fpel_txt, fmeso_txt, sim_obj, outcome_txt, rec_bool) {
  p_end <- finalParams(sim_obj)
  c_ssb <- round(unname(getSSB(p_end)["Cod"] / 1e9), 0)
  t_yld <- round(sum(getYield(p_end)) / 1e9, 0)
  list(
    name           = name,
    fcod           = fcod_txt,
    fpel           = fpel_txt,
    fmeso          = fmeso_txt,
    cod_ssb_kt     = c_ssb,
    cod_pct        = round(100 * c_ssb / precollapse_ssb, 0),
    total_yield_kt = t_yld,
    outcome        = outcome_txt,
    recovered      = rec_bool
  )
}

web_data <- list(
  hero = list(
    collapse_pct      = paste0("-", collapse_pct_val, "%"),
    forage_surge      = paste0(surge_ratio_val, "&times;"),
    larval_mort_ratio = paste0(mort_ratio_val, "&times;"),
    recovery_window   = "14&ndash;18 yrs"
  ),
  policy_start_year = 2026L,
  reference = list(
    precollapse_cod_ssb = precollapse_ssb,
    allee_threshold_ssb = round(res$allee_ssb_f15, 1)
  ),
  presets = list(
    status_quo        = list(fp = 0.15, fm = 0.05, fcod = 0.00, taper = FALSE),
    pelagic_rebalance = list(fp = 0.55, fm = 0.05, fcod = 0.00, taper = FALSE),
    ecosystem_pulse   = list(fp = 0.65, fm = 0.45, fcod = 0.15, taper = TRUE),
    premature_reopen  = list(fp = 0.15, fm = 0.05, fcod = 0.25, taper = FALSE)
  ),
  history = list(
    years         = as.integer(hist_years_keep),
    cod_ssb       = round(as.numeric(h_ssb[, "Cod"]), 1),
    forage_bio    = round(as.numeric(rowSums(h_bio[, c("Capelin", "Herring", "Mackerel")])), 1),
    sculpin_bio   = round(as.numeric(h_bio[, "Sculpin"]), 1),
    cod_yield     = round(as.numeric(h_yld[, "Cod"]), 1),
    pelagic_yield = round(as.numeric(rowSums(h_yld[, c("Capelin", "Herring", "Mackerel")])), 1),
    other_yield   = round(as.numeric(rowSums(h_yld[, c("Sculpin", "Plaice")])), 1)
  ),
  proj_years     = as.integer(keep_yrs),
  grid_axes      = list(fp = fp_vals, fm = fm_vals, fcod = fcod_vals),
  grid_sustained = grid_sustained,
  grid_taper     = grid_taper,
  size_diagnostics = list(
    log10_w    = round(log10(w_vec), 4),
    moratorium = extract_spectra_diag(params_low_f0),
    cultivated = extract_spectra_diag(params_high_f0)
  ),
  bifurcation = list(
    vs_fcod = list(
      f_up       = round(df_cod_up$scan_value, 4),
      ssb_up     = round(df_cod_up$ssb_kt, 1),
      f_down     = round(rev(df_cod_down$scan_value), 4),
      ssb_down   = round(rev(df_cod_down$ssb_kt), 1),
      collapse_f = round(res$collapse_f_cod, 3)
    ),
    vs_fpelagic = list(
      f_up       = round(df_pel_up$scan_value, 4),
      ssb_up     = round(df_pel_up$ssb_kt, 1),
      f_down     = round(rev(df_pel_down$scan_value), 4),
      ssb_down   = round(rev(df_pel_down$ssb_kt), 1),
      recovery_f = round(res$recovery_f_pel, 3)
    )
  ),
  policy_table = list(
    make_policy_row("1. Status Quo Moratorium",         "0.00 /yr",           "0.15 /yr",           "0.05 /yr",           scenarios$status_quo,        "Trapped in Low-Cod State", FALSE),
    make_policy_row("2. Premature Cod Reopening",       "0.00 -> 0.25 /yr",   "0.15 /yr",           "0.05 /yr",           scenarios$premature_reopen,  "Severely Depleted",        FALSE),
    make_policy_row("3. Mesopredator Harvest Only",     "0.00 /yr",           "0.15 /yr",           "0.45 /yr",           scenarios$mesopred_only,     "Partial / Insufficient",   FALSE),
    make_policy_row("4. Pelagic Forage Harvest",        "0.00 /yr",           "0.55 /yr",           "0.05 /yr",           scenarios$pelagic_rebalance, "Rebuilds to Cultivated",   TRUE),
    make_policy_row("5. Combined Harvest (Sustained)",  "0.00 /yr",           "0.55 /yr",           "0.35 /yr",           scenarios$combined_sustain,  "Rapid Unfished Recovery",  TRUE),
    make_policy_row("6. Adaptive Pulse + Cod Fishery",  "0.00 -> 0.15 /yr",   "0.65 -> 0.15 /yr",   "0.45 -> 0.05 /yr",   scenarios$ecosystem_pulse,   "Self-Sustaining Fishery",  TRUE)
  )
)

json_str <- toJSON(web_data, auto_unbox = TRUE, null = "null", na = "null", digits = 4)
writeLines(json_str, "docs/data/mizer_scenarios.json")
writeLines(paste0("window.MIZER_MODEL_DATA = ", json_str, ";"), "docs/model_data.js")
cat("Successfully generated all 6 figures and exported docs/model_data.js!\n")
