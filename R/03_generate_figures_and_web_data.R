# ==============================================================================
# R/03_generate_figures_and_web_data.R
#
# Generates publication-grade multi-panel figures (fig1-fig6, PDF + PNG) and
# exports exact precomputed mizer scenario grids (method = "tr_bdf2", dt = 0.05)
# with verified long-term asymptotic outcomes for the educational web explorer
# (docs/model_data.js and docs/data/mizer_scenarios.json).
# ==============================================================================

library(mizer)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(parallel)
source("R/00_model_functions.R")

dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("docs/data",     showWarnings = FALSE, recursive = TRUE)
dir.create("models/logs",   showWarnings = FALSE, recursive = TRUE)
sink("models/logs/03_generate_figures_and_web_data.log", split = TRUE)
NC <- min(8L, detectCores())

params_f15     <- readParams("models/cod_f15_params.rds")
params_high_f0 <- readParams("models/cod_high_f0_params.rds")
params_low_f0  <- readParams("models/cod_low_f0_params.rds")
params_t61     <- readParams("models/cod_t61_params.rds")
sim_collapse   <- readSim("models/sim_collapse_closure.rds")
res            <- readRDS("models/analysis_results.rds")

sim_status_quo     <- readSim("models/sims/sim_status_quo.rds")
sim_pelagic_reb    <- readSim("models/sims/sim_pelagic_rebalance.rds")
sim_joint_taper15  <- readSim("models/sims/sim_joint_taper_f15.rds")
sim_spare_cap_base <- readSim("models/sims/sim_spare_cap_base.rds")
sim_joint_sust15   <- readSim("models/sims/sim_joint_sustain_f15.rds")
sim_mesopred_only  <- readSim("models/sims/sim_mesopred_only.rds")
sim_premature      <- readSim("models/sims/sim_premature_reopen.rds")

theme_pub <- theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 10.8),
    plot.subtitle    = element_text(size = 8.9, colour = "#334e68"),
    legend.key.size  = unit(0.40, "cm"),
    legend.text      = element_text(size = 8.8),
    legend.title     = element_text(size = 9.2, face = "bold")
  )

save_pub_fig <- function(plot_obj, stem, width = 10.5, height = 5.6) {
  ggsave(file.path("paper/figures", paste0(stem, ".pdf")), plot_obj,
         width = width, height = height, device = "pdf")
  ggsave(file.path("paper/figures", paste0(stem, ".png")), plot_obj,
         width = width, height = height, dpi = 300)
}

kt_log_labels <- function(x) {
  vals <- x / 1e9
  ifelse(vals >= 1, format(round(vals, 1), trim = TRUE, scientific = FALSE),
         format(signif(vals, 2), trim = TRUE, scientific = FALSE))
}

# ------------------------------------------------------------------------------
# Figure 1: Collapse Ramp (t = 0..26) & 34-Year Directed-Cod Closure (t = 27..61)
#           Consistent use_cutoff = TRUE, kt axes, relative years, single legend
# ------------------------------------------------------------------------------
cat("1. Generating Figure 1 (Collapse & Closure Trajectory in Relative Years)...\n")
flush.console()

p1_bio <- suppressMessages(
  plotBiomass(sim_collapse, use_cutoff = TRUE, highlight = "Cod") +
    scale_x_continuous(name = "Simulation time t (years)") +
    scale_y_log10(name = "Cutoff biomass (kt)", labels = kt_log_labels) +
    geom_vline(xintercept = 27, linetype = "dashed", colour = "#b91c1c", linewidth = 0.6) +
    annotate("text", x = 28.2, y = 2.5e10,
             label = "Closure at t = 27 yr\n(E_cod = 0)",
             hjust = 0, vjust = 0, size = 3.0, colour = "#b91c1c", fontface = "bold") +
    labs(title = "A. Community Cutoff Biomass Trajectories (t = 0 to 61 yr)",
         subtitle = "Overfishing ramp (t = 0-26) collapses Cod; forage and sculpin release persists through 34-yr closure",
         colour = "Species", linetype = "Species", linewidth = "Species", size = "Species") +
    theme_pub
)

p1_yld <- suppressMessages(
  plotYieldGear(sim_collapse) +
    scale_x_continuous(name = "Simulation time t (years)") +
    scale_y_continuous(name = "Annual yield (kt / yr)",
                       labels = function(x) format(round(x / 1e9, 0), trim = TRUE)) +
    geom_vline(xintercept = 27, linetype = "dashed", colour = "#b91c1c", linewidth = 0.5) +
    guides(colour = "none", linewidth = "none", size = "none") +
    labs(title = "B. Annual Fishery Yield by Fleet (t = 0 to 61 yr)",
         subtitle = "Directed Cod_Trawl effort ceases at t = 27 yr while baseline pelagic, sculpin, and flatfish fleets continue",
         linetype = "Fleet") +
    theme_pub +
    theme(legend.key.width = unit(0.95, "cm"))
)

fig1 <- (p1_bio / p1_yld) + plot_layout(heights = c(1.08, 1.0))
save_pub_fig(fig1, "fig1_collapse_moratorium", width = 10.6, height = 7.8)

# ------------------------------------------------------------------------------
# Figure 2: Alternative Stable Biomass Size Spectra at E_cod = 0
# ------------------------------------------------------------------------------
cat("2. Generating Figure 2 (Alternative Stable Size Spectra)...\n")
flush.console()

p2_spec <- suppressMessages(
  plotSpectra2(
    params_high_f0, params_low_f0,
    name1 = "Cultivated (High-Cod)", name2 = "Depensatory (Low-Cod)",
    biomass = TRUE, per_log_size = TRUE, resource = TRUE
  ) +
    scale_x_log10(name = "Body weight w (g)") +
    labs(title = "A. Coexisting Equilibrium Biomass Spectra (E_cod = 0)",
         subtitle = "Identical parameters and baseline non-cod fishing; solid = Cultivated, dashed = Depensatory",
         y = "Biomass density per log weight (g)",
         colour = "Component") +
    theme_pub +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.key.width = unit(0.85, "cm"))
)

p2_rel <- suppressMessages(
  plotSpectraRelative(params_high_f0, params_low_f0) +
    scale_x_log10(name = "Body weight w (g)") +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    guides(colour = "none") +
    labs(title = "B. Symmetric Relative Spectral Difference",
         subtitle = "2 (N_low - N_high) / (N_low + N_high) across individual body weight",
         y = "Symmetric relative difference") +
    theme_pub
)

fig2 <- p2_spec + p2_rel + plot_layout(ncol = 2)
save_pub_fig(fig2, "fig2_size_spectra_cascade", width = 11.8, height = 5.6)

# ------------------------------------------------------------------------------
# Figure 3: Predator-Resolved Mortality, Feeding Ratio, Cumulative Hazard, & Diet
# ------------------------------------------------------------------------------
cat("3. Generating Figure 3 (Predator-Resolved Bottleneck & Cumulative Hazard)...\n")
flush.console()

w_vec  <- res$w_vec
dw_vec <- res$dw_vec
nursery_plot_mask <- w_vec >= 1e-3 & w_vec <= 100

pred_species <- c("Mackerel", "Herring", "Sculpin", "Capelin", "Plaice")
pm_low_sub   <- res$pm_pred_low[pred_species, nursery_plot_mask, drop = FALSE]
df_pred_mort <- do.call(rbind, lapply(pred_species, function(sp_nm) {
  data.frame(
    w        = w_vec[nursery_plot_mask],
    mort     = as.numeric(pm_low_sub[sp_nm, ]),
    source   = paste0(sp_nm, " (low-cod)"),
    stringsAsFactors = FALSE
  )
}))
df_tot_lines <- rbind(
  data.frame(w = w_vec[nursery_plot_mask],
             mort = colSums(res$pm_pred_low[, nursery_plot_mask]),
             regime = "Total predation (Depensatory)"),
  data.frame(w = w_vec[nursery_plot_mask],
             mort = colSums(res$pm_pred_high[, nursery_plot_mask]),
             regime = "Total predation (Cultivated)")
)

pred_cols <- c(
  "Mackerel (low-cod)" = "#b25806",
  "Herring (low-cod)"  = "#d95f02",
  "Sculpin (low-cod)"  = "#7c3aed",
  "Capelin (low-cod)"  = "#e69f00",
  "Plaice (low-cod)"   = "#15805c"
)

p3_mort <- ggplot() +
  annotate("rect", xmin = 0.001, xmax = 15, ymin = -Inf, ymax = Inf,
           fill = "#f59e0b", alpha = 0.10) +
  geom_area(data = df_pred_mort, aes(x = w, y = mort, fill = source), alpha = 0.82) +
  geom_line(data = subset(df_tot_lines, regime == "Total predation (Cultivated)"),
            aes(x = w, y = mort, linetype = regime), colour = "#0f172a", linewidth = 0.95) +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100),
                labels = c("1 mg", "10 mg", "0.1 g", "1 g", "10 g", "100 g")) +
  scale_fill_manual(values = pred_cols, name = "Predator (depensatory)") +
  scale_linetype_manual(values = c("Total predation (Cultivated)" = "dashed"), name = NULL) +
  labs(title = "A. Predator-Resolved Mortality on Early-Life Cod",
       subtitle = "Stacked predator contributions in depensatory state vs. cultivated total (dashed)",
       x = "Cod body weight w (g)",
       y = "Predation mortality (1 / yr)") +
  theme_pub +
  theme(legend.key.width = unit(0.85, "cm"))

p3_feed <- suppressMessages(
  plot2(
    getFeedingLevel(params_low_f0), getFeedingLevel(params_high_f0),
    "Depensatory (Low-Cod)", "Cultivated (High-Cod)",
    species = "Cod"
  ) +
    scale_x_log10(name = "Cod body weight w (g)") +
    annotate("rect", xmin = 0.001, xmax = 15, ymin = -Inf, ymax = Inf,
             fill = "#f59e0b", alpha = 0.10) +
    guides(colour = "none") +
    labs(title = "B. Size-Resolved Cod Feeding Level f(w)",
         subtitle = "1.4-15.0% lower below 15 g (plankton competition), but higher above 200 g",
         y = "Feeding level f(w)",
         linetype = "Equilibrium") +
    theme_pub +
    theme(legend.key.width = unit(0.90, "cm"))
)

# Cumulative hazard curves H(w) up to maturity (850 g)
cum_H_curve <- function(p_obj, eff_v, label) {
  p_ev <- p_obj; initial_effort(p_ev) <- eff_v
  g_c <- unclass(getEGrowth(p_ev))["Cod", ]
  m_c <- unclass(getMort(p_ev))["Cod", ]
  mk  <- w_vec <= 850
  H_w <- cumsum((m_c[mk] / g_c[mk]) * dw_vec[mk])
  data.frame(w = w_vec[mk], H = H_w, state = label, stringsAsFactors = FALSE)
}
df_H <- rbind(
  cum_H_curve(params_high_f0, EFFORT_MOR,  "Cultivated (E_cod = 0)"),
  cum_H_curve(params_f15,     EFFORT_HIST, "Baseline Fishery (E_cod = 0.15)"),
  cum_H_curve(params_low_f0,  EFFORT_MOR,  "Depensatory (E_cod = 0)")
)

p3_haz <- ggplot(df_H, aes(x = w, y = H, colour = state, linetype = state)) +
  annotate("rect", xmin = 0.001, xmax = 15, ymin = -Inf, ymax = Inf,
           fill = "#f59e0b", alpha = 0.10) +
  geom_line(linewidth = 1.1) +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 850),
                labels = c("1 mg", "10 mg", "0.1 g", "1 g", "10 g", "100 g", "850 g")) +
  scale_colour_manual(values = c(
    "Cultivated (E_cod = 0)"          = "#15805c",
    "Baseline Fishery (E_cod = 0.15)" = "#1d5fa7",
    "Depensatory (E_cod = 0)"         = "#c53929"
  ), name = "Equilibrium") +
  scale_linetype_manual(values = c(
    "Cultivated (E_cod = 0)"          = "solid",
    "Baseline Fishery (E_cod = 0.15)" = "dashed",
    "Depensatory (E_cod = 0)"         = "solid"
  ), name = "Equilibrium") +
  labs(title = "C. Cumulative Cohort Hazard Integral up to Maturity",
       subtitle = "H(w) = integral of (mu / g) dw; survival to 15 g drops 14-fold (83.1% to 6.0%)",
       x = "Cod body weight w (g)",
       y = "Cumulative hazard H(w)") +
  theme_pub +
  theme(legend.key.width = unit(0.90, "cm"))

p3_diet <- suppressMessages(
  plotDiet(params_f15, species = "Cod") +
    scale_x_log10(name = "Cod body weight w (g)") +
    labs(title = "D. Ontogenetic Diet Shift of Cod (Baseline E_cod = 0.15)",
         subtitle = "Transition from shared resource (< 10 g) to forage fishes and sculpin (> 200 g)",
         y = "Diet proportion") +
    theme_pub
)

fig3 <- (p3_mort + p3_feed) / (p3_haz + p3_diet)
save_pub_fig(fig3, "fig3_trophic_mechanisms", width = 12.2, height = 8.8)

# ------------------------------------------------------------------------------
# Figure 4: Hysteresis Scans, Newton Unstable Saddles, & Escape-Time Scaling
# ------------------------------------------------------------------------------
cat("4. Generating Figure 4 (Hysteresis Loops, Saddles, & Escape Times)...\n")
flush.console()

merge_fine_scan <- function(coarse_df, fine_df, eff_col) {
  f_min <- min(fine_df[[eff_col]])
  f_max <- max(fine_df[[eff_col]])
  keep_coarse <- coarse_df[coarse_df$effort_val < f_min - 1e-9 | coarse_df$effort_val > f_max + 1e-9,
                           c("effort_val", "cod_ssb_kt")]
  fine_std <- data.frame(effort_val = fine_df[[eff_col]], cod_ssb_kt = fine_df$cod_ssb_kt)
  out <- rbind(keep_coarse, fine_std)
  out[order(out$effort_val), , drop = FALSE]
}

cod_up_merged <- merge_fine_scan(res$scan_cod_up_df, res$fine_cod_df, "E_cod")
pel_up_merged <- merge_fine_scan(res$scan_pel_up_df, res$fine_pel_df, "E_pel")

df_cod_both <- rbind(
  transform(cod_up_merged, branch = "Forward: increasing E_cod (from high state)"),
  transform(res$scan_cod_down_df[, c("effort_val", "cod_ssb_kt")],
            branch = "Backward: decreasing E_cod (from low state)")
)
sad_pts <- data.frame(
  E_cod      = c(0.00, 0.15),
  cod_ssb_kt = res$saddle_summary_df$saddle_ssb_kt
)

p4_cod <- ggplot(df_cod_both, aes(x = effort_val, y = cod_ssb_kt, colour = branch, linetype = branch)) +
  annotate("rect", xmin = res$fine_cod_bracket[1], xmax = res$fine_cod_bracket[2],
           ymin = -Inf, ymax = Inf, fill = "#b91c1c", alpha = 0.18) +
  geom_line(linewidth = 1.05) +
  geom_point(data = sad_pts, aes(x = E_cod, y = cod_ssb_kt),
             inherit.aes = FALSE, shape = 23, size = 3.2, stroke = 1.0,
             fill = "#fef08a", colour = "#854d0e") +
  annotate("text", x = res$fine_cod_bracket[2] + 0.010, y = 850,
           label = sprintf("Transition bracket\n[%.3f, %.3f] /yr",
                           res$fine_cod_bracket[1], res$fine_cod_bracket[2]),
           hjust = 0, size = 2.85, colour = "#b91c1c", fontface = "bold") +
  annotate("text", x = 0.015, y = 395,
           label = "Newton unstable\nsaddles (diamonds)",
           hjust = 0, size = 2.75, colour = "#854d0e", fontface = "bold") +
  scale_colour_manual(values = c(
    "Forward: increasing E_cod (from high state)" = "#1d5fa7",
    "Backward: decreasing E_cod (from low state)" = "#c53929"
  )) +
  scale_linetype_manual(values = c(
    "Forward: increasing E_cod (from high state)" = "solid",
    "Backward: decreasing E_cod (from low state)" = "dashed"
  )) +
  labs(title = "A. Hysteresis vs. Directed Cod Effort E_cod",
       subtitle = "Shaded: 0.005 /yr bracket; diamonds: Newton saddles",
       x = "Directed cod effort E_cod (1 / yr)",
       y = "Settled Cod SSB (kt)", colour = NULL, linetype = NULL) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.width = unit(0.95, "cm"))

df_pel_both <- rbind(
  transform(pel_up_merged, branch = "Forward: increasing E_pel (from low state)"),
  transform(res$scan_pel_down_df[, c("effort_val", "cod_ssb_kt")],
            branch = "Backward: decreasing E_pel (from high state)")
)

p4_pel <- ggplot(df_pel_both, aes(x = effort_val, y = cod_ssb_kt, colour = branch, linetype = branch)) +
  annotate("rect", xmin = res$fine_pel_bracket[1], xmax = res$fine_pel_bracket[2],
           ymin = -Inf, ymax = Inf, fill = "#15805c", alpha = 0.18) +
  geom_line(linewidth = 1.05) +
  annotate("text", x = res$fine_pel_bracket[2] + 0.015, y = 260,
           label = sprintf("Escape bracket\n[%.3f, %.3f] /yr",
                           res$fine_pel_bracket[1], res$fine_pel_bracket[2]),
           hjust = 0, size = 2.85, colour = "#15805c", fontface = "bold") +
  scale_colour_manual(values = c(
    "Forward: increasing E_pel (from low state)"  = "#d97706",
    "Backward: decreasing E_pel (from high state)"= "#15805c"
  )) +
  scale_linetype_manual(values = c(
    "Forward: increasing E_pel (from low state)"  = "solid",
    "Backward: decreasing E_pel (from high state)"= "dashed"
  )) +
  labs(title = "B. Escape Hysteresis vs. Pelagic Effort (E_cod = 0)",
       subtitle = "High-cod state persists when E_pel returns to baseline 0.15 /yr",
       x = "Pelagic effort E_pel (1 / yr)",
       y = "Settled Cod SSB (kt)", colour = NULL, linetype = NULL) +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.width = unit(0.95, "cm"))

esc_valid <- subset(res$escape_time_df, !is.na(t_to_500kt_yr))
p4_esc <- ggplot(esc_valid, aes(x = E_pel, y = t_to_500kt_yr)) +
  geom_line(colour = "#1d5fa7", linewidth = 1.05) +
  geom_point(colour = "#1d5fa7", size = 2.4) +
  geom_hline(yintercept = 40, linetype = "dotted", colour = "#64748b") +
  annotate("text", x = 0.64, y = 55, label = "40-yr simulation window",
           hjust = 1, size = 2.8, colour = "#475569") +
  labs(title = "C. Direct Escape Time from Low State",
       subtitle = "Years at constant E_pel (E_cod = 0) to reach 500 kt Cod SSB",
       x = "Pelagic effort E_pel (1 / yr)",
       y = "Time to 500 kt Cod SSB (yr)") +
  theme_pub

fig4 <- p4_cod + p4_pel + p4_esc + plot_layout(ncol = 3, widths = c(1.05, 1.05, 0.9))
save_pub_fig(fig4, "fig4_hysteresis_bifurcation", width = 13.4, height = 5.3)

# ------------------------------------------------------------------------------
# Figure 5: Matched Rebuilding Comparisons & Phase-I Bycatch Thresholds
# ------------------------------------------------------------------------------
cat("5. Generating Figure 5 (Matched Policies, Capelin Sparing, HCR, & Bycatch)...\n")
flush.console()

p5_sq <- suppressMessages(
  plotBiomass(sim_status_quo, use_cutoff = TRUE, highlight = "Cod") +
    scale_x_continuous(name = "Intervention time tau (years)") +
    scale_y_log10(name = "Cutoff biomass (kt)", labels = kt_log_labels) +
    labs(title = "A. Status-Quo Closure (E_cod = 0)",
         subtitle = "Community remains locked in low-cod basin (~11.4 kt SSB at tau = 40 yr)") +
    theme_pub + theme(legend.position = "none")
)

p5_sp <- suppressMessages(
  plotBiomass(sim_spare_cap_base, use_cutoff = TRUE, highlight = "Cod") +
    scale_x_continuous(name = "Intervention time tau (years)") +
    scale_y_log10(name = "Cutoff biomass (kt)", labels = kt_log_labels) +
    geom_vline(xintercept = 20, linetype = "dashed", colour = "#1d5fa7", linewidth = 0.55) +
    annotate("text", x = 20.8, y = 3.0e10,
             label = "tau = 20 yr: taper +\nreopen E_cod = 0.15",
             hjust = 0, vjust = 0, size = 2.8, colour = "#1d5fa7", fontface = "bold") +
    labs(title = "B. Capelin-Sparing Pulse + Taper & Reopening",
         subtitle = "Harvesting Herring/Mackerel (0.65) + Sculpin (0.45) while keeping Capelin at 0.15",
         colour = "Species", linetype = "Species", linewidth = "Species", size = "Species") +
    theme_pub
)

extract_cod_ssb_df <- function(sim_obj, label) {
  sm <- unclass(getSSB(sim_obj))
  data.frame(
    tau     = as.numeric(rownames(sm)),
    cod_ssb = as.numeric(sm[, "Cod"]) / 1e9,
    policy  = label,
    stringsAsFactors = FALSE
  )
}
df_comp <- rbind(
  extract_cod_ssb_df(sim_status_quo,     "1. Status-quo closure (-> 15.8 kt)"),
  extract_cod_ssb_df(sim_premature,      "2. Unassisted reopen (-> 0.7 kt)"),
  extract_cod_ssb_df(sim_mesopred_only,  "3. Sculpin-only sustained (-> 1013 kt)"),
  extract_cod_ssb_df(sim_joint_sust15,   "4. Joint sustained + E_cod=0.15 (-> 665 kt)"),
  extract_cod_ssb_df(sim_joint_taper15,  "5. Joint tapered + E_cod=0.15 (-> 776 kt)"),
  extract_cod_ssb_df(sim_spare_cap_base, "6. Capelin-sparing tapered (-> 776 kt)"),
  data.frame(tau = 0:40, cod_ssb = res$hcr_det_traces$joint_ssb,
             policy = "7. State-dependent HCR (sw=23 -> 776 kt)", stringsAsFactors = FALSE)
)

p5_comp <- ggplot(df_comp, aes(x = tau, y = cod_ssb, colour = policy, linetype = policy)) +
  geom_vline(xintercept = 20, linetype = "dotted", colour = "grey55") +
  geom_line(linewidth = 1.0) +
  scale_colour_manual(values = c(
    "1. Status-quo closure (-> 15.8 kt)"          = "#64748b",
    "2. Unassisted reopen (-> 0.7 kt)"            = "#c53929",
    "3. Sculpin-only sustained (-> 1013 kt)"      = "#7c3aed",
    "4. Joint sustained + E_cod=0.15 (-> 665 kt)" = "#0d9488",
    "5. Joint tapered + E_cod=0.15 (-> 776 kt)"   = "#1d5fa7",
    "6. Capelin-sparing tapered (-> 776 kt)"      = "#15805c",
    "7. State-dependent HCR (sw=23 -> 776 kt)"    = "#b45309"
  )) +
  scale_linetype_manual(values = c(
    "1. Status-quo closure (-> 15.8 kt)"          = "solid",
    "2. Unassisted reopen (-> 0.7 kt)"            = "dashed",
    "3. Sculpin-only sustained (-> 1013 kt)"      = "dotdash",
    "4. Joint sustained + E_cod=0.15 (-> 665 kt)" = "dashed",
    "5. Joint tapered + E_cod=0.15 (-> 776 kt)"   = "solid",
    "6. Capelin-sparing tapered (-> 776 kt)"      = "solid",
    "7. State-dependent HCR (sw=23 -> 776 kt)"    = "twodash"
  )) +
  labs(title = "C. Matched Policy Cod SSB Trajectories (tau = 0-40 yr)",
       subtitle = "Legend arrows indicate verified tau = 340 yr asymptotic equilibrium fate",
       x = "Intervention time tau (years)",
       y = "Cod SSB (kt)", colour = "Policy (-> asymptotic SSB)", linetype = "Policy (-> asymptotic SSB)") +
  theme_pub +
  theme(legend.key.width = unit(1.0, "cm"))

byc_juv <- subset(res$bycatch_sweep_df,
                  policy_mode %in% c("Joint_Juvenile_Bycatch", "CapelinSpare_Juv_Bycatch"))
byc_long <- rbind(
  data.frame(effort = byc_juv$bycatch_effort, ssb = byc_juv$t40_cod_ssb_kt,
             series = ifelse(byc_juv$policy_mode == "Joint_Juvenile_Bycatch",
                             "Joint pulse (tau = 40 yr)", "Capelin-sparing (tau = 40 yr)")),
  data.frame(effort = byc_juv$bycatch_effort, ssb = byc_juv$asymp_cod_ssb_kt,
             series = ifelse(byc_juv$policy_mode == "Joint_Juvenile_Bycatch",
                             "Joint pulse (asymptotic)", "Capelin-sparing (asymptotic)"))
)

p5_byc <- ggplot(byc_long, aes(x = effort, y = ssb, colour = series, linetype = series)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.1) +
  scale_colour_manual(values = c(
    "Joint pulse (asymptotic)"       = "#1d5fa7",
    "Joint pulse (tau = 40 yr)"      = "#1d5fa7",
    "Capelin-sparing (asymptotic)"   = "#15805c",
    "Capelin-sparing (tau = 40 yr)"  = "#15805c"
  )) +
  scale_linetype_manual(values = c(
    "Joint pulse (asymptotic)"       = "solid",
    "Joint pulse (tau = 40 yr)"      = "dashed",
    "Capelin-sparing (asymptotic)"   = "solid",
    "Capelin-sparing (tau = 40 yr)"  = "dashed"
  )) +
  labs(title = "D. Sensitivity to Phase-I Juvenile Cod Bycatch",
       subtitle = "Joint pulse fails above 0.005 /yr (note misleading tau=40 overshoot at 0.01 /yr)",
       x = "Phase-I juvenile cod bycatch effort E_byc (1 / yr)",
       y = "Cod SSB (kt)", colour = NULL, linetype = NULL) +
  theme_pub +
  theme(legend.key.width = unit(1.0, "cm"))

fig5 <- (p5_sq + p5_sp) / (p5_comp + p5_byc)
save_pub_fig(fig5, "fig5_recovery_strategies", width = 12.8, height = 8.8)

# ------------------------------------------------------------------------------
# 6. Unified-Solver 3D Scenario Grid (method = "tr_bdf2", dt = 0.05) for
#    Figure 6 AND Exact Grid-Snapped Educational Web Explorer
# ------------------------------------------------------------------------------
cat("6. Computing canonical 3D scenario grid (tr_bdf2, dt = 0.05) with asymptotic verification...\n")
flush.console()

fp_vals   <- c(0.00, 0.15, 0.35, 0.50, 0.55, 0.65, 0.80)
fm_vals   <- c(0.00, 0.05, 0.20, 0.35, 0.45, 0.60)
fcod_vals <- c(0.00, 0.10, 0.15, 0.25, 0.40)

keep_tau     <- seq(0, 40, by = 2)
keep_idx     <- match(as.character(keep_tau), as.character(POLICY_YEARS))
nursery_mask <- w_vec >= 0.05 & w_vec <= 5.0

eval_one_combo <- function(idx_row) {
  i <- idx_row$i; j <- idx_row$j; k <- idx_row$k
  fp <- fp_vals[i]; fm <- fm_vals[j]; fcod <- fcod_vals[k]

  extract_node <- function(taper_flag) {
    fp2 <- if (taper_flag) 0.15 else fp
    fm2 <- if (taper_flag) 0.05 else fm
    sim40 <- run_policy_4gear(
      params_t61, fp1 = fp, fm1 = fm, fc1 = 0.00,
      fp2 = fp2, fm2 = fm2, fc2 = fcod,
      split_yr = 20L, yrs = POLICY_YEARS,
      dt = CANONICAL_DT, t_save = 1, method = CANONICAL_METHOD
    )
    ssb <- unclass(getSSB(sim40)) / 1e9
    bio <- unclass(getBiomass(sim40, use_cutoff = TRUE)) / 1e9
    yld <- unclass(getYield(sim40)) / 1e9
    nr  <- nrow(ssb)

    cum_tot <- sum(0.5 * (rowSums(yld[-nr, , drop = FALSE]) + rowSums(yld[-1, , drop = FALSE])))
    cum_cod <- sum(0.5 * (yld[-nr, "Cod"] + yld[-1, "Cod"]))

    # Verify post-policy asymptotic attractor under constant final effort
    eff_end <- c(Pelagic = fp2, Mesopredator = fm2, Flatfish = 0.15, Cod_Trawl = fcod)
    p_end40 <- finalParams(sim40)
    # Use tr_bdf2 dt=0.10 for autonomous constant-effort 300-yr asymptotic extension
    sim_as  <- project(p_end40, effort = eff_end, t_max = 300, dt = 0.10,
                       t_save = 150, method = CANONICAL_METHOD, progress_bar = FALSE)
    p_asymp <- finalParams(sim_as)
    asymp_ssb <- unname(getSSB(p_asymp)["Cod"] / 1e9)
    asymp_yld <- sum(getYield(p_asymp)) / 1e9
    lm_end    <- mean(unclass(getPredMort(p_end40))["Cod", nursery_mask])

    list(
      cod_ssb         = round(as.numeric(ssb[keep_idx, "Cod"]), 1),
      forage_bio      = round(as.numeric(rowSums(bio[keep_idx, c("Capelin", "Herring", "Mackerel")])), 1),
      capelin_bio     = round(as.numeric(bio[keep_idx, "Capelin"]), 1),
      sculpin_bio     = round(as.numeric(bio[keep_idx, "Sculpin"]), 1),
      cod_yield       = round(as.numeric(yld[keep_idx, "Cod"]), 1),
      pelagic_yield   = round(as.numeric(rowSums(yld[keep_idx, c("Capelin", "Herring", "Mackerel")])), 1),
      other_yield     = round(as.numeric(rowSums(yld[keep_idx, c("Sculpin", "Plaice")])), 1),
      t20_cod_ssb     = round(unname(ssb["20", "Cod"]), 1),
      t40_cod_ssb     = round(unname(ssb[nr, "Cod"]), 1),
      asymp_cod_ssb   = round(asymp_ssb, 1),
      asymp_recovered = asymp_ssb > 350,
      asymp_tot_yield = round(asymp_yld, 1),
      cum40_tot_yield = round(cum_tot, 1),
      cum40_cod_yield = round(cum_cod, 1),
      min_capelin_bio = round(min(bio[, "Capelin"]), 1),
      larval_mort_end = round(lm_end, 3)
    )
  }

  list(
    key = sprintf("%d_%d_%d", i - 1L, j - 1L, k - 1L),
    fp = fp, fm = fm, fcod = fcod,
    sus = extract_node(FALSE),
    tap = extract_node(TRUE)
  )
}

combo_df <- expand.grid(i = seq_along(fp_vals), j = seq_along(fm_vals), k = seq_along(fcod_vals))
scen_json_path <- "docs/data/mizer_scenarios.json"
reuse_grid <- file.exists(scen_json_path) &&
  !identical(tolower(Sys.getenv("RECOMPUTE_WEB_GRID")), "true") &&
  file.mtime(scen_json_path) >= file.mtime("models/cod_t61_params.rds")

grid_sustained <- list()
grid_taper     <- list()
surface_rows   <- list()

if (reuse_grid) {
  cat("   Reusing up-to-date canonical tr_bdf2 scenario grid from docs/data/mizer_scenarios.json...\n")
  cached_web <- fromJSON(scen_json_path, simplifyVector = TRUE)
  grid_sustained <- cached_web$grid_sustained
  grid_taper     <- cached_web$grid_taper
  for (r_idx in seq_len(nrow(combo_df))) {
    i <- combo_df$i[r_idx]; j <- combo_df$j[r_idx]; k <- combo_df$k[r_idx]
    fp <- fp_vals[i]; fm <- fm_vals[j]; fcod <- fcod_vals[k]
    if (abs(fcod - 0.15) < 1e-9) {
      ky <- sprintf("%d_%d_%d", i - 1L, j - 1L, k - 1L)
      tp <- grid_taper[[ky]]
      n_last <- length(tp$cod_ssb)
      surface_rows[[length(surface_rows) + 1L]] <- data.frame(
        fp              = fp,
        fm              = fm,
        t40_cod_ssb_kt  = tp$t40_cod_ssb,
        asymp_cod_kt    = tp$asymp_cod_ssb,
        asymp_recovered = tp$asymp_recovered,
        t40_tot_yld_kt  = tp$cod_yield[n_last] + tp$pelagic_yield[n_last] + tp$other_yield[n_last],
        cum40_tot_mt    = tp$cum40_tot_yield / 1000,
        stringsAsFactors = FALSE
      )
    }
  }
} else {
  combo_list <- split(combo_df, seq_len(nrow(combo_df)))
  grid_results <- mclapply(combo_list, eval_one_combo, mc.cores = NC)
  for (el in grid_results) {
    grid_sustained[[el$key]] <- el$sus
    grid_taper[[el$key]]     <- el$tap
    if (abs(el$fcod - 0.15) < 1e-9) {
      n_last <- length(el$tap$cod_ssb)
      surface_rows[[length(surface_rows) + 1L]] <- data.frame(
        fp              = el$fp,
        fm              = el$fm,
        t40_cod_ssb_kt  = el$tap$t40_cod_ssb,
        asymp_cod_kt    = el$tap$asymp_cod_ssb,
        asymp_recovered = el$tap$asymp_recovered,
        t40_tot_yld_kt  = el$tap$cod_yield[n_last] + el$tap$pelagic_yield[n_last] + el$tap$other_yield[n_last],
        cum40_tot_mt    = el$tap$cum40_tot_yield / 1000,
        stringsAsFactors = FALSE
      )
    }
  }
}

# STRICT CROSS-ARTIFACT CONSISTENCY ASSERTION (Referee Comment 4):
# Verify that the exported grid node for joint_taper_f15 (fp=0.65 -> idx 5, fm=0.45 -> idx 4, fcod=0.15 -> idx 2)
# matches R/02_analyse_bistability_and_scenarios.R identically!
flagship_grid <- grid_taper[["5_4_2"]]
flagship_r02  <- subset(res$policy_comparison_df, id == "joint_taper_f15")
stopifnot(abs(flagship_grid$t40_cod_ssb - round(flagship_r02$t40_cod_ssb_kt, 1)) < 1e-6)
stopifnot(abs(flagship_grid$cum40_tot_yield - round(flagship_r02$cum40_tot_yld_kt, 1)) < 1e-6)

# ------------------------------------------------------------------------------
# Figure 6: 2D Policy Response Surface (tau=40 SSB + Asymptotic Fate, Yield, & Cumulative Catch)
# ------------------------------------------------------------------------------
cat("7. Generating Figure 6 (Unified-Solver Policy Response & Trade-off Surface)...\n")
flush.console()

surf_df <- do.call(rbind, surface_rows)
surf_df$fp_fac <- factor(sprintf("%.2f", surf_df$fp), levels = sprintf("%.2f", fp_vals))
surf_df$fm_fac <- factor(sprintf("%.2f", surf_df$fm), levels = sprintf("%.2f", fm_vals))
surf_df$cell_lbl <- ifelse(surf_df$asymp_recovered,
                           paste0(round(surf_df$t40_cod_ssb_kt), "*"),
                           as.character(round(surf_df$t40_cod_ssb_kt)))

p6_ssb <- ggplot(surf_df, aes(x = fp_fac, y = fm_fac, fill = t40_cod_ssb_kt)) +
  geom_tile(aes(colour = asymp_recovered), linewidth = 0.75) +
  geom_text(aes(label = cell_lbl), size = 3.0, fontface = "bold",
            colour = ifelse(surf_df$t40_cod_ssb_kt > 320, "white", "#10243e")) +
  scale_fill_gradient2(low = "#fcebe9", mid = "#93c5fd", high = "#1d5fa7", midpoint = 300,
                       name = "tau=40 SSB\n(kt)") +
  scale_colour_manual(values = c(`TRUE` = "#0f172a", `FALSE` = "white"), guide = "none") +
  labs(title = "A. Cod SSB at tau = 40 yr (* = Rebuilds to 776 kt)",
       subtitle = "Tapered at tau = 20 (E_cod = 0.15); border/* = verified upper basin",
       x = "Phase-I pelagic effort E_pel (1 / yr)",
       y = "Phase-I sculpin effort E_meso (1 / yr)") +
  theme_pub

p6_yld <- ggplot(surf_df, aes(x = fp_fac, y = fm_fac, fill = t40_tot_yld_kt)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(t40_tot_yld_kt)), size = 3.0, fontface = "bold",
            colour = ifelse(surf_df$t40_tot_yld_kt > 520, "white", "#10243e")) +
  scale_fill_gradient(low = "#fef3c7", high = "#b45309", name = "tau=40 Yield\n(kt / yr)") +
  labs(title = "B. Projected Annual Total Yield at tau = 40 yr",
       subtitle = "Lower endpoint yield reflects forage cropping by rebuilt cod",
       x = "Phase-I pelagic effort E_pel (1 / yr)",
       y = "Phase-I sculpin effort E_meso (1 / yr)") +
  theme_pub

p6_cum <- ggplot(surf_df, aes(x = fp_fac, y = fm_fac, fill = cum40_tot_mt)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", cum40_tot_mt)), size = 3.0, fontface = "bold",
            colour = ifelse(surf_df$cum40_tot_mt > 32, "white", "#10243e")) +
  scale_fill_gradient(low = "#ecfdf5", high = "#15805c", name = "Cum. Yield\n(Mt, 40 yr)") +
  labs(title = "C. Cumulative 40-Yr Multispecies Catch (Mt)",
       subtitle = "Integrated harvest across all fleets over tau in [0, 40] yr",
       x = "Phase-I pelagic effort E_pel (1 / yr)",
       y = "Phase-I sculpin effort E_meso (1 / yr)") +
  theme_pub

fig6 <- p6_ssb + p6_yld + p6_cum + plot_layout(ncol = 3)
save_pub_fig(fig6, "fig6_policy_tradeoffs", width = 14.5, height = 5.0)

# ------------------------------------------------------------------------------
# 7. Export Educational Explorer Data (docs/model_data.js & mizer_scenarios.json)
#    Using relative simulation years (-61..0 collapse+closure, 0..40 intervention)
# ------------------------------------------------------------------------------
cat("8. Exporting docs/model_data.js and docs/data/mizer_scenarios.json...\n")
flush.console()

hist_t_raw  <- COLLAPSE_YEARS # 0..61
hist_keep   <- c(seq(0, 60, by = 2), 61)
h_idx       <- match(as.character(hist_keep), rownames(unclass(getSSB(sim_collapse))))
rel_hist_yr <- hist_keep - 61L # -61..0 (tau = 0 is intervention start)

h_ssb <- unclass(getSSB(sim_collapse))[h_idx, , drop = FALSE] / 1e9
h_bio <- unclass(getBiomass(sim_collapse, use_cutoff = TRUE))[h_idx, , drop = FALSE] / 1e9
h_yld <- unclass(getYield(sim_collapse))[h_idx, , drop = FALSE] / 1e9

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

baseline_f15_ssb <- round(unname(getSSB(params_f15)["Cod"] / 1e9), 1)
high_f0_ssb      <- round(unname(getSSB(params_high_f0)["Cod"] / 1e9), 1)
low_f0_ssb       <- round(unname(getSSB(params_low_f0)["Cod"] / 1e9), 1)
forage_f15_bio   <- sum(getBiomass(params_f15, use_cutoff = TRUE)[1:4]) / 1e9
forage_low_bio   <- sum(getBiomass(params_low_f0, use_cutoff = TRUE)[1:4]) / 1e9
surge_ratio_val  <- round(forage_low_bio / forage_f15_bio, 1)
surv_drop_ratio  <- round(res$haz_high_15$surv_cont / res$haz_low_15$surv_cont, 1)

pol_df <- res$policy_comparison_df
web_policy_table <- lapply(seq_len(nrow(pol_df)), function(r) {
  rw <- pol_df[r, ]
  list(
    name             = rw$strategy,
    ecap_ph1         = sprintf("%.2f", rw$ecap_ph1),
    ehm_ph1          = sprintf("%.2f", rw$ehm_ph1),
    emeso_ph1        = sprintf("%.2f", rw$emeso_ph1),
    ecod_ph2         = sprintf("%.2f", rw$ecod_ph2),
    tapered          = ifelse(rw$tapered, "Yes (tau=20)", "No"),
    t20_ssb_kt       = round(rw$t20_cod_ssb_kt, 1),
    t40_ssb_kt       = round(rw$t40_cod_ssb_kt, 1),
    asymp_ssb_kt     = round(rw$asymp_cod_ssb_kt, 1),
    cum40_tot_yld_kt = round(rw$cum40_tot_yld_kt, 0),
    t40_tot_yld_kt   = round(rw$t40_tot_yld_kt, 1),
    asymp_tot_yld_kt = round(rw$asymp_tot_yld_kt, 1),
    min_cap_ph1_kt   = round(rw$min_cap_ph1_kt, 0),
    recovered        = isTRUE(rw$asymp_recovered),
    outcome          = ifelse(rw$asymp_recovered, "Upper Cultivated Basin", "Lower Depensatory Basin")
  )
})

web_data <- list(
  solver = list(
    method = CANONICAL_METHOD,
    dt     = CANONICAL_DT,
    note   = "Exact precomputed mizer trajectories (no off-grid cross-basin interpolation)"
  ),
  hero = list(
    high_low_ssb     = paste0(round(high_f0_ssb), " vs ", low_f0_ssb, " kt"),
    forage_surge     = paste0(surge_ratio_val, "&times;"),
    survival_drop    = paste0(surv_drop_ratio, "&times;"),
    bycatch_limit    = "&lt; 0.01 /yr"
  ),
  policy_start_year = 0L,
  closure_year      = -34L,
  reference = list(
    baseline_f15_cod_ssb = baseline_f15_ssb,
    high_f0_cod_ssb      = high_f0_ssb,
    low_f0_cod_ssb       = low_f0_ssb,
    saddle_f0_cod_ssb    = round(res$saddle_summary_df$saddle_ssb_kt[1], 1),
    saddle_f15_cod_ssb   = round(res$saddle_summary_df$saddle_ssb_kt[2], 1)
  ),
  presets = list(
    status_quo        = list(fp = 0.15, fm = 0.05, fcod = 0.00, taper = FALSE),
    pelagic_rebalance = list(fp = 0.55, fm = 0.05, fcod = 0.00, taper = FALSE),
    ecosystem_pulse   = list(fp = 0.65, fm = 0.45, fcod = 0.15, taper = TRUE),
    premature_reopen  = list(fp = 0.15, fm = 0.05, fcod = 0.25, taper = FALSE)
  ),
  history = list(
    years         = as.integer(rel_hist_yr),
    cod_ssb       = round(as.numeric(h_ssb[, "Cod"]), 1),
    forage_bio    = round(as.numeric(rowSums(h_bio[, c("Capelin", "Herring", "Mackerel")])), 1),
    sculpin_bio   = round(as.numeric(h_bio[, "Sculpin"]), 1),
    cod_yield     = round(as.numeric(h_yld[, "Cod"]), 1),
    pelagic_yield = round(as.numeric(rowSums(h_yld[, c("Capelin", "Herring", "Mackerel")])), 1),
    other_yield   = round(as.numeric(rowSums(h_yld[, c("Sculpin", "Plaice")])), 1)
  ),
  proj_years     = as.integer(keep_tau),
  grid_axes      = list(fp = fp_vals, fm = fm_vals, fcod = fcod_vals),
  grid_sustained = grid_sustained,
  grid_taper     = grid_taper,
  size_diagnostics = list(
    note       = "Static reference comparison of the two E_cod = 0 equilibrium MizerParams objects (not recalculated per policy slider)",
    log10_w    = round(log10(w_vec), 4),
    moratorium = extract_spectra_diag(params_low_f0),
    cultivated = extract_spectra_diag(params_high_f0)
  ),
  bifurcation = list(
    vs_fcod = list(
      f_up            = round(cod_up_merged$effort_val, 4),
      ssb_up          = round(cod_up_merged$cod_ssb_kt, 1),
      f_down          = round(rev(res$scan_cod_down_df$effort_val), 4),
      ssb_down        = round(rev(res$scan_cod_down_df$cod_ssb_kt), 1),
      bracket_lo      = round(res$fine_cod_bracket[1], 3),
      bracket_hi      = round(res$fine_cod_bracket[2], 3)
    ),
    vs_fpelagic = list(
      f_up            = round(pel_up_merged$effort_val, 4),
      ssb_up          = round(pel_up_merged$cod_ssb_kt, 1),
      f_down          = round(rev(res$scan_pel_down_df$effort_val), 4),
      ssb_down        = round(rev(res$scan_pel_down_df$cod_ssb_kt), 1),
      bracket_lo      = round(res$fine_pel_bracket[1], 3),
      bracket_hi      = round(res$fine_pel_bracket[2], 3)
    )
  ),
  policy_table = web_policy_table
)

json_str <- toJSON(web_data, auto_unbox = TRUE, digits = 4, null = "null")
writeLines(json_str, "docs/data/mizer_scenarios.json")
writeLines(paste0("window.MIZER_MODEL_DATA = ", json_str, ";\n"), "docs/model_data.js")
cat("Successfully generated figures 1-6 and exported canonical web explorer data.\n")
if (sink.number() > 0) sink()
