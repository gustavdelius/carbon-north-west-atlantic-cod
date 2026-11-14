# ==============================================================================
# R/00_model_functions.R
#
# Custom mizer recruitment rate function incorporating sheltered inshore coastal
# spawning refuges (e.g. Smith Sound, Trinity Bay, and Gilbert Bay in
# Newfoundland; Rose 2007, Lilly et al. 2013) alongside offshore Beverton-Holt
# density-dependent recruitment.
# ==============================================================================

#' Beverton-Holt recruitment with inshore coastal spawning refuge
#'
#' Adds a small density-independent inshore refuge recruitment flux
#' `species_params$r_refuge` to standard `mizer::BevertonHoltRDD()`.
#' Registered on the model via `setRateFunction(params, "RDD", "RefugeBevertonHoltRDD")`.
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
