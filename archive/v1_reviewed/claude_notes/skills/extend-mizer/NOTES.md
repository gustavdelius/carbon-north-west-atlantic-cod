# Project Notes — `extend-mizer`

- **Custom refuge recruitment rate (`R/00_model_functions.R`):** `RefugeBevertonHoltRDD(rdi, species_params, params, t, ...)` wraps `BevertonHoltRDD(rdi, species_params, ...)` and adds `species_params$r_refuge` (representing sheltered coastal/fjord spawning refuges such as Smith Sound and Gilbert Bay decoupled from offshore pelagic predation).
- **Reload requirement:** Always `source("R/00_model_functions.R")` before calling `readParams()` or `project()` on saved `models/*.rds` objects so `RefugeBevertonHoltRDD` is present in the global environment (`getRateFunction(params, "RDD")`).
