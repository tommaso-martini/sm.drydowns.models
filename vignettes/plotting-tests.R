# 1) Profilo “pulito” 0–30 cm, 3 layer × 0.10 m (overlay t0 vs t1) -------------

# 3 layers spanning 0–0.30 m
depths <- c(0.05, 0.15, 0.25)   # layer centers [m]
dz     <- rep(0.10, 3)          # layer thickness [m]

# theta at two times (t0, t1) — t1 is drier
theta_t0 <- c(0.26, 0.24, 0.22)
theta_t1 <- c(0.24, 0.22, 0.20)

# Plot with "segments" (classic stairs)
plot_profile_stairs(theta_t1, depths, dz = dz, col = "black", lwd = 2,
                    xlim = c(0.18, 0.30), style = "segments",
                    main = "Stair profile (0–30 m) — t1 vs t0")
plot_profile_stairs(theta_t0, depths, dz = dz, col = "blue",  lwd = 2,
                    add = TRUE, style = "segments")
legend("bottomright", c("t1", "t0"), col = c("black","blue"), lwd = 2, bty = "n")

# 2) Layer irregolari, senza passare dz (la funzione li ricostruisce) ----------

# Irregular layer centers (e.g., 0–0.6 m domain)
depths <- c(0.06, 0.18, 0.32, 0.50)

# No dz: function infers top/bottom from midpoints between centers
theta_t0 <- c(0.30, 0.28, 0.26, 0.24)
theta_t1 <- c(0.27, 0.26, 0.24, 0.22)

plot_profile_stairs(theta_t1, depths, col = "black", lwd = 2,
                    xlim = c(0.20, 0.32), style = "segments",
                    main = "Irregular layers (dz inferred) — t1 vs t0")
plot_profile_stairs(theta_t0, depths, col = "steelblue", lwd = 2,
                    add = TRUE, style = "segments")
legend("bottomright", c("t1", "t0"), col = c("black","steelblue"), lwd = 2, bty = "n")

# 3) Rettangoli pieni con alpha (utile per evidenziare differenze) -------------
depths <- c(0.05, 0.15, 0.25)
dz     <- rep(0.10, 3)

theta_t0 <- c(0.25, 0.23, 0.21)
theta_t1 <- c(0.22, 0.205, 0.195)

# First draw t0 as light blue rectangles
plot_profile_stairs(theta_t0, depths, dz = dz,
                    style = "rect", col = "blue", alpha = 0.20, border = NA,
                    xlim = c(0.18, 0.27), main = "Rectangles with alpha — t0 (blue) vs t1 (red)")
# Then overlay t1 as light red rectangles + a dark spine
plot_profile_stairs(theta_t1, depths, dz = dz,
                    style = "rect", col = "red", alpha = 0.20, border = NA, add = TRUE)
plot_profile_stairs(theta_t1, depths, dz = dz,
                    style = "segments", col = "red", add = TRUE, lwd = 2)
legend("bottomright", c("t1 (red)", "t0 (blue)"),
       col = c("red","blue"), lwd = 2, bty = "n")

# 4) Overlay di più profili (timeline breve) -----------------------------------

depths <- c(0.05, 0.15, 0.25)
dz     <- rep(0.10, 3)

# Build a small synthetic timeline: 5 time steps drying down
theta_series <- rbind(
  c(0.26, 0.24, 0.22),
  c(0.255, 0.235, 0.215),
  c(0.248, 0.229, 0.210),
  c(0.242, 0.223, 0.205),
  c(0.236, 0.218, 0.200)
)  # [time x layer]; we'll plot by columns, so just transpose or index properly.

# --- Synthetic example for plot_daily_profiles() ------------------------------

# Load required plotting function (assume it's already in your environment)
# source("R/plots.R")   # if you saved it in your package folder

# --- geometry (0–30 m, 3 layers) ---
depths <- c(0.05, 0.15, 0.25)
dz     <- rep(0.10, 3)

# --- build fake res object similar to musec_simulate() output ---

n_days  <- 10
steps   <- n_days * 24   # hourly steps
n_layers <- length(depths)

# Create a synthetic drying pattern (theta decreases with time)
theta_mm <- sapply(1:steps, function(t) {
  base <- c(0.30, 0.28, 0.26)
  decay <- 0.001 * t       # linear drying rate
  pmax(base - decay, 0.05) # avoid negative values
})
rownames(theta_mm) <- paste0("layer", 1:n_layers)

# Create synthetic fluxes (mm/step)
evap_mm   <- rep(0.10, steps) * (1 + sin(seq(0, 4*pi, length.out=steps)) * 0.2)
transp_mm <- rep(0.05, steps) * (1 + cos(seq(0, 4*pi, length.out=steps)) * 0.1)
drain_mm  <- rep(0.02, steps)

# Build list like musec_simulate() returns
res <- list(
  series = list(
    theta_mm         = theta_mm,
    evap_mm          = evap_mm,
    transp_mm        = transp_mm,
    drain_bottom_mm  = drain_mm
  ),
  meta = list(structure = list(depths = depths, dz = dz))
)

# --- synthetic rain and ET0 vectors (hourly) ---
set.seed(42)
rain_vec <- rbinom(steps, size = 1, prob = 0.15) * runif(steps, 0, 1) * 2   # some random rainfall
ET0_vec  <- rep(3/24, steps) * (1 + 0.3*sin(seq(0, 2*pi, length.out = steps))) # ~3 mm/day sinusoidal

# --- run the wrapper (plots to device) ---
plot_daily_profiles(
  res, depths, dz,
  rain = rain_vec, ET0 = ET0_vec,
  file = NULL,                # set e.g. "profiles_demo.pdf" to save
  days = 2:10,                # days to plot
  xlim = c(0, 0.35),
  style = "segments",         # or "rect"
  add_areas = TRUE, alpha = 0.30,
  main_prefix = "Synthetic demo 0–30 m"
)