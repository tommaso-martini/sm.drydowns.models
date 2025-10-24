# ============================================================================
# DES MODEL (mono-layer daily -> optional sub-daily within-day integration)
# ============================================================================
# - Rainy days: simplified water balance (theta updated by rain + overflow; reset DES memory)
# - Dry days : sub-steps with DES-like evaporation, Feddes transpiration, and drainage
# - Drainage: by default uses a power-law K(theta); if analytical=TRUE, calls update_theta_analytical()
# - ET partition: E = ET0 * exp(-k_ext * LAI),  T = ET0 - E
# - Units:
#     * Ks in m/day
#     * rain_vec, ET_vec in mm/day
#     * dt in days (e.g., 1/24)
#     * dz, root_depth in m
#     * theta in m3/m3
#
# !!!NOTE: a precipitation vector is among the parameters: this is a leftover from
# a previous implementation involving a simplified infiltration mechanism (the
# simple water balance equation). The default is a vector of zeroes, as we are
# dealing with drydowns simulation.
# ============================================================================
#' DES drydown simulator (mono-layer) with daily rain/ET, LAI-based partition, and power-law drainage
#'
#' @param Ks         saturated hydraulic conductivity [m/day]
#' @param n          van Genuchten n parameter [-] (also used to shape K(\theta) power via p = 1/n)
#' @param theta_init initial water content [m3/m3]
#' @param ET_vec     daily potential evapotranspiration [mm/day] (required)
#' @param rain_vec   daily rainfall [mm/day] (default all zeros if omitted)
#' @param LAI_vec    daily leaf area index [-]; if NULL, uses \code{LAI_default}
#' @param k_ext      light extinction coefficient for ET partition [-], default 0.378
#' @param LAI_default scalar LAI used if \code{LAI_vec} is NULL (default 1)
#' @param theta_r    residual water content [m3/m3]
#' @param theta_s    saturated water content [m3/m3]
#' @param theta_fc   field capacity [m3/m3]
#' @param theta_star optimal water content (no stress) [m3/m3]
#' @param theta_wp   wilting point [m3/m3]
#' @param psi_b      bubbling pressure [m]
#' @param root_depth rooting depth [m] (default 0.60)
#' @param dt         sub-step length for dry days [days], e.g. 1/24 (default)
#' @param dz         layer thickness [m] (default 0.05)
#' @param desorp     use DES-style cumulative evaporation on dry days (default TRUE)
#' @param analytical if TRUE use \code{update_theta_analytical()} for drainage on dry days;
#'                   otherwise use the power-law K(\theta) (default FALSE)
#' @param debug      verbose console messages (default FALSE)
#' @param hourly     return sub-daily expanded output (true within-day evolution) (default TRUE)
#' @param rain_split placeholder for future rain distribution options (default "first")
#'
#' @return A list with:
#' \itemize{
#'   \item \code{series}: per-step (expanded) and daily time series:
#'     \describe{
#'       \item{\code{theta}}{sub-daily \theta (length = n_day * n_substeps) if \code{hourly=TRUE}, else daily}
#'       \item{\code{evap_mm}, \code{transp_mm}, \code{drain_mm}, \code{overflow_mm}}{sub-daily (or daily) fluxes in mm}
#'       \item{\code{theta_daily}, \code{evap_daily_mm}, \code{transp_daily_mm}, \code{drain_daily_mm}, \code{overflow_daily_mm}}{daily series}
#'     }
#'   \item \code{meta}: units, flags, structure parameters.
#' }
#' @export
des_simulate <- function(Ks,
                         n,
                         theta_init,
                         ET_vec,
                         rain_vec = NULL,
                         LAI_vec  = NULL,
                         k_ext = 0.378,
                         LAI_default = 1,
                         theta_r,
                         theta_s,
                         theta_fc,
                         theta_star,
                         theta_wp,
                         psi_b,
                         root_depth = 0.60,
                         dt = 1/24,
                         dz = 0.05,
                         desorp = TRUE,
                         analytical = FALSE,
                         debug  = FALSE,
                         hourly = TRUE,
                         rain_split = c("first", "uniform")) {
  
  rain_split <- match.arg(rain_split)
  
  # --- Basic checks ---
  stopifnot(is.finite(Ks), is.finite(n), is.finite(theta_init),
            is.finite(theta_r), is.finite(theta_s),
            is.finite(theta_fc), is.finite(theta_star), is.finite(theta_wp),
            is.finite(psi_b), is.finite(root_depth), is.finite(dt), is.finite(dz))
  
  if (!is.numeric(ET_vec) || length(ET_vec) < 1)
    stop("ET_vec must be a numeric vector of daily ET0 [mm/day].")
  
  n_day <- length(ET_vec)
  
  if (is.null(rain_vec)) {
    rain_vec <- rep(0, n_day)
  } else if (length(rain_vec) != n_day) {
    stop("rain_vec and ET_vec must have the same length (daily).")
  }
  
  if (is.null(LAI_vec)) {
    LAI_vec <- rep(LAI_default, n_day)
  } else if (length(LAI_vec) != n_day) {
    stop("LAI_vec must have the same length as ET_vec (daily).")
  }
  
  if (dt <= 0 || dt > 1)
    stop("dt must be in days and in (0,1]. For hourly, use dt = 1/24.")
  
  if (!exists("desorptivity", mode = "function"))
    stop("desorptivity() not found. Please provide it in your utils (see documentation).")
  if (analytical && !exists("update_theta_analytical", mode = "function"))
    stop("analytical=TRUE but update_theta_analytical() is not available in utils.")
  if (!exists("K_theta", mode = "function"))
    stop("K_theta() not found. Please provide it in hydraulics.R and export it.")
  
  # --- Pre-allocate daily outputs ---
  theta_daily_mean   <- numeric(n_day)  # daily mean of theta
  theta_end_daily    <- numeric(n_day)  # end-of-day theta (last sub-step)
  
  E_daily_mm         <- numeric(n_day)
  D_daily_mm         <- numeric(n_day)
  T_daily_mm         <- numeric(n_day)
  overflow_daily_mm  <- numeric(n_day)
  
  # --- Initial state & DES memory ---
  theta_now <- theta_init
  E_cum <- 1e-6  # m; keep > 0 to avoid division by zero
  
  # --- Substep size in hours (needed because K_theta returns m/h) ---
  sub_steps <- as.integer(round(1 / dt))
  if (sub_steps < 1) sub_steps <- 1
  dt_h <- dt * 24  # hours per sub-step
  
  # --- Traces for "hourly" output (true within-day evolution) ---
  theta_trace   <- if (hourly) numeric(0) else NULL
  evap_trace_mm <- if (hourly) numeric(0) else NULL
  transp_trace_mm <- if (hourly) numeric(0) else NULL
  drain_trace_mm  <- if (hourly) numeric(0) else NULL
  overflow_trace_mm <- if (hourly) numeric(0) else NULL
  
  # ======================= MAIN DAILY LOOP =======================
  for (day in seq_len(n_day)) {
    
    if (debug) cat(sprintf("\n===== Day %d =====\n", day))
    
    rain_mm <- max(0, rain_vec[day])   # mm/day
    ET0_mm  <- max(0, ET_vec[day])     # mm/day
    LAI_d   <- LAI_vec[day]
    
    if (rain_mm > 0) {
      # ---------- RAINY DAY (simplified) ----------
      if (debug) cat(" Rainy day - simplified water balance\n")
      
      # 1) add rain -> Δtheta = rain/1000 / dz
      theta_now <- theta_now + (rain_mm / 1000) / dz
      
      # 2) overflow if theta > theta_s
      overflow_m  <- pmax(0, (theta_now - theta_s) * dz)
      overflow_mm <- overflow_m * 1000
      theta_now   <- pmin(theta_now, theta_s)
      
      # 3) reset DES memory
      E_cum <- 1e-6
      
      # 4) daily accounting
      E_daily_mm[day]        <- 0
      T_daily_mm[day]        <- 0
      D_daily_mm[day]        <- 0
      overflow_daily_mm[day] <- overflow_mm
      
      # daily bookkeeping for rainy day (flat within the day)
      theta_end_daily[day]   <- theta_now         # end-of-day theta
      theta_daily_mean[day]  <- theta_now         # mean equals the flat value
      
      # 5) if hourly, replicate flat sub-steps for this rainy day
      if (hourly) {
        theta_trace         <- c(theta_trace, rep(theta_now, sub_steps))
        evap_trace_mm       <- c(evap_trace_mm, rep(0, sub_steps))
        transp_trace_mm     <- c(transp_trace_mm, rep(0, sub_steps))
        drain_trace_mm      <- c(drain_trace_mm, rep(0, sub_steps))
        overflow_trace_mm   <- c(overflow_trace_mm, rep(overflow_mm / sub_steps, sub_steps))
      }
      
      if (debug) {
        cat(sprintf("🌧️  Rain %.2f mm | theta_post %.4f | Overflow %.2f mm\n",
                    rain_mm, theta_now, overflow_mm))
      }
      
    } else {
      # ---------- DRY DAY (sub-step integration) ----------
      if (debug) cat("Dry day - sub-step integration\n")
      
      theta_day_sum <- 0  # accumulator for daily mean theta
      
      # Partition ET0 at the day level; then distribute across sub-steps
      E_day_req_m <- (ET0_mm * exp(-k_ext * LAI_d)) / 1000  # m/day
      T_day_req_m <- (ET0_mm / 1000) - E_day_req_m          # m/day
      E_sub_req_m <- E_day_req_m / sub_steps
      T_sub_req_m <- T_day_req_m / sub_steps
      
      # daily accumulators (mm/day)
      E_day_acc_mm <- 0
      D_day_acc_mm <- 0
      T_day_acc_mm <- 0
      
      # convert Ks to m/h for K_theta()
      Ks_h <- Ks / 24
      
      for (s in seq_len(sub_steps)) {
        theta_pre <- theta_now
        
        # (1) Drainage (above field capacity) -> use K_theta (m/h)
        dtheta_drain <- 0
        if (theta_pre > theta_fc) {
          if (analytical) {
            # analytical solver expects Ks in m/day and dt in days -> consistent as-is
            theta_after <- update_theta_analytical(theta_pre, dt, Ks, dz,
                                                   theta_r, theta_s, n)
            dtheta_drain <- max(0, theta_pre - theta_after)
          } else {
            # K_theta returns m/h; use hours in the step (dt_h)
            Ki_h <- K_theta(theta_pre, Ks = Ks_h, n = n,
                            theta_s = theta_s, theta_r = theta_r, theta_fc = theta_fc)
            dtheta_drain <- min(max(Ki_h, 0) * (dt_h / dz), theta_pre - theta_fc)
            dtheta_drain <- max(0, dtheta_drain)
          }
        }
        D_day_acc_mm <- D_day_acc_mm + dtheta_drain * dz * 1000
        
        # (2) Evaporation (DES cumulative limiter)
        dtheta_E <- 0
        if (desorp && E_sub_req_m > 0) {
          Theta <- (theta_pre - theta_r) / (theta_s - theta_r)
          Theta <- min(max(Theta, 0), 1)
          
          S_val <- desorptivity(Ks = Ks_h, n = n,
                                theta_s = theta_s, theta_r = theta_r,
                                psi_b = psi_b, Theta = Theta)   # m/sqrt(day)
          
          e_star <- (S_val^2) / (2 * E_cum)    # m
          
          if (e_star <= E_sub_req_m) {
            E_cum_old <- E_cum
            E_cum <- E_cum * sqrt(1 + (S_val / E_cum)^2)
            E_take_m <- max(0, E_cum - E_cum_old)
          } else {
            E_take_m <- E_sub_req_m
            E_cum    <- E_cum + E_take_m
          }
          
          avail_after_drain <- max(theta_pre - dtheta_drain - theta_r, 0)
          dtheta_E <- min(E_take_m / dz, avail_after_drain)
        }
        E_day_acc_mm <- E_day_acc_mm + dtheta_E * dz * 1000
        
        # (3) Transpiration (Feddes + top-layer fraction)
        dtheta_T <- 0
        if (T_sub_req_m > 0) {
          alpha    <- feddes_stress(theta_pre, theta_star, theta_wp)
          frac_top <- min(1, max(0, dz / root_depth)) # fraction of top layer wrt root zone
          T_i_m    <- T_sub_req_m * alpha * frac_top  # m in this layer
          dtheta_T <- min(T_i_m / dz, max(theta_pre - theta_r, 0))
        }
        T_day_acc_mm <- T_day_acc_mm + dtheta_T * dz * 1000
        
        # (4) Update theta
        theta_now <- max(theta_r, theta_pre - dtheta_drain - dtheta_E - dtheta_T)
        
        # accumulate for daily mean
        theta_day_sum <- theta_day_sum + theta_now
        
        # (5) Store traces if hourly
        if (hourly) {
          theta_trace       <- c(theta_trace, theta_now)
          evap_trace_mm     <- c(evap_trace_mm, dtheta_E * dz * 1000)
          transp_trace_mm   <- c(transp_trace_mm, dtheta_T * dz * 1000)
          drain_trace_mm    <- c(drain_trace_mm, dtheta_drain * dz * 1000)
          overflow_trace_mm <- c(overflow_trace_mm, 0)
        }
      }
      
      # daily bookkeeping (END and MEAN)
      theta_end_daily[day]   <- theta_now                  # end-of-day theta
      theta_daily_mean[day]  <- theta_day_sum / sub_steps  # daily mean theta
      
      E_daily_mm[day]        <- E_day_acc_mm
      D_daily_mm[day]        <- D_day_acc_mm
      T_daily_mm[day]        <- T_day_acc_mm
      overflow_daily_mm[day] <- 0
      
      if (debug) {
        cat(sprintf("Evap %.2f mm | Transp %.2f mm | Drain %.2f mm | theta_post %.4f\n",
                    E_day_acc_mm, T_day_acc_mm, D_day_acc_mm, theta_now))
      }
    }
  } # end for day
  
  # --- Build output ---
  if (hourly) {
    series <- list(
      theta              = theta_trace,
      evap_mm            = evap_trace_mm,
      transp_mm          = transp_trace_mm,
      drain_mm           = drain_trace_mm,
      overflow_mm        = overflow_trace_mm,
      
      # Daily (now both)
      theta_daily        = theta_daily_mean,   # NOTE: daily *mean* (retro-compat)
      theta_end_daily    = theta_end_daily,    # NEW: end-of-day theta (last sub-step)
      
      evap_daily_mm      = E_daily_mm,
      transp_daily_mm    = T_daily_mm,
      drain_daily_mm     = D_daily_mm,
      overflow_daily_mm  = overflow_daily_mm,
      
      day_index          = rep(seq_len(n_day), each = sub_steps)
    )
  } else {
    series <- list(
      theta              = theta_daily_mean,   # daily mean if not expanding
      theta_end_daily    = theta_end_daily,    # NEW
      
      evap_mm            = E_daily_mm,
      transp_mm          = T_daily_mm,
      drain_mm           = D_daily_mm,
      overflow_mm        = overflow_daily_mm
    )
  }
  
  meta <- list(
    units = list(Ks = "m/day", rain = "mm/day", ET = "mm/day",
                 fluxes = "mm", theta = "m3/m3"),
    params = list(k_ext = k_ext, LAI_default = LAI_default,
                  dt_days = dt, dz_m = dz, root_depth_m = root_depth,
                  analytical = analytical, desorp = desorp),
    notes = "Daily forcing with within-day sub-steps on dry days; rainy days simplified."
  )
  
  list(series = series, meta = meta)
}
