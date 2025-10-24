#' muSEC simulator (mono- or multi-layer)
#'
#' Simulates muSEC with SEC evaporation at the surface layer, drainage (exponential or power-law),
#' and transpiration with Feddes stress, for mono- and multi-layer configurations.
#'
#' @param theta_init numeric scalar (mono-layer) or numeric vector (multi-layer, top→bottom) in m3/m3
#' @param structure list:
#'   - mono: list(dz = <thickness [m]>)
#'   - multi: list(depths = <centres [m]>, dz = <thicknesses [m]>)
#' @param Ks,n,theta_s,theta_r,theta_fc,theta_wp,theta_star hydraulic and stress parameters
#' @param ET numeric [mm/step], optional (alternative to E+T)
#' @param E,T numeric [mm/step], optional (if already partitioned)
#' @param LAI numeric scalar/vector; if NULL uses LAI_default
#' @param k_ext light-extinction coefficient (default 0.378)
#' @param LAI_default default LAI when missing (default 1)
#' @param timestep "hourly" or "daily" (if daily → hourly expansion with sinusoid)
#' @param hourly_out if TRUE and mono-layer, also returns daily means
#' @param exp_drain if TRUE uses exponential K calibrated at s_fc; else power-law K_theta()
#' @param roots_depth max rooting depth [m] for multi-layer; default = profile bottom
#' @param debug verbose diagnostics
#' @return list(mode = "mono"/"multi", series = ..., meta = ...)
#' @export
musec_simulate <- function(
    theta_init,
    structure,
    Ks, n, theta_s, theta_r, theta_fc, theta_wp, theta_star,
    ET = NULL, E = NULL, T = NULL, LAI = NULL,
    k_ext = 0.378, LAI_default = 1,
    timestep = c("hourly","daily"),
    hourly_out = FALSE,
    exp_drain = TRUE,
    roots_depth = NULL,
    debug = FALSE
){
  timestep <- match.arg(timestep)
  
  # ---- determine mono vs multi ----
  is_multi <- length(theta_init) > 1L
  if (is_multi) {
    stopifnot(is.list(structure), !is.null(structure$depths), !is.null(structure$dz))
    depths <- as.numeric(structure$depths)
    dz     <- as.numeric(structure$dz)
    stopifnot(length(depths) == length(theta_init), length(dz) == length(theta_init))
    n_layer <- length(theta_init)
  } else {
    stopifnot(is.list(structure), !is.null(structure$dz))
    dz_loc <- as.numeric(structure$dz)
    n_layer <- 1
    stopifnot(is.finite(dz_loc), dz_loc > 0)
  }
  
  # ---- build forcing (E, T) in mm/step ----
  if (is.null(E) || is.null(T)) {
    if (is.null(ET)) stop("Provide ET (or, alternatively, E and T).")
    if (is.null(LAI)) LAI <- LAI_default
    if (length(LAI) == 1L) LAI <- rep(LAI, length(ET))
    stopifnot(length(LAI) == length(ET))
    E <- pmax(ET * exp(-k_ext * LAI), 0)
    T <- pmax(ET - E, 0)
  } else {
    stopifnot(length(E) == length(T))
  }
  
  # ---- daily → hourly expansion (if needed) ----
  if (timestep == "daily") {
    E <- expand_daily_to_hourly_sinusoidale(E)
    T <- expand_daily_to_hourly_sinusoidale(T)
  }
  steps <- length(E)
  stopifnot(length(T) == steps)
  
  mm_to_m <- function(mm) mm / 1000
  
  if (is_multi) {
    # ===================== MULTI-LAYER =====================
    theta_profile <- as.numeric(theta_init)
    
    sec_state <- sec_create_state(
      theta0  = theta_profile[1],
      Ks      = Ks,
      nVG     = 1/n,
      L       = dz[1],
      theta_s = theta_s,
      theta_r = theta_r,
      flags   = list(stage_two = TRUE, continue_one_after_two = TRUE, legacy_E2 = TRUE)
    )
    
    theta_mat       <- matrix(NA_real_, nrow = n_layer, ncol = steps)
    evap_mm_vec     <- numeric(steps)
    transp_mm_vec   <- numeric(steps)
    drain_layer_mm  <- matrix(0, nrow = n_layer, ncol = steps)
    drain_bottom_mm <- numeric(steps)
    
    theta_avg_5  <- theta_avg_10 <- theta_avg_15 <- theta_avg_20 <- theta_avg_30 <- theta_avg_60 <- numeric(steps)
    
    if (is.null(roots_depth)) roots_depth <- max(depths) + dz[n_layer]/2
    idx_ET <- which(depths <= roots_depth)
    w_ET   <- if (length(idx_ET)) rep(1/length(idx_ET), length(idx_ET)) else numeric(0)
    
    for (t in seq_len(steps)) {
      if (debug) cat(sprintf("\n── timestep %d ─────────────────────────\n", t))
      theta_pre <- theta_profile
      
      # SEC synchronisation (surface layer)
      sec_state$SatT <- (theta_pre[1] - theta_r) / (theta_s - theta_r)
      sec_state$SatT <- min(max(sec_state$SatT, 0), 1)
      
      # 1) drainage per layer
      dtheta_drain  <- numeric(n_layer)
      vol_drain_m   <- numeric(n_layer)
      for (i in seq_len(n_layer)) {
        if (theta_pre[i] <= theta_fc) next
        
        if (exp_drain) {
          s    <- (theta_pre[i] - theta_r) / (theta_s - theta_r)
          s_fc <- (theta_fc     - theta_r) / (theta_s - theta_r)
          beta <- 1/n
          denom <- exp(beta*(1 - s_fc)) - 1
          Ki <- if (denom <= 0) 0 else Ks * (exp(beta*(s - s_fc)) - 1) / denom
        } else {
          Ki <- K_theta(theta_pre[i], Ks, n, theta_s, theta_r, theta_fc) # power-law
        }
        dtheta_drain[i] <- min(max(Ki, 0) * (1 / dz[i]), theta_pre[i] - theta_fc)
        dtheta_drain[i] <- max(dtheta_drain[i], 0)
        vol_drain_m[i]  <- dtheta_drain[i] * dz[i]
        drain_layer_mm[i, t] <- vol_drain_m[i] * 1000
      }
      
      # 2) evaporation (top layer, SEC) + 3) transpiration (distributed)
      E_m <- mm_to_m(E[t])  # m per hour
      T_m <- mm_to_m(T[t])  # m per hour
      
      dtheta_E <- dtheta_T <- numeric(n_layer)
      
      if (E_m > 0) {
        sec_step  <- sec_step_state(sec_state, E0_m = E_m, dt = 1)  # dt=1 hour
        sec_state <- sec_step$state
        theta_surf <- sec_state$SatT * (theta_s - theta_r) + theta_r
        dtheta_E[1] <- max(0, min(theta_pre[1] - theta_surf, theta_pre[1] - theta_r))
      }
      
      if (T_m > 0 && length(idx_ET)) {
        for (k in seq_along(idx_ET)) {
          i <- idx_ET[k]
          alpha <- feddes_stress(theta_pre[i], theta_star, theta_wp)
          T_i   <- (T_m * w_ET[k]) / dz[i] * alpha
          dtheta_T[i] <- min(T_i, max(theta_pre[i] - theta_r, 0))
          dtheta_T[i] <- max(dtheta_T[i], 0)
        }
      }
      
      # 4) update \theta and propagate drainage downwards
      theta_post <- pmax(theta_pre - dtheta_E - dtheta_T - dtheta_drain, theta_r)
      
      residuo <- 0
      for (i in seq_len(n_layer)) {
        residuo <- residuo + vol_drain_m[i]
        if (residuo <= 0) next
        if (i < n_layer) {
          for (j in seq.int(i + 1, n_layer)) {
            space <- (theta_s - theta_post[j]) * dz[j]
            take  <- min(space, residuo)
            theta_post[j] <- theta_post[j] + take / dz[j]
            residuo <- residuo - take
            if (residuo <= 0) break
          }
        }
      }
      drain_bottom_mm[t] <- max(residuo, 0) * 1000
      
      # commit
      theta_profile       <- theta_post
      theta_mat[, t]      <- theta_post
      evap_mm_vec[t]      <- sum(dtheta_E * dz) * 1000
      transp_mm_vec[t]    <- sum(dtheta_T * dz) * 1000
      
      theta_avg_5[t]  <- theta_medio_profondita(theta_profile, depths, dz, 0.05)
      theta_avg_10[t] <- theta_medio_profondita(theta_profile, depths, dz, 0.10)
      theta_avg_15[t] <- theta_medio_profondita(theta_profile, depths, dz, 0.15)
      theta_avg_20[t] <- theta_medio_profondita(theta_profile, depths, dz, 0.20)
      theta_avg_30[t] <- theta_medio_profondita(theta_profile, depths, dz, 0.30)
      theta_avg_60[t] <- theta_medio_profondita(theta_profile, depths, dz, 0.60)
    }
    
    series <- list(
      theta_mm          = theta_mat,
      evap_mm           = evap_mm_vec,
      transp_mm         = transp_mm_vec,
      drain_layer_mm    = drain_layer_mm,
      drain_bottom_mm   = drain_bottom_mm,
      theta_avg_5       = theta_avg_5,
      theta_avg_10      = theta_avg_10,
      theta_avg_15      = theta_avg_15,
      theta_avg_20      = theta_avg_20,
      theta_avg_30      = theta_avg_30,
      theta_avg_60      = theta_avg_60
    )
    
    return(list(
      mode   = "multi",
      series = series,
      meta   = list(units = list(E="mm/step", T="mm/step", theta="m3/m3"),
                    exp_drain=exp_drain, k_ext=k_ext, LAI_default=LAI_default,
                    structure=list(depths=depths, dz=dz))
    ))
  }
  
  # ===================== MONO-LAYER =====================
  theta_post <- as.numeric(theta_init)
  theta_vec  <- numeric(steps)
  
  sec_state <- sec_create_state(
    theta0  = theta_post,
    Ks      = Ks,
    nVG     = 1/n,
    L       = dz_loc,
    theta_s = theta_s,
    theta_r = theta_r,
    flags   = list(stage_two = FALSE, continue_one_after_two = TRUE, legacy_E2 = FALSE)
  )
  
  # --- Initialize outputs ---
  theta_mat       <- matrix(NA_real_, n_layer, steps)  # [layers x steps]
  evap_mm_vec     <- numeric(steps)                     # E per step [mm per timestep]
  transp_mm_vec   <- numeric(steps)                     # T per step [mm per timestep]
  drain_bottom_mm <- numeric(steps)                     # percolazione profonda [mm per timestep]
  
  for (t in seq_len(steps)) {
    if (debug) cat(sprintf("\n── timestep %d ─────────────────────────\n", t))
    theta_pre <- theta_post
    
    # drainage above field capacity
    dtheta_drain <- 0
    if (theta_pre >= theta_fc) {
      if (exp_drain) {
        s    <- (theta_pre - theta_r) / (theta_s - theta_r)
        s_fc <- (theta_fc  - theta_r) / (theta_s - theta_r)
        beta <- 1/n
        denom <- exp(beta*(1 - s_fc)) - 1
        Ki <- if (denom <= 0) 0 else Ks * (exp(beta*(s - s_fc)) - 1) / denom
      } else {
        Ki <- K_theta(theta_pre, Ks, n, theta_s, theta_r, theta_fc)
      }
      dtheta_drain <- min(max(Ki, 0) * (1 / dz_loc), theta_pre - theta_fc)
      dtheta_drain <- max(dtheta_drain, 0)
    }
    
    # evaporation (SEC) and transpiration (Feddes)
    E_m <- mm_to_m(E[t])
    T_m <- mm_to_m(T[t])
    
    sec_state$SatT <- (theta_pre - theta_r) / (theta_s - theta_r)
    sec_state$SatT <- min(max(sec_state$SatT, 0), 1)
    
    dtheta_E <- 0
    if (E_m > 0) {
      sec_step  <- sec_step_state(sec_state, E0_m = E_m, dt = 1)  # dt=1 hour
      sec_state <- sec_step$state
      theta_surf <- sec_state$SatT * (theta_s - theta_r) + theta_r
      dtheta_E <- max(0, min(theta_pre - theta_surf, theta_pre - theta_r))
    }
    
    if (dz_loc <= 0.6) T_m <- T_m * dz_loc / 0.6
    
    dtheta_T <- 0
    if (T_m > 0) {
      alpha  <- feddes_stress(theta_pre, theta_star, theta_wp)
      T_i    <- (T_m * alpha) / dz_loc
      dtheta_T <- min(max(T_i, 0), max(theta_pre - theta_r, 0))
    }
    
    theta_post <- pmax(theta_pre - dtheta_E - dtheta_T - dtheta_drain, theta_r)
    theta_mat[1, t] <- theta_post
    theta_vec[t] <- theta_post
    
    # commit per-step
    evap_mm_vec[t]       <- dtheta_E * dz_loc * 1000
    transp_mm_vec[t]     <- dtheta_T * dz_loc * 1000
    drain_bottom_mm[t]   <- (dtheta_drain * dz_loc) * 1000
  
  }
  
  theta_daily <- NULL
  if (hourly_out) {
    n_days <- floor(steps/24)
    if (n_days >= 1) theta_daily <- tapply(theta_vec, rep(seq_len(n_days), each=24), mean)
  }
  
  series = list(theta = theta_vec, 
                theta_daily = theta_daily,
                theta_mm          = theta_mat,
                evap_mm           = evap_mm_vec,
                transp_mm         = transp_mm_vec,
                drain_bottom_mm   = drain_bottom_mm)
  
  
  list(
    mode   = "mono",
    series = series,
    meta   = list(units = list(E="mm/step", T="mm/step", theta="m3/m3"),
                  exp_drain=exp_drain, k_ext=k_ext, LAI_default=LAI_default,
                  structure=list(dz = dz_loc))
  )
}
