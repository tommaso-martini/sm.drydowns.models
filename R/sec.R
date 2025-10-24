sec_step_state <- function(state, E0_m, dt) {
  with(state, {
    kappa <- Ks; gamma <- nVG; deltL <- L * (theta_s - theta_r); grad <- 1
    beta <- gamma
    K <- Ks * 1/(exp(beta) - 1) * (exp(beta * (SatT)) - 1)
    evap_I <- E0_m * K * grad / (E0_m + K * grad)   # Stage I
    E_due_0 <- min(1/24/1000, E0_m/2) * dt
    y <- (E_due_0 * E0_m * dt) / (E0_m * dt - E_due_0)
    SatC <- log(y/kappa * dt * (exp(beta) - 1) + 1) / beta
    
    stage <- if (SatT > SatC) "I" else "II"
    stage <- if (stage2_reached) "II" else "I"
    
    evap_rate <- evap_I
    
    if (stage_two && stage == "II" && SatT > 0) {
      entering_now <- (tII == 0)
      if (entering_now) {
        E2_0   <- if (legacy_E2) min(1/24/1000, E0_m/2) else evap_I
        xi     <- XI0_MM/1000
        theta_eff <- SatT*(theta_s-theta_r) + theta_r
        dtheta <- theta_eff - theta_r/2
        tII    <- dt
      } else {
        tII    <- tII + dt
      }
      evap_rate <- flux_Es_II_B2(tII, dt, E2_0, xi, dtheta)
      
      if(entering_now)
        evap_rate <- E2_0
      
      if (!continue_one_after_two) stop_evap_after_stage2 <- TRUE
    }
    
    if (!stage_two && stage == "II") {
      if (!continue_one_after_two) {
        evap_rate <- 0
      } else {
        evap_rate <- evap_I
      }
    }
    
    cumE_mm <- cumE_mm + evap_rate * 1000 * dt
    SatT <- max(SatT - ((evap_rate) * dt) / deltL, 0)
    
    if (!stage2_reached && SatT <= SatC) {
      stage2_reached <- TRUE
      if (!stage_two && !continue_one_after_two) stop_evap_after_stage2 <- TRUE
    }
    
    state_out <- state
    state_out$SatT    <- SatT
    state_out$cumE_mm <- cumE_mm
    state_out$tII     <- tII
    state_out$E2_0    <- E2_0
    state_out$xi      <- xi
    state_out$dtheta  <- dtheta
    state_out$stage2_reached <- stage2_reached
    
    list(evap_mm = evap_rate * 1000 * dt, SatT = SatT, state = state_out)
  })
}

sec_create_state <- function(theta0, Ks, nVG, L, theta_s, theta_r,
                             flags = list(stage_two = TRUE,
                                          continue_one_after_two = TRUE,
                                          legacy_E2 = FALSE),
                             XI0_MM = 10) {
  list(
    SatT   = (theta0 - theta_r) / (theta_s - theta_r),
    Sat0   = NA_real_,
    cumE_mm = 0, cumQ_mm = 0,
    tII    = 0L,
    E_cum  = 1e-6,
    E2_0   = NA_real_,
    xi     = NA_real_,
    dtheta = NA_real_,
    
    stage2_reached          = FALSE,
    continue_one_after_two  = flags$continue_one_after_two,
    
    Ks = Ks, nVG = nVG, L = L,
    theta_s = theta_s, theta_r = theta_r,
    stage_two = flags$stage_two,
    legacy_E2 = flags$legacy_E2,
    XI0_MM    = XI0_MM
  )
}

# Stage-II B2 kernels (yours)
flux_Es_II_B2 <- function(t, dt, E2_0, xi, de_theta) cum_Es_II_B2(t,E2_0,xi,de_theta) - cum_Es_II_B2(t-dt,E2_0,xi,de_theta)

cum_Es_II_B2  <- function(t, E2_0, xi, de_theta) if (t <= 0) 0 else sqrt(xi)/sqrt(xi + 2*E2_0*t/de_theta)*(2*E2_0*t + de_theta*xi)
