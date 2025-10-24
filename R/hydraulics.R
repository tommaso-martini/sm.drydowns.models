#' Power-law unsaturated hydraulic conductivity function K(θ)
#'
#' Computes the relative hydraulic conductivity around field capacity using
#' a power-law relation between soil saturation and conductivity.
#' The relationship is defined as:
#' \deqn{
#'   K(\theta) = K_s \left(\frac{S - S_{fc}}{1 - S_{fc}}\right)^{p_K},
#' }
#' where \eqn{S = (\theta - \theta_r) / (\theta_s - \theta_r)} and
#' \eqn{S_{fc}} is the relative saturation at field capacity.
#'
#' This simplified power-law form provides a continuous transition between
#' field capacity and saturation, avoiding the exponential tail of the
#' full van Genuchten–Mualem model. It is used in both `muSEC` and `DES`
#' simulations for drainage flux estimation.
#'
#' @param theta current volumetric water content [m3/m3]
#' @param Ks saturated hydraulic conductivity [m/h]
#' @param n van Genuchten n parameter [-]
#' @param theta_s saturated water content [m3/m3]
#' @param theta_r residual water content [m3/m3]
#' @param theta_fc field capacity water content [m3/m3]
#' @param pK optional exponent for the power-law (default \code{1/n})
#'
#' @return hydraulic conductivity [m/h]
#' @export
#'
#' @examples
#' K_theta(theta = 0.25, Ks = 0.5, n = 3, theta_s = 0.4, theta_r = 0.05, theta_fc = 0.2)
#'
K_theta <- function(theta, Ks, n, theta_s, theta_r, theta_fc, pK = 1/n) {
  S     <- (theta - theta_r) / (theta_s - theta_r)
  S_fc  <- (theta_fc - theta_r) / (theta_s - theta_r)
  S     <- max(min(S, 1), 0)
  S_fc  <- max(min(S_fc, 1), 0)
  S_eff <- max(S - S_fc, 0) / max(1 - S_fc, .Machine$double.eps)
  Ks * (S_eff ^ max(pK, 1e-6))
}

#' Feddes stress function
#'
#' Calculates water stress factor for transpiration
#' 
#' 1 above theta_star, 0 below theta_wp, linear in between
#'
#' @param theta current water content [m3/m3]
#' @param theta_star optimal water content (no stress) [m3/m3]
#' @param theta_wp wilting point [m3/m3]
#' @return stress factor alpha (0-1) [-]
feddes_stress <- function(theta, theta_star, theta_wp) {
  if (theta >= theta_star) return(1)
  if (theta <= theta_wp)   return(0)
  (theta - theta_wp) / (theta_star - theta_wp)
}

#' Analytical drainage update
#'
#' Updates theta using analytical solution for drainage above field capacity
#'
#' @param theta0 initial water content [m3/m3]
#' @param dt time step [d]
#' @param Ks saturated hydraulic conductivity [m/d]
#' @param dz layer thickness [m]
#' @param theta_r residual water content [m3/m3]
#' @param theta_s saturated water content [m3/m3]
#' @param n van Genuchten n parameter [-]
#' @return updated water content [m3/m3]
update_theta_analytical <- function(theta0, dt, Ks, dz, theta_r, theta_s, n) {
  alpha <- (n - 1) / n
  C <- Ks / (dz * (theta_s - theta_r)^(1/n))
  y0 <- theta0 - theta_r
  y1 <- max(0, y0^alpha - alpha * C * dt)^(1/alpha)
  return(y1 + theta_r)
}

#' Desorptivity calculation
#'
#' Calculates desorptivity for stage-1 evaporation
#'
#' @param Ks saturated hydraulic conductivity [m/d]
#' @param n van Genuchten n parameter [-]
#' @param theta_s saturated water content [m3/m3]
#' @param theta_r residual water content [m3/m3]
#' @param psi_b bubbling pressure [m]
#' @param Theta relative saturation [-]
#' @return desorptivity S [m/sqrt(d)]
desorptivity <- function(Ks, n, theta_s, theta_r, psi_b, Theta) {
  if (!is.finite(Ks) || !is.finite(psi_b))
    stop("Ks and psi_b must be finite (check units: m/day and m respectively).")
  eps <- 1/n
  num <- 16 * (eps - 3) * Ks * psi_b * (theta_s - theta_r)
  denom <- 3 * (eps + 3) * (eps + 5)
  fact <- Theta^((eps + 5)/4)
  valore <- (num / denom)^(1/2)
  
  valore * fact
}

