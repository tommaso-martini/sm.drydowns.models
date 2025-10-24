#' Expand daily totals to hourly with a sinusoidal shape
#'
#' Multiplies each daily value by a normalized 24-hour sinusoidal profile so that
#' the sum over 24 hours equals the original daily input.
#'
#' @param daily_values Numeric vector of daily totals (e.g., mm/day).
#' @return Numeric vector of length 24 * length(daily_values) in the same units, per hour.
#' @export
expand_daily_to_hourly_sinusoidale <- function(daily_values) {
  profile <- sinusoidal_profile()
  unlist(lapply(daily_values, function(x) x * profile))
}

#' Normalized sinusoidal hourly profile (24 values)
#'
#' @param hours Number of hourly slots (default 24).
#' @return Numeric vector summing to 1.
#' @keywords internal
sinusoidal_profile <- function(hours = 24) {
  t <- seq(0, 2 * pi, length.out = hours)
  prof <- (sin(t - pi/2) + 1) / 2
  prof / sum(prof)
}

#' Build multi-layer geometry from target depth and two thicknesses
#'
#' @param Z_t Total target depth (m), e.g. 0.60.
#' @param dz_top Thickness (m) for the first `n_top_layers` near-surface layers.
#' @param dz_rest Thickness (m) repeated below until reaching `Z_t`.
#' @param n_top_layers Integer, number of near-surface layers with `dz_top`.
#' @return A list with `depths` (centers, m) and `dz` (thicknesses, m).
#' @export
compute_dz_depths <- function(Z_t, dz_top, dz_rest, n_top_layers){
  depths <- seq(dz_top / 2, by = dz_top, length.out = n_top_layers)
  z_curr <- n_top_layers * dz_top
  while ((z_curr + dz_rest) <= Z_t + 1e-6) {
    center <- z_curr + dz_rest / 2
    depths <- c(depths, center)
    z_curr <- z_curr + dz_rest
  }
  dz <- c(rep(dz_top, n_top_layers), rep(dz_rest, length(depths) - n_top_layers))
  list(dz = dz, depths = depths)
}

#' Infer a theta profile from cumulative depth-averaged targets
#'
#' @param depths Vector of layer centers (m).
#' @param dz Vector of layer thicknesses (m), same length as `depths`.
#' @param avg_targets Named numeric like c("0.30"=0.28, "0.60"=0.30) with depths in meters.
#' @return Numeric vector `theta` per layer.
#' @export
infer_profile_from_averages <- function(depths, dz, avg_targets) {
  target_depths <- as.numeric(names(avg_targets))
  cum_means     <- as.numeric(avg_targets)
  inc_means <- cum_means
  if (length(cum_means) >= 2) {
    for (i in 2:length(cum_means)) {
      z0 <- target_depths[i - 1]; z1 <- target_depths[i]
      inc_means[i] <- (cum_means[i] * z1 - cum_means[i - 1] * z0) / (z1 - z0)
    }
  }
  theta_out <- numeric(length(depths))
  for (k in seq_along(depths)) {
    z <- depths[k]
    j <- which(z <= target_depths)[1]
    if (is.na(j)) j <- length(target_depths)
    theta_out[k] <- inc_means[j]
  }
  theta_out
}

#' Depth-averaged theta up to z_max
#'
#' @param theta Layer thetas (m3/m3).
#' @param depths Layer centers (m).
#' @param dz Layer thicknesses (m).
#' @param z_max Upper depth (m) for averaging (e.g., 0.30).
#' @return Scalar depth-averaged theta in [0, z_max].
#' @export
theta_medio_profondita <- function(theta, depths, dz, z_max) {
  top    <- depths - dz/2
  bottom <- depths + dz/2
  frac <- pmax(0, pmin(bottom, z_max) - top) / dz
  w <- frac * dz
  sum(theta * w) / sum(w)
}