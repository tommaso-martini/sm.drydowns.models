#' Stair-step soil moisture profile (theta vs depth)
#'
#' Draws a layer-wise “stair” profile of volumetric water content \eqn{\theta}
#' versus depth, using either filled rectangles or classic segment outlines.
#' If layer thicknesses (`dz`) are not provided, top/bottom boundaries are
#' inferred from layer centers (`depths`) via midpoints.
#'
#' @param theta Numeric vector of volumetric water content per layer
#'   (m^3 m^-3). Length must equal `length(depths)`.
#' @param depths Numeric vector of layer centers (m), topmost near 0.
#' @param col Color for profile drawing. Default `"black"`.
#' @param lwd Line width for outlines/spines. Default `2`.
#' @param add Logical; if `TRUE`, add to existing plot. Default `FALSE`.
#' @param xlim Numeric length-2; x-axis limits for \eqn{\theta}. Default `c(0, 0.40)`.
#' @param dz Optional numeric vector of layer thicknesses (m); if `NULL`,
#'   the function infers boundaries from `depths`.
#' @param border Rectangle border color when `style="rect"`. Default `NA`.
#' @param alpha Fill transparency (0–1) for `style="rect"`. Use `NA` to disable
#'   transparency (i.e., solid `col`). Default `0.15`.
#' @param style Either `"rect"` (filled bars + spine) or `"segments"`
#'   (classic stair outline). Default `"rect"`.
#' @param xlab,ylab Axis labels. Defaults to expression(theta) and `"Depth (m)"`.
#' @param ... Further graphical arguments passed to `plot()` when `add=FALSE`.
#'
#' @details
#' If `dz` is missing, the function constructs layer top/bottom boundaries
#' by taking midpoints between adjacent `depths`; the domain is clamped to
#' start at 0 m depth for the top layer. When `dz` is supplied, boundaries are
#' simply `depths ± dz/2`.
#'
#' @return Invisibly returns a list with two numeric vectors:
#'   \item{top}{Layer top depths (m).}
#'   \item{bottom}{Layer bottom depths (m).}
#'
#' @examples
#' # Three layers in the top 30 cm
#' depths <- c(0.05, 0.15, 0.25)
#' dz     <- rep(0.10, 3)
#' theta_t0 <- c(0.26, 0.24, 0.22)
#' theta_t1 <- c(0.23, 0.215, 0.205)
#'
#' # Segments style: t1 baseline (black) with t0 overlay (blue)
#' plot_profile_stairs(theta_t1, depths, dz = dz,
#'                     style = "segments", xlim = c(0.18, 0.30),
#'                     main = "Stair profile — t1 vs t0")
#' plot_profile_stairs(theta_t0, depths, dz = dz, style = "segments",
#'                     add = TRUE, col = "blue")
#'
#' # Rectangles with transparency
#' plot_profile_stairs(theta_t0, depths, dz = dz,
#'                     style = "rect", col = "blue", alpha = 0.25,
#'                     xlim = c(0.18, 0.30), main = "Rectangular fill")
#' plot_profile_stairs(theta_t1, depths, dz = dz,
#'                     style = "segments", add = TRUE, col = "black")
#'
#' @seealso [plot_daily_profiles()] for automated day-by-day overlays with titles.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics plot rect segments legend polygon points lines abline matplot
#' @importFrom grDevices adjustcolor rgb
#' @importFrom utils head tail
#' @export
plot_profile_stairs <- function(theta, depths,
                                col = "black", lwd = 2, add = FALSE,
                                xlim = c(0, 0.40),
                                dz = NULL, border = NA, alpha = 0.15,
                                style = c("rect","segments"),
                                xlab = expression(theta), ylab = "Depth (m)", ...) {
  stopifnot(is.numeric(theta), is.numeric(depths))
  stopifnot(length(theta) == length(depths))
  style <- match.arg(style)
  
  # --- 1) build layer boundaries and thicknesses ---
  if (is.null(dz)) {
    # Boundaries from centers: midpoints between adjacent centers; symmetric extension at ends
    zc <- depths
    if (length(zc) == 1L) {
      # single layer: thickness unknown; assume dz = 2*zc (top=0)
      top <- max(0, zc - zc)          # 0
      bot <- zc + zc                   # 2*zc
    } else {
      mid <- (head(zc, -1) + tail(zc, -1)) / 2
      top <- c(max(0, 2*zc[1] - mid[1]), mid)  # extrapolate top; enforce >= 0
      top[1] <- 0
      bot <- c(mid, 2*zc[length(zc)] - mid[length(mid)])
    }
    # thickness per layer
    dz_loc <- bot - top
    # top/bottom per layer
    top_layer <- top
    bottom_layer <- bot
  } else {
    stopifnot(length(dz) == length(depths))
    dz_loc <- dz
    top_layer    <- depths - dz_loc/2
    bottom_layer <- depths + dz_loc/2
    # enforce top >= 0
    top_layer[1] <- max(0, top_layer[1])
  }
  
  # --- 2) init plot window ---
  y_limits <- rev(range(c(top_layer, bottom_layer), finite = TRUE))
  if (!add) {
    plot(NULL, xlim = xlim, ylim = y_limits, xlab = xlab, ylab = ylab, ...)
  }
  
  # --- 3) draw profile ---
  if (style == "rect") {
    # semi-transparent rectangles for clarity
    col_fill <- if (is.na(alpha)) col else grDevices::adjustcolor(col, alpha.f = alpha)
    for (i in seq_along(theta)) {
      rect(xleft = min(xlim[1], theta[i]), xright = theta[i],
           ybottom = bottom_layer[i], ytop = top_layer[i],
           col = col_fill, border = border, lwd = lwd)
      # draw the vertical spine (optional)
      segments(x0 = theta[i], x1 = theta[i], y0 = top_layer[i], y1 = bottom_layer[i],
               col = col, lwd = lwd)
      if (i < length(theta)) {
        # horizontal connector on the boundary between layers
        segments(x0 = theta[i], x1 = theta[i+1],
                 y0 = bottom_layer[i], y1 = bottom_layer[i],
                 col = col, lwd = lwd)
      }
    }
  } else { # "segments" — classic stair outline (your original intent)
    for (i in seq_along(theta)) {
      segments(x0 = theta[i], x1 = theta[i],
               y0 = top_layer[i], y1 = bottom_layer[i],
               col = col, lwd = lwd)
      if (i < length(theta)) {
        segments(x0 = theta[i], x1 = theta[i+1],
                 y0 = bottom_layer[i], y1 = bottom_layer[i],
                 col = col, lwd = lwd)
      }
    }
  }
  
  invisible(list(top = top_layer, bottom = bottom_layer))
}

#' Plot daily soil moisture profiles (t vs t-1) with optional rain/ET titles
#'
#' @param res   List returned by `musec_simulate()`. Must contain
#'              `res$series$theta_mm` [layers x steps] and `res$series$evap_mm` [steps].
#' @param depths Numeric vector of layer centers [m], length = n_layers.
#' @param dz     Numeric vector of layer thicknesses [m], length = n_layers.
#' @param rain   Optional numeric vector [steps] with rainfall per time step [mm].
#' @param ET0    Optional numeric vector [steps] with ET0 per time step [mm].
#' @param file   Optional PDF file path. If non-NULL, plots are saved into a PDF.
#'               If NULL, plots are sent to the current device.
#' @param days   Optional integer vector of day indices to plot. Defaults to 2:n_days
#'               (day 1 is skipped because it needs t-1 to compare).
#' @param xlim   X-axis limits for theta (default c(0, 0.40)).
#' @param style  "segments" (classic stair outline) or "rect" (filled bars + spine).
#' @param add_areas Logical. If TRUE, when daily rain==0, highlight per-layer water loss
#'               area between t and t-1 (semi-transparent).
#' @param alpha  Alpha for filled areas when add_areas=TRUE and style="rect".
#' @param main_prefix Character prefix used in plot titles.
#'
#' @return (Invisibly) a list with meta info (n_days, plotted days). Primarily called for its side-effect (plots).
#' @export
plot_daily_profiles <- function(res, depths, dz,
                                rain = NULL, ET0 = NULL,
                                file = NULL, days = NULL,
                                xlim = c(0, 0.40),
                                style = c("segments", "rect"),
                                add_areas = FALSE, alpha = 0.30,
                                main_prefix = "Daily profile") {
  style <- match.arg(style)
  
  # ---- Basic checks ----
  stopifnot(is.list(res),
            is.numeric(depths), is.numeric(dz),
            length(depths) == length(dz))
  
  if (is.null(res$series$theta_mm))
    stop("`res$series$theta_mm` not found. Did you run musec_simulate() with standard output?")
  if (is.null(res$series$evap_mm))
    warning("`res$series$evap_mm` not found. Titles will omit E_cum.")
  
  theta_mat <- res$series$theta_mm   # [layers x steps]
  if (!is.matrix(theta_mat)) stop("`res$series$theta_mm` must be a matrix [layers x steps].")
  if (nrow(theta_mat) != length(depths))
    stop("nrow(res$series$theta_mm) != length(depths).")
  if (!is.numeric(theta_mat)) storage.mode(theta_mat) <- "numeric"
  
  steps  <- ncol(theta_mat)
  n_days <- floor(steps / 24L)
  if (n_days < 2L) stop("Not enough steps to form >= 2 full days.")
  
  # Time series (optional) checks
  have_rain <- !is.null(rain)
  have_ET0  <- !is.null(ET0)
  
  if (have_rain) {
    if (!is.numeric(rain) || length(rain) != steps)
      stop("`rain` must be numeric vector of length == steps.")
  }
  if (have_ET0) {
    if (!is.numeric(ET0) || length(ET0) != steps)
      stop("`ET0` must be numeric vector of length == steps.")
  }
  
  # Evap cumulative (for title)
  E_series <- res$series$evap_mm
  have_E   <- is.numeric(E_series) && length(E_series) == steps
  E_cum    <- if (have_E) cumsum(E_series) else rep(NA_real_, steps)
  
  # Days to plot: default 2:n_days (needs t-1 to compare)
  if (is.null(days)) {
    days <- 2:n_days
  } else {
    days <- as.integer(days)
    days <- days[days >= 2L & days <= n_days]
    if (length(days) == 0L) stop("No valid day indices to plot (must be in 2..n_days).")
  }
  
  # Prepare device
  if (!is.null(file)) {
    grDevices::pdf(file, width = 7, height = 9)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  # ---- Loop days ----
  for (d in days) {
    t0 <- (d - 1L) * 24L     # last hour of day d-1
    t1 <- d * 24L            # last hour of day d
    idx <- (t0 + 1L):t1      # the 24 hours of day d
    
    # Daily sums (if series provided)
    rain_day <- if (have_rain) sum(rain[idx], na.rm = TRUE) else NA_real_
    ET0_day  <- if (have_ET0)  sum(ET0[idx],  na.rm = TRUE) else NA_real_
    Ecum_t1  <- if (have_E)    E_cum[t1] else NA_real_
    
    # Build title
    title_parts <- c(sprintf("%s %d", main_prefix, d))
    if (have_rain) title_parts <- c(title_parts, sprintf("P: %.1f mm", rain_day))
    if (have_ET0)  title_parts <- c(title_parts, sprintf("ET0: %.1f mm", ET0_day))
    if (have_E)    title_parts <- c(title_parts, sprintf("E_cum: %.1f mm", Ecum_t1))
    main_title <- paste(title_parts, collapse = "  —  ")
    
    # Extract profiles
    theta_t0 <- theta_mat[, t0]
    theta_t1 <- theta_mat[, t1]
    
    # Draw base profile (t1) and overlay (t0)
    bb <- plot_profile_stairs(theta_t1, depths, dz = dz,
                              col = "black", lwd = 2,
                              xlim = xlim, style = style,
                              main = main_title)
    plot_profile_stairs(theta_t0, depths, dz = dz,
                        col = "blue", lwd = 2,
                        add = TRUE, style = style)
    
    # Optional water-loss areas (only when no rain that day)
    if (add_areas && have_rain && isTRUE(rain_day == 0)) {
      top    <- bb$top
      bottom <- bb$bottom
      for (i in seq_along(depths)) {
        if (!is.na(theta_t1[i]) && !is.na(theta_t0[i]) && theta_t1[i] < theta_t0[i]) {
          polygon(
            x = c(theta_t1[i], theta_t0[i], theta_t0[i], theta_t1[i]),
            y = c(top[i],      top[i],      bottom[i],   bottom[i]),
            col = grDevices::rgb(0, 0, 1, alpha), border = NA
          )
        }
      }
    }
    
    legend("bottomright",
           legend = c(sprintf("t=%d (end day %d)", t1, d),
                      sprintf("t=%d (end day %d)", t0, d-1)),
           col = c("black", "blue"), lwd = 2, lty = 1, bty = "n")
  }
  
  invisible(list(n_days = n_days, plotted = days))
}