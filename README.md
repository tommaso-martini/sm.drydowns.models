# sm.drydowns.models — DES and muSEC Models for Soil Moisture Drydowns

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

This R package provides the **reference implementation** of two simplified hydrological models for soil moisture drydowns:

- **DES** — the *Desorptivity-based* model;
- **muSEC** — the *multilayer Surface Evaporative Capacitor* extension accounting for vertical variability and transpiration.

The package enables reproducible ecohydrological simulations, sensitivity analyses, and teaching applications.  
It reproduces model behaviour from the associated preprint:

📄 Martini, T., Gentile, A., Gisolo, D., Olivero, A., & Ferraris, S. (2025).  
*From near surface to root zone soil water losses: a new model validated with field TDR and remotely sensed data.*  
**Preprint:** [https://doi.org/10.2139/ssrn.5648370](https://doi.org/10.2139/ssrn.5648370)

---

## 🚀 Installation

Install locally (while developing):

```r
# install.packages("devtools")  # if not already installed
devtools::install_local(".")
library(sm.drydowns.models)
```

Once the repository is public, you can install directly from GitHub:

```r
devtools::install_github("tommasomartini/sm.drydowns.models")
```

---

## 💧 Quick examples

### 1️⃣ DES (mono-layer drydown)

```r
res_des <- des_simulate(
  Ks = 1.0, n = 1/5, theta_init = 0.30,
  ET_vec = rep(4, 5), rain_vec = rep(0, 5),
  theta_s = 0.40, theta_r = 0.005,
  theta_fc = 0.17, theta_wp = 0.01, theta_star = 0.16,
  psi_b = 0.08, dt = 1/24, dz = 0.10
)

plot(res_des$series$theta, type = "l",
     xlab = "Hour (sub-steps)", ylab = expression(theta),
     main = "DES drydown (constant ET = 4 mm/day)")
```

### 2️⃣ muSEC (multi-layer)

```r
theta0 <- 0.30
dz <- 0.10
H  <- 24 * 5
ET <- rep(4/24, H)

res_musec <- musec_simulate(
  theta_init = theta0,
  structure  = list(dz = dz),
  Ks = 5/24, n = 1/11,
  theta_s = 0.40, theta_r = 0.005,
  theta_fc = 0.17, theta_wp = 0.01, theta_star = 0.16,
  ET = ET, timestep = "hourly",
  exp_drain = TRUE, hourly_out = TRUE
)

plot(res_musec$series$theta, type = "l",
     xlab = "Hour", ylab = expression(theta),
     main = "muSEC mono-layer (constant ET)")
```

For full, commented examples (ET partitioning, vertical profiles, visualization):

```r
vignette("musec-examples-annotated", package = "sm.drydowns.models")
```

---

## 📖 Citation

If you use this package, please cite:

> Martini, T. & Olivero, A. (2025).  
> *sm.drydowns.models: DES and muSEC Models for Soil Moisture Drydowns (R package).*  
> Version 1.0.0. MIT License.  
> [https://doi.org/10.2139/ssrn.5648370](https://doi.org/10.2139/ssrn.5648370)

---

© 2025 Tommaso Martini & Aurora Olivero —  
Interuniversity Department of Regional and Urban Studies and Planning (DIST), University of Turin  
Licensed under the [MIT License](LICENSE).
