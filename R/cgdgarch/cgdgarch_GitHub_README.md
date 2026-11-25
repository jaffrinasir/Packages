# cgdgarch — Sparse ARCH/GARCH estimation in R

**Penalized conditional-volatility modelling and automatic lag selection with coordinate gradient descent.**

`cgdgarch` is an R package for estimating univariate ARCH and GARCH models using **coordinate gradient descent (CGD)** and sparsity-inducing penalties. It is designed for analysts who want to start from a reasonably broad lag specification and determine which ARCH and GARCH coefficients should remain active, without manually comparing every possible subset.

**Release in this repository:** `0.1.14` · **License:** GPL-3 · **R:** 4.1.0 or newer · **Compiled code:** Rcpp / C++17.

> **Source and release archive:** This repository contains the unpacked package source *and* an exact, byte-for-byte copy of the original `cgdgarch_0.1.14.tar.gz` release, renamed **`cgdgarch.tar.gz`**. Renaming the archive does not change the package name or version; installation will still report `cgdgarch` 0.1.14. The repository documentation may be more detailed than the README embedded inside the preserved release archive.

## Why cgdgarch?

Traditional GARCH workflows require a candidate order to be specified in advance. `cgdgarch` focuses on **simultaneous parameter estimation and sparse lag selection**, with tools to choose a regularization parameter, inspect coefficient paths, assess conditional-variance dynamics, and calculate Gaussian conditional risk measures.

The package emphasizes:

- **Six penalty families:** adaptive LASSO, SCAD, MCP, arctan, Laplace, and composite L1–Gaussian (CLG).
- **Rcpp acceleration:** compiled routines for the conditional-variance recursion, Gaussian quasi-likelihood, score and diagonal curvature calculations; a pure-R fallback is available through `use_cpp = FALSE` in the estimator.
- **Regularization tuning:** automatically generated dyadic/logarithmic lambda grids, information-criterion paths, and time-ordered cross-validation.
- **Model inspection:** fitted coefficients, active lags, information criteria, coefficient solution paths, fitted/filtered variance traces, and standardized-error dependence summaries.
- **Empirical comparison:** an optional helper to compare fixed-order, non-penalized fits with `rugarch` and `fGarch`.
- **Risk tools:** Gaussian conditional Value-at-Risk (VaR) and exception-based calibration/backtesting summaries.
- **Reproducibility:** simulation tools, tests, a workflow vignette, bundled R Journal materials, and an optional interactive Shiny dashboard.

## Statistical model and estimation

For a series \(r_t\), the package uses a constant-mean univariate specification,

\[
 r_t = \mu + \varepsilon_t, \qquad
 \varepsilon_t = \sigma_t z_t,
\]

with ARCH/GARCH conditional variance

\[
 \sigma_t^2 = \omega + \sum_{i=1}^{q}\alpha_i\varepsilon_{t-i}^2
                     + \sum_{j=1}^{s}\beta_j\sigma_{t-j}^2.
\]

The intercept \(\omega\) is positive and the lag coefficients are constrained to be nonnegative. A penalty on the lag coefficients promotes sparse estimates; the variance intercept is not penalized. The current fitting procedure uses **Gaussian quasi-likelihood**, so the Gaussian distribution should not be interpreted as having been verified for every application.

The main call is `cgd_garch(y, q, s, penalty, lambda, ...)`. To fit a pure ARCH model, use `cgd_arch(y, q, ...)` (equivalent to `s = 0`).

## Installation

### From the included release archive (recommended for exact version 0.1.14)

Open R or RStudio in the directory containing `cgdgarch.tar.gz`, then run:

```r
install.packages("cgdgarch.tar.gz", repos = NULL, type = "source")
library(cgdgarch)
packageVersion("cgdgarch")  # 0.1.14
```

Or provide its full path, for example:

```r
pkg <- "C:/Users/mjaf8/Documents/Github/Packages/R/cgdgarch R package/cgdgarch.tar.gz"
install.packages(pkg, repos = NULL, type = "source")
```

Since the package includes C++17 code, building from source requires an appropriate compiler toolchain (for example, a compatible Rtools installation for Windows) and the `Rcpp` dependency. See [Windows and other installation instructions](docs/INSTALLATION.md).

### From the unpacked source tree

Inside the repository folder, run:

```r
install.packages(".", repos = NULL, type = "source")
```

If calling R from a shell, the equivalent is `R CMD INSTALL "."`. The unpacked source makes functions, tests, help pages, and C++ routines inspectable and version-controllable.

## Five-minute example

```r
library(cgdgarch)

# Simulate a sparse GARCH(3,1) process: the second ARCH lag is zero.
set.seed(123)
sim <- simulate_garch(
  n = 1000,
  omega = 0.02,
  alpha = c(0.12, 0, 0.06),
  beta = 0.78,
  return_details = TRUE
)
y <- sim$y

# Fit a regularization path and let BIC select lambda.
pen <- penalty_mcp(gamma = 2.5)
ic <- cgd_lambda_ic_path(
  y = y, q = 3, s = 1,
  penalty = pen,
  lambda_method = "log",
  nlambda = 40,
  lambda_min_ratio = 1e-4,
  fit_unpenalized = FALSE,
  max_iter = 200
)

# Inspect the BIC-selected model and fitted coefficients.
ic$best
fit <- ic$fits[[ic$best$index]]
coef(fit)
summary(fit)

# Information-criterion profile.
plot(ic)

# Conditional VaR forecast (Gaussian quantile assumption).
cgd_var(fit, level = c(0.95, 0.99), n.ahead = 5)
```

**Interpreting the output:** zero lag coefficients indicate positions that are inactive in the reported solution. The BIC-selected tuning parameter is identified by `ic$best`. Results from a single simulated realization are examples of package use, **not** evidence of an estimator's repeated-sample selection probability.

## Penalty families

| Penalty | Constructor | Main idea |
|---|---|---|
| Adaptive LASSO | `penalty_alasso(weights = NULL, gamma = 2)` | Weighted L1 shrinkage and lag selection. |
| SCAD | `penalty_scad(a = 3.7)` | Nonconvex penalty with less shrinkage of sufficiently large coefficients. |
| MCP | `penalty_mcp(gamma = 3)` | Nonconvex penalty whose slope decreases with coefficient magnitude. |
| Arctan | `penalty_arctan(scale = 0.1)` | Bounded nonconvex penalty with a positive transition scale. |
| Laplace | `penalty_laplace(scale = 0.5)` | Bounded exponential-type penalty with a positive transition scale. |
| Composite L1–Gaussian | `penalty_clg(alpha = 0.7, kappa = 2)` | L1 component plus smooth Gaussian-type shrinkage. |

The **arctan and Laplace** constructors both take `scale`, not the older `gamma`/`eps` names. For CLG, a zero L1 mixture weight makes the penalty smooth rather than an exact subset-selection penalty.

## Primary functions

| Area | Functions | Purpose |
|---|---|---|
| Fit | `cgd_garch()`, `cgd_arch()` | Penalized or unpenalized estimation for a specified `lambda`. |
| Tuning | `lambda_grid_cgd()` | Create dyadic or logarithmic regularization grids. |
| Model selection | `cgd_lambda_ic_path()` | Compute IC paths and BIC-selected fits. |
| Validation | `cv_cgd_garch()` | Expanding-window cross-validation. |
| Paths | `solution_path()` | Estimate and inspect parameter paths across lambda values. |
| Simulation | `simulate_garch()` | Generate ARCH/GARCH series; optionally return true latent variance. |
| Diagnostics | `cgd_diagnostic_summary()`, `cgd_error_dependence()` | Summarize standardized-error checks and ACF/PACF data. |
| Variance plots | `cgd_variance_trace()`, `plot_variance_trace()` | Compare stored fitted and explicitly filtered variance traces. |
| Benchmarking | `cgd_compare_garch()` | Optional comparison to `rugarch` and `fGarch`. |
| Risk | `cgd_var()`, `cgd_var_backtest()` | Gaussian conditional VaR, coverage and independence tests. |
| Interactive display | `cgd_shiny_app()` | Optional Shiny dashboard. |

Additional standard R methods include `coef()`, `fitted()`, `residuals()`, `predict()`, `plot()`, `print()`, and `summary()` for supported fitted objects.

## Common workflows

### 1. Compare the regularization paths

```r
sp <- solution_path(
  y, q = 3, s = 1,
  penalty = penalty_scad(),
  nlambda = 40,
  lambda_method = "log",
  fit_unpenalized = FALSE,
  max_iter = 200
)
plot(sp, x_axis = "log_lambda", include_intercept = FALSE)
```

Lambda values and coefficient matrices are also available directly as `sp$lambda` and `sp$coefficients`. If a plot requires customization or a particular package version has a plotting-method issue, the stored values can be plotted using ordinary base R graphics without refitting the models.

### 2. Examine diagnostics

```r
dg <- cgd_diagnostic_summary(fit)
print(dg)

dep <- cgd_error_dependence(fit, lag.max = 24)
plot_error_acf(dep)
plot_error_pacf(dep)

trace <- cgd_variance_trace(fit)
head(trace)
plot_variance_trace(trace)
```

The diagnostic p-values and approximate uncertainty screens are **not** formal post-selection inference. Standardized errors/innovations can retain dependence or depart from normality even if the GARCH fit appears plausible.

### 3. Use an empirical data series

```r
# Built into base R; no external download required.
dax <- 100 * diff(log(EuStockMarkets[, "DAX"]))

empirical <- cgd_lambda_ic_path(
  dax, q = 2, s = 2,
  penalty = penalty_arctan(scale = 0.1),
  nlambda = 40,
  lambda_method = "log",
  fit_unpenalized = FALSE
)
selected <- empirical$fits[[empirical$best$index]]
coef(selected)
```

The candidate order is a modelling choice; an ACF/PACF of squared returns can motivate a low-order candidate but does not uniquely identify the ARCH and GARCH orders.

### 4. Compare against other R packages (optional)

```r
# Install rugarch/fGarch separately if those comparisons are needed.
cmp <- cgd_compare_garch(dax, q = 2, s = 2, strict = FALSE)
print(cmp)
```

A fixed-order coefficient may be zero because of a nonnegativity boundary or reporting tolerance, **not** because a fixed-order benchmark performed penalized selection. Exact equality across engines is not expected because initialization and optimizer conventions differ.

### 5. Evaluate Gaussian VaR

```r
v <- cgd_var(selected, level = c(0.95, 0.99), n.ahead = 5)
print(v)

# This is an in-sample diagnostic, not a forecasting backtest.
bt <- cgd_var_backtest(selected, level = 0.99)
print(bt)
```

For an out-of-sample VaR evaluation, provide holdout `realized` returns and aligned one-step-ahead `sigma` forecasts to `cgd_var_backtest()`. Gaussian-tail VaR needs particular caution if standardized innovations have heavy tails.

## Limitations and appropriate interpretation

- The estimator is **univariate** and based on **Gaussian quasi-likelihood**; Student-t and skewed innovations, asymmetric GARCH extensions, and multivariate GARCH are not provided by this release.
- ARCH/GARCH coefficients are constrained to be nonnegative, but the optimizer does **not** impose \(\sum_i\alpha_i+\sum_j\beta_j<1\) as a hard stationarity constraint. Check persistence if covariance stationarity is required.
- The optimization problem is **nonconvex** for several supported penalties. Fitted solutions should not be treated as guaranteed global optima.
- Fixed-order estimates may land at zero because of parameter boundaries; tiny positive estimates may also be printed as zero at limited decimal precision.
- Information-criterion and diagnostic summaries are useful model-development tools, not proofs that the chosen GARCH order or conditional distribution is correct.
- `cgd_var_backtest(fit)` performs **in-sample** calibration by default, not an out-of-sample predictive assessment.

## Source layout

```text
cgdgarch R package/
├── DESCRIPTION              Package metadata, version and dependencies
├── NAMESPACE                Registered exported functions and methods
├── R/                       Estimation, penalties, tuning, diagnostics, risk
├── src/                     Rcpp/C++17 implementation
├── man/                     Function help pages
├── tests/                   Automated tests
├── vignettes/               Worked workflow
├── inst/                    Shiny dashboard, citation and article material
├── NEWS.md                  Package release history
├── README.md                This project overview
├── cgdgarch.tar.gz          Exact original 0.1.14 archive, renamed
└── docs/                    Additional usage, installation and developer guides
```

See [Installation](docs/INSTALLATION.md), [User guide](docs/USER_GUIDE.md), [Development and testing](docs/DEVELOPMENT.md), and [Release history](NEWS.md).

## Research background and citation

The adaptive-LASSO CGD methodology underlying the package is described in:

> Nasir, M. J. M., Khan, R. N., Nair, G., & Nur, D. (2025). Coordinate gradient descent algorithm in adaptive LASSO for pure ARCH and pure GARCH models. *Computational Statistics*, **40**, 3527–3561. https://doi.org/10.1007/s00180-025-01642-1

For a machine-readable citation bundled with the installed package, use:

```r
citation("cgdgarch")
```

## Maintainer and contributors

**Maintainer:** Muhammad Jaffri Mohd Nasir — [jaffri.mn@umk.edu.my](mailto:jaffri.mn@umk.edu.my)

**Contributors:** Gopalan Nair, Darfiana Nur, and Ahmad Hakiim Jamaluddin. See `DESCRIPTION` for the package's formal author/contributor roles.

## License

This package declares the **GNU General Public License, version 3 (GPL-3)** in `DESCRIPTION`. Review the licensing obligations before redistribution or incorporation into another project.
