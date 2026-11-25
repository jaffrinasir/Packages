# Packages

A collection of open-source research and teaching packages in **R, Python,
and related computational tools**.

## Featured R package: `cgdgarch`

[`cgdgarch`](R/cgdgarch/) implements **penalized ARCH/GARCH estimation
using coordinate gradient descent (CGD)** for sparse volatility modelling
in R. It supports penalized lag selection, tuning and solution-path
visualization, diagnostics, and Gaussian conditional risk calculations.

**Latest documented version: `0.1.14`** · **GPL-3** ·
**R ≥ 4.1.0** · **Rcpp / C++17**

### What's new in `cgdgarch` 0.1.14?

- Better information-criterion plots: legend colours and line types agree
  with AIC, BIC, Hannan–Quinn, and optional Shibata curves.
- Improved coefficient solution-path plots with `lambda` or
  `log(lambda)` on the horizontal axis.
- BIC-selected lambda reference lines and the option to hide the
  unpenalized intercept from coefficient-path displays.
- Coefficient matrices stored by information-criterion path objects,
  avoiding the need to refit models for path inspection.
- Expanded documentation and tests for plotting and path inspection.

These are **version 0.1.14 features**, not a claim of a newer release.
See the package [changelog](R/cgdgarch/NEWS.md) for the full history.

### Explore and install

| Resource | Location |
|---|---|
| Package overview, examples and API guide | [R/cgdgarch/README.md](R/cgdgarch/README.md) |
| Source code and `DESCRIPTION` | [R/cgdgarch/](R/cgdgarch/) |
| Release history | [R/cgdgarch/NEWS.md](R/cgdgarch/NEWS.md) |
| Source release archive | [R/cgdgarch/cgdgarch.tar.gz](R/cgdgarch/cgdgarch.tar.gz) |
| Detailed installation guide | [R/cgdgarch/docs/INSTALLATION.md](R/cgdgarch/docs/INSTALLATION.md) |

To install the included source archive locally (with an appropriate
R/C++ toolchain):

```r
install.packages(
  "R/cgdgarch/cgdgarch.tar.gz",
  repos = NULL,
  type = "source"
)
library(cgdgarch)
packageVersion("cgdgarch")
```

The path above is relative to the repository root. Alternatively, supply
an absolute file path to the archive.

### Main capabilities

`cgdgarch` currently offers adaptive LASSO, SCAD, MCP, arctan, Laplace,
and composite L1–Gaussian penalties; automatic lambda grids and
information-criterion paths; time-ordered cross-validation; simulation
and standardized-innovation diagnostics; comparisons with `rugarch` and
`fGarch` (optional); and Gaussian conditional Value-at-Risk tools.

Its estimation scope is **univariate Gaussian quasi-likelihood**.
Non-Gaussian innovations and multivariate GARCH extensions are not part
of version 0.1.14.

## Maintainer

**Muhammad Jaffri Mohd Nasir** — Universiti Malaysia Kelantan  
[jaffri.mn@umk.edu.my](mailto:jaffri.mn@umk.edu.my)

Formal package authorship and contributor roles are in
[`R/cgdgarch/DESCRIPTION`](R/cgdgarch/DESCRIPTION).
