# <PACKAGE NAME> Manifold-aware Bayesian Additive Tree Models for Density Estimation on Complex Domains

<!-- badges: start -->
<!--
[![R-CMD-check](https://github.com/<ORG>/<REPO>/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/<ORG>/<REPO>/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/<PKGNAME>)](https://CRAN.R-project.org/package=<PKGNAME>)
-->
<!-- badges: end -->

This repository contains an R package that implements the methods in:

> Diaz-Ray, I., Sang, H., Hu, G., & Lu, L. (2025). *Nonparametric Density Estimation on Complex Domains using Manifold-Aware Bayesian Additive Tree Models*. **Computational Statistics & Data Analysis**, 108335.

The package provides Bayesian nonparametric density estimation tools for data supported on **complex / non-Euclidean domains** (e.g., constrained regions, manifolds/surfaces, networks, and other geometries) using **manifold-aware Bayesian additive tree models**.

---

## Installation

### From GitHub (development version)

```r
# install.packages("remotes")
remotes::install_github("<ORG>/<REPO>")
