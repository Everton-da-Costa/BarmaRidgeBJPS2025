# BarmaRidgeBJPS2025

[![Status](https://img.shields.io/badge/Status-Submitted-lightgrey.svg)](https://projecteuclid.org/journals/brazilian-journal-of-probability-and-statistics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the R package and associated data for the scientific article:

**"Numerical stability enhancements in beta autoregressive moving average model estimation"** by Cribari-Neto, F., Costa, E., and Fonseca, R. V.
Submitted to the *Brazilian Journal of Probability and Statistics*. 

---

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [✨ Key Features](#-key-features)
- [📂 Repository Structure](#-repository-structure)
- [🛠️ Installation](#️-installation)
- [🚀 Getting Started & Example (Vignette)](#-getting-started--example-vignette)
- [🎓 Citation](#-citation)
- [🤝 Contributing](#-contributing)
- [📄 License](#-license)
- [📬 Contact](#-contact)

---

## 🎯 Project Motivation

The Beta Autoregressive Moving Average ($\beta$ARMA) model is a powerful tool for analyzing time series data bounded between 0 and 1. However, standard parameter estimation using conditional maximum likelihood can suffer from **numerical instability**, which often occurs when the log-likelihood function has flat regions, leading to convergence failures or unreliable estimates.

This R package implements a ridge penalization scheme that adds a penalty term to the log-likelihood function, enhancing its curvature and promoting stable convergence.

The effectiveness of these methods is demonstrated through a detailed case study on modeling the relative humidity in Brasília, Brazil. 

---

## ✨ Key Features

This package provides a robust toolkit for stable BARMA model estimation.

* **Ridge-Penalized BARMA Model:** The core `barma()` function is enhanced with a `penalty` argument to apply the ridge penalization scheme.
* **Core Estimation Engine:** The mathematical foundation is implemented in a series of functions for computing the log-likelihood (`loglik_*`), score vector (`score_vector_*`), and information matrix (`inf_matrix_*`), with variants for both standard and ridge-penalized estimation.
* **Vignette as a Case Study:** A detailed vignette (`relative_humidity_brasilia.Rmd`) serves as a practical guide and portfolio piece, demonstrating how to diagnose and solve numerical instability in a real-world application.

---

## 📂 Repository Structure

The repository is structured as a standard R package for clarity and reproducibility.

```plaintext
.
├── R/                  # Source code for all R functions.
├── data/               # Processed data included in the package (.rda).
├── data-raw/           # Raw data and scripts used to process it.
├── man/                # R package documentation files for functions.
├── reports/            # Pre-rendered, static vignettes for archival.
├── vignettes/          # Detailed tutorial and case study (.Rmd).
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

---

## 🛠️ Installation
This research compendium can be installed as an R package directly from GitHub. This is the recommended method as it handles all dependencies automatically.

First, ensure you have the `remotes` package. If not, install it from CRAN:

```R
if (!require("remotes")) {
  install.packages("remotes")
}
```

Then, install the package from GitHub:

```R
remotes::install_github("everton-da-costa/BarmaRidgeBJPS2025", 
                        dependencies = TRUE,
                        build_vignettes = TRUE)
```

**Prerequisites**

This package requires the following external R packages. You can run the command below to ensure all dependencies are installed on your system before proceeding.

```R
install.packages(c("doMC", "doRNG", "foreach", "lbfgs", "Rdpack", "dplyr", "ggplot2", "gridExtra", "zoo"))
```

**Last Tested Environment**
The scripts were last successfully tested on:

* **R version:** 4.4.2
* **Platform:** x86_64-pc-linux-gnu (64-bit)

---

## 🚀 Getting Started & Example (Vignette)

The best way to understand and replicate the analysis is through the package vignette, which provides a detailed, narrated code example.

**1. List Available Vignettes**
After installation, you can see all available vignettes with the following command:

```R
# Lists all tutorials for this package
browseVignettes("BarmaRidgeBJPS2025")
```

**2. Open the Vignette**
The main vignette showcases the practical application of the package's methods.

* `simulated_ts_example`: (Numerical Demonstration) Reproduces the simulation study from the paper, showing how standard CMLE fails and PCMLE succeeds in a controlled environment where the true parameters are known.

* `relative_humidity_brasilia`: (Empirical Application) An end-to-end project demonstrating how to solve numerical instability when modeling the relative humidity in Brasília. It covers the comparison between CMLE, PCMLE, and bootstrap-based estimates.

You can open the vignette directly from your R console to view the full analysis and code.

```R
# Open the numerical example with simulated data
vignette("simulated_ts_example", package = "BarmaRidgeBJPS2025")
```

```R
# Open the empirical application with Brasília humidity data
vignette("relative_humidity_brasilia", package = "BarmaRidgeBJPS2025")
```

📄 Accessing Archival Reports

To ensure long-term reproducibility and provide a static record of the output, we have made pre-rendered HTML versions of the vignettes available in the reports/ directory of this repository.

These present the output generated by the code on a specific date. You can download these files directly from GitHub and open them in any web browser to see the final report without needing to install the package or run the R code. For complete context, each report concludes with a "Reproducibility" section detailing the exact versions of R, the operating system, and all packages used to generate the results.

---

## 🎓 Citation

If you use this code or data in your research, please cite the original article:

```bibtex
@article{CribariNeto+Costa+Fonseca_2025,
  title     = {Numerical stability enhancements in beta autoregressive moving average model estimation},
  author    = {Cribari-Neto, F. and Costa, E. and Fonseca, R. V.},
  journal   = {Brazilian Journal of Probability and Statistics},
  year      = {2025}
}
```

## 🤝 Contributing
Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## 📄 License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## 📬 Contact
For questions, suggestions, or issues related to the code, please contact:

Everton da Costa  
📧 everto.cost@gmail.com 
