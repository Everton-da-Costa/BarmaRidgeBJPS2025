# BarmaRidgeBJPS2025

[![Status](https://img.shields.io/badge/Status-Submitted-lightgrey.svg)](https://projecteuclid.org/journals/brazilian-journal-of-probability-and-statistics)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/Everton-da-Costa/BarmaRidgeBJPS2025/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Everton-da-Costa/BarmaRidgeBJPS2025/actions/workflows/R-CMD-check.yaml)

This repository contains the R package and associated data for the scientific article:

**"Numerical stability enhancements in beta autoregressive moving average model estimation"** by Cribari-Neto, F., Costa, E., and Fonseca, R. V.
Submitted to the *Brazilian Journal of Probability and Statistics*. 

---

## 📚 Table of Contents

- [🎯 Project Motivation](#-project-motivation)
- [🏆 Publication & Journal Quality](#-publication--journal-quality)
- [🛠️ Key Skills Demonstrated](#️-key-skills-demonstrated)
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

The Beta Autoregressive Moving Average (BARMA) model is a powerful tool for analyzing time series data bounded between 0 and 1. However, standard parameter estimation using conditional maximum likelihood can suffer from **numerical instability**.  This often occurs when the log-likelihood function has flat regions, leading to convergence failures or unreliable, implausible estimates. 

This project introduces and implements a ridge penalization scheme that adds a penalty term to the log-likelihood function, enhancing its curvature and promoting stable convergence.

The effectiveness of these methods is demonstrated through a detailed case study on modeling the relative humidity in Brasília, Brazil. 

---

## 🏆 Publication & Journal Quality

This research has been submitted to the **Brazilian Journal of Probability and Statistics**, the official publication of the **Brazilian Statistical Association (ABE)** and supported by the **Institute of Mathematical Statistics (IMS)**.

The journal is a respected venue for methodological advances in statistics. Its key bibliometric indicators (2023 data) include:

* **SJR (SCImago Journal Rank):** 0.251
* **H-Index:** 23
* **CiteScore:** 1.2
* **Impact Factor:** 0.55

<a href="https://www.scimagojr.com/journalsearch.php?q=19900192736&amp;tip=sid&amp;exact=no" title="SCImago Journal &amp; Country Rank"><img border="0" src="https://www.scimagojr.com/journal_img.php?id=19900192736" alt="SCImago Journal &amp; Country Rank" /></a>

*(Note: Metrics reflect the journal's standing in the field of Statistics and Probability)*

---

## 🛠️ Key Skills Demonstrated

This project showcases a range of advanced data science and statistical engineering skills:

* **Statistical Algorithm Design:** Implementing a **Ridge Penalization** scheme from scratch to solve non-convergence in Maximum Likelihood Estimation.
* **R Package Development:** Creating a structured, documented, and installable R package.
* **High-Performance Computing:** Utilizing the `doMC` package for parallel processing in simulation studies (Linux environments).
* **Reproducible Research:** Authoring detailed R Markdown vignettes that serve as reproducible case studies.
* **Numerical Optimization:** Working with advanced optimization routines (`lbfgs`) to handle complex likelihood surfaces.

---

## ✨ Key Features

This package provides an toolkit for stable BARMA model estimation.

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
├── vignettes/          # Detailed tutorial and case study (.Rmd).
├── reports/            # Pre-rendered HTML reports.
├── man/                # R package documentation files for functions.
├── CODE_OF_CONDUCT.md  # Contributor Code of Conduct.
├── DESCRIPTION         # Package metadata and dependencies.
├── NAMESPACE           # Manages the package's namespace.
├── LICENSE             # MIT License file.
└── README.md           # This file.
```

---

## 🛠️ Local Installation Guide

### 1. Install Dependencies
Before installing this package, please ensure the required dependencies are installed from CRAN:

```R
install.packages(c("doMC", "doRNG", "foreach", "lbfgs", "Rdpack", "dplyr", "ggplot2", "gridExtra", "zoo"))
```

**⚠️ Compatibility Note for Windows Users:**
This package uses doMC for parallel computing. **doMC is not available on Windows.** Windows users may install the package but might encounter issues running the specific parallelized vignettes unless they use a Linux subsystem (WSL) or modify the code to use doParallel.

### 2. Install the Package
Save the provided `.tar.gz` file to your computer. Run the following command in R, replacing `path/to/file/` with the actual location of the file:

```R
install.packages("path/to/file/BarmaRidgeBJPS2025_1.0.0.tar.gz", 
                 repos = NULL, 
                 type = "source")
```
*(Note: Change the version number in the filename if it differs from `1.0.0`)*

**⚠️ Note on Installation Time:**
The installation process builds the package vignettes, which involve running **bootstrap procedures** to reproduce the examples from the article. 
**This process may take several minutes.**
Please note that while the original simulation study utilized high-performance parallel computing (8+ cores), the package vignettes are restricted to run on **1 core** to ensure compatibility across all systems, even though the code supports parallelization.


---

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

* `relative_humidity_brasilia`: (Portfolio Case Study) An end-to-end data science project demonstrating how to diagnose and solve numerical instability when modeling relative humidity. It covers the comparison between standard MLE, penalized MLE (PMLE), and bootstrap-based estimates.

You can open the vignette directly from your R console to view the full analysis and code.

```R
# Open the main application case study
vignette("relative_humidity_brasilia", package = "BarmaRidgeBJPS2025")
```

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
