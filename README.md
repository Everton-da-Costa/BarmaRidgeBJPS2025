# BarmaRidgeBJPS2025

[![DOI](https://img.shields.io/badge/DOI-10.1214/25--BJPS645-blue.svg)](https://doi.org/10.1214/25-BJPS645)
[![GitHub Downloads](https://img.shields.io/github/downloads/Everton-da-Costa/BarmaRidgeBJPS2025/total.svg?logo=github&color=blue)](https://github.com/Everton-da-Costa/BarmaRidgeBJPS2025/releases)
[![R-CMD-check](https://github.com/Everton-da-Costa/BarmaRidgeBJPS2025/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Everton-da-Costa/BarmaRidgeBJPS2025/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the R package and associated data for the scientific article:

**"Numerical stability enhancements in beta autoregressive moving average model estimation"** by Cribari-Neto, F., Costa, E., and Fonseca, R. V. Published in the ***Brazilian Journal of Probability and Statistics*** (2025), Vol. 39, No. 4, pp. 410-437.

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

The Beta Autoregressive Moving Average (BARMA) model is a powerful tool for analyzing time series data bounded between 0 and 1. However, standard estimation via conditional maximum likelihood can suffer from **numerical instability**, such as convergence failures or implausible estimates. 

This project implements a **ridge penalization** scheme to enhance the log-likelihood curvature, ensuring stable convergence. For cases where penalization alone is insufficient, the package also provides a **bootstrap-based estimation strategy**.

---

## 🏆 Publication & Journal Quality

This research was published in the **Brazilian Journal of Probability and Statistics (BJPS)**, the official journal of the Brazilian Statistical Association (ABE). The journal's standing speaks to the rigor of the underlying research, as reflected in its key metrics:

* **Impact Factor:** 0.55
* **CiteScore:** 1.2
* **SJR:** 0.251
* **H-Index:** 23

[![SCImago Journal & Country Rank](https://www.scimagojr.com/journal_img.php?id=19900192736)](https://www.scimagojr.com/journalsearch.php?q=19900192736&tip=sid&exact=no)

*(Note: Metrics reflect the journal's standing in the field of Statistics and Probability)*

---

## 🛠️ Key Skills Demonstrated

This project showcases a range of advanced data science and statistical engineering skills:

* **Statistical Algorithm Design:** Implementing a **Ridge Penalization** scheme from scratch to solve non-convergence in Conditional Maximum Likelihood Estimation.
* **R Package Development:** Creating a structured, documented, and installable R package.
* **High-Performance Computing:** Utilizing the `doMC` package for parallel processing in simulation studies (Linux environments).
* **Reproducible Research:** Authoring detailed R Markdown vignettes that serve as reproducible case studies.
* **Numerical Optimization:** Working with advanced optimization routines (`lbfgs`) to handle complex likelihood surfaces.

---

## ✨ Key Features

* **Ridge-Penalized BARMA:** Core estimation with analytical gradients for enhanced stability.
* **Bootstrap Estimation:** Non-parametric block bootstrap for reliable estimates in extreme scenarios.
* **Reproducible Research:** Comprehensive vignettes replicating the relative humidity (Brasília) and reservoir volume (Itaparica) case studies.

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
@article{Cribari_Costa_Fonseca_2025,
  author  = {Cribari-Neto, F. and Costa, E. and Fonseca, R. V.},
  title   = {Numerical stability enhancements in beta autoregressive moving average model estimation},
  journal = {Brazilian Journal of Probability and Statistics},
  year    = {2025},
  volume  = {39},
  number  = {4},
  pages   = {410--437},
  doi     = {10.1214/25-BJPS645}
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
