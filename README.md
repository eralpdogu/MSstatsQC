[![Bioconductor-devel Build Status](http://bioconductor.org/shields/build/devel/bioc/MSstatsQC.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/MSstatsQC/)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)

# MSstatsQC: Longitudinal System Suitability Monitoring & Quality Control

**MSstatsQC** is an open-source R package and Shiny application designed to streamline quality control (QC) and system suitability testing (SST) for targeted proteomic experiments.

By bridging the gap between traditional statistical process control (SPC) and modern **Machine Learning (ML)**, MSstatsQC offers a comprehensive framework for real-time monitoring of mass spectrometric assays.

### 🚀 Key Features

* **Advanced Statistical Process Control:** Implements simultaneous and time-weighted control charts alongside change point analysis to monitor mean and variability.
* **Next-Gen Machine Learning Integration:**
    * Built on the powerful **H2O.ai** engine.
    * Utilizes a novel ML process that uses SPC metrics as input features.
    * Includes a built-in **synthetic data generation** function to train robust models even with limited historical data.
* **Comprehensive Metric Tracking:** Automatically tracks critical metrics including Total Peak Area, Retention Time, FWHM, and Peak Asymmetry.
* **Early Detection:** Significantly improves the ability to detect chromatographic and instrumental failures early, reducing the cost of control and preventing data loss.

### ⚠️ System Requirements & Installation

This package leverages **h2o** for its machine learning capabilities.

1.  **Java Requirement:** You must have **Java (JDK 8 or higher)** installed on your system to run the ML modules.
2.  **Installation:**
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MSstatsQC")
```
# Developers version (GitHub) installation:
```r
BiocManager::install("eralpdogu/MSstatsQC")
```
## Quick Start
Here is a simple example to create a decision map:

```r
library(MSstatsQC)
data <- MSstatsQC::S9Site54
data <- DataProcess(data)
DecisionMap(S9Site54,
    method = "XmR", peptideThresholdRed = 0.25, peptideThresholdYellow = 0.10,
    L = 1, U = 20, type = "mean", title = "Decision map", listMean = NULL, listSD = NULL
)
```
