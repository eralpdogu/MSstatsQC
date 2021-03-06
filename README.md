# Longitudinal system suitability monitoring and quality control for targeted proteomic experiments
Statistical process control (SPC) is a general and well-established method of quality control (QC) which can be used to monitor and improve the quality of a process such as LC MS/MS. `MSstatsQC` is an open-source R package and Shiny web application for statistical analysis and monitoring of quality control and system suitability testing (SST) samples produced by spectrometry-based proteomic experiments. Our framework termed `MSstatsQC` is available through www.msstats.org/msstatsqc. It uses SPC tools to track metrics including total peak area, retention time, full width at half maximum (FWHM) and peak asymmetry for proteomic experiments. We introduce simultaneous and time weighted control charts and change point analysis to monitor mean and variability of metrics. Proposed longitudinal monitoring approach significantly improves the ability of real time monitoring, early detection and prevention of chromatographic and instrumental problems of mass spectrometric assays, thereby, reducing cost of control and failure.

Installation

To install this package, start R and enter:
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")