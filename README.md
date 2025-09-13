# Wu et al. 2025 Proteomics Data Analysis

This repository contains the R scripts used for the initial processing and statistical analysis of proteomics data, supporting the findings in the upcoming publication: *NUDT5 regulates purine metabolism and thiopurine sensitivity by interacting with PPAT*.

---

## Overview

The scripts in this repository were used to:

* Clean and normalize raw mass spectrometry data for total protein, phosphoprotein, and phosphopeptide levels.
* Perform statistical tests, including t-tests and two-way ANOVA, to identify proteins and peptides with significant changes across different experimental conditions.
* Categorize proteins into groups based on their response to treatment and genotype.
* Generate the processed data tables and preliminary plots that served as the basis for further visualization and analysis in the paper.

---

## Repository Contents

### `data/`

This folder contains the raw input data files for the analysis.

* `pathway_list.rds`: An R data file containing a list of pathways for enrichment analysis.
* `PCF-ZW-7426--05--2024_total_PD1.xlsx`: Raw data for total protein.
* `PCF-ZW-7426--05--2024_phosphoprotein_PD1.xlsx`: Raw data for phosphoprotein.
* `PCF-ZW-7426--05--2024_phosphopeptide_PD1.xlsx`: Raw data for phosphopeptide.

### `script/`

This folder contains the R scripts used to perform the analysis.

* `data_cleaning.R`: Handles data import, normalization, statistical comparisons, and grouping of proteins and peptides.
* `plotting.R`: Generates visualizations, including scatter plots and heatmaps.
* `HyperGeometricTest.R`: Script for performing hypergeometric tests, for pathway enrichment analysis.
* `beautify_df.R`: A utility script for formatting data frames for better readability or export.

---

## Usage

The primary script to run is `data_cleaning.R`. It depends on the data files located in the `data/` folder and saves intermediate files in the `work/` directory. The `plotting.R` script then reads these intermediate files to generate the visualizations.
