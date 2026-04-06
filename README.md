# Rethinking global hotspots of high-priority viral zoonoses to guide spillover risk monitoring

[![Preprint](https://img.shields.io/badge/Publication-Paper-D95C5C)](https://www.researchsquare.com/article/rs-7076613/v1)
[![Zenodo](https://img.shields.io/badge/Zenodo-Data-11557c)](https://doi.org/10.5281/zenodo.16778284)

## Overview

This repository contains the R scripts used for predicting and mapping wild mammal hosts for WHO high-priority viral groups. The pipeline includes assembling the datasets, running the analyses, and mapping the hotspots and trends of exposure pressures.

### Target Viral Groups

The scripts iterate through the target viral groups. The abbreviations used across the file names and scripts are:

| Code | Viral Group |
| :--- | :--- |
| **`bcov`** | Betacoronavirus |
| **`hpv`** | Henipavirus |
| **`ebv`** | Ebolavirus & Marburgvirus |
| **`phv`** | Phlebovirus |
| **`ohv`** | Orthonairovirus |
| **`mmv`** | Mammarenavirus |
| **`flv`** | Flavivirus |

-----

## Repository Structure

### 1\. Spatial Overlaps & Pseudo-Negatives

  * **`1a_blueprint_pn_overlaps_ebv.R`** to **`1f_blueprint_pn_overlaps_phv.R`**
    Computes the spatial overlaps between known hosts and non-hosts based on species ranges. This step is critical for identifying and sampling valid "pseudo-negative" species for the machine learning pipeline.

### 2\. Dataset Assembly

  * **`2a_dataset_assem_ebv.R`** to **`2f_dataset_assem_phv.R`**
    Assembles the final feature datasets for each viral group. This includes computing and assigning specific instance weights for the modeling phase, strictly accounting for high-evidence hosts, low-evidence hosts, and the generated pseudo-negatives.

### 3\. Modelling Pipeline

  * **`3a_modelling_ebv.R`** to **`3f_modelling_phv.R`**
    These scripts run the hyperparameter tuning, model training, and validation of the ensemble models using a nested cross-validation. They ultimately output the predicted probabilities for both in-sample and out-of-sample mammal species.

### 4\. Hotspot Mapping

  * **`4a_observedhotspots.R`**
    Maps the global richness of observed hosts.
  * **`4b_predictedhotspots.R`**
    Maps the global richness of observed and predicted hosts based on the model outputs.

### 5\. Figures

  * **`5_PredictedHostsFigures.R`**
    Reproduces the main text and supplementary figures *(Spatial maps are generated in their respective scripts).*

### 6\. Cumulative Hotspots

  * **`6_cumulative_hotspots.R`**
    Aggregates the viral hazard maps and identifies the cumulative hotspots.

### 7\. Human Exposure Trends

  * **`7_population_deforestation_trends.R`**
    Extracts and analyzes population growth and deforestation rates within hotspot areas. Reproduces **Figure 4** .

-----

## Software Requirements

R version and core packages:

  * **R:** version `[4.5.2]`
  * **mlr3:** version `[1.3.0]` *(and the mlr3 ecosystem: mlr3tuning, mlr3learners)*
  * **tidyverse:** version `[2.0.0]`
  * **terra:** version `[1.8-86]`
