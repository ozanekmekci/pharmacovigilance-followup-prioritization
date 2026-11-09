# Pharmacovigilance Follow-up Prioritization

Machine-learning project focused on prioritizing pharmacovigilance (PV) follow-up cases when specialist review capacity is limited.

## Project overview

The project uses FDA FAERS-derived reference data as an empirical calibration source for a reproducible synthetic dataset. The final modeling dataset contains **15,256 synthetic cases**, **22 candidate predictors**, two synthetic identifiers, and one binary target indicating whether serious outcome information is recorded during later observed follow-up.

The analysis is framed as a **case-prioritization / ranking problem**, not as a causal clinical-risk model.

## Workflow

1. **Synthetic-data generation**  
   `R/01_generate_synthetic_data.R` uses the reference dataset to generate the final synthetic dataset with a fixed random seed (`set.seed(42)`).

2. **Synthetic-data validation**  
   `R/02_validate_synthetic_data.R` checks structural integrity, distributions, dependencies, target behavior, and predictive plausibility. The final synthetic dataset passed all 38 structural checks.

3. **Modeling and prioritization**  
   `R/03_modeling.R` fits an ordinary logistic-regression model and ranks cases by predicted probability for capacity-constrained follow-up.

4. **Capacity-sensitivity analysis**  
   Model performance is evaluated at approximately 10%, 20%, and 30% review capacity to examine the trade-off between case coverage and review efficiency.

## Key results

- ROC-AUC: **0.654**
- PR-AUC: **0.091**
- At ~20% review capacity, the model captured **34.2%** of positive cases.
- Lift at ~20% review capacity: **1.71x** relative to random case selection.

These results describe performance in a synthetic-data simulation and should not be interpreted as demonstrated real-world clinical or operational benefit.

## Repository structure

```text
.
├── Data/
│   ├── Reference/
│   │   └── faers_reference_data_extended_v5.csv
│   └── Synthetic/
│       └── synthetic_pv_data_v2.csv
├── R/
│   ├── 01_generate_synthetic_data.R
│   ├── 02_validate_synthetic_data.R
│   └── 03_modeling.R
├── outputs/
│   └── modeling/
├── Report/
│   └── project.pdf
└── Machine Learning Project.Rproj
```

## Reproducing the analysis

Open the project from the repository root and run the scripts in this order:

```r
source("R/01_generate_synthetic_data.R")
source("R/02_validate_synthetic_data.R")
source("R/03_modeling.R")
```

The scripts use project-relative paths. Required R packages for the generation/modeling workflow include `tidyverse`, `here`, and `tidymodels`.

## Data note

The reference dataset is used as an empirical calibration and validation source. Complete real FAERS records are not copied directly into the synthetic dataset. The generated dataset is intended for an academic machine-learning exercise and does not provide a privacy guarantee or a clinically validated risk score.

## Report

The final project report is available at `Report/project.pdf`.
