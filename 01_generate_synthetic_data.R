# =============================================================================
# Synthetic Dataset v2 Generator
# Machine Learning for Pharmacovigilance Follow-up Prioritization
# =============================================================================
#
# Uses Reference Data v5 as the empirical calibration source.
# Complete real FAERS rows are not copied.
#
# Expected folders:
# Machine Learning Project/
#   R/01_generate_synthetic_data.R
#   Data/Reference/faers_reference_data_extended_v5.csv
#   Data/Synthetic/
# =============================================================================


# 1. Setup 

library(tidyverse)
library(here)

set.seed(42)
n_synthetic <- 15256

reference_path <- here(
  "Data",
  "Reference",
  "faers_reference_data_extended_v5.csv"
)

# Keep the historical v1 output untouched while v2 is being validated.
output_path <- here(
  "Data",
  "Synthetic",
  "synthetic_pv_data_v2.csv"
)

if (!file.exists(reference_path)) {
  stop(
    "faers_reference_data_extended_v5.csv was not found in Data/Reference/"
  )
}


# 2. Read Reference Data v5 

reference_data <- read_csv(
  reference_path,
  show_col_types = FALSE
)

required_columns <- c(
  "caseid",
  "primaryid",
  "age_years",
  "age_known",
  "event_date_known",
  "reporter_type",
  "patient_sex",
  "num_drugs",
  "has_concomitant_drug",
  "suspect_dose_known",
  "num_reactions",
  "suspect_indication_known",
  "suspect_start_date_known",
  "reaction_medication_error",
  "reaction_systemic_symptoms",
  "reaction_respiratory",
  "reaction_dermatological",
  "reaction_gastrointestinal",
  "reaction_treatment_ineffectiveness",
  "indication_known",
  "indication_oncology",
  "indication_respiratory",
  "indication_autoimmune",
  "indication_neurological",
  "new_serious_info_followup"
)

stopifnot(
  nrow(reference_data) == n_synthetic,
  all(required_columns %in% names(reference_data))
)


# 3. Create new synthetic identifiers

# IDs are new and are not used to generate predictors or the target.

synthetic_caseid <- 91000000 + seq_len(n_synthetic)
synthetic_primaryid <- synthetic_caseid * 10 + 2


# 4. Report-context / information-availability block 

# Sample these variables together to preserve their empirical combinations.

context_frequency <- reference_data |>
  count(
    reporter_type,
    patient_sex,
    age_known,
    event_date_known,
    suspect_dose_known,
    suspect_indication_known,
    suspect_start_date_known,
    indication_known,
    name = "frequency"
  )

synthetic_context <- context_frequency |>
  slice_sample(
    n = n_synthetic,
    replace = TRUE,
    weight_by = frequency
  ) |>
  select(-frequency)


# 5. Age

# Unknown age stays NA.
# Known age is sampled from the empirical distribution of known ages.

age_frequency <- reference_data |>
  filter(age_known == 1) |>
  count(age_years, name = "frequency")

synthetic_age <- rep(NA_real_, n_synthetic)

n_known_age <- sum(synthetic_context$age_known == 1)

synthetic_age[synthetic_context$age_known == 1] <- sample(
  age_frequency$age_years,
  size = n_known_age,
  replace = TRUE,
  prob = age_frequency$frequency
)


# Temporary age bands used only for clinical-profile sampling.

synthetic_age_band <- case_when(
  is.na(synthetic_age) ~ "Unknown",
  synthetic_age < 45 ~ "<45",
  synthetic_age >= 45 & synthetic_age < 65 ~ "45-64",
  synthetic_age >= 65 ~ "65+"
)

reference_with_age_band <- reference_data |>
  mutate(
    age_band = case_when(
      is.na(age_years) ~ "Unknown",
      age_years < 45 ~ "<45",
      age_years >= 45 & age_years < 65 ~ "45-64",
      age_years >= 65 ~ "65+"
    )
  )


# 6. Medication-complexity block 

# num_drugs and has_concomitant_drug are sampled together.

medication_frequency <- reference_data |>
  count(
    num_drugs,
    has_concomitant_drug,
    name = "frequency"
  )

synthetic_medication <- medication_frequency |>
  slice_sample(
    n = n_synthetic,
    replace = TRUE,
    weight_by = frequency
  ) |>
  select(-frequency)


# 7. Clinical-profile block 

# Full age conditioning remains active:
#   indication_known x age_band
#
# Reporter dependence is preserved partially through a 50/50 blend:
#   50% reporter-aware sampling
#   50% broader age/indication sampling
#
# All clinical variables in a sampled profile stay together.

clinical_profile_columns <- c(
  "num_reactions",
  "reaction_medication_error",
  "reaction_systemic_symptoms",
  "reaction_respiratory",
  "reaction_dermatological",
  "reaction_gastrointestinal",
  "reaction_treatment_ineffectiveness",
  "indication_oncology",
  "indication_respiratory",
  "indication_autoimmune",
  "indication_neurological"
)


# Broad empirical pool:
# indication_known x age_band x clinical profile

clinical_frequency_broad <- reference_with_age_band |>
  count(
    indication_known,
    age_band,
    num_reactions,
    reaction_medication_error,
    reaction_systemic_symptoms,
    reaction_respiratory,
    reaction_dermatological,
    reaction_gastrointestinal,
    reaction_treatment_ineffectiveness,
    indication_oncology,
    indication_respiratory,
    indication_autoimmune,
    indication_neurological,
    name = "frequency"
  )


# Reporter-aware empirical pool:
# indication_known x age_band x reporter_type x clinical profile

clinical_frequency_reporter <- reference_with_age_band |>
  count(
    indication_known,
    age_band,
    reporter_type,
    num_reactions,
    reaction_medication_error,
    reaction_systemic_symptoms,
    reaction_respiratory,
    reaction_dermatological,
    reaction_gastrointestinal,
    reaction_treatment_ineffectiveness,
    indication_oncology,
    indication_respiratory,
    indication_autoimmune,
    indication_neurological,
    name = "frequency"
  )


# For each synthetic case:
# 1 = reporter-aware pool
# 0 = broader age/indication pool
#
# This helper is not exported.

reporter_aware_branch <- rbinom(
  n = n_synthetic,
  size = 1,
  prob = 0.50
)


# Store sampled clinical profiles here before returning them
# to the original synthetic row order.

synthetic_clinical_list <- list()
list_index <- 1


for (indication_status in c(0, 1)) {

  for (band in c(
    "Unknown",
    "<45",
    "45-64",
    "65+"
  )) {

    stratum_rows <- which(
      synthetic_context$indication_known == indication_status &
        synthetic_age_band == band
    )

    if (length(stratum_rows) == 0) {
      next
    }


    # The broader pool is also the fallback for reporter-specific strata.

    broad_frequency <- clinical_frequency_broad |>
      filter(
        indication_known == indication_status,
        age_band == band
      )

    stopifnot(
      nrow(broad_frequency) > 0
    )


    # Broad branch 

    broad_rows <- stratum_rows[
      reporter_aware_branch[stratum_rows] == 0
    ]

    if (length(broad_rows) > 0) {

      sampled_broad <- broad_frequency |>
        slice_sample(
          n = length(broad_rows),
          replace = TRUE,
          weight_by = frequency
        ) |>
        select(
          all_of(clinical_profile_columns)
        ) |>
        mutate(
          row_id = broad_rows,
          .before = 1
        )

      synthetic_clinical_list[[list_index]] <- sampled_broad
      list_index <- list_index + 1
    }


    # Reporter-aware branch 

    reporter_rows <- stratum_rows[
      reporter_aware_branch[stratum_rows] == 1
    ]

    if (length(reporter_rows) > 0) {

      reporter_values <- unique(
        synthetic_context$reporter_type[reporter_rows]
      )

      for (reporter in reporter_values) {

        rows_for_reporter <- reporter_rows[
          synthetic_context$reporter_type[reporter_rows] == reporter
        ]

        reporter_frequency <- clinical_frequency_reporter |>
          filter(
            indication_known == indication_status,
            age_band == band,
            reporter_type == reporter
          )


        # If the reporter-specific empirical stratum is empty,
        # fall back to the broader age/indication pool.

        if (nrow(reporter_frequency) == 0) {

          sampling_frequency <- broad_frequency

        } else {

          sampling_frequency <- reporter_frequency
        }


        sampled_reporter <- sampling_frequency |>
          slice_sample(
            n = length(rows_for_reporter),
            replace = TRUE,
            weight_by = frequency
          ) |>
          select(
            all_of(clinical_profile_columns)
          ) |>
          mutate(
            row_id = rows_for_reporter,
            .before = 1
          )

        synthetic_clinical_list[[list_index]] <- sampled_reporter
        list_index <- list_index + 1
      }
    }
  }
}


# Combine all sampled profiles and restore the original row order.

synthetic_clinical <- bind_rows(
  synthetic_clinical_list
) |>
  arrange(row_id) |>
  select(-row_id)


stopifnot(
  nrow(synthetic_clinical) == n_synthetic
)


# 8. Synthetic v2 target-generation score 

# This is a synthetic generation score, not a clinical risk score.
# The score is used only to generate the target and is not exported.

# Age component: maximum 2 points.

age_score <- case_when(
  is.na(synthetic_age) ~ 0,
  synthetic_age < 45 ~ 0,
  synthetic_age >= 45 & synthetic_age < 65 ~ 1,
  synthetic_age >= 65 ~ 2
)

# Medication component: maximum 2 points.

medication_score <- case_when(
  synthetic_medication$num_drugs == 1 ~ 0,
  synthetic_medication$num_drugs >= 2 &
    synthetic_medication$num_drugs <= 4 ~ 1,
  synthetic_medication$num_drugs >= 5 ~ 2
)

# Reporter component: maximum 1 point.

reporter_score <- ifelse(
  synthetic_context$reporter_type %in% c(
    "Consumer",
    "Physician"
  ),
  1,
  0
)

# Indication-family component: maximum 1 point.
#
# Unknown indication must stay NA in the dataset.
# coalesce(..., 0) is used only here so unknown indication gives zero points.

indication_score <- ifelse(
  coalesce(
    synthetic_clinical$indication_oncology,
    0
  ) == 1 |
    coalesce(
      synthetic_clinical$indication_respiratory,
      0
    ) == 1,
  1,
  0
)

# Reaction-family component: maximum 1 point.

reaction_score <- ifelse(
  synthetic_clinical$reaction_systemic_symptoms == 1 |
    synthetic_clinical$reaction_gastrointestinal == 1,
  1,
  0
)

generation_score <- (
  age_score +
    medication_score +
    reporter_score +
    indication_score +
    reaction_score
)

stopifnot(
  all(generation_score >= 0),
  all(generation_score <= 7)
)


# 9. Propensity groups and stochastic target

propensity_group <- case_when(
  generation_score == 0 ~ "Very Low",
  generation_score == 1 ~ "Low",
  generation_score == 2 ~ "Moderate",
  generation_score == 3 ~ "Elevated",
  generation_score >= 4 ~ "High"
)

target_probability <- case_when(
  propensity_group == "Very Low" ~ 0.015,
  propensity_group == "Low" ~ 0.04,
  propensity_group == "Moderate" ~ 0.065,
  propensity_group == "Elevated" ~ 0.11,
  propensity_group == "High" ~ 0.13
)

synthetic_target <- rbinom(
  n = n_synthetic,
  size = 1,
  prob = target_probability
)


# 10. Assemble the final Synthetic Dataset v2 

synthetic_data <- tibble(
  caseid = synthetic_caseid,
  primaryid = synthetic_primaryid,
  age_years = synthetic_age,
  age_known = synthetic_context$age_known,
  event_date_known = synthetic_context$event_date_known,
  reporter_type = synthetic_context$reporter_type,
  patient_sex = synthetic_context$patient_sex,
  num_drugs = synthetic_medication$num_drugs,
  has_concomitant_drug = synthetic_medication$has_concomitant_drug,
  suspect_dose_known = synthetic_context$suspect_dose_known,
  num_reactions = synthetic_clinical$num_reactions,
  suspect_indication_known = synthetic_context$suspect_indication_known,
  suspect_start_date_known = synthetic_context$suspect_start_date_known,
  reaction_medication_error =
    synthetic_clinical$reaction_medication_error,
  reaction_systemic_symptoms =
    synthetic_clinical$reaction_systemic_symptoms,
  reaction_respiratory =
    synthetic_clinical$reaction_respiratory,
  reaction_dermatological =
    synthetic_clinical$reaction_dermatological,
  reaction_gastrointestinal =
    synthetic_clinical$reaction_gastrointestinal,
  reaction_treatment_ineffectiveness =
    synthetic_clinical$reaction_treatment_ineffectiveness,
  indication_known = synthetic_context$indication_known,
  indication_oncology =
    synthetic_clinical$indication_oncology,
  indication_respiratory =
    synthetic_clinical$indication_respiratory,
  indication_autoimmune =
    synthetic_clinical$indication_autoimmune,
  indication_neurological =
    synthetic_clinical$indication_neurological,
  new_serious_info_followup = synthetic_target
)


# 11. Essential structural checks 

# Detailed realism checks belong in the separate validation workflow.

reaction_columns <- c(
  "reaction_medication_error",
  "reaction_systemic_symptoms",
  "reaction_respiratory",
  "reaction_dermatological",
  "reaction_gastrointestinal",
  "reaction_treatment_ineffectiveness"
)

indication_columns <- c(
  "indication_oncology",
  "indication_respiratory",
  "indication_autoimmune",
  "indication_neurological"
)

# Basic dataset checks.

stopifnot(
  nrow(synthetic_data) == n_synthetic,
  ncol(synthetic_data) == 25,
  identical(names(synthetic_data), required_columns),
  !anyDuplicated(synthetic_data$caseid),
  !anyDuplicated(synthetic_data$primaryid),
  !any(synthetic_data$caseid %in% reference_data$caseid),
  !any(synthetic_data$primaryid %in% reference_data$primaryid),
  all(synthetic_data$num_drugs >= 1),
  all(synthetic_data$num_reactions >= 1),
  all(synthetic_data$new_serious_info_followup %in% c(0, 1))
)

# Age rules.

stopifnot(
  all(
    is.na(
      synthetic_data$age_years[
        synthetic_data$age_known == 0
      ]
    )
  ),
  all(
    !is.na(
      synthetic_data$age_years[
        synthetic_data$age_known == 1
      ]
    )
  )
)

# Medication rule.

stopifnot(
  !any(
    synthetic_data$num_drugs == 1 &
      synthetic_data$has_concomitant_drug == 1
  )
)

# indication_known = 1 should only occur in valid context combinations.

stopifnot(
  !any(
    synthetic_data$indication_known == 1 &
      synthetic_data$suspect_indication_known == 0
  )
)

# Indication-category missingness rules.

unknown_indication_data <- synthetic_data |>
  filter(indication_known == 0) |>
  select(all_of(indication_columns))

known_indication_data <- synthetic_data |>
  filter(indication_known == 1) |>
  select(all_of(indication_columns))

stopifnot(
  all(is.na(as.matrix(unknown_indication_data))),
  all(!is.na(as.matrix(known_indication_data))),
  all(as.matrix(known_indication_data) %in% c(0, 1))
)

# Reaction-category count cannot exceed num_reactions.

reaction_category_count <- rowSums(
  synthetic_data |>
    select(all_of(reaction_columns))
)

stopifnot(
  all(
    reaction_category_count <= synthetic_data$num_reactions
  )
)


# 12. Save Synthetic Dataset v2 

dir.create(
  here("Data", "Synthetic"),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  synthetic_data,
  output_path,
  na = ""
)


# 13. Short generation summary 

cat("\nSynthetic Dataset v2 generation complete.\n")
cat("Rows:", nrow(synthetic_data), "\n")
cat("Columns:", ncol(synthetic_data), "\n")

cat(
  "Reference v5 target prevalence:",
  round(
    mean(reference_data$new_serious_info_followup) * 100,
    2
  ),
  "%\n"
)

cat(
  "Synthetic v2 target prevalence:",
  round(
    mean(synthetic_data$new_serious_info_followup) * 100,
    2
  ),
  "%\n"
)

cat("Saved to:", output_path, "\n")

# synthetic_data is the final dataset.
# synthetic_context, synthetic_medication and synthetic_clinical are
# intermediate generation objects.
