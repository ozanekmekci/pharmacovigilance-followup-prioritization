
# Final modeling - Synthetic Dataset v2
# Machine Learning for Pharmacovigilance Follow-up Prioritization


# 1. Packages

library(tidyverse)
library(tidymodels)
library(here)


# 2. Load the frozen Synthetic Dataset v2

synthetic_path <- here("Data", "Synthetic", "synthetic_pv_data_v2.csv")
output_dir <- here("outputs", "modeling")

synthetic_data <- read_csv(
  synthetic_path,
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
  nrow(synthetic_data) == 15256,
  ncol(synthetic_data) == 25,
  all(required_columns %in% names(synthetic_data)),
  all(synthetic_data$new_serious_info_followup %in% c(0, 1))
)

synthetic_data <- synthetic_data |>
  select(all_of(required_columns)) |>
  mutate(
    new_serious_info_followup = factor(
      new_serious_info_followup,
      levels = c(0, 1),
      labels = c("0", "1")
    ),
    reporter_type = factor(reporter_type),
    patient_sex = factor(patient_sex)
  )


# 3. Train/test split

set.seed(123)

model_split <- initial_split(
  synthetic_data,
  prop = 0.80,
  strata = new_serious_info_followup
)

train_data <- training(model_split)
test_data <- testing(model_split)

split_summary <- tibble(
  partition = c("Full", "Training", "Test"),
  rows = c(nrow(synthetic_data), nrow(train_data), nrow(test_data)),
  positives = c(
    sum(synthetic_data$new_serious_info_followup == "1"),
    sum(train_data$new_serious_info_followup == "1"),
    sum(test_data$new_serious_info_followup == "1")
  ),
  positive_rate = c(
    mean(synthetic_data$new_serious_info_followup == "1"),
    mean(train_data$new_serious_info_followup == "1"),
    mean(test_data$new_serious_info_followup == "1")
  )
)

print(split_summary)


# 4. Preprocessing 

# Structural indication NAs are represented as 0 for modeling only.
# indication_known remains in the model so unknown indication is still identifiable.

indication_columns <- c(
  "indication_oncology",
  "indication_respiratory",
  "indication_autoimmune",
  "indication_neurological"
)

train_data <- train_data |>
  mutate(across(all_of(indication_columns), ~ replace_na(.x, 0)))

test_data <- test_data |>
  mutate(across(all_of(indication_columns), ~ replace_na(.x, 0)))

model_recipe <- recipe(
  new_serious_info_followup ~ .,
  data = train_data
) |>
  update_role(caseid, primaryid, new_role = "id") |>
  step_impute_median(age_years) |>
  step_dummy(all_nominal_predictors()) |>
  step_zv(all_predictors())


# 5. Logistic regression 

logistic_spec <- logistic_reg() |>
  set_engine("glm") |>
  set_mode("classification")

logistic_workflow <- workflow() |>
  add_recipe(model_recipe) |>
  add_model(logistic_spec)

logistic_fit <- fit(
  logistic_workflow,
  data = train_data
)


# 6. Test-set predictions

prob_predictions <- predict(
  logistic_fit,
  new_data = test_data,
  type = "prob"
)

test_predictions <- test_data |>
  select(caseid, primaryid, new_serious_info_followup) |>
  bind_cols(prob_predictions)

ranked_predictions <- test_predictions |>
  arrange(desc(.pred_1), caseid, primaryid) |>
  mutate(rank = row_number())


# 7. Overall evaluation 

roc_auc_value <- roc_auc(
  test_predictions,
  truth = new_serious_info_followup,
  .pred_1,
  event_level = "second"
) |>
  pull(.estimate)

pr_auc_value <- pr_auc(
  test_predictions,
  truth = new_serious_info_followup,
  .pred_1,
  event_level = "second"
) |>
  pull(.estimate)

test_prevalence <- mean(
  test_predictions$new_serious_info_followup == "1"
)

total_true_positives <- sum(
  test_predictions$new_serious_info_followup == "1"
)

metrics_table <- tibble(
  roc_auc = roc_auc_value,
  pr_auc = pr_auc_value,
  test_positive_prevalence = test_prevalence
)

print(metrics_table)


# 8. Capacity-sensitivity analysis

# Meaningful extension:
# evaluate the same final test ranking at fixed 10%, 20%, and 30% review capacity.

capacity_table <- tibble(
  capacity = c(0.10, 0.20, 0.30)
) |>
  mutate(
    total_test_cases = nrow(ranked_predictions),
    reviewed_cases = ceiling(capacity * total_test_cases),
    actual_capacity = reviewed_cases / total_test_cases,
    total_true_positives = total_true_positives
  ) |>
  rowwise() |>
  mutate(
    positives_captured = sum(
      ranked_predictions$rank <= reviewed_cases &
        ranked_predictions$new_serious_info_followup == "1"
    ),
    recall_at_k = positives_captured / total_true_positives,
    precision_at_k = positives_captured / reviewed_cases,
    lift_at_k = precision_at_k / test_prevalence,
    random_expected_recall = actual_capacity,
    random_expected_precision = test_prevalence,
    random_lift = 1
  ) |>
  ungroup()

print(capacity_table)


# 9. Coefficient table 

# Coefficients describe associations within the synthetic dataset only.

logistic_engine <- extract_fit_engine(logistic_fit)

coefficient_table <- broom::tidy(logistic_engine) |>
  transmute(
    term,
    coefficient = estimate,
    odds_ratio = exp(estimate)
  )

print(coefficient_table)


# 10. Evaluation curves 

roc_curve_data <- roc_curve(
  test_predictions,
  truth = new_serious_info_followup,
  .pred_1,
  event_level = "second"
)

roc_plot <- ggplot(
  roc_curve_data,
  aes(x = 1 - specificity, y = sensitivity)
) +
  geom_path() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "Logistic regression ROC curve",
    subtitle = paste("ROC-AUC =", round(roc_auc_value, 3)),
    x = "False-positive rate",
    y = "True-positive rate"
  )

pr_curve_data <- pr_curve(
  test_predictions,
  truth = new_serious_info_followup,
  .pred_1,
  event_level = "second"
)

pr_plot <- ggplot(
  pr_curve_data,
  aes(x = recall, y = precision)
) +
  geom_path() +
  geom_hline(
    yintercept = test_prevalence,
    linetype = "dashed"
  ) +
  labs(
    title = "Logistic regression precision-recall curve",
    subtitle = paste("PR-AUC =", round(pr_auc_value, 3)),
    x = "Recall",
    y = "Precision"
  )


# Cumulative capture curve for the capacity-sensitivity extension.
capture_curve_data <- ranked_predictions |>
  mutate(
    review_capacity = rank / nrow(ranked_predictions),
    positives_captured_cumulative = cumsum(
      new_serious_info_followup == "1"
    ),
    recall_cumulative = positives_captured_cumulative / total_true_positives
  )

capacity_plot <- ggplot(
  capture_curve_data,
  aes(x = review_capacity, y = recall_cumulative)
) +
  geom_line() +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  geom_point(
    data = capacity_table,
    aes(x = actual_capacity, y = recall_at_k),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    title = "Positive capture by review capacity",
    subtitle = "Dashed line represents expected capture under random ranking",
    x = "Share of test cases reviewed",
    y = "Share of positives captured"
  )


# 11. Save outputs 

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  metrics_table,
  file.path(output_dir, "logistic_metrics.csv")
)

write_csv(
  capacity_table,
  file.path(output_dir, "capacity_sensitivity.csv")
)

write_csv(
  coefficient_table,
  file.path(output_dir, "logistic_coefficients.csv")
)

write_csv(
  ranked_predictions,
  file.path(output_dir, "test_predictions.csv")
)

ggsave(
  file.path(output_dir, "roc_curve.png"),
  roc_plot,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_dir, "pr_curve.png"),
  pr_plot,
  width = 7,
  height = 5
)

ggsave(
  file.path(output_dir, "capacity_sensitivity.png"),
  capacity_plot,
  width = 7,
  height = 5
)

cat("\nFinal Synthetic v2 modeling complete.\n")
cat("Outputs saved to:", output_dir, "\n")
