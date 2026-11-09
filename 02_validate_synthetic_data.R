# =============================================================================
# Synthetic Dataset v2 validation - FAERS reference v5 versus frozen v2 candidate
# Run from project root: Rscript R/02_validate_synthetic_data.R
# Base R only. Does not regenerate or modify the synthetic candidate CSV.
# =============================================================================

# 1. Project paths and inputs --------------------------------------------------
find_root <- function() {
  candidates <- c(getwd(), dirname(getwd()), dirname(dirname(getwd())))
  for (p in candidates) {
    if (file.exists(file.path(p, "Data", "Reference", "faers_reference_data_extended_v5.csv")) &&
        file.exists(file.path(p, "Data", "Synthetic", "synthetic_pv_data_v2.csv"))) {
      return(normalizePath(p))
    }
  }
  stop("Run this script from the project root (or R/). Both input CSVs are required.")
}
root <- find_root()
out_dir <- file.path(root, "outputs", "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
reference <- read.csv(file.path(root, "Data", "Reference", "faers_reference_data_extended_v5.csv"),
                      stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE)
synthetic <- read.csv(file.path(root, "Data", "Synthetic", "synthetic_pv_data_v2.csv"),
                      stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE)

id_cols <- c("caseid", "primaryid")
target <- "new_serious_info_followup"
context_cols <- c("reporter_type", "patient_sex", "age_known", "event_date_known",
                  "suspect_dose_known", "suspect_indication_known",
                  "suspect_start_date_known", "indication_known")
reaction_cols <- c("reaction_medication_error", "reaction_systemic_symptoms",
                   "reaction_respiratory", "reaction_dermatological",
                   "reaction_gastrointestinal", "reaction_treatment_ineffectiveness")
indication_cols <- c("indication_oncology", "indication_respiratory",
                     "indication_autoimmune", "indication_neurological")
continuous_cols <- c("age_years", "num_drugs", "num_reactions")
nominal_cols <- c("reporter_type", "patient_sex")
expected_cols <- c("caseid", "primaryid", "age_years", "age_known", "event_date_known",
                   "reporter_type", "patient_sex", "num_drugs", "has_concomitant_drug",
                   "suspect_dose_known", "num_reactions", "new_serious_info_followup",
                   "suspect_indication_known", "suspect_start_date_known", reaction_cols[1:3],
                   "indication_known", indication_cols[1:3], reaction_cols[4:6],
                   indication_cols[4])
all_predictors <- setdiff(expected_cols, c(id_cols, target))
binary_cols <- setdiff(all_predictors, c(continuous_cols, nominal_cols))

save_output <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, na = "")
}
rate <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
proportion <- function(v, value) if (!length(v)) NA_real_ else mean(v == value, na.rm = TRUE)
fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits))

# 2. Layer 1: hard structural checks ------------------------------------------
if (!all(expected_cols %in% names(reference)) || !all(expected_cols %in% names(synthetic))) {
  stop("Missing required input columns; cannot run complete validation.")
}
checks <- data.frame(check = character(), result = character(), violations = integer(),
                     status = character(), stringsAsFactors = FALSE)
add_check <- function(check, violations, result) {
  checks <<- rbind(checks,
                   data.frame(check = check, result = as.character(result),
                              violations = as.integer(violations),
                              status = if (violations == 0) "PASS" else "FAIL"))
}
add_check("exact_row_count", abs(nrow(synthetic) - 15256), nrow(synthetic))
add_check("exact_column_count", abs(ncol(synthetic) - 25), ncol(synthetic))
add_check("required_schema", length(setdiff(expected_cols, names(synthetic))) +
            length(setdiff(names(synthetic), expected_cols)),
          "25 required names; no extras")
add_check("unique_caseid", sum(duplicated(synthetic$caseid)) + sum(is.na(synthetic$caseid)),
          length(unique(synthetic$caseid)))
add_check("unique_primaryid", sum(duplicated(synthetic$primaryid)) + sum(is.na(synthetic$primaryid)),
          length(unique(synthetic$primaryid)))
for (nm in id_cols) {
  add_check(paste0("no_reference_id_overlap_", nm),
            sum(as.character(synthetic[[nm]]) %in% as.character(reference[[nm]])),
            "Synthetic ID versus Reference ID")
}
for (nm in nominal_cols) {
  invalid <- is.na(synthetic[[nm]]) | !(synthetic[[nm]] %in% unique(reference[[nm]]))
  add_check(paste0("valid_categories_", nm), sum(invalid),
            paste(sort(unique(synthetic[[nm]])), collapse = " | "))
}
for (nm in binary_cols) {
  val <- synthetic[[nm]]
  invalid <- !is.na(val) & !(val %in% c(0, 1))
  # Only clinical indication indicators are allowed structural missingness.
  missing_disallowed <- is.na(val) & !(nm %in% indication_cols)
  add_check(paste0("binary_0_1_", nm), sum(invalid | missing_disallowed),
            paste0("invalid/nonstructural missing=", sum(invalid | missing_disallowed)))
}
add_check("target_binary_nonmissing", sum(is.na(synthetic[[target]]) |
                                              !(synthetic[[target]] %in% c(0, 1))),
          paste(sort(unique(synthetic[[target]])), collapse = ","))
for (nm in c("num_drugs", "num_reactions")) {
  invalid <- is.na(synthetic[[nm]]) | !is.finite(synthetic[[nm]]) |
    synthetic[[nm]] < 1 | synthetic[[nm]] != floor(synthetic[[nm]])
  add_check(paste0("positive_integer_", nm), sum(invalid),
            paste0("min=", min(synthetic[[nm]], na.rm = TRUE)))
}
add_check("age_known_0_implies_age_NA",
          sum(synthetic$age_known == 0 & !is.na(synthetic$age_years), na.rm = TRUE),
          "Structural age missingness")
add_check("age_known_1_implies_age_observed",
          sum(synthetic$age_known == 1 & is.na(synthetic$age_years), na.rm = TRUE),
          "Known ages observed")
add_check("observed_ages_finite",
          sum(!is.na(synthetic$age_years) & !is.finite(synthetic$age_years)),
          "All observed ages numeric")
add_check("single_drug_no_concomitant",
          sum(synthetic$num_drugs == 1 & synthetic$has_concomitant_drug != 0, na.rm = TRUE),
          "num_drugs=1 -> has_concomitant_drug=0")
add_check("known_indication_requires_suspect_indication",
          sum(synthetic$indication_known == 1 & synthetic$suspect_indication_known != 1,
              na.rm = TRUE), "Accepted context implication")
context_key <- function(d) apply(d[context_cols], 1, paste, collapse = "||")
add_check("report_context_joint_combinations_valid",
          sum(!context_key(synthetic) %in% unique(context_key(reference))),
          "Every synthetic context combination occurs in reference")
unknown <- synthetic$indication_known == 0
known <- synthetic$indication_known == 1
add_check("unknown_indication_all_four_NA",
          sum(rowSums(!is.na(synthetic[unknown, indication_cols, drop = FALSE])) > 0),
          "All four fields must be true NA")
add_check("known_indication_all_four_binary",
          sum(rowSums(is.na(as.matrix(synthetic[known, indication_cols, drop = FALSE])) |
              matrix(!(as.matrix(synthetic[known, indication_cols, drop = FALSE]) %in% c(0, 1)),
                     nrow = sum(known), ncol = length(indication_cols))) > 0),
          "Known indication fields observed as 0/1")
react_count <- rowSums(synthetic[reaction_cols] == 1, na.rm = TRUE)
add_check("reaction_category_count_not_above_num_reactions",
          sum(react_count > synthetic$num_reactions, na.rm = TRUE),
          "Flagged reactions <= number of reaction terms")
save_output(checks, "structural_checks.csv")
structural_pass <- all(checks$status == "PASS")
if (!structural_pass) {
  stop("STRUCTURAL FAILURE: inspect structural_checks.csv. Synthetic v2 must not be frozen.")
}

# 3. Layer 2: all 22 predictor marginal distributions -------------------------
marginals <- data.frame(variable = character(), kind = character(), statistic = character(),
                        denominator = character(), reference = numeric(), synthetic = numeric(),
                        difference = numeric(), unit = character())
add_marginal <- function(v, kind, stat, denominator, ref, syn, unit) {
  marginals <<- rbind(marginals, data.frame(variable = v, kind = kind,
    statistic = as.character(stat), denominator = denominator,
    reference = ref, synthetic = syn, difference = syn - ref, unit = unit))
}
for (nm in binary_cols) {
  refx <- reference[[nm]]; synx <- synthetic[[nm]]
  denom <- if (nm %in% indication_cols) "indication_known=1 (nonmissing)" else "all cases"
  add_marginal(nm, "binary", "prevalence_1", denom,
               100 * rate(refx == 1), 100 * rate(synx == 1), "percentage points")
  add_marginal(nm, "binary", "missingness", "all cases",
               100 * mean(is.na(refx)), 100 * mean(is.na(synx)), "percentage points")
}
for (nm in nominal_cols) {
  for (category in sort(unique(c(reference[[nm]], synthetic[[nm]])))) {
    add_marginal(nm, "categorical", category, "all cases",
                 100 * proportion(reference[[nm]], category),
                 100 * proportion(synthetic[[nm]], category), "percentage points")
  }
}
for (nm in continuous_cols) {
  summary_fun <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE), Q25 = unname(quantile(x, .25, na.rm = TRUE)),
      Q75 = unname(quantile(x, .75, na.rm = TRUE)),
      Q90 = unname(quantile(x, .90, na.rm = TRUE)),
      Q95 = unname(quantile(x, .95, na.rm = TRUE)),
      Q99 = unname(quantile(x, .99, na.rm = TRUE)))
  rr <- summary_fun(reference[[nm]]); ss <- summary_fun(synthetic[[nm]])
  for (stat in names(rr)) add_marginal(nm, "continuous_or_count", stat,
     "nonmissing observations", unname(rr[stat]), unname(ss[stat]), "original units")
  add_marginal(nm, "continuous_or_count", "missingness", "all cases",
     100 * mean(is.na(reference[[nm]])), 100 * mean(is.na(synthetic[[nm]])), "percentage points")
}
stopifnot(setequal(unique(marginals$variable), all_predictors))
save_output(marginals, "marginal_comparison.csv")

# 4. Layer 3: within-block and secondary dependency diagnostics ---------------
deps <- data.frame(block = character(), priority = character(), relationship = character(),
                   stratum = character(), metric = character(), reference = numeric(),
                   synthetic = numeric(), difference = numeric(), unit = character())
add_dep <- function(block, priority, relationship, stratum, metric, r, s, unit = "percentage points") {
  deps <<- rbind(deps, data.frame(block = block, priority = priority,
     relationship = relationship, stratum = as.character(stratum), metric = metric,
     reference = r, synthetic = s, difference = s - r, unit = unit))
}
phi <- function(d, a, b) {
  x <- d[[a]]; y <- d[[b]]; ok <- complete.cases(x, y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
context_binary <- setdiff(context_cols, nominal_cols)
for (i in seq_len(length(context_binary) - 1)) {
  for (j in (i+1):length(context_binary)) {
    a <- context_binary[i]; b <- context_binary[j]
    add_dep("report_context", "intended", paste(a, b, sep = " x "), "all cases",
            "phi_correlation", phi(reference, a, b), phi(synthetic, a, b), "correlation")
  }
}
for (category_var in nominal_cols) {
  for (category in sort(unique(reference[[category_var]]))) {
    for (binary in context_binary) {
      sr <- reference[reference[[category_var]] == category, binary]
      ss <- synthetic[synthetic[[category_var]] == category, binary]
      add_dep("report_context", "intended", paste(category_var, binary, sep = " x "),
              category, "binary_prevalence", 100 * rate(sr == 1), 100 * rate(ss == 1))
    }
  }
}
# Medication complexity: every observed count as well as interpretable bands.
med_band <- function(d) ifelse(d$num_drugs == 1, "1", ifelse(d$num_drugs <= 4, "2-4", "5+"))
for (group in c("1", "2-4", "5+")) {
  r <- reference$has_concomitant_drug[med_band(reference) == group]
  s <- synthetic$has_concomitant_drug[med_band(synthetic) == group]
  add_dep("medication", "intended", "num_drugs x has_concomitant_drug", group,
          "concomitant_prevalence", 100 * rate(r == 1), 100 * rate(s == 1))
}
for (drug_n in sort(unique(c(reference$num_drugs, synthetic$num_drugs)))) {
  rr <- reference$has_concomitant_drug[reference$num_drugs == drug_n]
  ss <- synthetic$has_concomitant_drug[synthetic$num_drugs == drug_n]
  if (length(rr) >= 25 && length(ss) >= 25)
    add_dep("medication", "intended", "num_drugs x has_concomitant_drug",
            drug_n, "count_specific_concomitant_prevalence", 100*rate(rr==1),100*rate(ss==1))
}
react_band <- function(d) ifelse(d$num_reactions == 1, "1",
                         ifelse(d$num_reactions == 2, "2",
                         ifelse(d$num_reactions <= 4, "3-4", "5+")))
for (group in c("1", "2", "3-4", "5+")) {
  for (nm in reaction_cols) {
    rr <- reference[[nm]][react_band(reference) == group]
    ss <- synthetic[[nm]][react_band(synthetic) == group]
    add_dep("clinical_profile", "intended", paste("num_reactions x", nm),
            group, "reaction_prevalence", 100*rate(rr==1), 100*rate(ss==1))
  }
}
for (i in 1:(length(reaction_cols) - 1)) {
  for (j in (i + 1):length(reaction_cols)) {
    a <- reaction_cols[i]; b <- reaction_cols[j]
    add_dep("clinical_profile", "intended", paste(a,b,sep = " x "), "all cases",
            "joint_positive_prevalence", 100*mean(reference[[a]]==1 & reference[[b]]==1),
            100*mean(synthetic[[a]]==1 & synthetic[[b]]==1))
  }
}
clinical_counts <- function(d, cols) rowSums(d[cols] == 1, na.rm = TRUE)
for (k in 0:6) add_dep("clinical_profile", "intended", "number of reaction categories",
    k, "count_distribution", 100*mean(clinical_counts(reference,reaction_cols)==k),
    100*mean(clinical_counts(synthetic,reaction_cols)==k))
add_dep("clinical_profile", "intended", "reaction categories / num_reactions", "all cases",
   "mean_ratio", mean(clinical_counts(reference,reaction_cols)/reference$num_reactions),
   mean(clinical_counts(synthetic,reaction_cols)/synthetic$num_reactions), "ratio")
ind_pattern <- function(d) apply(d[indication_cols], 1, function(x)
   if (all(is.na(x))) "UNKNOWN" else paste(x, collapse = ""))
rp <- ind_pattern(reference); sp <- ind_pattern(synthetic)
for (pattern in sort(unique(c(rp,sp)))) {
  # Known-indication denominator for patterns; unknown included separately.
  if (pattern != "UNKNOWN") {
    r <- rp[reference$indication_known == 1]; s <- sp[synthetic$indication_known == 1]
    add_dep("clinical_profile", "intended", "indication four-bit pattern",pattern,
            "pattern_prevalence_among_known", 100*mean(r==pattern),100*mean(s==pattern))
  }
}
for (k in 0:4) {
  rr <- clinical_counts(reference[reference$indication_known==1,,drop=FALSE], indication_cols)
  ss <- clinical_counts(synthetic[synthetic$indication_known==1,,drop=FALSE], indication_cols)
  add_dep("clinical_profile", "intended", "indication category count", k,
          "count_prevalence_among_known",100*mean(rr==k),100*mean(ss==k))
}
for (rcol in reaction_cols) for (icol in indication_cols) for (flag in c(0,1)) {
   ra <- reference[[rcol]][!is.na(reference[[icol]]) & reference[[icol]]==flag]
   sa <- synthetic[[rcol]][!is.na(synthetic[[icol]]) & synthetic[[icol]]==flag]
   add_dep("clinical_profile", "intended", paste(rcol, icol,sep=" x "),
           paste0("indication=",flag),"reaction_prevalence_among_known",
           100*rate(ra==1),100*rate(sa==1))
}
# Intentionally simplified cross-block relationships - not hard failures.
age_band <- function(d) ifelse(is.na(d$age_years), "Unknown",
                       ifelse(d$age_years < 45,"<45",ifelse(d$age_years<65,"45-64","65+")))
for (nm in c(reaction_cols, indication_cols)) {
  for (axis in c("age_band", "drug_band", "reporter_type", "patient_sex")) {
    rgroup <- switch(axis, age_band=age_band(reference),drug_band=med_band(reference),
                     reporter_type=reference$reporter_type,patient_sex=reference$patient_sex)
    sgroup <- switch(axis, age_band=age_band(synthetic),drug_band=med_band(synthetic),
                     reporter_type=synthetic$reporter_type,patient_sex=synthetic$patient_sex)
    for (g in sort(unique(c(rgroup, sgroup)))) {
      rv <- reference[[nm]][rgroup==g]; sv <- synthetic[[nm]][sgroup==g]
      add_dep("cross_block", "intentionally_simplified", paste(axis,nm,sep=" x "),g,
              "clinical_prevalence_nonmissing",100*rate(rv==1),100*rate(sv==1))
    }
  }
}
save_output(deps, "dependency_comparison.csv")

# 5. Independently recreate accepted generation score in BOTH datasets -------
score_data <- function(d) {
  age <- ifelse(is.na(d$age_years), 0,
                ifelse(d$age_years >= 65, 2, ifelse(d$age_years >= 45, 1, 0)))
  medication <- ifelse(d$num_drugs >= 5, 2, ifelse(d$num_drugs >= 2, 1, 0))
  reporter <- as.integer(d$reporter_type %in% c("Consumer", "Physician"))
  # Unknown indications remain NA in d; score alone treats them as zero.
  indication <- as.integer((!is.na(d$indication_oncology) & d$indication_oncology==1) |
                           (!is.na(d$indication_respiratory) & d$indication_respiratory==1))
  reaction <- as.integer(d$reaction_systemic_symptoms==1 |
                         d$reaction_gastrointestinal==1)
  total <- age + medication + reporter + indication + reaction
  group <- c("Very Low", "Low", "Moderate", "Elevated")[pmin(total+1,4)]
  group[total>=4] <- "High"
  stopifnot(all(total>=0 & total<=7))
  data.frame(age=age, medication=medication, reporter=reporter,
             indication=indication, reaction=reaction, total=total, group=group,
             y=d[[target]])
}
rscore <- score_data(reference); sscore <- score_data(synthetic)

score_dist <- data.frame(scale=character(), category=character(), reference_n=integer(),
  synthetic_n=integer(), reference_pct=numeric(), synthetic_pct=numeric(), difference_pp=numeric())
for (sc in c(as.character(0:7), "Very Low", "Low", "Moderate", "Elevated", "High")) {
  if (sc %in% as.character(0:7)) {
    rn <- sum(rscore$total==as.integer(sc)); sn <- sum(sscore$total==as.integer(sc)); kind <- "raw_0_to_7"
  } else {
    rn <- sum(rscore$group==sc); sn <- sum(sscore$group==sc); kind <- "propensity_group"
  }
  score_dist <- rbind(score_dist, data.frame(scale=kind,category=sc,
    reference_n=rn,synthetic_n=sn,reference_pct=100*rn/nrow(reference),
    synthetic_pct=100*sn/nrow(synthetic),
    difference_pp=100*(sn/nrow(synthetic)-rn/nrow(reference))))
}
save_output(score_dist, "score_distribution.csv")

target_diag <- data.frame(dimension=character(), category=character(),
   reference_n=integer(), synthetic_n=integer(), reference_positives=integer(),
   synthetic_positives=integer(), reference_target_pct=numeric(),
   synthetic_target_pct=numeric(), difference_pp=numeric())
add_target <- function(dim, cat, rmask, smask) {
  rn <- sum(rmask); sn <- sum(smask)
  ry <- sum(rscore$y[rmask]==1); sy <- sum(sscore$y[smask]==1)
  target_diag <<- rbind(target_diag, data.frame(dimension=dim,category=as.character(cat),
      reference_n=rn,synthetic_n=sn,reference_positives=ry,synthetic_positives=sy,
      reference_target_pct=if(rn) 100*ry/rn else NA_real_,
      synthetic_target_pct=if(sn) 100*sy/sn else NA_real_,
      difference_pp=if(rn && sn) 100*(sy/sn - ry/rn) else NA_real_))
}
add_target("overall","all",rep(TRUE,nrow(reference)),rep(TRUE,nrow(synthetic)))
for (v in 0:7) add_target("raw_score",v,rscore$total==v,sscore$total==v)
for (v in c("Very Low","Low","Moderate","Elevated","High"))
  add_target("propensity_group",v,rscore$group==v,sscore$group==v)
for (nm in c("age","medication","reporter","indication","reaction")) {
  for (v in sort(unique(c(rscore[[nm]], sscore[[nm]]))))
    add_target(paste0("score_component_",nm),v,rscore[[nm]]==v,sscore[[nm]]==v)
}
save_output(target_diag,"target_diagnostics.csv")

clinical_target <- data.frame(predictor=character(), predictor_value=integer(),
  denominator=character(), reference_n=integer(), synthetic_n=integer(),
  reference_target_pct=numeric(), synthetic_target_pct=numeric(), difference_pp=numeric())
for (nm in c(reaction_cols, indication_cols)) for (v in c(0,1)) {
  rm <- !is.na(reference[[nm]]) & reference[[nm]]==v
  sm <- !is.na(synthetic[[nm]]) & synthetic[[nm]]==v
  rv <- 100*rate(reference[[target]][rm]==1)
  sv <- 100*rate(synthetic[[target]][sm]==1)
  clinical_target <- rbind(clinical_target, data.frame(predictor=nm, predictor_value=v,
    denominator=if(nm %in% indication_cols) "indication_known=1" else "all cases",
    reference_n=sum(rm), synthetic_n=sum(sm), reference_target_pct=rv,
    synthetic_target_pct=sv, difference_pp=sv-rv))
}
save_output(clinical_target,"clinical_target_comparison.csv")

# 6. Layer 5: score-family leave-one-out ---------------------------------------
auc_rank <- function(scores, outcomes) {
  positives <- outcomes == 1; n1 <- sum(positives); n0 <- sum(!positives)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(scores, ties.method = "average")[positives]) - n1*(n1+1)/2) / (n1*n0)
}
sensitivity <- data.frame(removed_family=character(), reference_auc=numeric(),
    synthetic_auc=numeric(), reference_auc_loss=numeric(), synthetic_auc_loss=numeric(),
    reference_score_target_correlation=numeric(),
    synthetic_score_target_correlation=numeric())
full_r_auc <- auc_rank(rscore$total,rscore$y); full_s_auc <- auc_rank(sscore$total,sscore$y)
for (family in c("none","age","medication","reporter","indication","reaction")) {
  rr <- if(family=="none") rscore$total else rscore$total-rscore[[family]]
  ss <- if(family=="none") sscore$total else sscore$total-sscore[[family]]
  ra <- auc_rank(rr,rscore$y); sa <- auc_rank(ss,sscore$y)
  sensitivity <- rbind(sensitivity,data.frame(removed_family=family,
      reference_auc=ra,synthetic_auc=sa, reference_auc_loss=full_r_auc-ra,
      synthetic_auc_loss=full_s_auc-sa,
      reference_score_target_correlation=cor(rr,rscore$y),
      synthetic_score_target_correlation=cor(ss,sscore$y)))
}
save_output(sensitivity,"score_family_sensitivity.csv")

# 7. Layer 6: logistic sanity-check 
# Identifiers, generated score/group and future information never enter X.
# Structural indication NA becomes 0 FOR MODEL INPUT ONLY. indication_known stays.
# Age median is computed only from the train partition; age_known remains included.
set.seed(123)
positive <- which(synthetic[[target]]==1)
negative <- which(synthetic[[target]]==0)
train_idx <- sort(c(sample(positive, floor(.8*length(positive))),
                    sample(negative, floor(.8*length(negative)))))
test_idx <- setdiff(seq_len(nrow(synthetic)),train_idx)
train <- synthetic[train_idx,,drop=FALSE]; test <- synthetic[test_idx,,drop=FALSE]
train_median_age <- median(train$age_years,na.rm=TRUE)
prepare <- function(d, training_levels) {
  x <- d[all_predictors]
  x$age_years[is.na(x$age_years)] <- train_median_age
  for (nm in indication_cols) x[[nm]][is.na(x[[nm]])] <- 0
  for (nm in nominal_cols) {
    unexpected <- !x[[nm]] %in% training_levels[[nm]]
    if (any(unexpected)) stop("Unseen test level in ", nm, "; review categorical encoding.")
    x[[nm]] <- factor(x[[nm]],levels=training_levels[[nm]])
  }
  x
}
training_levels <- lapply(train[nominal_cols], function(x) sort(unique(x)))
train_x <- prepare(train,training_levels)
test_x <- prepare(test,training_levels)
# Standard glm uses treatment contrasts and trains on training data alone.
train_model <- cbind(data.frame(y=train[[target]]),train_x)
logit <- suppressWarnings(glm(y ~ .,data=train_model,family=binomial()))
if (!isTRUE(logit$converged)) warning("Validation-only logistic regression did not converge")
pred <- as.numeric(predict(logit,newdata=test_x,type="response"))
if (any(!is.finite(pred))) stop("Non-finite logistic predictions")
ytest <- test[[target]]
roc_auc <- auc_rank(pred,ytest)
# Trapezoidal area under the precision-recall curve (group tied scores).
order_idx <- order(pred,decreasing=TRUE)
y_sorted <- ytest[order_idx]; p_sorted <- pred[order_idx]
last_idx <- which(!duplicated(p_sorted,fromLast=TRUE))
tp <- cumsum(y_sorted)[last_idx]; retrieved <- last_idx
recall <- c(0,tp/sum(ytest==1))
precision <- c(1,tp/retrieved)
pr_auc <- sum(diff(recall) * (head(precision,-1)+tail(precision,-1))/2)
review_n <- ceiling(.2*length(ytest))
recall20 <- sum(y_sorted[seq_len(review_n)]==1)/sum(ytest==1)
model_check <- data.frame(metric=c("train_n","test_n","test_prevalence_pct",
  "ROC_AUC","PR_AUC_trapezoid","Recall_at_20pct","converged","train_age_median"),
  value=c(nrow(train),nrow(test),100*mean(ytest==1),roc_auc,pr_auc,recall20,
          as.integer(logit$converged),train_median_age))
save_output(model_check,"predictive_plausibility.csv")

# 8. Evidence-based concise status table --------------------------------------
# Numerical differences are diagnostics: only structural rules have automatic hard-fail cutoffs. 
summary_rows <- data.frame(category=c(
  "Structural integrity","Marginal realism","Intended within-block dependencies",
  "Secondary cross-block dependencies","Target prevalence","Target propensity gradient",
  "Score/group distribution","Individual clinical predictor-target relationships",
  "Score-family balance","Predictive plausibility","Overall recommendation"),
  status=c(if(structural_pass) "PASS" else "FAIL",rep("REVIEW",9),
           if(structural_pass) "REVIEW" else "FAIL"),
  evidence=c(
    paste(nrow(checks),"hard checks;",sum(checks$status=="FAIL"),"failed"),
    paste(length(unique(marginals$variable)),"predictors compared"),
    "Report context, medication, clinical profile in dependency_comparison.csv",
    "Intentionally unconditioned relationships; interpret, do not hard-fail",
    sprintf("Reference %.2f%%; Synthetic %.2f%%",100*mean(reference[[target]]),100*mean(synthetic[[target]])),
    "Target prevalence by all raw scores, groups and score families",
    "Distribution of scores 0-7 and five propensity groups",
    "Reaction and indication indicators, 0/1; known indications only",
    "AUC change after removing each of five families",
    sprintf("Held-out ROC-AUC %.3f; PR-AUC %.3f; Recall@20%% %.3f",roc_auc,pr_auc,recall20),
    if(structural_pass) "Review diagnostics before human approval to freeze" else "Structural failure blocks freeze"))
save_output(summary_rows,"validation_summary.csv")
cat("\nSynthetic v2 validation completed. Outputs:", out_dir,"\n")
cat("Structural checks:",sum(checks$status=="PASS"),"PASS /",
    sum(checks$status=="FAIL"),"FAIL\n")
cat("Target prevalence: reference",fmt(100*mean(reference[[target]])),"%; synthetic",
    fmt(100*mean(synthetic[[target]])),"%\n")
cat("Held-out: ROC-AUC",fmt(roc_auc),"PR-AUC",fmt(pr_auc),
    "Recall@20%",fmt(recall20),"\n")
if (!structural_pass) stop("Cannot freeze Synthetic v2: hard structural failure detected.")
