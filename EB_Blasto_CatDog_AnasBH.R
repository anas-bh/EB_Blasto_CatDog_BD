# =============================================================================
# Title: Molecular Epidemiology of Enterocytozoon bieneusi and Blastocystis
#        in Cats and Dogs in Dhaka and Gazipur, Bangladesh
# Authors: Anas Bin Harun, Abdullah Al Bayazid, Md. Farhan Hasan,
#          Ainun Nahar, Maksuda Akter Mily, Jannatul Bakia Tamanna,
#          Joynti Saha, Md. Ashiqur Rahman, S. H. M. Faruk Siddiki,
#          Md. Morshedur Rahman, Md Robiul Karim
# Journal: Food and Waterborne Parasitology
# Description: R script for statistical analysis including chi-square tests,
#              Firth penalized logistic regression, ROC curve analysis, and
#              Random Forest variable importance.
# Developed by: Anas Bin Harun, DVM, MS
# =============================================================================


# =============================================================================
# SECTION 1: LOAD REQUIRED PACKAGES
# =============================================================================
# Install any missing packages before running:
# install.packages(c("readxl", "logistf", "pROC", "randomForest"))

options(error = NULL)   # Ensure errors print as messages and do not interrupt execution

library(readxl)       # Reading Excel data files
library(logistf)      # Firth penalized likelihood logistic regression
library(pROC)         # ROC curve analysis and AUC calculation
library(randomForest) # Random Forest classification and variable importance


# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================
# Place CatDog_AnasBH.xlsx in the same folder as this script, or update the
# path below to the full file location on your system.

cat_data <- read_excel("CatDog_AnasBH.xlsx", sheet = "Cat")
dog_data <- read_excel("CatDog_AnasBH.xlsx", sheet = "Dog")


# =============================================================================
# SECTION 3: VARIABLE OVERVIEW
# =============================================================================
# cat sheet columns:
#   loc     — sampling location (Dhaka, Gazipur)
#   age     — age group (<1 years, >=1 year)
#   sex     — sex (Female, Male)
#   breed   — breed type (Local, Cross)
#   source  — animal source (FnF, Shelter/Breeder, Stray)
#   deworm  — deworming history (Yes, No)
#   cohab   — cohabitation with other animals (Yes, No)
#   out     — outdoor access (Indoor, Outdoor)
#   food    — food type (Meat-based, Mixed)
#   blasto  — Blastocystis outcome (0 = negative, 1 = positive)
#   eb      — E. bieneusi outcome (0 = negative, 1 = positive)
#   hunt_excluded — hunting variable (excluded from models; see Methods)
#
# dog sheet columns:
#   loc     — sampling location (Dhaka, Gazipur)
#   age     — age group (<1 years, 1-2 years, >2 years)
#   sex     — sex (Female, Male)
#   breed   — breed type (Cross, Local)
#   deworm  — deworming history (Yes, No)
#   cohab   — cohabitation with other animals (Yes, No)
#   out     — outdoor access (Yes, No)
#   food    — food type (Meat-based (cooked/raw), Mixed)
#   water   — water source (Filtered water, Tap water)
#   blasto  — Blastocystis outcome (0 = negative, 1 = positive)
#   eb      — E. bieneusi outcome (0 = negative, 1 = positive)
#   hunt_excluded — hunting variable (excluded from models; see Methods)


# =============================================================================
# SECTION 4: CONVERT TO FACTORS
# =============================================================================

cat_factor_vars <- c("loc","age","sex","breed","source","deworm","cohab","out","food")
cat_data[cat_factor_vars] <- lapply(cat_data[cat_factor_vars], factor)
cat_data$eb     <- as.integer(cat_data$eb)
cat_data$blasto <- as.integer(cat_data$blasto)

dog_factor_vars <- c("loc","age","sex","breed","deworm","cohab","out","food","water")
dog_data[dog_factor_vars] <- lapply(dog_data[dog_factor_vars], factor)
dog_data$eb     <- as.integer(dog_data$eb)
dog_data$blasto <- as.integer(dog_data$blasto)

# Sample sizes
cat("Sample sizes:\n")
cat("  Cats: n =", nrow(cat_data),
    "| E. bieneusi positives =", sum(cat_data$eb),
    "| Blastocystis positives =", sum(cat_data$blasto), "\n")
cat("  Dogs: n =", nrow(dog_data),
    "| E. bieneusi positives =", sum(dog_data$eb),
    "| Blastocystis positives =", sum(dog_data$blasto), "\n")


# =============================================================================
# SECTION 5: CHI-SQUARE TESTS (UNIVARIABLE ANALYSIS)
# =============================================================================

run_chi_sq <- function(df, outcome_var, predictor_vars, label) {
  cat(sprintf("\n--- Chi-square univariable analysis: %s ---\n", label))
  cat(sprintf("%-20s %10s %5s %10s\n", "Variable", "Chi-sq", "df", "p-value"))
  for (v in predictor_vars) {
    tbl  <- table(df[[v]], df[[outcome_var]])
    test <- tryCatch(chisq.test(tbl, correct = FALSE),
                     warning = function(w) suppressWarnings(chisq.test(tbl, correct = FALSE)))
    sig  <- ifelse(test$p.value < 0.05, " *", "")
    cat(sprintf("%-20s %10.3f %5g %10.4f%s\n",
                v, as.numeric(test$statistic),
                test$parameter, test$p.value, sig))
  }
}

cat_preds <- c("loc","age","sex","breed","source","deworm","cohab","out","food")
dog_preds <- c("loc","age","sex","breed","deworm","cohab","out","food","water")

cat("\n=== CHI-SQUARE RESULTS — CATS ===\n")
run_chi_sq(cat_data, "eb",     cat_preds, "Cat E. bieneusi")
run_chi_sq(cat_data, "blasto", cat_preds, "Cat Blastocystis")

cat("\n=== CHI-SQUARE RESULTS — DOGS ===\n")
run_chi_sq(dog_data, "eb",     dog_preds, "Dog E. bieneusi")
run_chi_sq(dog_data, "blasto", dog_preds, "Dog Blastocystis")


# =============================================================================
# SECTION 6: FIRTH PENALIZED LOGISTIC REGRESSION (MULTIVARIABLE ANALYSIS)
# =============================================================================
# Firth penalized likelihood logistic regression was used in place of standard
# maximum likelihood logistic regression due to sparse positive counts across
# subgroups and quasi-complete separation in certain predictor categories.
# Reference: Heinze G, Schemper M. Stat Med. 2002;21(16):2409-2419.
#
# Reference categories were selected as the group with the lowest infection
# prevalence per variable.

run_firth <- function(data, outcome_var, predictor_vars, ref_levels, label) {
  cat(sprintf("\n--- Firth logistic regression: %s ---\n", label))

  for (v in predictor_vars) {
    data[[v]] <- factor(data[[v]])
    if (!is.null(ref_levels[[v]])) {
      data[[v]] <- relevel(data[[v]], ref = ref_levels[[v]])
    }
  }

  formula_str <- paste(outcome_var, "~", paste(predictor_vars, collapse = " + "))

  result <- tryCatch({
    model   <- logistf(as.formula(formula_str), data = data)
    or_vals <- exp(coef(model))
    ci_lo   <- exp(model$ci.lower)
    ci_hi   <- exp(model$ci.upper)
    p_vals  <- model$prob
    data.frame(
      Variable = names(or_vals),
      OR       = round(or_vals, 3),
      CI_Lower = round(ci_lo,   3),
      CI_Upper = round(ci_hi,   3),
      P_value  = round(p_vals,  4),
      row.names = NULL
    )
  }, error = function(e) {
    cat(sprintf("  Model could not be fitted: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (!is.null(result)) print(result)
  return(invisible(result))
}

# Reference categories — lowest EB prevalence per variable
cat_ref_eb <- list(
  loc    = "Dhaka",
  age    = ">=1 year",
  sex    = "Male",
  breed  = "Cross",
  source = "Stray",
  deworm = "Yes",
  cohab  = "No",
  out    = "Indoor",
  food   = "Mixed"
)

cat_ref_blasto <- list(
  loc    = "Dhaka",
  age    = ">=1 year",
  sex    = "Male",
  breed  = "Cross",
  source = "Stray",
  deworm = "Yes",
  cohab  = "No",
  out    = "Indoor",
  food   = "Mixed"
)

dog_ref_eb <- list(
  loc    = "Dhaka",
  age    = "1-2 years",
  sex    = "Male",
  breed  = "Local",
  deworm = "Yes",
  cohab  = "No",
  out    = "Yes",
  food   = "Meat-based (cooked/raw)",
  water  = "Filtered water"
)

dog_ref_blasto <- list(
  loc    = "Dhaka",
  age    = "1-2 years",
  sex    = "Male",
  breed  = "Local",
  deworm = "Yes",
  cohab  = "No",
  out    = "Yes",
  food   = "Meat-based (cooked/raw)",
  water  = "Filtered water"
)

cat("\n=== FIRTH REGRESSION RESULTS — CATS ===\n")
cat_eb_res     <- run_firth(cat_data, "eb",     cat_preds, cat_ref_eb,     "Cat E. bieneusi")
cat_blasto_res <- run_firth(cat_data, "blasto", cat_preds, cat_ref_blasto, "Cat Blastocystis")

cat("\n=== FIRTH REGRESSION RESULTS — DOGS ===\n")
dog_eb_res     <- run_firth(dog_data, "eb",     dog_preds, dog_ref_eb,     "Dog E. bieneusi")
dog_blasto_res <- run_firth(dog_data, "blasto", dog_preds, dog_ref_blasto, "Dog Blastocystis")


# =============================================================================
# SECTION 7: ROC CURVE ANALYSIS (SUPPLEMENTARY FIGURES S1-S2)
# =============================================================================
# Results are presented as exploratory analyses only. Given the limited number
# of positive cases across all models, AUC values should not be interpreted
# as confirmatory evidence of predictive performance.

cat_data$eb_f     <- factor(cat_data$eb,     levels=c(0,1), labels=c("Negative","Positive"))
cat_data$blasto_f <- factor(cat_data$blasto, levels=c(0,1), labels=c("Negative","Positive"))
dog_data$eb_f     <- factor(dog_data$eb,     levels=c(0,1), labels=c("Negative","Positive"))
dog_data$blasto_f <- factor(dog_data$blasto, levels=c(0,1), labels=c("Negative","Positive"))

run_roc <- function(df, outcome_f, predictors, label) {
  tryCatch({
    formula_str <- paste(outcome_f, "~", paste(predictors, collapse = "+"))
    model   <- glm(as.formula(formula_str), data = df, family = binomial)
    probs   <- predict(model, type = "response")
    roc_obj <- roc(df[[outcome_f]], probs, quiet = TRUE)
    auc_val <- auc(roc_obj)
    plot(roc_obj,
         main = paste("ROC Curve -", label),
         sub  = paste("AUC =", round(auc_val, 3),
                      "| n positives =", sum(df[[sub("_f","",outcome_f)]])),
         col  = "#185FA5", lwd = 2)
    abline(a=0, b=1, lty=2, col="gray60")
    cat(sprintf("AUC (%s): %.4f\n", label, auc_val))
    return(invisible(roc_obj))
  }, error = function(e) {
    cat(sprintf("  ROC could not be computed (%s): %s\n", label, conditionMessage(e)))
    return(NULL)
  })
}

cat("\n=== ROC CURVE ANALYSIS ===\n")

# Figure S1: Cats
par(mfrow = c(1, 2))
roc_cat_eb     <- run_roc(cat_data, "eb_f",     cat_preds, "Cat E. bieneusi")
roc_cat_blasto <- run_roc(cat_data, "blasto_f", cat_preds, "Cat Blastocystis")

# Figure S2: Dogs
par(mfrow = c(1, 2))
roc_dog_eb     <- run_roc(dog_data, "eb_f",     dog_preds, "Dog E. bieneusi")
roc_dog_blasto <- run_roc(dog_data, "blasto_f", dog_preds, "Dog Blastocystis")

par(mfrow = c(1, 1))


# =============================================================================
# SECTION 8: RANDOM FOREST - VARIABLE IMPORTANCE (SUPPLEMENTARY FIGURES S1-S2)
# =============================================================================
# Results are presented as exploratory analyses only.

run_rf <- function(df, outcome_f, predictors, label, ntree=500, seed=123) {
  tryCatch({
    set.seed(seed)
    formula_str <- paste(outcome_f, "~", paste(predictors, collapse = "+"))
    rf_model <- randomForest(as.formula(formula_str),
                             data = df, importance = TRUE, ntree = ntree)
    varImpPlot(rf_model,
               main = paste("Variable Importance -", label),
               sub  = paste("n positives =",
                            sum(df[[sub("_f","",outcome_f)]])))
    cat(sprintf("\nVariable importance (%s):\n", label))
    print(importance(rf_model))
    return(invisible(rf_model))
  }, error = function(e) {
    cat(sprintf("  Random Forest could not be fitted (%s): %s\n", label, conditionMessage(e)))
    return(NULL)
  })
}

cat("\n=== RANDOM FOREST — CATS ===\n")
rf_cat_eb     <- run_rf(cat_data, "eb_f",     cat_preds, "Cat E. bieneusi")
rf_cat_blasto <- run_rf(cat_data, "blasto_f", cat_preds, "Cat Blastocystis")

cat("\n=== RANDOM FOREST — DOGS ===\n")
rf_dog_eb     <- run_rf(dog_data, "eb_f",     dog_preds, "Dog E. bieneusi")
rf_dog_blasto <- run_rf(dog_data, "blasto_f", dog_preds, "Dog Blastocystis")


# =============================================================================
# SECTION 9: SESSION INFO
# =============================================================================
cat("\n=== SESSION INFORMATION ===\n")
sessionInfo()

# =============================================================================
# END OF SCRIPT
# =============================================================================
