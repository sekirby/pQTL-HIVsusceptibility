##### data & library loading #####
packages <- c("dplyr", "tidyr", "stats", "data.table","purrr", "emmeans", "car", "MASS", "broom", "rstatix")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


data <- read.csv("dummy_data.csv")

# factors
data$HIV <- as.factor(data$HIV)
data$rs143864957 <- as.factor(data$rs143864957)
data$rs2920282 <- as.factor(data$rs2920282)
data$rs8106317 <- as.factor(data$rs8106317)
data$rs147944114 <- as.factor(data$rs147944114)
data$CST <- as.factor(data$CST)
data$CONTRA <- as.factor(data$CONTRA)
data$studyArm <- as.factor(data$studyArm)
data$condom_coded <- as.factor(data$condom_coded)
data$SEX_7DAYS <- as.numeric(data$SEX_7DAYS)
data$HSV_E <- as.factor(data$HSV_E)
data$SITE <- as.factor(data$SITE)
data$TOTAL.BAC.N <- scale(data$TOTAL.BAC)

# filters
N589 <- data %>% filter(CST != "NA")
N580 <- N589 %>% filter(inflammation != "NA")

##### inflammation scoring #####
cytokines <- c("IL_1B", "IL_6", "IP_10", "MCP_1", "MIP_1A", "MIP_1B", "TNF_A", "IL_8", "IL_1A")

# Calculate the upper quartile (75th percentile) for each cytokine
upper_quartiles <- apply(N580[cytokines], 2, quantile, probs = 0.75, na.rm = TRUE)

# Create binary columns
for (cytokine in cytokines) {
  colname <- paste0(cytokine, "_high")
  N580[[colname]] <- as.integer(N580[[cytokine]] > upper_quartiles[cytokine])
}

# Calculate inflammation score
high_cols <- paste0(cytokines, "_high")
N580$inflammation_score <- rowSums(N580[high_cols], na.rm = TRUE)

N580$inflammation_binary <- ifelse(N580$inflammation_score < 3, 0L, 1L)
N580$inflammation_binary <- as.factor(N580$inflammation_binary)


##### SECTION HEADER: "The genital proteome is genetically regulated and linked to HIV acquisition" #####
### HIV acqusition odds #####
vars <- c("rs2920282","rs8106317","rs147944114","rs143864957")
prots <- c("A2ML1","PSCA","GSTO1","CD177")

run_models <- function(preds, group_name) {
  
  results <- data.frame()
  
  for (i in preds) {
    
    HIV_model <- glm(
      HIV ~ N589[[i]] + C1 + AGE + as.factor(CONTRA) + CST + SITE + studyArm + condom_coded + SEX_7DAYS + HSV_E,
      family = binomial,
      data = N589
    )
    
    summary_model <- summary(HIV_model)$coefficients
    
    coeff <- summary_model[2, "Estimate"]
    std_err <- summary_model[2, "Std. Error"]
    
    co_lower <- coeff - 1.96 * std_err
    co_upper <- coeff + 1.96 * std_err
    
    or <- exp(coeff)
    ci_lower <- exp(co_lower)
    ci_upper <- exp(co_upper)
    
    p_value <- summary_model[2, "Pr(>|z|)"]
    
    results <- rbind(
      results,
      data.frame(
        Group = group_name,
        Variable = i,
        COEFF = round(coeff, 4),
        STDERR = round(std_err, 4),
        CO_Lower = round(co_lower, 4),
        CO_Upper = round(co_upper, 4),
        OR = round(or, 4),
        CI_Lower = round(ci_lower, 4),
        CI_Upper = round(ci_upper, 4),
        p_value = p_value
      )
    )
  }
  
  # FDR correction
  results$FDR <- p.adjust(results$p_value, method = "fdr")
  
  results
}

var_results <- run_models(vars, "Variant")
prot_results <- run_models(prots, "Protein")

all_HIV_results <- rbind(var_results, prot_results)

print("Predictors of HIV acquisition - binomial logistic regressions")
print(all_HIV_results, row.names = FALSE)

### does HIV risk differ between study arm?
interaction_model <- glm(
  HIV ~ rs143864957 * studyArm +
    C1 + AGE + as.factor(CONTRA) + CST + SITE +
    condom_coded + SEX_7DAYS + HSV_E,
  family = binomial,
  data = N589
)

summary(interaction_model)


### predicted probability #####
N589$A2ML1_col <- N589[["A2ML1"]]

# Fit model
HIV_model_A2ML1 <- glm(
  HIV ~ A2ML1_col + C1 + AGE + as.factor(CONTRA) + CST + SITE +
    studyArm + condom_coded + SEX_7DAYS + HSV_E,
  family = binomial,
  data = N589
)

# Prediction data
pred_data <- data.frame(
  A2ML1_col = seq(
    from = min(N589$A2ML1_col, na.rm = TRUE),
    to = max(N589$A2ML1_col, na.rm = TRUE),
    length.out = 100
  ),
  C1 = mean(N589$C1, na.rm = TRUE),
  AGE = mean(N589$AGE, na.rm = TRUE),
  CONTRA = names(which.max(table(N589$CONTRA))),
  CST = names(which.max(table(N589$CST))),
  SITE = names(which.max(table(N589$SITE))),
  studyArm = 1,
  condom_coded = 1,
  SEX_7DAYS = mean(N589$SEX_7DAYS, na.rm = TRUE),
  HSV_E = names(which.max(table(N589$HSV_E)))
)

# Ensure factor levels match model data
pred_data$CONTRA <- factor(pred_data$CONTRA, levels = levels(as.factor(N589$CONTRA)))
pred_data$CST <- factor(pred_data$CST, levels = levels(N589$CST))
pred_data$SITE <- factor(pred_data$SITE, levels = levels(N589$SITE))
pred_data$studyArm <- factor(pred_data$studyArm, levels = levels(N589$studyArm))
pred_data$HSV_E <- factor(pred_data$HSV_E, levels = levels(N589$HSV_E))
pred_data$condom_coded <- factor(pred_data$condom_coded, levels = levels(N589$condom_coded))


# Predictions
predictions <- predict(
  HIV_model_A2ML1,
  newdata = pred_data,
  type = "link",
  se.fit = TRUE
)

# Convert to probabilities
pred_data$probability <- plogis(predictions$fit)
pred_data$lower <- plogis(predictions$fit - 1.96 * predictions$se.fit)
pred_data$upper <- plogis(predictions$fit + 1.96 * predictions$se.fit)


##### SECTION HEADER: "A2ML1 is a predictor of genital mucosal inflammation" #####
### inflammation odds #####
HIV_vars <- c("rs143864957", "A2ML1")
all_inflam_results <- data.frame()

for (i in HIV_vars) {
  # Fit
  inflam_model <- glm(inflammation_binary ~ N580[[i]] + C1 + AGE + as.factor(CONTRA) + CST + SITE,
                      family = binomial, data = N580)
  
  # Extract
  summary_model <- summary(inflam_model)$coefficients
  coeff <- summary_model[2, "Estimate"]
  std_err <- summary_model[2, "Std. Error"]
  
  # Compute confidence intervals and odds ratio
  co_lower <- coeff - 1.96 * std_err
  co_upper <- coeff + 1.96 * std_err
  or <- exp(coeff)
  ci_lower <- exp(co_lower)
  ci_upper <- exp(co_upper)
  p_value <- summary_model[2, "Pr(>|z|)"]
  
  # Store
  inflam_result <- data.frame(
    Variable = i,
    COEFF = round(coeff, 4),
    STDERR = round(std_err, 4),
    CO_Lower = round(co_lower, 4),
    CO_Upper = round(co_upper, 4),
    OR = round(or, 4),
    CI_Lower = round(ci_lower, 4),
    CI_Upper = round(ci_upper, 4),
    p_value = p_value
  )
  
  # Combine
  all_inflam_results <- rbind(all_inflam_results, inflam_result)
}

print("Predictors of inflammation - binomial logistic regressions")
print(all_inflam_results, row.names = FALSE)

# Reshape N580 to long format with only the cytokines of interest
data_long <- N580 %>%
  dplyr::select(all_of(c(cytokines))) %>%  
  pivot_longer(cols = all_of(cytokines),
               names_to = "Cytokine",
               values_to = "Measurement")

# Calculate Q1, Median, Q3
quartiles <- data_long %>%
  group_by(Cytokine) %>%
  summarise(
    Q1 = quantile(Measurement, 0.25, na.rm = TRUE),
    Median = median(Measurement, na.rm = TRUE),
    Q3 = quantile(Measurement, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(quartiles)

### A2ML1 x longitudinal inflammatory changes #####
# see script "a2ml1_longitudinal_inflammation.R"


##### SECTION HEADER: "The rs143864957-A2ML1 association does not interact with biological or epidemiological strata" #####
### A2ML1 ~ rs143864957*covariate models #####
# null model
m0 <- lm(A2ML1 ~ rs143864957 + CST + C1 + AGE + as.factor(CONTRA) + SITE, data = N589)

# CST interaction
mCST <- lm(A2ML1 ~ rs143864957*CST + C1 + AGE + as.factor(CONTRA) + SITE, data = N589)
summary(mCST)
anova(m0, mCST)

# contraceptive interaction
mCONTRA <- lm(A2ML1 ~ rs143864957*as.factor(CONTRA) + CST + C1 + AGE + SITE, data = N589)
summary(mCONTRA)
anova(m0, mCONTRA)

# site interaction
mSITE <- lm(A2ML1 ~ rs143864957*SITE + CST + C1 + AGE + as.factor(CONTRA), data = N589)
summary(mSITE)
anova(m0, mSITE)


##### SECTION HEADER: "Vaginal microbiome composition and bacterial spectral count are independently associated with A2ML1 abundance" #####
### A2ML1 ~ proteome covariates #####
m.full <- lm(
  A2ML1 ~ CST + as.factor(CONTRA) + SITE + C1 + AGE,
  data = N589
)

# CST
m.noCST <- update(m.full, . ~ . - CST)
p_CST <- anova(m.noCST, m.full)$`Pr(>F)`[2]

# CONTRA
m.noCONTRA <- update(m.full, . ~ . - as.factor(CONTRA))
p_CONTRA <- anova(m.noCONTRA, m.full)$`Pr(>F)`[2]

# SITE
m.noSITE <- update(m.full, . ~ . - SITE)
p_SITE <- anova(m.noSITE, m.full)$`Pr(>F)`[2]

pvals <- c(
  CST = p_CST,
  CONTRA = p_CONTRA,
  SITE = p_SITE
)

pvals
p.adjust(pvals, method = "BH")

### Post hoc pairwise CST comparisons ######
m.full <- lm(
  A2ML1 ~ CST + as.factor(CONTRA) + SITE + C1 + AGE,
  data = N589
)

summary(m.full)

emm <- emmeans(m.full, ~ CST)

unad <- pairs(emm, adjust = 'none')
unad_df <- as.data.frame(unad)
unad_df[, c("contrast", "estimate", "p.value")]
pairs(emm, adjust = "none") |> confint()

adj <- pairs(emm, adjust = 'BH')
adj_df <- as.data.frame(adj)
adj_df[, c("contrast", "estimate", "p.value")]

# Print median A2ML1 by CST to console
medians <- aggregate(A2ML1 ~ CST, data = N589, FUN = median, na.rm = TRUE)
print(medians)

### A2ML1 x Bacterial spectral count #####
mTOTALBAC <- lm(A2ML1 ~ TOTAL.BAC + C1 + AGE + as.factor(CONTRA) + SITE, data = N589)
summary(mTOTALBAC)
confint(mTOTALBAC)

# colinearity
cor.test(N589$TOTAL.BAC, N589$hcount, method = "spearman")

mVIF <- lm(A2ML1 ~ TOTAL.BAC + hcount, data = N589)
vif(mVIF)

### A2ML1 ~ CST + TOTAL.BAC independence ######
mINDEP <- lm(A2ML1 ~ CST + TOTAL.BAC + SITE + as.factor(CONTRA) + C1 + AGE, data = N589)
summary(mINDEP)
confint(mINDEP)

# partial F-test
m.cst <- lm(
  A2ML1 ~ CST + C1 + AGE + as.factor(CONTRA) + SITE,
  data = N589
)

m.cst.bac <- lm(
  A2ML1 ~ CST + TOTAL.BAC + C1 + AGE + as.factor(CONTRA) + SITE,
  data = N589
)

anova(m.cst, m.cst.bac)

## does cst add information beyond bacterial load
m.bac <- lm(
  A2ML1 ~ TOTAL.BAC + C1 + AGE + as.factor(CONTRA) + SITE,
  data = N589
)

anova(m.bac, m.cst.bac)

##### SECTION HEADER: "Microbiome characteristics exhibit nominal interaction with the association between A2ML1 and inflammation" #####
### A2ML1 ~ CST + TOTAL.BAC + inflammation #####
m0 <- lm(
  A2ML1 ~ CST + TOTAL.BAC + SITE + as.factor(CONTRA) + C1 + AGE,
  data = N580
)

m1 <- lm(
  A2ML1 ~ CST + TOTAL.BAC + inflammation_score +
    SITE + as.factor(CONTRA) + C1 + AGE,
  data = N580
)

summary(m0)$coefficients
summary(m1)$coefficients

confint(m0)
confint(m1)

# Percent attenuation after adjustment
# CST-D
cst_beta0 <- coef(m0)["CSTD"]
cst_beta1 <- coef(m1)["CSTD"]

cst_attenuation <- 100 * (
  abs(cst_beta0) - abs(cst_beta1)
) / abs(cst_beta0)

# TOTAL.BAC
bac_beta0 <- coef(m0)["TOTAL.BAC"]
bac_beta1 <- coef(m1)["TOTAL.BAC"]

bac_attenuation <- 100 * (
  abs(bac_beta0) - abs(bac_beta1)
) / abs(bac_beta0)

cst_attenuation
bac_attenuation

### A2ML1 ~ inflammation * TOTAL.BAC #####
totalbac_nb <- glm.nb(
  inflammation_score ~ A2ML1 * TOTAL.BAC.N +
    C1 + CONTRA + SITE + AGE,
  data = N580
)

broom::tidy(
  totalbac_nb,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::filter(
    term %in% c(
      "A2ML1",
      "TOTAL.BAC.N",
      "A2ML1:TOTAL.BAC.N"
    )
  )

### exploratory breakpoint analysis #####
# fit at a fixed threshold
fit_threshold_model <- function(data, threshold) {
  
  data <- data %>%
    mutate(
      high_bacterial = TOTAL.BAC > threshold
    )
  
  fit <- glm.nb(
    inflammation_score ~
      A2ML1 * high_bacterial +
      C1 + CONTRA + SITE + AGE,
    data = data
  )
  
  list(
    threshold = threshold,
    fit = fit,
    logLik = as.numeric(logLik(fit))
  )
}

# search for best threshold
find_best_threshold <- function(
    data,
    lower_quantile = 0.15,
    upper_quantile = 0.85,
    step = 0.01
) {
  
  candidate_thresholds <- quantile(
    data$TOTAL.BAC,
    probs = seq(
      lower_quantile,
      upper_quantile,
      by = step
    ),
    na.rm = TRUE
  )
  
  candidate_thresholds <- unique(
    as.numeric(candidate_thresholds)
  )
  
  results <- lapply(
    candidate_thresholds,
    function(th)
      fit_threshold_model(
        data,
        threshold = th
      )
  )
  
  logLik_values <- sapply(
    results,
    function(x) x$logLik
  )
  
  best_index <- which.max(logLik_values)
  
  list(
    threshold =
      results[[best_index]]$threshold,
    
    fit =
      results[[best_index]]$fit,
    
    logLik =
      results[[best_index]]$logLik,
    
    thresholds =
      candidate_thresholds,
    
    logLik_profile =
      logLik_values
  )
}

# observed change point statistic
change_point_stat <- function(data) {
  
  null_fit <- glm.nb(
    inflammation_score ~
      A2ML1 + C1 + CONTRA + SITE + AGE,
    data = data
  )
  
  best_fit <- find_best_threshold(data)
  
  LR <- 2 * (
    best_fit$logLik -
      as.numeric(logLik(null_fit))
  )
  
  list(
    LR = LR,
    threshold = best_fit$threshold,
    best_fit = best_fit,
    null_fit = null_fit
  )
}

# simulate under H0
simulate_null_nb <- function(data, null_fit) {
  
  sim_data <- data
  
  sim_data$inflammation_score <- rnbinom(
    n = nrow(data),
    mu = fitted(null_fit),
    size = null_fit$theta
  )
  
  sim_data
}

# bootstrap entire threshold search
supLR_bootstrap <- function(
    data,
    B = 1000,
    seed = 123
) {
  
  set.seed(seed)
  
  observed <- change_point_stat(data)
  
  observed_LR <- observed$LR
  
  observed_threshold <- observed$threshold
  
  null_fit <- observed$null_fit
  
  LR_boot <- rep(NA, B)
  
  threshold_boot <- rep(NA, B)
  
  for (b in seq_len(B)) {
    
    sim_data <- simulate_null_nb(
      data,
      null_fit
    )
    
    boot_result <- tryCatch(
      change_point_stat(sim_data),
      error = function(e) NULL
    )
    
    if (!is.null(boot_result)) {
      
      LR_boot[b] <- boot_result$LR
      
      threshold_boot[b] <- boot_result$threshold
    }
  }
  
  LR_boot <- na.omit(LR_boot)
  
  threshold_boot <- na.omit(threshold_boot)
  
  p_value <- mean(
    LR_boot >= observed_LR
  )
  
  threshold_ci <- quantile(
    threshold_boot,
    c(0.025, 0.975)
  )
  
  list(
    observed_threshold =
      observed_threshold,
    
    observed_LR =
      observed_LR,
    
    p_value =
      p_value,
    
    threshold_ci =
      threshold_ci,
    
    threshold_boot =
      threshold_boot,
    
    LR_boot =
      LR_boot,
    
    observed =
      observed
  )
}

## fit and add features
# ~30 minutes
cp_test <- supLR_bootstrap(
  N580,
  B = 1000
)

discovered_threshold <-
  cp_test$observed_threshold

## reporting
cat(
  "Threshold:",
  round(cp_test$observed_threshold, 3),
  "\n"
)

cat(
  "95% CI:",
  round(cp_test$threshold_ci[1], 3),
  "-",
  round(cp_test$threshold_ci[2], 3),
  "\n"
)

cat(
  "Bootstrap p-value:",
  signif(cp_test$p_value, 3),
  "\n"
)
