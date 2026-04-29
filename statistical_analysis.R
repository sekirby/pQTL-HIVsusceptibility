##### data & library loading #####
packages <- c("dplyr", "tidyr", "stats", "data.table", "glue", "ggplot2", "purrr", "patchwork", "mediation")

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
data$SITE <- as.factor(data$SITE)
data$TOTAL.BAC.N <- scale(as.numeric(data$TOTAL.BAC))

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

##### HIV acqusition odds #####
vars <- c("rs2920282","rs8106317","rs147944114","rs143864957")
prots <- c("A2ML1","PSCA","GSTO1","CD177")
cats <- c(vars, prots)

all_HIV_results <- data.frame()

for (i in cats) {
  # Fit
  HIV_model <- glm(HIV ~ N589[[i]] + C1 + AGE + as.factor(CONTRA) + CST + SITE,
                   family = binomial, data = N589)
  
  # Extract
  summary_model <- summary(HIV_model)$coefficients
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
  HIV_result <- data.frame(
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
  all_HIV_results <- rbind(all_HIV_results, HIV_result)
}

print("Predictors of HIV acquisition - binomial logistic regressions")
print(all_HIV_results, row.names = FALSE)

##### predicted probability #####
N589$A2ML1_col <- N589[["A2ML1"]]

# Fit the model
HIV_model_A2ML1 <- glm(HIV ~ A2ML1_col + C1 + AGE + as.factor(CONTRA) + CST + SITE,
                       family = binomial, data = N589)

# Create prediction data
pred_data <- data.frame(
  A2ML1_col = seq(from = min(N589$A2ML1_col, na.rm = TRUE),
                  to = max(N589$A2ML1_col, na.rm = TRUE),
                  length.out = 100),
  C1 = mean(N589$C1, na.rm = TRUE),
  AGE = mean(N589$AGE, na.rm = TRUE),
  CONTRA = names(which.max(table(N589$CONTRA))),
  CST = names(which.max(table(N589$CST))),
  SITE = names(which.max(table(N589$SITE)))
)

# Generate predictions
predictions <- predict(HIV_model_A2ML1, 
                       newdata = pred_data, 
                       type = "link",
                       se.fit = TRUE)

# calculate CI
pred_data$probability <- plogis(predictions$fit)
pred_data$lower <- plogis(predictions$fit - 1.96 * predictions$se.fit)
pred_data$upper <- plogis(predictions$fit + 1.96 * predictions$se.fit)


##### inflammation odds #####
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
  
  # confidence intervals and odds ratio
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
               values_to = "Measurement") %>%
  mutate(
    Cytokine_clean = Cytokine %>%
      gsub("_", "-", .) %>%             # Replace underscore with dash
      gsub("A$", "\u03B1", .) %>%       # Replace terminal A with Greek alpha
      gsub("B$", "\u03B2", .)           # Replace terminal B with Greek beta
  )

# Calculate Q1, Median, Q3
quartiles <- data_long %>%
  group_by(Cytokine_clean) %>%
  summarise(
    Q1 = quantile(Measurement, 0.25, na.rm = TRUE),
    Median = median(Measurement, na.rm = TRUE),
    Q3 = quantile(Measurement, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

print(quartiles)

##### A2ML1 bacterial load/inflammation "functional relevance" #####
functional_relevance <- N580 %>%
  filter(!is.na(CST)) %>%
  split(.$CST) %>%
  map_dfr(~{
    # Skip if A2ML1 is NA or has no variation
    if (all(is.na(.x$A2ML1)) || length(unique(na.omit(.x$A2ML1))) < 2) {
      return(NULL)
    }
    
    tryCatch({
      # Logi
      inflammation_model <- glm(inflammation_binary ~ A2ML1 + C1 + AGE + CONTRA + SITE, family = binomial, data = .x)
      inflammation_summary <- summary(inflammation_model)$coefficients
      inflammation_beta <- inflammation_summary[2, 1]
      inflammation_se <- inflammation_summary[2, 2]
      inflammation_p <- inflammation_summary[2, 4]
      inflammation_ci_low <- inflammation_beta - 1.96 * inflammation_se
      inflammation_ci_high <- inflammation_beta + 1.96 * inflammation_se
      
      # Lin
      bacterial_model <- lm(TOTAL.BAC.N ~ A2ML1 + C1 + AGE + CONTRA + SITE, data = .x)
      bacterial_summary <- summary(bacterial_model)$coefficients
      bacterial_beta <- bacterial_summary[2, 1]
      bacterial_se <- bacterial_summary[2, 2]
      bacterial_p <- bacterial_summary[2, 4]
      bacterial_ci_low <- bacterial_beta - 1.96 * bacterial_se
      bacterial_ci_high <- bacterial_beta + 1.96 * bacterial_se
      
      tibble(
        CST = unique(.x$CST),
        
        # Inflammation model
        inflammation_beta = inflammation_beta,
        inflammation_se = inflammation_se,
        inflammation_CI_low = inflammation_ci_low,
        inflammation_CI_high = inflammation_ci_high,
        inflammation_p = inflammation_p,
        
        # Bacterial model
        bacterial_beta = bacterial_beta,
        bacterial_se = bacterial_se,
        bacterial_CI_low = bacterial_ci_low,
        bacterial_CI_high = bacterial_ci_high,
        bacterial_p = bacterial_p,
        
        n = nrow(.x)
      )
    }, error = function(e) {
      message("Skipping CST = ", unique(.x$CST), " due to error: ", e$message)
      NULL
    })
  })

print(functional_relevance)

# Reshape to long
plot_data <- functional_relevance %>%
  dplyr::select(
    CST,
    inflammation_beta, inflammation_CI_low, inflammation_CI_high,
    bacterial_beta, bacterial_CI_low, bacterial_CI_high
  ) %>%
  pivot_longer(
    cols = -CST,
    names_to = c("Outcome", ".value"),
    names_pattern = "(.*)_(beta|CI_low|CI_high)"
  ) %>%
  mutate(
    Outcome = case_when(
      Outcome == "inflammation" ~ "Inflammation (logistic)",
      Outcome == "bacterial" ~ "Bacterial Load (linear)",
      TRUE ~ Outcome
    ),
    CST = factor(CST,
                 levels = c("A", "B", "C", "D"),
                 labels = c("CST-LC", "CST-LI", "CST-GV", "CST-PM")
    )
  )

##### A2ML1 mediation by CST ######
set.seed(1234)
full_snp_levels <- levels(factor(N580$rs143864957))
full_contra_levels <- levels(factor(N580$CONTRA))

# inflammation outcome
full_snp_levels <- levels(factor(N580$rs143864957))
full_contra_levels <- levels(factor(N580$CONTRA))

boot_fallback <- list()

mediation_inflam <- N580 %>%
  split(.$CST) %>%
  imap(~{
    relevant_vars <- c("rs143864957", "A2ML1", "C1", "inflammation_binary", "AGE", "CONTRA", "SITE")
    
    # Filter and drop unused levels
    .x <- .x %>%
      filter(complete.cases(across(all_of(relevant_vars)))) %>%
      mutate(
        rs143864957 = droplevels(factor(rs143864957)),
        CONTRA = droplevels(factor(CONTRA))
      )
    
    # Check levels
    if (length(levels(.x$rs143864957)) < 2) {
      message(glue("CST {.y}: rs143864957 has less than 2 levels after filtering — skipping group"))
      return(NULL)
    }
    
    if (length(levels(.x$CONTRA)) < 2) {
      message(glue("CST {.y}: CONTRA has less than 2 levels after filtering — skipping group"))
      return(NULL)
    }
    
    if (nrow(.x) < 50) {
      message(glue("CST {.y}: fewer than 50 observations after filtering — skipping group"))
      return(NULL)
    }
    
    result <- list()
    
    # Reduced model with bootstrap
    result$reduced_bootstrap <- tryCatch({
      mediate(
        model.m = lm(A2ML1 ~ rs143864957 + C1 + AGE, data = .x),
        model.y = glm(inflammation_binary ~ rs143864957 + A2ML1 + C1 + AGE,
                      data = .x, family = binomial),
        treat = "rs143864957",
        mediator = "A2ML1",
        boot = TRUE,
        sims = 2000
      )
    }, error = function(e) {
      message(glue("Reduced bootstrap model failed for group {.y}: {e$message}"))
      boot_fallback[[.y]] <<- "Reduced model bootstrap failed"
      return(NULL)
    })
    
    # Full model without bootstrap
    result$full_noboot <- tryCatch({
      mediate(
        model.m = lm(A2ML1 ~ rs143864957 + C1 + AGE + CONTRA + SITE, data = .x),
        model.y = glm(inflammation_binary ~ rs143864957 + A2ML1 + C1 + AGE + CONTRA + SITE,
                      data = .x, family = binomial),
        treat = "rs143864957",
        mediator = "A2ML1",
        boot = FALSE
      )
    }, error = function(e) {
      message(glue("Full model (no boot) failed for group {.y}: {e$message}"))
      boot_fallback[[.y]] <<- "Full model no boot failed"
      return(NULL)
    })
    
    return(result)
  })

# Format results
mediation_inflam_results <- imap_dfr(
  mediation_inflam,
  ~{
    reduced <- .x$reduced_bootstrap
    full <- .x$full_noboot
    
    reduced_summary <- if (!is.null(reduced)) summary(reduced) else NULL
    full_summary <- if (!is.null(full)) summary(full) else NULL
    
    tibble(
      CST = .y,
      
      # Reduced model
      reduced_acme_estimate = if (!is.null(reduced)) reduced$d0 else NA_real_,
      reduced_acme_ci_lower = if (!is.null(reduced_summary)) reduced_summary$d0.ci[1] else NA_real_,
      reduced_acme_ci_upper = if (!is.null(reduced_summary)) reduced_summary$d0.ci[2] else NA_real_,
      reduced_acme_p = if (!is.null(reduced)) reduced$d0.p else NA_real_,
      
      reduced_ade_estimate = if (!is.null(reduced)) reduced$z0 else NA_real_,
      reduced_ade_ci_lower = if (!is.null(reduced_summary)) reduced_summary$z0.ci[1] else NA_real_,
      reduced_ade_ci_upper = if (!is.null(reduced_summary)) reduced_summary$z0.ci[2] else NA_real_,
      reduced_ade_p = if (!is.null(reduced)) reduced$z0.p else NA_real_,
      
      reduced_total_effect_estimate = if (!is.null(reduced)) reduced$tau.coef else NA_real_,
      reduced_total_effect_ci_lower = if (!is.null(reduced_summary)) reduced_summary$tau.ci[1] else NA_real_,
      reduced_total_effect_ci_upper = if (!is.null(reduced_summary)) reduced_summary$tau.ci[2] else NA_real_,
      reduced_total_effect_p = if (!is.null(reduced)) reduced$tau.p else NA_real_,
      
      reduced_prop_med_estimate = if (!is.null(reduced)) reduced$n0 else NA_real_,
      reduced_prop_med_ci_lower = if (!is.null(reduced_summary)) reduced_summary$n0.ci[1] else NA_real_,
      reduced_prop_med_ci_upper = if (!is.null(reduced_summary)) reduced_summary$n0.ci[2] else NA_real_,
      reduced_prop_med_p = if (!is.null(reduced)) reduced$n0.p else NA_real_,
      
      # Full model
      full_acme_estimate = if (!is.null(full)) full$d0 else NA_real_,
      full_acme_p = if (!is.null(full)) full$d0.p else NA_real_,
      
      full_ade_estimate = if (!is.null(full)) full$z0 else NA_real_,
      full_ade_p = if (!is.null(full)) full$z0.p else NA_real_,
      
      full_total_effect_estimate = if (!is.null(full)) full$tau.coef else NA_real_,
      full_total_effect_p = if (!is.null(full)) full$tau.p else NA_real_,
      
      full_prop_med_estimate = if (!is.null(full)) full$n0 else NA_real_,
      full_prop_med_p = if (!is.null(full)) full$n0.p else NA_real_
    )
  }
)

print(mediation_inflam_results)

##### A2ML1-inflammation-bacterial load interaction by CST #####
packages <- c("MASS", "sandwich", "lmtest", "AER", "broom")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Fit
tolerance_nb <- glm.nb(inflammation_score ~ A2ML1 * TOTAL.BAC.N + C1 + CONTRA + SITE + AGE, data = N580)

tolerance_results <- broom::tidy(tolerance_nb, conf.int = TRUE)

tolerance_results <- broom::tidy(tolerance_nb, conf.int = TRUE) %>%
  mutate(
    IRR = exp(estimate),
    IRR_conf.low = exp(conf.low),
    IRR_conf.high = exp(conf.high)
  ) %>%
  mutate(
    IRR = round(IRR, 5),
    IRR_conf.low = round(IRR_conf.low, 6),
    IRR_conf.high = round(IRR_conf.high, 6)
  ) %>%
  dplyr::select(
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    IRR,
    IRR_conf.low,
    IRR_conf.high,
    p.value
  )

print(tolerance_results)

##### threshold analysis #####
two_phase_model_negbin <- function(data, initial_threshold = NULL, plot_results = TRUE) {
  if (is.null(initial_threshold)) {
    initial_threshold <- median(data$TOTAL.BAC, na.rm = TRUE)
  }
  
  # Construct design matrix
  design_mat <- model.matrix(~ C1 + CONTRA + SITE + AGE, data = data)[, -1]  # drop intercept
  
  # Number of covariates
  p <- ncol(design_mat)
  
  # Objective function
  two_phase_fn <- function(params, data, design_mat) {
    threshold <- params[1]
    
    alpha1 <- params[2]
    beta1  <- params[3]
    alpha2 <- params[4]
    beta2  <- params[5]
    
    covariate_coefs <- params[6:(5 + p)]
    
    below_threshold <- data$TOTAL.BAC <= threshold
    
    eta_main <- ifelse(below_threshold,
                       alpha1 + beta1 * data$A2ML1,
                       alpha2 + beta2 * data$A2ML1)
    
    eta_covs <- as.numeric(design_mat %*% covariate_coefs)  # shared covariate effects
    
    eta <- eta_main + eta_covs
    
    mu <- exp(eta)
    y <- data$inflammation_score
    
    # Estimate theta from NB model
    theta <- glm.nb(inflammation_score ~ A2ML1 + C1 + CONTRA + SITE + AGE, data = data)$theta
    
    dev_residuals <- 2 * (lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
                            theta * log(theta) + y * log(mu) -
                            (y + theta) * log(mu + theta))
    
    return(-sum(dev_residuals, na.rm = TRUE))
  }
  
  # Initial parameter guesses
  initial_params <- c(
    initial_threshold,
    0,  # alpha1
    0,  # beta1
    0,  # alpha2
    0,  # beta2
    rep(0, p)  # covariate coefficients
  )
  
  # optimization
  fit <- optim(
    par = initial_params,
    fn = two_phase_fn,
    data = data,
    design_mat = design_mat,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Extract fitted parameters
  fitted_params <- fit$par
  names(fitted_params) <- c("threshold", "alpha1", "beta1", "alpha2", "beta2", colnames(design_mat))
  
  # plotting
  if (plot_results) {
    threshold <- fitted_params["threshold"]
    below_threshold <- data$TOTAL.BAC <= threshold
    
    eta_main <- ifelse(below_threshold,
                       fitted_params["alpha1"] + fitted_params["beta1"] * data$A2ML1,
                       fitted_params["alpha2"] + fitted_params["beta2"] * data$A2ML1)
    
    eta_covs <- as.numeric(design_mat %*% fitted_params[6:length(fitted_params)])
    
    eta <- eta_main + eta_covs
    predicted <- exp(eta)
    
    plot(data$A2ML1, data$inflammation_score,
         main = "Two-Phase Negative Binomial Fit",
         xlab = "A2ML1", ylab = "Inflammation Score", pch = 19, col = "grey")
    points(data$A2ML1, predicted, col = "blue", pch = 16)
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
  }
  
  return(list(
    fitted_params = fitted_params,
    logLik = -fit$value,
    convergence = fit$convergence
  ))
}

## Fit Model and Add Features
two_phase_result <- two_phase_model_negbin(N580, plot_results = F)
print(paste("Two-Phase Threshold:", round(two_phase_result$fitted_params["threshold"], 3)))


discovered_threshold <- two_phase_result$fitted_params["threshold"]

data <- N580 %>%
  mutate(
    high_bacterial = TOTAL.BAC > discovered_threshold,
    bacterial_above_threshold = ifelse(high_bacterial, TOTAL.BAC - discovered_threshold, 0)
  )


## Likelihood Ratio Test
below_threshold <- filter(data, TOTAL.BAC <= discovered_threshold)
above_threshold <- filter(data, TOTAL.BAC > discovered_threshold)

model_below_nb <- glm.nb(inflammation_score ~ A2ML1 + C1 + CONTRA + SITE + AGE, data = below_threshold)
model_above_nb <- glm.nb(inflammation_score ~ A2ML1 + C1 + CONTRA + SITE + AGE, data = above_threshold)
model_pooled_nb <- glm.nb(inflammation_score ~ A2ML1 + C1 + CONTRA + SITE + AGE, data = data)

lr_stat <- -2 * (logLik(model_pooled_nb) - (logLik(model_below_nb) + logLik(model_above_nb)))
lr_p <- pchisq(as.numeric(lr_stat), df = 2, lower.tail = FALSE)

print(paste("Likelihood Ratio test statistic:", round(lr_stat, 3)))
print(paste("LR test p-value:", round(lr_p, 4)))

## Interaction Model and Confidence Intervals
interaction_model_nb <- glm.nb(inflammation_score ~ A2ML1 * high_bacterial + C1 + CONTRA + SITE + AGE, data = data)
summary(interaction_model_nb)

# Below threshold
beta_below <- coef(interaction_model_nb)["A2ML1"]
se_below <- summary(interaction_model_nb)$coefficients["A2ML1", "Std. Error"]
ci_below <- beta_below + c(-1.96, 1.96) * se_below

# Above threshold
beta_above <- beta_below + coef(interaction_model_nb)["A2ML1:high_bacterialTRUE"]

vcov_matrix <- vcov(interaction_model_nb)
se_above <- sqrt(vcov_matrix["A2ML1", "A2ML1"] +
                   vcov_matrix["A2ML1:high_bacterialTRUE", "A2ML1:high_bacterialTRUE"] +
                   2 * vcov_matrix["A2ML1", "A2ML1:high_bacterialTRUE"])

ci_above <- beta_above + c(-1.96, 1.96) * se_above

## Bootstrap Threshold CI
set.seed(123)

bootstrap_threshold <- function(data, n_iter = 1000) {
  bootstrap_thresholds <- numeric(n_iter)
  
  # Initialize progress tracking
  start_time <- Sys.time()
  cat(sprintf("Starting bootstrap with %d iterations...\n", n_iter))
  cat("Progress: [", rep(" ", 50), "]   0%", sep="")
  
  for (i in 1:n_iter) {
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    result <- tryCatch(two_phase_model_negbin(boot_data, plot_results = FALSE), 
                       error = function(e) list(threshold = NA))
    bootstrap_thresholds[i] <- result$fitted_params["threshold"]
    
    update_freq <- max(ceiling(n_iter * 0.02), 10)
    if (i %% update_freq == 0 || i == n_iter) {
      progress <- i / n_iter
      filled_bars <- round(50 * progress)
      
      # time estimates
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (i < n_iter) {
        est_total_time <- elapsed_time / progress
        time_remaining <- est_total_time - elapsed_time
        time_remaining_str <- if (time_remaining > 60) {
          sprintf("ETA: %.1f min", time_remaining / 60)
        } else {
          sprintf("ETA: %.0f sec", time_remaining)
        }
      } else {
        time_remaining_str <- sprintf("Done in %.1f sec", elapsed_time)
      }
      
      progress_bar <- paste0(c(rep("=", filled_bars), rep(" ", 50 - filled_bars)), collapse = "")
      
      cat("\rProgress: [", progress_bar, "] ", sprintf("%3.0f%% ", progress * 100), 
          time_remaining_str, sep="")

      flush.console()
    }
  }
  
  cat("\n") 
  successful_boots <- sum(!is.na(bootstrap_thresholds))
  cat(sprintf("Bootstrap completed: %d/%d iterations successful (%.1f%%)\n", 
              successful_boots, n_iter, 100 * successful_boots / n_iter))
  
  return(na.omit(bootstrap_thresholds))
}

bootstrap_thresholds <- bootstrap_threshold(data, n_iter = 1000)
n_iter <- 1000
threshold_ci <- quantile(bootstrap_thresholds, c(0.025, 0.975))

# Extract summary and variance-covariance matrix
summary_int <- summary(interaction_model_nb)
vcov_matrix <- vcov(interaction_model_nb)

# Below threshold: coefficient, SE, CI, p-value
beta_below <- coef(interaction_model_nb)["A2ML1"]
se_below <- summary_int$coefficients["A2ML1", "Std. Error"]
ci_below <- beta_below + c(-1.96, 1.96) * se_below
z_below <- beta_below / se_below
p_below <- 2 * pnorm(-abs(z_below))

# Above threshold: A2ML1 + interaction
interaction_coef <- coef(interaction_model_nb)["A2ML1:high_bacterialTRUE"]
beta_above <- beta_below + interaction_coef

# Standard error for linear combination
se_above <- sqrt(
  vcov_matrix["A2ML1", "A2ML1"] +
    vcov_matrix["A2ML1:high_bacterialTRUE", "A2ML1:high_bacterialTRUE"] +
    2 * vcov_matrix["A2ML1", "A2ML1:high_bacterialTRUE"]
)

ci_above <- beta_above + c(-1.96, 1.96) * se_above
z_above <- beta_above / se_above
p_above <- 2 * pnorm(-abs(z_above))

## COMPREHENSIVE SUMMARY TABLE

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("                    TWO-PHASE NEGATIVE BINOMIAL MODEL SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n\n")

# Threshold Results
cat("THRESHOLD ANALYSIS:\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(sprintf("%-35s: %8.3f\n", "Estimated Threshold", discovered_threshold))
cat(sprintf("%-35s: %8.3f to %8.3f\n", "Bootstrap 95% CI", threshold_ci[1], threshold_ci[2]))
cat(sprintf("%-35s: %8d (%.1f%% successful)\n", "Bootstrap Iterations", 
            n_iter, 100*length(bootstrap_thresholds)/n_iter))
cat(sprintf("%-35s: %8d (%.1f%%)\n", "Observations Below Threshold", 
            nrow(below_threshold), 100*nrow(below_threshold)/nrow(data)))
cat(sprintf("%-35s: %8d (%.1f%%)\n", "Observations Above Threshold", 
            nrow(above_threshold), 100*nrow(above_threshold)/nrow(data)))

cat("\nSTRUCTURAL BREAK TEST:\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(sprintf("%-35s: %8.3f\n", "Likelihood Ratio Statistic", lr_stat))
cat(sprintf("%-35s: %8.4f\n", "P-value", lr_p))
cat(sprintf("%-35s: %8s\n", "Significance", ifelse(lr_p < 0.05, "Yes***", 
                                                   ifelse(lr_p < 0.10, "Yes*", "No"))))

cat("\nA2ML1 EFFECT COEFFICIENTS:\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(sprintf("%-35s: %8.4f [%8.4f, %8.4f], p = %e\n", 
            "Below Threshold (95% CI)", beta_below, ci_below[1], ci_below[2], p_below))
cat(sprintf("%-35s: %8.4f [%8.4f, %8.4f], p = %e\n", 
            "Above Threshold (95% CI)", beta_above, ci_above[1], ci_above[2], p_above))

cat("\nINCIDENCE RATE RATIOS (IRR):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
irr_below <- exp(beta_below)
irr_ci_below <- exp(ci_below)
irr_above <- exp(beta_above)
irr_ci_above <- exp(ci_above)

cat(sprintf("%-35s: %8.4f [%8.4f, %8.4f]\n", "Below Threshold (95% CI)", 
            irr_below, irr_ci_below[1], irr_ci_below[2]))
cat(sprintf("%-35s: %8.4f [%8.4f, %8.4f]\n", "Above Threshold (95% CI)", 
            irr_above, irr_ci_above[1], irr_ci_above[2]))

cat("\nMODEL DIAGNOSTICS:\n")
cat("-" ,50, "\n")
cat(sprintf("%-35s: %8.3f\n", "Pooled Model AIC", AIC(model_pooled_nb)))
cat(sprintf("%-35s: %8.3f\n", "Below Threshold AIC", AIC(model_below_nb)))
cat(sprintf("%-35s: %8.3f\n", "Above Threshold AIC", AIC(model_above_nb)))
cat(sprintf("%-35s: %8.3f\n", "Combined Separate Models AIC", 
            AIC(model_below_nb) + AIC(model_above_nb)))
cat(sprintf("%-35s: %8.3f\n", "Pooled Model Theta", model_pooled_nb$theta))
cat(sprintf("%-35s: %8.3f\n", "Below Threshold Theta", model_below_nb$theta))
cat(sprintf("%-35s: %8.3f\n", "Above Threshold Theta", model_above_nb$theta))

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("Notes: *** p < 0.05, * p < 0.10\n")
cat("IRR > 1: Positive association; IRR < 1: Negative association\n")
cat(paste(rep("=", 80), collapse=""), "\n")


###### bacterial quartiles #####
bacterial_load_summary <- N580 %>%
  group_by(CST) %>%
  summarise(
    n                     = n(),
    mean_bacterial_load   = mean(TOTAL.BAC, na.rm = TRUE),
    median_bacterial_load = median(TOTAL.BAC, na.rm = TRUE),
    min_bacterial_load    = min(TOTAL.BAC, na.rm = TRUE),
    max_bacterial_load    = max(TOTAL.BAC, na.rm = TRUE),
    Q1_bacterial_load     = quantile(TOTAL.BAC, 0.25, na.rm = TRUE),
    Q3_bacterial_load     = quantile(TOTAL.BAC, 0.75, na.rm = TRUE)
  )

print(bacterial_load_summary)

