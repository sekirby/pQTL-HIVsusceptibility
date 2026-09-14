packages <- c("dplyr", "readxl", "broom")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


data <- read_excel("L:/NLHG/Sam_Kirby/Sam_Kirby_Local/cap004_a2ml1_cytochange/CAP004_longitudinalcytoXa2ml1.xlsx")
N580 <- data %>% filter(CST != "NA")

cytokines <- c("IL_1B","IL_6","IP_10","MCP_1","MIP_1A","MIP_1B","IL_8","IL_1A","TNF_A")

# calculate 75th percentile thresholds from baseline visits only
baseline_thresholds <- N580 %>%
  group_by(PID) %>%
  arrange(days_til_cytokine) %>%
  slice(1) %>%
  ungroup() %>%
  summarise(across(all_of(cytokines), ~quantile(., 0.75, na.rm = TRUE)))

# apply thresholds to all visits and compute inflammation score (0-9)
binary_cols <- mapply(
  function(col, thresh) as.integer(col > thresh),
  N580[cytokines],
  baseline_thresholds[cytokines]
)
colnames(binary_cols) <- paste0(cytokines, "_high")

data_scored <- N580 %>%
  bind_cols(as.data.frame(binary_cols)) %>%
  mutate(inflammation_score = rowSums(across(all_of(paste0(cytokines, "_high"))), na.rm = TRUE))

# create one row per person with baseline, change scores, and covariates
# does not matter which covariate pulled as they are repeated across visits (ie in long format)
data_wide_inflammation <- data_scored %>%
  group_by(PID) %>%
  arrange(days_til_cytokine) %>%
  summarise(
    A2ML1_base            = first(A2ML1_base),
    days_followup         = last(days_til_cytokine) - first(days_til_cytokine),
    inflammation_baseline = first(inflammation_score),
    inflammation_followup = last(inflammation_score),
    C1                    = first(C1),
    AGE                   = first(AGE),
    SITE                  = first(SITE),
    CST                   = first(CST),
    CONTRA                = first(CONTRA)
  ) %>%
  ungroup()

# model change in inflammation score with covariates
model_inflammation <- lm(
  inflammation_followup ~
    A2ML1_base +
    inflammation_baseline +
    days_followup +
    C1 +
    AGE +
    as.factor(SITE) +
    as.factor(CST) +
    as.factor(CONTRA),
  data = data_wide_inflammation
)

# print
summary(model_inflammation)
tidy(model_inflammation)

tidy(model_inflammation, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value)