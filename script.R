required_pkgs <- c("tidyverse", "survival", "survminer", "broom", "flexsurv")
installed <- rownames(installed.packages())
for (p in required_pkgs) if (!(p %in% installed)) install.packages(p)
library(tidyverse)
library(survival)
library(survminer)
library(broom)      # for tidy model output
library(flexsurv)

raw   <- read_csv("/home/cyberpvnk/Downloads/UTF-8breast_cancer_survival(1).csv",
                  col_types = cols())

# Convert dates & derive survival object
bc <- raw %>%
  mutate(date_surg = as.Date(Date_of_Surgery,  format = "%d-%b-%y"),
         date_last = as.Date(Date_of_Last_Visit, format = "%d-%b-%y"),
         time      = as.numeric(date_last - date_surg),      # in days
         status    = if_else(Patient_Status == "Dead", 1, 0),
         # Factor versions of key categoricals
         Tumour_Stage = factor(Tumour_Stage),
         Gender       = factor(Gender),
         Histology    = factor(Histology),
         ER_status    = factor(`ER status`),
         PR_status    = factor(`PR status`),
         HER2_status  = factor(`HER2 status`),
         Surgery_type = factor(Surgery_type)) %>%
  # Remove impossible or missing times
  filter(!is.na(time) & time >= 0)

surv_obj <- with(bc, Surv(time, status))


num_vars <- bc %>% select(where(is.numeric))
summary(num_vars)


num_vars %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ name, scales = "free_x") +
  theme_minimal()


fit_all <- survfit(surv_obj ~ 1, data = bc)



ggsurvplot(fit_all, conf.int = TRUE, risk.table = TRUE,
           title = "Overall survival – all stages")




# -----------------------------
# PART B – Focus on Tumour_Stage
# -----------------------------

# (a) Kaplan–Meier curves per stage
fit_stage <- survfit(surv_obj ~ Tumour_Stage, data = bc)
ggsurvplot(fit_stage, conf.int = FALSE, risk.table = TRUE,
           xlab = "Days since surgery", ylab = "Survival probability",
           legend.title = "Stage", legend.labs = levels(bc$Tumour_Stage),
           title = "Kaplan–Meier survival by tumour stage")

# (b) Quartiles & CIs per stage
quartiles <- map_df(levels(bc$Tumour_Stage), function(s) {
  tmp <- survfit(surv_obj ~ 1, data = bc %>% filter(Tumour_Stage == s))
  tibble(Stage = s,
         q25 = quantile(tmp, probs = 0.25)[[1]],
         q50 = quantile(tmp, probs = 0.50)[[1]],
         q75 = quantile(tmp, probs = 0.75)[[1]])
})
print(quartiles)


logrank <- survdiff(surv_obj ~ Tumour_Stage, data = bc)
print(logrank)



cox_stage <- coxph(surv_obj ~ strata(Tumour_Stage), data = bc)
base <- survfit(cox_stage)
ggsurvplot(base, data=bc, legend.title = "Stage (Cox)",
           legend.labs = levels(bc$Tumour_Stage),
           title = "Cox‑based survival estimates by stage")


-----------------------------
  # PART C – Multivariable modelling
  # -----------------------------

# (a) Cox proportional‑hazards model (stage forced in)
factor_vars <- c("Gender", "Histology", "ER_status", "PR_status", "HER2_status", "Surgery_type")
valid_factor_vars <- factor_vars[sapply(bc[factor_vars], \(x) nlevels(x) > 1)]
if (length(setdiff(factor_vars, valid_factor_vars)) > 0) {
  message("Removed variables with only one level: ", paste(setdiff(factor_vars, valid_factor_vars), collapse = ", "))
}

full_formula <- as.formula(paste(
  "surv_obj ~ Tumour_Stage + Age +",
  paste(valid_factor_vars, collapse = " + "),
  "+ Protein1 + Protein2 + Protein3 + Protein4"))

full_cox <- coxph(full_formula, data = bc)
step_cox <- MASS::stepAIC(full_cox, k = 2, trace = TRUE)

print(summary(step_cox))
print(broom::tidy(step_cox, exponentiate = TRUE, conf.int = TRUE))
print(cox.zph(step_cox, transform = "km"))   

# Hazard ratios & 95 % CIs
broom::tidy(step_cox, exponentiate = TRUE, conf.int = TRUE)

# Significance of Tumour_Stage (partial‑likelihood ratio test)
cox.zph(step_cox, transform = "km")   # also tests PH assumption

#to-do: parametric AFT model
