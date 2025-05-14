# Core survival-analysis tooling
library(survival)      # Surv(), survfit(), coxph(), survreg(), cox.zph(), …
library(survminer)     # gg-friendly plotting helpers
library(dplyr)         # cleaning / manipulation
library(lubridate)     # safer date handling
library(MASS)          # stepAIC() for stepwise variable selection



bc <- read.csv("/cloud/project/breast_cancer_survival.csv", stringsAsFactors = TRUE)

#Part A 

bc <- bc %>% 
  mutate(
    date_surg = dmy(Date_of_Surgery),          # lubridate::dmy()
    date_last = dmy(Date_of_Last_Visit),
    time      = as.numeric(date_last - date_surg),   # follow-up in *days*
    event     = ifelse(toupper(Patient_Status) == "DEAD", 1, 0)
  ) %>% 
  filter(!is.na(time))  # drop 17 rows with missing dates


# Continuous variables
summary(dplyr::select(bc, Age, starts_with("Protein"), time))


# Categorical variables
lapply(dplyr::select(bc, Gender, Tumour_Stage, Histology,
              `ER.status`, `PR.status`, `HER2.status`, Surgery_type), table)

# Censoring proportion
mean(bc$event == 0) 



#Part B: Kaplan Meier curves

surv_obj <- Surv(time = bc$time, event = bc$event)

km_stage <- survfit(surv_obj ~ Tumour_Stage, data = bc)

ggsurvplot(
  km_stage,
  conf.int     = TRUE,
  risk.table   = TRUE,
  censor.shape = "|",
  legend.title = "Tumour stage",
  xlab = "Days since surgery", ylab = "Survival probability"
)


q <- quantile(km_stage, probs = c(.25, .50, .75))



survdiff(surv_obj ~ Tumour_Stage, data = bc)         # default = log-rank
# For a Wilcoxon‐type early-event-focused test:
survdiff(surv_obj ~ Tumour_Stage, data = bc, rho = 1)




cox_stage <- coxph(surv_obj ~ Tumour_Stage, data = bc)

newdat <- data.frame(Tumour_Stage = levels(bc$Tumour_Stage))
cox_surv <- survfit(cox_stage, newdata = newdat)


ggsurvplot(
  cox_surv,
  data        = newdat,        
  conf.int    = TRUE,
  legend.title = "Tumour stage",
  xlab = "Days since surgery", ylab = "Adjusted survival"
)


#Part C

# Candidate covariates
bc <- bc |> dplyr::select(where(function(x)
  !(is.factor(x) && nlevels(x) == 1)))

full <- coxph(surv_obj ~ Tumour_Stage + Age + Gender +
                Protein1 + Protein2 + Protein3 + Protein4 +
                Histology  +
                `HER2.status` + Surgery_type,
              data = bc, x = TRUE)

# Keep Tumour_Stage regardless of p-value
scope <- list(lower = ~ Tumour_Stage,
              upper = formula(full)[-2])   # remove Surv(…) part

stepmod <- stepAIC(full, scope = scope, direction = "both", trace = FALSE)
summary(stepmod)


exp(cbind(HR = coef(stepmod), confint(stepmod)))



# Likelihood-ratio test by refitting without stage
no_stage <- update(stepmod, . ~ . - Tumour_Stage)
anova(no_stage, stepmod, test = "LRT")



cox.zph(stepmod)      # Schoenfeld global & per-covariate tests
ggcoxzph(cox.zph(stepmod))   # optional plots
ggsurvplot(
  survfit(stepmod),
  data = bc,             # <-- needed for ggsurvplot to find covariates
  fun  = "cloglog",      # log-minus-log transform
  xlab = "Log(days)", ylab = "log{-log S(t)}"
)



bc <- bc |> filter(time > 0 & !is.na(time))   # keep only positive times
surv_obj <- Surv(bc$time, bc$event)           # rebuild the Surv object


dists <- c("weibull", "lognormal", "loglogistic", "exponential")

fits <- lapply(dists, \(d)
               survreg(surv_obj ~ Tumour_Stage + Protein4, data = bc, dist = d)
)
AICs <- sapply(fits, AIC); names(AICs) <- dists; AICs
best <- fits[[which.min(AICs)]]
summary(best)


