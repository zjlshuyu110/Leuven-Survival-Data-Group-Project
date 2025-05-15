#----------------------------
# 1. Load the Data
#----------------------------
data <- read.csv("breast_cancer_survival.csv", header = TRUE, stringsAsFactors = FALSE)

#----------------------------
# 2. Convert Date Variables
#----------------------------
data$Date_of_Surgery    <- as.Date(data$Date_of_Surgery, format = "%d-%b-%y")
data$Date_of_Last_Visit <- as.Date(data$Date_of_Last_Visit, format = "%d-%b-%y")

#----------------------------
# 3. Basic Descriptive Statistics
#----------------------------
# Summary of all variables
summary(data)
sapply(data, function(x) sum(is.na(x)))

continuous_vars <- c("Age", "Protein1", "Protein2", "Protein3", "Protein4")
summary(data[ , continuous_vars])
library(ggplot2)
ggplot(data, aes(x = Age)) +
  geom_histogram(bins = 10, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Age", x = "Age", y = "Frequency")

# Histogram of Protein1 (repeat similarly for Protein2, Protein3, Protein4 as needed)
ggplot(data, aes(x = Protein1)) +
  geom_histogram(bins = 15, fill = "lightgreen", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Protein1", x = "Protein1", y = "Count")


table(data$Gender)
table(data$Tumour_Stage)
table(data$Histology)
# For variables like "ER status", use:
table(data$`ER status`)
table(data$`PR status`)
table(data$`HER2 status`)
table(data$Surgery_type)
# 5. Create Follow-Up Time and Handle Censoring
#----------------------------
# Compute follow-up time in days (if Date_of_Last_Visit is available)

data$Followup_Days <- as.numeric(data$Date_of_Last_Visit - data$Date_of_Surgery)
summary(data$Followup_Days)

# "Dead" signifies the event (1) and "Alive" is censored (0).
data$Event <- ifelse(data$Patient_Status == "Dead", 1, 0)
table(data$Patient_Status)
table(data$Event)


# Part B:
library(survival)
library(survminer)

# Fit Kaplan-Meier survival curves by Tumour_Stage.
km_fit <- survfit(Surv(Followup_Days, Event) ~ Tumour_Stage, data = data)

# Plot the survival curves with confidence intervals and a risk table.
ggsurvplot(km_fit,
           data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "Tumour Stage",
           legend.labs = levels(as.factor(data$Tumour_Stage)),
           xlab = "Follow-up Time (days)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curves by Tumour Stage")
#(b)
# Get the unique levels for Tumour_Stage.
stages <- levels(as.factor(data$Tumour_Stage))

# Loop over each stage to compute the quartiles.
for(stage in stages){
  cat("\n---\nTumour Stage:", stage, "\n")
  
  # Fit a Kaplan-Meier curve for the subset for this stage.
  km_stage <- survfit(Surv(Followup_Days, Event) ~ 1, data = data[data$Tumour_Stage == stage, ])
  
  # Compute the quantiles at the 25th, 50th, and 75th percentiles.
  q <- quantile(km_stage, probs = c(0.25, 0.5, 0.75))
  print(q)
}

#c:
# Perform the log-rank test
logrank_test <- survdiff(Surv(Followup_Days, Event) ~ Tumour_Stage, data = data)
print(logrank_test)

#d# Fit a Cox proportional hazards model
cox_fit <- coxph(Surv(Followup_Days, Event) ~ Tumour_Stage, data = data)
summary(cox_fit)

# Estimate the survival function from the Cox model.
cox_surv <- survfit(cox_fit)

# Create a new data frame covering all levels of Tumour_Stage
newdata <- data.frame(Tumour_Stage = levels(as.factor(data$Tumour_Stage)))

cox_surv <- survfit(cox_fit, newdata = newdata)

ggsurvplot(cox_surv,
           data = newdata,
           conf.int = TRUE,
           legend.title = "Tumour Stage",
           legend.labs = newdata$Tumour_Stage,
           xlab = "Follow-up Time (days)",
           ylab = "Survival Probability",
           title = "Cox Model Estimated Survival Curves by Tumour Stage")

# P.C
# Convert appropriate variables to factors
data$Tumour_Stage  <- as.factor(data$Tumour_Stage)
data$Histology     <- as.factor(data$Histology)
# Remove they are constant:
# data$ER.status     <- as.factor(data$ER.status)
# data$PR.status     <- as.factor(data$PR.status)
data$HER2.status   <- as.factor(data$HER2.status)
data$Surgery_type  <- as.factor(data$Surgery_type)
data$Gender        <- as.factor(data$Gender)


library(survival)

# Build the Cox proportional hazards model including all covariates.

cox_mod <- coxph(Surv(Followup_Days, Event) ~ 
                   Age + Protein1 + Protein2 + Protein3 + Protein4 +
                   Tumour_Stage + Histology + HER2.status + Surgery_type,
                 data = data)


summary(cox_mod)


cox_mod_reduced <- update(cox_mod, . ~ . - Tumour_Stage)
anova_test <- anova(cox_mod_reduced, cox_mod, test = "LRT")
print(anova_test)


cox_ph_test <- cox.zph(cox_mod)
print(cox_ph_test)

plot(cox_ph_test)

censoring_rate <- sum(data$Event == 0, na.rm = TRUE) / nrow(data)
print(paste("Censoring Rate:", round(censoring_rate * 100, 2), "%"))



#P.c.(B)
library(survival)

data <- data[!is.na(data$Followup_Days), ]

data <- data[data$Followup_Days > 0, ]


# Fit an AFT model with a Weibull distribution. (ER.status and PR.status are omitted.)
aft_model <- survreg(Surv(Followup_Days, Event) ~ 
                       Age + Protein1 + Protein2 + Protein3 + Protein4 +
                       Tumour_Stage + Histology + HER2.status + Surgery_type,
                     data = data, dist = "weibull")


summary(aft_model)



summary(aft_model)  # Get point estimates and standard errors
confint(aft_model)  # Get confidence intervals

exp(coef(aft_model))  # Exponentiate to get time ratios

exp(confint(aft_model))  # Confidence intervals for time ratios

summary(cox_mod)  # Cox model results
summary(aft_model)  # AFT model results


library(survminer)
km_fit <- survfit(Surv(Followup_Days, Event) ~ Tumour_Stage, data = data)

ggsurvplot(km_fit,
           data = data,
           conf.int = TRUE,
           pval = TRUE,
           risk.table = TRUE,
           legend.title = "Tumour Stage",
           xlab = "Follow-up Time (days)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curves by Tumour Stage")





