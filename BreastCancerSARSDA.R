library(survival)
library(ggplot2)
data <- read.csv("/Users/lingjong/Downloads/breast_cancer_survival.csv")

attach(data)

# We convert the following to appropriate date types:

data$Date_of_Surgery <- as.Date(data$Date_of_Surgery, format="%d-%b-%y")
data$Date_of_Last_Visit <- as.Date(data$Date_of_Last_Visit, format="%d-%b-%y")
data$Patient_Status <- factor(data$Patient_Status, levels = c("Alive", "Dead"))



data$Survival_Time <- as.numeric(data$Date_of_Last_Visit - data$Date_of_Surgery)
