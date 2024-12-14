# Load necessary libraries
library(dplyr)  # For data manipulation

# Define paths for input and output files
input_path <- "data/processed/agg_dummies.csv"
output_path <- "results/population_estimate.csv"

# Load the dataset
population_data <- read.csv(input_path)

# Filter data and create variables
population_data <- population_data %>%
  filter(!is.na(age)) %>%  # Remove rows with missing age
  mutate(
    income_per_capita = income / household_size,  # Example transformation
    high_income = ifelse(income_per_capita > 50000, 1, 0)  # Example binary variable
  )

# Run the population model (example: Poisson regression)
model <- glm(count ~ age + income_per_capita + gender, 
             data = population_data, family = poisson())

# Output summary of the model
summary(model)

# Save model predictions to a file
population_data <- population_data %>%
  mutate(predicted_count = predict(model, type = "response"))
write.csv(population_data, output_path, row.names = FALSE)

# Notes:
# - Replace 'count', 'age', 'income_per_capita', and 'gender' with actual variable names.
# - Adjust the model formula based on the methodology described in the paper.
