# Murray and Larry (2005)
# Variables
population_size <- as.numeric(readline('Size of the original population: ')) # Size of the population
z_score <- 1.96 # Z score for C.I = 0.95 (95%)
population_sd <- 0.5 # Standard Deviation of the population. Use 50% if not sure
margin_of_error <- 0.05 # Margin of error of 5%

# Calculate sample size
sample_size <- (z_score^2 * population_sd^2 * population_size) / (((margin_of_error^2) * (population_size - 1)) + (z_score^2 * population_sd^2))

# Print the result
cat('The sample size is ', sample_size)