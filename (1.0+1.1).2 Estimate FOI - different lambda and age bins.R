library(rstan)
library(ggplot2)
library(tidyr)
library(dplyr)
options(mc.cores=4)
library("RColorBrewer")

# Set random seed for reproducibility
set.seed(1230)

# Number of individuals
N <- 200 

# Simulate continuous ages between 0 and 50 years (this will be reused)
age_continuous <- runif(N, min=0, max=10) 

# List of true lambda (FOI) values to simulate
lambda_values <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8)

# Function to calculate the probability of seropositivity given the exact age
# X(a) = 1 -e^(-lmabda * a)
ProbSeropos <- function(a, lambda) {
  1 - exp(-lambda * a)
}

# Store the results
continuous_est <- list() # estimates from the continuous model
binned_est <- list() # # estimates from the binned model

# List to store the simulated survey designs for each lambda
survey_designs <- list()

################################################################################
# 1. Estimate FOI with continuous age data
################################################################################
# The continuous age model
cont_age_model <- stan_model("1.0 Estimate FOI - continuous age and serostatus.stan")

# Run the continuous model for each lambda values and store the results
for (true_lambda in lambda_values) {
  Y_true_serostatus <- numeric(N) 
  prob_seropos <- numeric(N)
  
  # Simulate serostatus y based on true lambda and continuous ages
  for (i in 1:N) {
    prob_seropos[i] <- ProbSeropos(age_continuous[i], true_lambda)
    # Y(a) ~ Bernoulli(X(a))
    Y_true_serostatus[i] <- rbinom(1, 1, prob_seropos[i])
  }
  
  # Create the data frame to store the simulated survey design and sort by age
  sim_survey_design <- tibble(
    age = age_continuous,
    prob = prob_seropos,
    Y = Y_true_serostatus
  ) %>% 
    arrange(age)
  
  # Store the survey design for the current lambda in the list
  survey_designs[[paste0("lambda=", true_lambda)]] <- sim_survey_design
  
  # Simulated data for Stan (continuous age model)
  cont_sim_data <- list(
    N = N,
    age = sim_survey_design$age,
    Y = sim_survey_design$Y
  )
  
  # Fit the continuous age model
  cont_fit <- sampling(
    cont_age_model,   
    data = cont_sim_data,           
    chains = 4,                
    iter = 1000            
  )
  
  # Store the estimates of lambda (continuous age model)
  lambda_est <- rstan::extract(cont_fit, "lambda")[[1]] 
  continuous_est[[paste0("lambda=", true_lambda)]] <- lambda_est
}


################################################################################
# Function to store summary statistics and estimates of lambda for continuous models
# Use of the function:
# - This function processes the results of the estimated lambdas for each true
#   lambda value from the continuous model. It calculates the mean and standard
#   deviation of the estimated lambdas and stores both the summary statistics and
#   the individual estimates.
#
# Input:
# - continuous_est: A list where each element contains the estimated lambda
#   values for different true lambda values.
#
# Output:
# - A list with two data frames:
#     - summary: Contains the summary statistics (mean and standard deviation) 
#       of the estimated lambda values for each true lambda.
#     - estimates: Contains the individual estimated lambda values for all 
#       iterations.
################################################################################
StoreContinuousResults <- function(continuous_est) {
  
  # store summary statistics (meand and sd)
  cont_summary_df <- data.frame()
  # store estimates of lambda
  cont_estimates_df <- data.frame()
  
  for (lambda_key in names(continuous_est)) {
    lambda_est <- continuous_est[[lambda_key]]
    # calculate the mean and sd from the estimated lambda
    temp_summary <- data.frame(model_type = "Continuous", 
                               lambda_key = lambda_key,
                               mean_lambda = mean(lambda_est),
                               sd_lambda = sd(lambda_est))
    # the estimated lambda for each individuals of each chain under different lambda (2000 *7)
    temp_df <- data.frame(lambda_est = lambda_est, 
                          model_type = "Continuous", 
                          lambda_key = lambda_key)
    
    cont_summary_df <- rbind(cont_summary_df, temp_summary)
    cont_estimates_df <- rbind(cont_estimates_df, temp_df)
  }
  
  return(list(summary = cont_summary_df, estimates = cont_estimates_df))
}

# Extract summary and estimates for continuous models
cont_results <- StoreContinuousResults(continuous_est)
# model type; lamda_key; mean; sd
cont_summary_df <- cont_results$summary
# lambda_est; model type; lambda_key
cont_estimates_df <- cont_results$estimates

# Extract true lambda values from the lambda_key (remove "lambda=")
cont_summary_df <- cont_summary_df %>%
  mutate(true_lambda = as.numeric(sub("lambda=", "", gsub("_.*", "", lambda_key))))


# Plot using ggplot2 and add lines for true lambda and estimated lambda
ggplot(cont_estimates_df, aes(x = lambda_est, fill = model_type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  facet_wrap(~ lambda_key, scales = "free", ncol = 2) +
  labs(title = "Lambda Estimates for Continuous Models", 
       x = "Estimated Lambda", y = "Frequency") +
  
  # Add red vertical lines for true lambda
  geom_vline(data = cont_summary_df, aes(xintercept = true_lambda), color = "red", linetype = "dashed") +
  
  # Add blue vertical lines for estimated mean lambda
  geom_vline(data = cont_summary_df, aes(xintercept = mean_lambda), color = "blue", linetype = "solid") +
  
  # Add text for summary statistics
  geom_text(data = cont_summary_df, 
            aes(label = paste("Mean:", round(mean_lambda, 3), "\nSD:", round(sd_lambda, 3)),
                x = Inf, y = Inf), 
            hjust = 1.1, vjust = 2, size = 3) + 
  scale_fill_manual(values = c("#CC9966", "grey"))


################################################################################
# 2. Estimate FOI with binned age data
################################################################################
# The binned age model
binned_age_model <- stan_model("1.1 Estimate FOI - binned age and serostaus.stan")

# List of bin widths to use for binned age model
bin_sizes <- c(1:5)

# Run the binned model for each lambda value and bin width
for (true_lambda in lambda_values) {
  # Retrieve the previously stored data for continuous models
  sim_survey_design <- survey_designs[[paste0("lambda=", true_lambda)]]
  
  for (bin_size in bin_sizes) {
    # Categorize the individual ages into bins
    bin_breaks <- seq(0, 80, by = bin_size) 
    age_bins <- cut(sim_survey_design$age, breaks = bin_breaks,  
                    include.lowest = TRUE, right = FALSE)
    age_lower <- floor(sim_survey_design$age / bin_size) * bin_size
    age_upper <- age_lower + bin_size
    
    # Create the data frame for binned ages
    binned_survey_design <- tibble(
      age_lower = age_lower,
      age_upper = age_upper,
      prob = sim_survey_design$prob,   
      Y = sim_survey_design$Y          
    ) %>%
      arrange(age_lower)
    
    # Simulated data for Stan (binned age model)
    bin_sim_data <- list(
      N = N,
      age_lower = binned_survey_design$age_lower,
      age_upper = binned_survey_design$age_upper,
      Y = binned_survey_design$Y
    )
    
    # Fit the binned age model
    bin_fit <- sampling(
      binned_age_model,   
      data = bin_sim_data,           
      chains = 4,                
      iter = 1000            
    )
    
    # Store the estimates of lambda (binned age model)
    lambda_est <- rstan::extract(bin_fit, "lambda")[[1]] 
    binned_est[[paste0("lambda=", true_lambda, "_bin=", bin_size)]] <- lambda_est
  }
}


# Store binned results
StoreBinnedResults <- function(binned_est) {
  # store summary statistics (mean and sd)
  bin_summary_df <- data.frame()
  # store estimates of lambda
  bin_estimates_df <- data.frame()
  
  for (lambda_bin_key in names(binned_est)) {
    lambda_est <- binned_est[[lambda_bin_key]]
    # calculate the mean and sd from the estimated lambda
    temp_summary <- data.frame(model_type = "Binned", 
                               lambda_key = lambda_bin_key,
                               mean_lambda = mean(lambda_est),
                               sd_lambda = sd(lambda_est))
    # the estimated lambda for each individual of each chain under different lambda
    temp_df <- data.frame(lambda_est = lambda_est, 
                          model_type = "Binned", 
                          lambda_key = lambda_bin_key)
    
    bin_summary_df <- rbind(bin_summary_df, temp_summary)
    bin_estimates_df <- rbind(bin_estimates_df, temp_df)
  }
  
  return(list(summary = bin_summary_df, estimates = bin_estimates_df))
}

# Extract summary and results for binned models
bin_results <- StoreBinnedResults(binned_est)
bin_summary_df <- bin_results$summary
bin_estimates_df <- bin_results$estimates

# Add corresponding columns to the data frame to combine and plot
bin_summary_df <- bin_summary_df %>%
  mutate(bin_width = as.numeric(sub(".*_bin=", "", lambda_key))) %>%
  mutate(true_lambda = as.numeric(sub("lambda=", "", gsub("_bin=.*", "", lambda_key)))) 

bin_estimates_df <- bin_estimates_df %>%
  mutate(bin_width = as.numeric(sub(".*_bin=", "", lambda_key))) %>%
  mutate(true_lambda = as.numeric(sub("lambda=", "", gsub("_bin=.*", "", lambda_key))))

cont_estimates_df <- cont_estimates_df %>%
  mutate(true_lambda = as.numeric(sub("lambda=", "", gsub("_.*", "", lambda_key)))) %>%
  mutate(bin_width = NA)  # bin_width is not applicable for continuous model

cont_summary_df <- cont_summary_df %>%
  mutate(bin_width = NA) 

# Combine continuous and binned results for comparison
combined_summary_df <- rbind(cont_summary_df, bin_summary_df)
combined_estimates_df <- rbind(cont_estimates_df, bin_estimates_df)

combined_estimates_df$model_type <- factor(combined_estimates_df$model_type, levels = c("Continuous", "Binned"))
combined_estimates_df$bin_width[combined_estimates_df$model_type == "Continuous"] <- 0
combined_summary_df$bin_width[combined_summary_df$model_type == "Continuous"] <- 0


ggplot(combined_estimates_df, aes(x = lambda_est, fill = model_type)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  
  # Separate by lambda and bin width using facet_grid
  facet_grid(bin_width ~ true_lambda, scales = "free_x") +
  
  labs(title = "Lambda Estimates: Continuous vs. Binned Models", 
       x = "Estimated Lambda", y = "Frequency") +
  
  # Add red vertical lines for true lambda
  geom_vline(data = combined_summary_df, aes(xintercept = true_lambda, color = "True Lambda"), linetype = "dashed") +
  
  # Add blue vertical lines for estimated mean lambda
  geom_vline(data = combined_summary_df, aes(xintercept = mean_lambda, color = "Estimated Lambda"), linetype = "solid") +
  
  # Add text for summary statistics
  geom_text(data = combined_summary_df, 
            aes(label = paste("Mean:", round(mean_lambda, 3), "\nSD:", round(sd_lambda, 3)),
                x = Inf, y = Inf), 
            hjust = 1.1, vjust = 2, size = 3) +
  
  # Set colors for the histograms
  scale_fill_manual(values = c("#666666", "#CC9966")) +
  
  # Set colors for the vertical lines (True Lambda and Estimated Lambda)
  scale_color_manual(values = c("True Lambda" = "red", 
                                "Estimated Lambda" = "blue"))  +
  
  # Add secondary Y-axis for bin width, but without the numbers or ticks
  #scale_y_continuous(sec.axis = sec_axis(~ ., name = "Bin Width", breaks = NULL)) +
  
  # Add legend title and theme adjustments
  labs(color = "Legend", fill = "Model Type") +
  theme_minimal() +
  theme(
    legend.position = "top"
  )

################################################################################
# Create a heatmap between true lambda and bin width
################################################################################
# Calculate the absolute error between the estimated and true lambda
combined_summary_df <- combined_summary_df %>%
  mutate(abs_error = abs(mean_lambda - true_lambda))

# Plot the heatmap with a single color and varying transparency
ggplot(combined_summary_df, aes(x = factor(true_lambda), 
                                y = factor(bin_width), fill = abs_error)) +
  geom_tile(color = "white", aes(alpha = abs_error)) + 
  
  scale_fill_gradient(low = "blue", high = "blue", name = "Absolute Error") + 
  scale_alpha_continuous(range = c(0.1, 1), name = "Error Transparency") +  
  
  # Add titles and axis labels
  labs(title = "Error in Estimated Lambda by True Lambda and Bin Width",
       x = "True Lambda",
       y = "Bin Width (years)") +
  theme_minimal() 





