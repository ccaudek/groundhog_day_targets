
library(cmdstanr)
library(posterior)
library(bayesplot)
library(rstan)
library(loo)
color_scheme_set("brightblue")

cmdstan_version()

set.seed(42)
n <- 100
x <- runif(n, 30, 40)
a <- 2
b <- 1.5
y <- a + b * x + rnorm(n, 0, 4)
plot(x, y)

stan_data <- list(
  N = n,
  x = x,
  y = y
)


mod <- cmdstan_model("R/stan_code/reg1.stan")

fm <- lm(y ~ x + (1 | ))
coef(fm)


# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 1000,
  iter_warmup = 500,
  chains = 4
)

# Summary of the fit
fit$summary()



### hierarchical model

# Set random seed for reproducibility
set.seed(1234)

# Number of groups
J <- 5

# Number of observations per group
n_per_group <- 100

# Total number of observations
N <- J * n_per_group

# Generate global intercept and mean slope across groups
alpha <- 5
mu <- 2

# Generate standard deviation of slopes across groups and residual standard deviation
tau <- 0.5
sigma <- 1

# Generate group-specific slopes
beta <- rnorm(J, mean = mu, sd = tau)

# Generate group identifiers and pre-allocate vectors for x and y
group <- rep(1:J, each = n_per_group)
x = rep(NA, N)
y = rep(NA, N)

# Generate data for each group
for (j in 1:J) {
  idx <- which(group == j)
  x[idx] <- runif(n_per_group, 0, 10)  # Predictor variable x in [0, 10]
  y[idx] <- rnorm(n_per_group, mean = alpha + beta[j] * x[idx], sd = sigma)
}

# Combine into a data frame
data_df <- data.frame(y = y, x = x, group = factor(group))

# Show the first few rows of the data
head(data_df)

mod <- cmdstan_model("R/stan_code/reg1.stan")


stan_data <- list(
  N = N, 
  J = J, 
  y = data_df$y, 
  x = data_df$x, 
  group = as.integer(data_df$group)
  )



# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 1000,
  iter_warmup = 500,
  chains = 4
)

# Summary of the fit
fit$summary()



# Tutorial ---------------------------------------------------------------------

# https://vasishth.github.io/bayescogsci/book/ch-introstan.html#ref-Stan2021

# target notation
stan_code <- "
data {
  int<lower = 1> N;  // Total number of trials
  vector[N] y;  // Score in each trial
}

parameters {
  real mu;
  real<lower = 0> sigma;
}

model {
  // Priors:
  target += normal_lpdf(mu | 0, 20);
  target += lognormal_lpdf(sigma | 3, 1);
  // Likelihood:
  for(i in 1:N)
    target += normal_lpdf(y[i] | mu, sigma);
}
"

# Sampling notation
stan_code <- "
data {
  int<lower = 1> N;  // Total number of trials
  vector[N] y;  // Score in each trial
}

parameters {
  real mu;
  real<lower = 0> sigma;
}

model {
  // Priors:
  mu ~ normal(0, 20);
  sigma ~ exponential(1.0/10);
  // Likelihood:
  y ~ normal(mu, sigma); 
}
"

stan_file <- write_stan_file(stan_code)
mod <- cmdstan_model(stan_file)


Y <- rnorm(n = 100, mean = 3, sd = 10)

stan_data <- list(
  N = length(Y), 
  y = Y
)

# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  chains = 4
)

fit$summary()
fit$summary(variables = c("mu", "sigma", "lp__"), "mean", "sd")

fit$summary("mu", pr_lt_half = ~ mean(. <= 0.0))

fit$summary(
  variables = NULL,
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(.0275, .975))
)

mcmc_hist(fit$draws("mu"))

stanfit <- rstan::read_stan_csv(fit$output_files())

mcmc_trace(stanfit, pars = c("mu", "sigma"))

mcmc_intervals(stanfit, pars = c("mu", "sigma"))

mcmc_areas(
  stanfit, 
  pars = c("mu", "sigma"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.95, # 95%
  point_est = "mean"
)


fit$diagnostic_summary()

# Let's go back to the APRL data.
tar_load(params_happiness_df)
glimpse(params_happiness_df)

# Get the data of only ema_number == 1

ema1_df <- params_happiness_df |> 
  dplyr::filter(ema_number == 1) |> 
  dplyr::select(user_id, is_reversal, trial, happiness)

ema1_df$subj_ix <- as.numeric(as.factor(as.character(ema1_df$user_id)))
ema1_df$trial_num <-  ema1_df$trial - 1 # trial starts from 

ema1_sorted_df <- ema1_df[order(ema1_df$subj_ix, ema1_df$trial_num), ]

N <- length(unique(ema1_sorted_df$user_id))
T_per_subject <- 30
T <- N * T_per_subject

# Create the list to be passed to cmdstan
stan_data <- list(
  N = N,
  T = T,
  subj_ix = ema1_sorted_df$subj_ix,
  trial_num = ema1_sorted_df$trial_num,
  Y = ema1_sorted_df$happiness
)




# Preparing the data for Stan in a more robust and efficient manner
prepare_stan_data <- function(df) {
  
  # Filter and select columns, sort by 'user_id' and 'trial'
  prepared_df <- df %>%
    filter(ema_number == 1) %>%
    select(user_id, trial, instant_mood) %>%
    arrange(user_id, trial)
  
  # Convert 'user_id' to a numerical index for Stan
  prepared_df$subj_ix <- as.integer(as.factor(prepared_df$user_id))
  
  # Make trial numbers zero-based
  prepared_df$trial_num <- prepared_df$trial - 1
  
  # Number of unique subjects and trials per subject
  N <- n_distinct(prepared_df$user_id)
  T_per_subject <- max(prepared_df$trial_num) + 1
  T <- N * T_per_subject
  
  # Create the list to be passed to cmdstan
  stan_data <- list(
    N = N,
    T = T,
    subj_ix = prepared_df$subj_ix,
    trial_num = prepared_df$trial_num,
    Y = prepared_df$instant_mood
  )
  
  return(stan_data)
}


# Identify unique user_ids
unique_user_ids <- unique(prl_df$user_id)

# Randomly sample 10 user_ids
random_user_ids <- sample(unique_user_ids, 10)

# Filter the original data frame to include only these 10 user_ids
filtered_df <- prl_df %>% 
  filter(user_id %in% random_user_ids)

stan_data <- prepare_stan_data(filtered_df)
stan_data <- prepare_stan_data(prl_df)


mod <- cmdstan_model("R/stan_code/M0a.stan")


# Fit the model to the data
fit <- mod$sample(
  data = stan_data,
  iter_sampling = 2000,
  iter_warmup = 1000,
  parallel_chains = 4,
  chains = 4,
  refresh = 500
)

fit$diagnostic_summary()

# In order to use this method you must compute and save the pointwise 
# log-likelihood in your Stan program.
loo_m0a <- fit$loo(cores = 4)
print(loo_m0a)

fit$summary()


fit_vb <- mod$variational(
  data = stan_data,
  seed = 123
)

fit$diagnostic_summary()
fit$summary()

summary_fit <- fit$summary()
# Display the column names of the summary data frame
colnames(summary_fit)
summary_fit$variable




