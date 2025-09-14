rm(list = ls())
if (dev.cur() != 1) dev.off()

if (!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if (!require(purrr)) install.packages("purrr"); library(purrr)
if (!require(readxl)) install.packages("readxl"); library(readxl)

###############################################################################
# Data Transformations
###############################################################################
source("code/helpers.R")

# --------------------------------- Currencies ---------------------------------
currency.df <- readRDS("data/currency_data/df_short.rds") |> 
  # scenario 1: removing all the currencies that joined the Euro in 1999
  removing.eu.currencies() |>
  
  # scenario 2: keeping all the currencies that joined the Euro in 1999 and
  # excluding the currencies whose start date is much later (tolerance)
  # keeping.eu.currencies(tolerance = 5) |>
  
  # removing all currencies with low volatility
  removing.low.sd.currencies() |> 
  
  # filtering the dataset to start from the latest first available observation across all currencies.
  filter_by_latest_start() |>
  
  # filtering the dataset to end at the earliest last available observation across all currencies.
  filter_by_earliest_end() 

# --------------------------------- Macro Data ---------------------------------

macro.df.m <- read_xlsx("data/macroeconomic_data/all_data.xlsx", sheet = "monthly") |>
  dplyr::select(-ends_with("100")) |> 
  # standardising variables
  mutate(across(
    .cols = -date,  
    .fns = ~ scale(., center = TRUE, scale = TRUE)[, 1]
  ))

# --------------------------------- Final Data ---------------------------------

# for each macroeconomic variable, trims leading and trailing NAs, then we 
# left-joins it (by date) to the currency dataset over the overlapping date range.
# In return, we a list of data frames, each containing the currency data and one macro variable.
df.list <- create_macro_currency_datasets(macro.df.m, currency.df)

###############################################################################
# Gibbs Sampler
###############################################################################

# Importing the gibbs sampler function
source("code/gibbs.sampler.R")

# Settings
K <- 3
S_bar <- 24
num_iter <- 500
burnin <- 100
seed <- 42

# Number of macroeconomic variables to run
macro_index <- 1:length(names(df.list))
selected_vars <- names(df.list)[macro_index]

results <- run.gibbs.sampler(
  list = df.list,
  macro_index = macro_index,
  S_bar = S_bar,
  K = K,
  num_iter = num_iter,
  burnin = burnin,
  seed = seed
)

# write_rds(results, paste0("data/results/K=", K, "&seed=", seed,"&tau_bv=0.01&macrostandardised&rho_g=fixed&lambda_v=fixed.RDS"))
# results <- read_rds(paste0("data/results/K=", K, "&seed=", seed,"&tau_bv=0.01&macrostandardised.RDS"))
