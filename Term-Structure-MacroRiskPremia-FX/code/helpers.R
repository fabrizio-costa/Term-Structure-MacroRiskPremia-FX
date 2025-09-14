if (!require(lubridate)) install.packages("lubridate"); library(lubridate)

###############################################################################
# Data Transformations
###############################################################################

# this function removes all the main currencies that joined the Euro in 1999
removing.eu.currencies <- function(.data){
  .data |> 
    dplyr::select(-c(ATS, BEF, FRF, DEM, IEP, ITL, PTE, ESP))
}

# this function removes some currencies which have low volatility 
removing.low.sd.currencies <- function(.data){
  .data |> 
    dplyr::select(-c(SAR,HKD, QAR, JOD))
}

# this function keeps the eu currencies that joined the Euro in 1999 and removes
# all the other currencies whose starting date is far from 1990 
keeping.eu.currencies <- function(df, euro_currencies = c("ATS", "BEF", "FRF", "DEM", "IEP", "ITL", "PTE", "ESP"), tolerance = 5) {
  currency_cols <- setdiff(names(df), "date")
  
  # first non-NA date for each currency
  start_dates <- df %>%
    pivot_longer(cols = all_of(currency_cols), names_to = "currency", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(currency) %>%
    summarise(start_date = min(date), .groups = "drop")
  
  # get the earliest start date among the EU currencies
  euro_start <- start_dates %>%
    filter(currency %in% euro_currencies) %>%
    summarise(min_start = max(start_date)) %>%
    pull(min_start)
  
  # get all currencies whose start date is within `tolerance` months of euro_start
  valid_currencies <- start_dates %>%
    filter(start_date <= euro_start %m+% months(tolerance)) %>%
    pull(currency)
  
  # filter the data to keep only those currencies (including euro ones)
  df_filtered <- df %>%
    dplyr::select(date, all_of(valid_currencies))
  
  return(df_filtered)
}

# this function filters the dataset to start from the latest first available observation
# across all currencies. Keeps only rows from the latest common start date onward 
# and drops any rows with NA values.

filter_by_latest_start <- function(df) {
  
  currency_cols <- setdiff(names(df), "date")
  
  # for each currency, find the first date with a non-NA
  start_dates <- df %>%
    pivot_longer(cols = all_of(currency_cols), names_to = "currency", values_to = "value") %>%
    group_by(currency) %>%
    filter(!is.na(value)) %>%
    summarise(start_date = min(date), .groups = "drop")
  
  # find the latest of all the starting dates
  latest_start <- max(start_dates$start_date)
  
  # filter to keep only dates from the latest start onward
  df_filtered <- df %>%
    filter(date >= latest_start) %>%
    drop_na()  
  
  return(df_filtered)
}

# The function filters the dataset to end at the earliest last available observation 
# across all currencies. Keeps only rows up to the earliest common end date. 
filter_by_earliest_end <- function(df) {
  currency_cols <- setdiff(names(df), "date")
  
  # Find the last date with non-NA value for each currency
  end_dates <- df %>%
    pivot_longer(cols = all_of(currency_cols), names_to = "currency", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(currency) %>%
    summarise(end_date = max(date), .groups = "drop")
  
  # Find the earliest of the end dates
  earliest_end <- min(end_dates$end_date)
  
  # Filter to keep only dates up to and including the earliest end date
  df_filtered <- df %>%
    filter(date <= earliest_end) %>%
    drop_na()
  
  return(df_filtered)
}


create_macro_currency_datasets <- function(macro_df, currency_df) {
  # All macro vars excluding the date
  macro_vars <- setdiff(names(macro_df), "date")
  
  # Function to trim leading and trailing NAs
  trim_na_ends <- function(x, date) {
    not_na <- !is.na(x)
    if (all(!not_na)) return(tibble(date = date, value = NA_real_))  # if all NA
    first <- which(not_na)[1]
    last  <- rev(which(not_na))[1]
    tibble(date = date[first:last], value = x[first:last])
  }
  
  # Create a named list to store results
  result_list <- vector("list", length(macro_vars))
  names(result_list) <- macro_vars
  
  for (var in macro_vars) {
    # Step 1: trim the macro variable
    trimmed <- trim_na_ends(macro_df[[var]], macro_df$date) %>%
      rename(!!var := value)
    
    # Step 2: calculate the common date range
    common_range <- range(trimmed$date)
    currency_trimmed <- currency_df %>%
      filter(date >= common_range[1], date <= common_range[2])
    
    # Step 3: join macro to currency
    joined <- currency_trimmed %>%
      left_join(trimmed, by = "date")
    
    # Step 4: store in list
    result_list[[var]] <- joined
  }
  
  return(result_list)
}

###############################################################################
# Gibbs sampler
###############################################################################

#' Runs the Gibbs sampler for a selected subset of macroeconomic variables.
#' For each macro variable in df.list[macro_index], extracts g_t and R_t and runs the sampler.
#' Returns a named list of sampler outputs.

run.gibbs.sampler <- function(list, macro_index, S_bar, K, num_iter, burnin, seed) {
  selected_vars <- names(df.list)[macro_index]
  result_list <- vector("list", length = length(selected_vars))
  names(result_list) <- selected_vars
  
  for (i in seq_along(selected_vars)) {
    var <- selected_vars[i]
    df <- df.list[[var]]
    
    message("Running Gibbs sampler for: ", var)
    
    g_t <- df[[var]]
    R <- df %>% dplyr::select(-date, -all_of(var)) %>% as.matrix() 
    R <- scale(R, center = TRUE, scale = TRUE)
    result_list[[i]] <- gibbs_sampler_full(
      g_t = g_t,
      R = R,
      S_bar = S_bar,
      K = K,
      num_iter = num_iter,
      burnin = burnin,
      seed = seed
    )
  }
  
  return(result_list)
}



