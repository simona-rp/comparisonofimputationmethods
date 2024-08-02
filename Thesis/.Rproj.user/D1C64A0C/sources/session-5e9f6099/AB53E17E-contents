# Install and load required packages----
if(is.na(packageDate("pacman"))) install.packages("pacman")
library(pacman)
p_load(skimr, tidyr, ggplot2, dplyr, mice, cheapr, bestNormalize, corrplot, tibble, Metrics, Hmisc, VIM, microbenchmark, MASS, wesanderson, bench, scales, ks)

# Import data from csv files----
wine_red <- read.csv(file = "winequality-red.csv", 
                     header = TRUE,
                     sep = ";",
                     quote = "\"")

wine_red['type'] = 0

wine_white <- read.csv(file = "winequality-white.csv", 
                       header = TRUE,
                       sep = ";",
                       quote = "\"")
wine_white['type'] = 1

## Merge the two datasets into one dataset and save as RDS file
wine <- rbind(wine_red, wine_white)
saveRDS(wine, "wine.rds")
data <- readRDS("wine.rds")

# Data analysis----
## Data size and classes----
tibble(data)

## Classes in the data
classes <- data %>%
  mutate(label = case_when(type == 0 ~ "Red wine",
                           type == 1 ~ "White wine")) %>%
  ggplot(aes(label)) +
  geom_bar(alpha = 0.7, fill="black", color="white") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.25, "lines"),
    text = element_text(size = 20)) +
  xlab("") +
  ylab("Frequency")

data <- data %>%
  filter(type == 1) %>%
  dplyr::select(-type)
saveRDS(data, "data.rds")
data <- readRDS("data.rds")

## Descriptive analysis----
skimr::skim(data)

### Skewness and kurtosis
data_desc <- psych::describe(data, type = 2)

### Histograms and boxplots
###Prepare data
data_rename <- function(data){
  data_renamed <- data %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = recode(variable, 
                             "fixed.acidity" = "Fixed Acidity (g/L)",
                             "volatile.acidity" = "Volatile Acidity (g/L)",
                             "citric.acid" = "Citric Acid (g/L)",
                             "residual.sugar" = "Residual Sugar (g/L)",
                             "chlorides" = "Chlorides (g/L)",
                             "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                             "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                             "density" = "Density (g/ml)",
                             "pH" = "pH",
                             "sulphates" = "Sulphates (g/L)",
                             "alcohol" = "Alcohol (% vol)",
                             "quality" = "Quality"))
  return(data_renamed)
}

data_plot <- data_rename(data)

hist <- data_plot %>%
  ggplot2::ggplot(aes(value)) +
  facet_wrap(~ variable, scales = "free") +
  geom_histogram(bins = 15, alpha = 0.7, fill="black", color="white") +
  theme_classic() +
  theme(legend.position="none",
        panel.spacing = unit(0.25, "lines"),
        text = element_text(size = 15)) +
  xlab("") +
  ylab("Frequency")


boxplot <- data_plot %>%
  ggplot2::ggplot(aes(variable, value)) +
  geom_boxplot() +
  facet_wrap(~ variable, scales = "free") +
  theme_classic() +
  theme(legend.position="none",
        panel.spacing = unit(0.25, "lines"),
        text = element_text(size = 15)) +
  xlab("") +
  ylab("")


### Normality
nor_results <- sapply(data, function(x) {
  shapiro_result <- shapiro.test(x)
  p_value <- shapiro_result$p.value
  normality_indicator <- ifelse(p_value <= 0.05, "not normal", "normal")
  c(p_value, normality_indicator)
})

variables_unit <- c("Fixed Acidity (g/L)","Volatile Acidity (g/L)", "Citric Acid (g/L)","Residual Sugar (g/L)","Chlorides (g/L)","Free Sulfur Dioxide (mg/L)","Total Sulfur Dioxide (mg/L)","Density (g/ml)","pH","Sulphates (g/L)","Alcohol (% vol)","Quality")
nor_results_df <- as.data.frame(t(nor_results))
colnames(nor_results_df) <- c("p_value", "normality_indicator")
nor_results_df$variable <- variables_unit
nor_results_df <- nor_results_df[,c("variable", "p_value", "normality_indicator")]
as_tibble(nor_results_df)

### Outliers
z_scores <- function(data) {

  # Calculate Z-scores
  z_scores <- as.data.frame(scale(data))
  
  # Count Z-scores > 3 or < -3 for each variable
  extreme_counts <- z_scores %>%
    summarise(across(everything(), ~ sum(. > 3 | . < -3, na.rm = TRUE)))
  
  list(
    z_scores = z_scores,
    extreme_counts = extreme_counts
  )
}
outliers_mean <- z_scores(data)

z_scores_adj <- function(data) {

  # Function to compute the modified Z-score 
  z_score_adj <- function(x) {
    med <- median(x, na.rm = TRUE)
    mad <- median(abs(x - med), na.rm = TRUE)
    z <- 0.6745 * (x - med) / mad
    return(z)
  }
  
  # Apply the modified Z-score function to each column in the dataframe
  z_scores_adj <- reframe(data, across(everything(), z_score_adj))
  # Count Z-scores > 3 or < -3 for each variable
  extreme_counts <- reframe(z_scores_adj, 
                            across(everything(), ~ sum(. > 3.5 | . < -3.5, na.rm = TRUE)))
  
  list(
    z_scores = z_scores,
    extreme_counts = extreme_counts
  )
}

outliers_median <- z_scores_adj(data)

outliers <- tibble(variables_unit,
                   unlist(outliers_mean$extreme_counts),
                   unlist(outliers_median$extreme_counts))
                                                                                 
names(outliers) <- c("variable", "outliers_mean", "outliers_median")
outliers

### Duplicates
duplicates <- sum(duplicated(data), na.rm = TRUE) #937 duplicates


# Missing data generation----
## Default specifications of ampute()
data_mar_def <- ampute(data, 0.3, mech ="MAR")
def_pattern <- md.pattern(data_mar_def$amp, 
                          rotate.names = TRUE)

## Set seed for reproducibility 
set.seed(12345)

## Random pattern generation
patterns <- matrix(sample(c(0, 1), 12 * ncol(data), 
                          replace = TRUE), 
                   nrow = 12, 
                   ncol = ncol(data), 
                   byrow = TRUE)

## Frequencies of each pattern
frequencies <- rnorm(12, 100, 5) 
frequencies <- frequencies / sum(frequencies)
data_mar <- ampute(data, 0.3, 
                   mech = "MAR", 
                   patterns = patterns, 
                   freq = frequencies)
md.pattern(data_mar$amp, rotate.names = TRUE)

## Amount of missing data
mar_cells <- num_na(data_mar$amp) / unlisted_length(data) # ca. 14%
mar_cases <- sum(row_any_na(data_mar$amp)) / nrow(data) # ca. 30%

## Per column and row
mar_col <- col_na_counts(data_mar$amp)
mal_row <- row_na_counts(data_mar$amp)

## Histograms
hist_mar <- as_tibble(mal_row) %>%
  ggplot2::ggplot(aes(value)) +
  geom_histogram(binwidth = 1, alpha = 0.7, fill="black", color="white") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.25, "lines"),
    text = element_text(size = 20)) +
  xlab("Missing values per case") +
  ylab("Frequency") +
  scale_x_continuous(breaks = 0:ncol(data_mar$amp))

## Transform data to a tibble
data_mar <- data_mar$amp  
data_mar <- as_tibble(data_mar,.name_repair = "minimal")
colnames(data_mar) <- colnames(data)
saveRDS(data_mar, "data_mar.rds")
data_mar <- readRDS("data_mar.rds")

# Software overview----
## Package use
package_use <- data.frame(package = c("Hmisc", "mice", "mi", "VIM", "yaImpute"),
                          downloads = c(180929, 48515, 19894, 12339, 3012))

## Plot
packages <- ggplot(package_use, aes(x = package, y = downloads)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.25, "lines"),
    text = element_text(size = 20)) +
  xlab("Package") +
  ylab("Number of downloads per month")

# Data preparation----
## Normality approximation with log, Box-Cox and Yeo-Johnson-Transformation
### Log10 Transformation
data_mar_log10 <- log10(data_mar+1)

### Data preparation for Box-Cox and Yeo-Johnson
matrix <- as.matrix(data_mar)
names(matrix) <- NULL

### Box-Cox-Transformation
data_mar_trans_bc <- bestNormalize::boxcox(matrix+1)
data_mar_bc <- data_mar_trans_bc$x.t
data_mar_trans_bc$lambda

### Yeo-Johnson-Transformation
data_mar_trans_yj <- yeojohnson(matrix, eps = 0.001, standardize = TRUE)
data_mar_yj <- data_mar_trans_yj$x.t
data_mar_trans_yj$lambda

### Transform results to a tibble
data_mar_bc <- as_tibble(data_mar_bc, .name_repair = "minimal")
colnames(data_mar_bc) <- colnames(data_mar)

data_mar_yj <- as_tibble(data_mar_yj, .name_repair = "minimal")
colnames(data_mar_yj) <- colnames(data_mar)
saveRDS(data_mar_yj, "data_mar_yj.rds")
data_mar_yj <- readRDS("data_mar_yj.rds")

## Histograms to check transformed distributions
### Data preparation
data_mar_log10_long <- data_rename(data_mar_log10)
data_mar_bc_long <- data_rename(data_mar_bc)
data_mar_yj_long <- data_rename(data_mar_yj)
data_original <- data_rename(data)

data_mar_log10_plot <- cbind(data_original, data_mar_log10_long[,2])
colnames(data_mar_log10_plot) <- c("variable", "original", "log10")
data_mar_log10_plot <- data_mar_log10_plot %>%
  pivot_longer(-variable, names_to = "method", values_to = "value")

data_mar_bc_plot <- cbind(data_original, data_mar_bc_long[,2])
colnames(data_mar_bc_plot) <- c("variable", "original", "bc")
data_mar_bc_plot <- data_mar_bc_plot %>%
  pivot_longer(-variable, names_to = "method", values_to = "value")

data_mar_yj_plot <- cbind(data_original, data_mar_yj_long[,2])
colnames(data_mar_yj_plot) <- c("variable", "original", "yj")
data_mar_yj_plot <- data_mar_yj_plot %>%
  pivot_longer(-variable, names_to = "method", values_to = "value")

hist_log10 <- data_mar_log10_plot %>%
  ggplot2::ggplot(aes(x = value, fill = method)) +
  geom_histogram(bins = 15, alpha = 1, color = "black") +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual("Data",values = c("log10" = "#D69C4E", "original" = "#046C9A"), labels = c("log10" = "Log10 transformed", "original" = "Original")) +
  theme_classic() +
  theme(legend.position = "top",
        panel.spacing = unit(0.25, "lines"),
        text = element_text(size = 15)) +
  xlab("") +
  ylab("Frequency")

hist_bc <- data_mar_bc_plot %>%
  ggplot2::ggplot(aes(x = value, fill = method)) +
  geom_histogram(bins = 15, alpha = 1, color = "black") +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual("Data",values = c("bc" = "#D69C4E", "original" = "#046C9A"), labels = c("bc" = "Box-Cox transformed", "original" = "Original")) +
  theme_classic() +
  theme(legend.position = "top",
        panel.spacing = unit(0.25, "lines"),
        text = element_text(size = 15)) +
  xlab("") +
  ylab("Frequency")

hist_yj <- data_mar_yj_plot %>%
  ggplot2::ggplot(aes(x = value, fill = method)) +
  geom_histogram(bins = 15, alpha = 1, color = "black") +
  facet_wrap(~ variable, scales = "free") +
  scale_fill_manual("Data",values = c("yj" = "#D69C4E", "original" = "#046C9A"), labels = c("yj" = "Yeo-Johnson transformed", "original" = "Original")) +
  theme_classic() +
  theme(legend.position = "top",
        panel.spacing = unit(0.25, "lines"),
        text = element_text(size = 15)) +
  xlab("") +
  ylab("Frequency")

## Check skewness after transformations
data_mar_bc_skew <- psych::describe(data_mar_bc, type = 2)
data_mar_yj_skew <- psych::describe(data_mar_yj, type = 2)
data_mar_log10_skew <- psych::describe(data_mar_log10, type = 2)
data_mar_skew <- psych::describe(data_mar, type = 2)
data_skew <- psych::describe(data, type = 2)
skewness <- tibble(variable = variables_unit, skew_original = data_skew$skew, skew_mar = data_mar_skew$skew, skew_log10 = data_mar_log10_skew$skew, skew_bc = data_mar_bc_skew$skew, skew_yj = data_mar_yj_skew$skew)

## Check normality after transformations
nor_results_mar_bc <- sapply(data_mar_bc, function(x) {
  shapiro_result <- shapiro.test(x)
  p_value <- shapiro_result$p.value
  normality_indicator <- ifelse(p_value <= 0.05, "not normal", "normal")
  c(p_value, normality_indicator)
})

nor_results_mar_yj <- sapply(data_mar_yj, function(x) {
  shapiro_result <- shapiro.test(x)
  p_value <- shapiro_result$p.value
  normality_indicator <- ifelse(p_value <= 0.05, "not normal", "normal")
  c(p_value, normality_indicator)
})

nor_results_mar_log10 <- sapply(data_mar_log10, function(x) {
  shapiro_result <- shapiro.test(x)
  p_value <- shapiro_result$p.value
  normality_indicator <- ifelse(p_value <= 0.05, "not normal", "normal")
  c(p_value, normality_indicator)
})

nor_results_mar_df <- as.data.frame(cbind(t(nor_results_mar_bc), t(nor_results_mar_yj), t(nor_results_mar_log10)))
colnames(nor_results_mar_df) <- c("p_value_bc", "normality_indicator_bc","p_value_yj", "normality_indicator_yj", "p_value_log10", "normality_indicator_log10")
nor_results_mar_df$variable <- variables_unit
nor_results_mar_df <- nor_results_mar_df[,c("variable", "p_value_bc", "normality_indicator_bc","p_value_yj", "normality_indicator_yj", "p_value_log10", "normality_indicator_log10")]
as_tibble(nor_results_mar_df)

## Format final dataset with MAR
data_mar_yj <- as_tibble(data_mar_yj)

# Predictor analysis----
##Correlations
### Original data
cor_matrix <- cor(data)
cor_plot <- corrplot(cor_matrix, method = "color", 
                     type = "upper", tl.col = "black",
                     addCoef.col = "black", number.cex = 0.6)
### MAR data
cor_matrix_mar <- cor(data_mar_yj, use = "pairwise.complete.obs")
cor_plot_mar <- corrplot(cor_matrix_mar, method = "color", 
                     type = "upper", tl.col = "black", 
                     addCoef.col = "black", number.cex = 0.6)

## Influx and outflux
flux <- flux(data_mar_yj)
flux <- rownames_to_column(flux, var = "variable")
flux_overview <- tibble(variable = variables_unit,
                        pobs = flux[,2],
                        influx = flux[,3],
                        outflux = flux[,4])
fluxplot <- fluxplot(data_mar_yj, main = "Influx-outflux pattern", cex.lab = 1.5, cex = 1.3)

# Imputations----

## Simulation preparation
set.seed(0987)  # Setting a seed for reproducibility
first_seed <- 12345 # First seed from the section Missing Data generation
seeds <- sample(10000:100000, 49) # Generating 49 different seeds
seeds_sim <- c(first_seed, seeds) # Seeds for generating 50 different data sets with missing values

## Simulation of 50 different datasets with MAR missing data and a Yeo-Johnson transformation
multi_ampute_yj <- function(data, num_runs = 50) {

  # Initialize vectors to store results
  mar_cells_list <- numeric(num_runs)
  mar_cases_list <- numeric(num_runs)
  mar_col_list <- numeric(num_runs)
  mar_row_list <- numeric(num_runs)
  
  # Create a list to store the datasets
  datasets_list <- vector("list", num_runs)
  trans_list <- vector("list", num_runs)
  
  for (i in 1:num_runs) {
    set.seed(seeds_sim[i])
    
    ## Apply Yeo-Johnson transformation 
    matrix_sim <- as.matrix(data)
    colnames(matrix_sim) <- colnames(data)  # Ensure the column names match the original data
    yj_objects <- yeojohnson(matrix_sim, eps = 0.001, standardize = TRUE)
    matrix_sim <- yj_objects$x.t
    colnames(matrix_sim) <- colnames(data)  # Ensure the column names are retained after transformation
    
    ## Random pattern generation
    patterns_sim <- matrix(sample(c(0, 1), 12 * ncol(data), replace = TRUE), 
                       nrow = 12, ncol = ncol(data), byrow = TRUE)
    
    ## Frequencies of each pattern
    frequencies_sim <- rnorm(12, 100, 5)
    frequencies_sim <- frequencies_sim / sum(frequencies_sim)
    
    ## Ampute data
    data_mar_sim <- ampute(matrix_sim, 0.3, mech = "MAR", patterns = patterns_sim, freq = frequencies_sim)
    
    ## Store the dataset in the list
    datasets_list[[i]] <- as.data.frame(data_mar_sim$amp)
    trans_list[[i]] <- yj_objects
    
    ## Amount of missing data
    mar_cells_sim <- num_na(data_mar_sim$amp) / unlisted_length(data)
    mar_cases_sim <- sum(row_any_na(data_mar_sim$amp)) / nrow(data)
    
    mar_cells_list[i] <- mar_cells_sim
    mar_cases_list[i] <- mar_cases_sim
    
    ## Per column and row
    mar_col_sim <- col_na_counts(data_mar_sim$amp)
    mar_row_sim <- row_na_counts(data_mar_sim$amp)
    
    mar_col_list[i] <- mar_col_sim
    mar_row_list[i] <- mar_row_sim
  }
  
  # Calculate the means
  mean_mar_cells <- mean(mar_cells_list)
  mean_mar_cases <- mean(mar_cases_list)
  mean_mar_col <- mean(mar_col_list)
  mean_mar_row <- mean(mar_row_list)
  
  return(list(mean_mar_cells = mean_mar_cells, 
              mean_mar_cases = mean_mar_cases,
              mean_mar_col = mean_mar_col,
              mean_mar_row = mean_mar_row,
              datasets = datasets_list,
              trans = trans_list))
}

## Create datasets----
### With MAR missing data and a Yeo Johnson transformation
datasets_amputed_yj <- multi_ampute_yj(data, num_runs = 50)
saveRDS(datasets_amputed_yj, "datasets_amputed_yj.rds")
datasets_amputed_yj <- readRDS("datasets_amputed_yj.rds")
#### Check missingness amount
datasets_amputed_yj$mean_mar_cells #ca 15% 
datasets_amputed_yj$mean_mar_cases #ca 30%
#### Save datasets with missing values
datasets_mar_yj <- datasets_amputed_yj$datasets
saveRDS(datasets_mar_yj, "datasets_mar_yj.rds")
datasets_mar_yj <- readRDS("datasets_mar_yj.rds")
### Transform original data with Yeo-Johnson to be comparable
data_yj <- as.matrix(data)
data_yj <- yeojohnson(data_yj, eps = 0.001, standardize = TRUE)$x.t
colnames(data_yj) <- colnames(data)
saveRDS(data_yj, "data_yj.rds")
data_yj <- readRDS("data_yj.rds")

# Functions----
## Function to perform different imputations on the 50 MAR datasets 
imputation <- function(datasets_list, original_data, method) {
  
  # Initialize lists to store imputed datasets, MAE, RMSE, means, std dev, skewness
  imputed_datasets_list <- vector("list", length(datasets_list))
  mae_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  rmse_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  mean_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  sd_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  skewness_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  
  for (i in 1:length(datasets_list)) {
    data_mar <- datasets_list[[i]]
    
    ## Apply the specified imputation method using ifelse
    if (method == "mean") {
      mean_imp <- mice(data_mar, method = "mean", m = 1, maxit = 1, printFlag = FALSE)
      imputed_data <- complete(mean_imp)
    } else if (method == "median") {
      imputed_data <- data_mar %>% dplyr::mutate_all(~ Hmisc::impute(., fun = median))
    } else if (method == "zero") {
      imputed_data <- data_mar %>% dplyr::mutate_all(~ Hmisc::impute(., fun = function(x) 0))
    } else if (method == "regression") {
      methods <- make.method(data_mar)
      methods['quality'] <- 'polr'
      methods[names(data_mar)[names(data_mar) != 'quality']] <- 'norm.predict'
      data_mar$quality <- factor(data_mar$quality, ordered = TRUE)
      reg_imp <- mice(data_mar, method = methods, m = 1, maxit = 1, printFlag = FALSE)
      imputed_data <- complete(reg_imp)
    } else if (method == "stochastic regression") {
      methods <- make.method(data_mar)
      methods['quality'] <- 'polr'
      methods[names(data_mar)[names(data_mar) != 'quality']] <- 'norm.nob'
      data_mar$quality <- factor(data_mar$quality, ordered = TRUE)
      sreg_imp <- mice(data_mar, method = methods, m = 1, maxit = 1, printFlag = FALSE)
      imputed_data <- complete(sreg_imp)
    } else if (method == "knn") {
      imputed_data <- kNN(data_mar, 
                          variable = colnames(data_mar),
                          k = 5,
                          dist_var = colnames(data_mar),
                          trace = FALSE,
                          imp_var = FALSE,
                          addRF = FALSE,
                          onlyRF = FALSE,
                          addRandom = FALSE,
                          useImputedDist = FALSE)
    } else {
      stop("Unknown imputation method")
    }
    
    ## Convert the ordinal variable 'quality' to numeric
    imputed_data$quality <- as.numeric(as.character(imputed_data$quality))
    
    ## Store the imputed dataset in the list
    imputed_datasets_list[[i]] <- imputed_data
    
    ## Calculate MAE and RMSE for each variable
    for (j in 1:ncol(original_data)) {
      missing_indices <- which(is.na(data_mar[[j]]))
      mae_list[i, j] <- mae(original_data[missing_indices, j], imputed_data[missing_indices, j])
      rmse_list[i, j] <- rmse(original_data[missing_indices, j], imputed_data[missing_indices, j])
      }
    
    ## Calculate mean, standard deviation, and skewness for each variable
    desc_stats <- psych::describe(imputed_data, type = 2)
    for (j in 1:ncol(original_data)) {
      mean_list[i, ] <- desc_stats$mean
      sd_list[i, ] <- desc_stats$sd
      skewness_list[i, ] <- desc_stats$skew
    }
    
  }
  
  # Calculate average metrics per feature
  average_mae <- colMeans(mae_list, na.rm = TRUE)
  average_rmse <- colMeans(rmse_list, na.rm = TRUE)
  average_mean <- colMeans(mean_list, na.rm = TRUE)
  average_sd <- colMeans(sd_list, na.rm = TRUE)
  average_skewness <- colMeans(skewness_list, na.rm = TRUE)
  
  ## Create histograms for each variable using facet_wrap
  imputed_data_long <- bind_rows(imputed_datasets_list) %>% 
    pivot_longer(cols = everything(), names_to = "name", values_to = "value")
  
  histogram <- ggplot(data = imputed_data_long, aes(x = value)) +
    geom_histogram(alpha = 0.7, fill = "black", color = "white") +
    facet_wrap(~ name, scales = "free") +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.25, "lines"),
      text = element_text(size = 14)
    ) +
    xlab("") +
    ylab("Frequency")
  
  return(list(imputed_datasets = imputed_datasets_list, 
              imputed_datasets_long = imputed_data_long,
              average_mae = average_mae, 
              average_rmse = average_rmse,
              average_mean = average_mean,
              average_sd = average_sd,
              average_skewness = average_skewness,
              histogram = histogram
              ))
}

## Function to reverse the Yeo-Johnson transformation of the imputed datasets if needed
reverse_yj <- function(imputed_datasets_list, transformations_list) {
  # Create list to store reverse transformed datasets
  reverse_imputed_datasets <- vector("list", length(imputed_datasets_list))
  
  # Loop through each combination of datasets
  for (i in 1:length(imputed_datasets_list)) {
    ## Extract each imputed dataset and each yeojohnson object
    imputed_data <- imputed_datasets_list[[i]]
    trans_data <- transformations_list[[i]]
    
    ## Reverse transform
    rev_imputed_data <- predict(trans_data, newdata = as.matrix(imputed_data), inverse = TRUE)
    
    ## Transform reverse imputed data to dataframe and round quality to be discrete
    rev_imputed_data <- as.data.frame(rev_imputed_data)
    colnames(rev_imputed_data) <- colnames(imputed_data)
    #rev_imputed_data <- rev_imputed_data %>%
      #mutate(quality = round(quality))
    
    reverse_imputed_datasets[[i]] <- rev_imputed_data
  }
  return(list(reverse_imputed_datasets = reverse_imputed_datasets))
}

## Function to calculate mean regression coefficients and coverage rates
coverage <- function(datasets_list, original_data) {

  # Fit the original model to get the number of coefficients
  fit_original <- polr(factor(quality, ordered = TRUE) ~ ., data = original_data, method = "logistic", Hess = TRUE)
  original_coefs <- c(coef(fit_original), fit_original$zeta)
  num_coefs <- length(original_coefs)
  
  # For storing coverage rate data
  coverage_counts <- matrix(0, nrow = length(datasets_list), ncol = num_coefs+1)
  coefficients_list <- matrix(NA, nrow = length(datasets_list), ncol = num_coefs+1)
  
  for (i in 1:length(datasets_list)) {
    imputed_data <- datasets_list[[i]]
    
    ## Fit linear model on imputed data
    fit_imputed <- polr(factor(quality, ordered = TRUE) ~ ., data = imputed_data, method = "logistic", Hess = TRUE)
    
    ## Extract coefficients and confidence intervals
    imputed_coefs <- c(coef(fit_imputed), fit_imputed$zeta)
    num_imputed_coefs <- length(imputed_coefs)
    coefficients_list[i, 1:num_imputed_coefs] <- imputed_coefs
    
    ## Get the standard errors from the imputed model
    imputed_se <- summary(fit_imputed)$coefficients[, "Std. Error"]
    
    ## Calculate Wald confidence intervals
    z_critical <- qnorm(0.975)
    lower_ci <- imputed_coefs - z_critical * imputed_se
    upper_ci <- imputed_coefs + z_critical * imputed_se

    ## Count if original coefficient is within the confidence interval of imputed model
    for (k in 1:length(original_coefs)) {
      if (original_coefs[k] >= lower_ci[k] && original_coefs[k] <= upper_ci[k]) {
        coverage_counts[i, k] <- 1
      } else {
        coverage_counts[i, k] <- 0
      }
     }
    
  }
  
  coverage_rates <- colMeans(coverage_counts)
  mean_coefficients <- colMeans(coefficients_list, na.rm = TRUE)
  
  ## Combine coverage rates with coefficient names
  original_coefs <- c(original_coefs, rep(NA, num_imputed_coefs - num_coefs))
  
  num_intercepts <- length(fit_imputed$zeta)
  coefficient_names <- c(names(coef(fit_original)), paste0("Intercept", 1:num_intercepts))
  coverage_rates_named <- tibble(variable = coefficient_names,
                                original_coefficients = original_coefs,
                                mean_coefficients = mean_coefficients[1:num_imputed_coefs],
                                coverage_rate = coverage_rates[1:num_imputed_coefs])
  
  return(list(coverage_rates = coverage_rates_named,
              coefficients = coefficients_list))
}

## Function to create convergence plots
convergence_plot <- function(data, method, n) {
  if (method == "mean") {
    # Transform the array into a long format
    chain_mean_imp_df <- as.data.frame(as.table(data))
    names(chain_mean_imp_df) <- c("variable", "iteration", "Imputation", "mean") 
    # Rename Chain to Imputation
    chain_mean_imp_df$Imputation <- gsub("Chain ", "Imputation ", chain_mean_imp_df$Imputation)
    # Transform iteration to factor to preserve order
    chain_mean_imp_df$iteration <- factor(chain_mean_imp_df$iteration, ordered = TRUE,  levels = c(1:n))
    # Rename variables
    chain_mean_imp_df <- chain_mean_imp_df %>%
      mutate(variable = recode(variable, 
                               "fixed.acidity" = "Fixed Acidity (g/L)",
                               "volatile.acidity" = "Volatile Acidity (g/L)",
                               "citric.acid" = "Citric Acid (g/L)",
                               "residual.sugar" = "Residual Sugar (g/L)",
                               "chlorides" = "Chlorides (g/L)",
                               "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                               "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                               "density" = "Density (g/ml)",
                               "pH" = "pH",
                               "sulphates" = "Sulphates (g/L)",
                               "alcohol" = "Alcohol (% vol)"))
    
    # Create convergence plot
    mean_plot <- ggplot(chain_mean_imp_df, 
                        aes(x = iteration, y = mean, group = interaction(variable, Imputation), 
                            color = Imputation)) +
      geom_line() +
      facet_wrap(~ variable, scales = "free_y") + 
      scale_color_manual(" ", values = wes_palette(n = 5, name = "Darjeeling2")) +
      theme_minimal() +
      xlab("Iteration") +
      ylab("Iteration mean") +
      theme(legend.position = "top",
            panel.spacing = unit(0.25, "lines"),
            text = element_text(size = 15))
    
    return(list(mean = mean_plot))
    
  } else if (method == "var") {
    # Transform the array into a long format
    chain_var_imp_df <- as.data.frame(as.table(data))
    names(chain_var_imp_df) <- c("variable", "iteration", "Imputation", "variance") 
    # Rename Chain to Imputation
    chain_var_imp_df$Imputation <- gsub("Chain ", "Imputation ", chain_var_imp_df$Imputation)
    # Transform iteration to factor to preserve order
    chain_var_imp_df$iteration <- factor(chain_var_imp_df$iteration, ordered = TRUE,  levels = c(1,2,3,4,5,6,7,8,9,10))
    # Rename variables
    chain_var_imp_df <- chain_var_imp_df %>%
      mutate(variable = recode(variable, 
                               "fixed.acidity" = "Fixed Acidity (g/L)",
                               "volatile.acidity" = "Volatile Acidity (g/L)",
                               "citric.acid" = "Citric Acid (g/L)",
                               "residual.sugar" = "Residual Sugar (g/L)",
                               "chlorides" = "Chlorides (g/L)",
                               "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                               "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                               "density" = "Density (g/ml)",
                               "pH" = "pH",
                               "sulphates" = "Sulphates (g/L)",
                               "alcohol" = "Alcohol (% vol)"))
    # Create convergence plot
    var_plot <- ggplot(chain_var_imp_df, 
                       aes(x = iteration, y = variance, group = interaction(variable, Imputation), 
                           color = Imputation)) +
      geom_line() +
      facet_wrap(~ variable, scales = "free_y") + 
      scale_color_manual(" ", values = wes_palette(n = 5, name = "Darjeeling2")) +
      theme_minimal() +
      xlab("Iteration") +
      ylab("Iteration variance") +
      theme(legend.position = "top",
            panel.spacing = unit(0.25, "lines"),
            text = element_text(size = 15))
    
    return(list(var = var_plot))
  }
  else {
    "Unknown method"
  }
}

## Function to perform multiple imputation
multiple_imputation <- function(datasets_list, original_data, initialisation_list, method) {
  
  # Fit the original model to get the number of coefficients
  fit_original <- polr(factor(quality, ordered = TRUE) ~ ., data = original_data, method = "logistic", Hess = TRUE)
  original_coefs <- coef(fit_original)
  original_zetas <- fit_original$zeta
  original_coefs <- c(original_coefs, original_zetas)
  num_coefs <- length(original_coefs)

  # For storing coverage rate data
  coverage_counts <- matrix(0, nrow = length(datasets_list), ncol = num_coefs)
  coefficients_list <- matrix(NA, nrow = length(datasets_list), ncol = num_coefs)
  
  # For storing multiple imputation
  total_variance_list <- matrix(NA, nrow = length(datasets_list), ncol = num_coefs)
  ri_variance_list <- matrix(NA, nrow = length(datasets_list), ncol = num_coefs)
  fmi_list <- matrix(NA, nrow = length(datasets_list), ncol = num_coefs)
  
  # Initialize lists to store imputed datasets, MAE, RMSE, means, std dev, skewness
  imputed_datasets_list <- vector("list", length(datasets_list))
  mae_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  rmse_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  mean_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  sd_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  skewness_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  
  for (i in 1:length(datasets_list)) {
    data_mar <- datasets_list[[i]]
    data_init <- initialisation_list[[i]]
    
    ## Apply the specified imputation method using ifelse
    if (method == "mice") {
      ## Specify methods
      methods <- make.method(data_mar)
      methods['quality'] <- 'polr'
      methods[names(data_mar)[names(data_mar) != 'quality']] <- 'norm.nob'
      data_mar$quality <- factor(data_mar$quality, ordered = TRUE)
      
      ## Create mids object with 5 imputed datasets
      data_imp <- mice(data_mar, method = methods, m = 5, maxit = 10, visitSequence = "monotone", data_init = data_init, printFlag = FALSE)
      
      ## Impute data
      imputed_data <- complete(data_imp)
      
      
    } else if (method == "pmm") {
      ## Create mids object with 5 imputed datasets
      data_imp <- mice(data_mar, method = "pmm", m = 5, maxit = 10, visitSequence = "monotone", data_init = data_init, printFlag = FALSE)
      
      ## Impute data
      imputed_data <- complete(data_imp)
    }
    else {
      stop("Unknown imputation method")
    }
    
    ## Analyse mids object with polr and pool results
    analysis <- with(data_imp, polr(as.factor(quality) ~ ., data = complete(data_imp), method = "logistic", Hess = TRUE))
    ### Combine the results across all imputations and create a summary
    results1 <- pool(analysis)
    results2 <- summary(results1)
    ### Extract relevant metrics 
    estimates <- results1$pooled$estimate # Combined estimates
    total_variance <- results1$pooled$t # Total variance = ubar + b
    ri_variance <- results1$pooled$riv # Relative increase in variance due to imputation
    fmi <- results1$pooled$fmi # Fraction if missing information
    se <- results2$std.error
    
    ## Convert the ordinal variable 'quality' to numeric
    imputed_data$quality <- as.numeric(as.character(imputed_data$quality))
    
    ## Store the results
    imputed_datasets_list[[i]] <- imputed_data
    coefficients_list[i,] <- estimates
    total_variance_list[i,] <- total_variance
    ri_variance_list[i,] <- ri_variance
    fmi_list[i,] <- fmi
    
    ## Calculate average multiple imputation metrics
    average_total_variance <- colMeans(total_variance_list, na.rm = TRUE)
    average_ri_variance <- colMeans(ri_variance_list, na.rm = TRUE)
    average_fmi <- colMeans(fmi_list, na.rm = TRUE)
    
    ## Calculate Wald confidence intervals
    z_critical <- qnorm(0.975)
    lower_ci <- estimates - z_critical * se
    upper_ci <- estimates + z_critical * se
    
    ## Count if original coefficient is within the confidence interval of imputed model
    for (k in 1:length(original_coefs)) {
      if (original_coefs[k] >= lower_ci[k] && original_coefs[k] <= upper_ci[k]) {
        coverage_counts[i, k] <- 1
      } else {
        coverage_counts[i, k] <- 0
      }
    }
    
    ## Calcuate coverage rates and mean coefficients
    coverage_rates <- colMeans(coverage_counts)
    mean_coefficients <- colMeans(coefficients_list, na.rm = TRUE)
    
    ## Combine coverage rates with coefficient names
    quality_names <- c("level 1", "level 2", "level 3", "level 4", "level 5", "level 6") # Define names for the ordinal variable levesl
    names_tail <- tail(seq_along(names(original_coefs)), 6) # Specifiy names to be changed
    names(original_coefs)[names_tail] <- quality_names # Change the names
    
    coefficient_names <- names(original_coefs)
    coverage_rates_named <- tibble(variable = coefficient_names,
                                   original_coefficients = original_coefs,
                                   mean_coefficients = mean_coefficients,
                                   coverage_rates = coverage_rates)
    
    ## Calculate MAE and RMSE for each variable
    for (j in 1:ncol(original_data)) {
      missing_indices <- which(is.na(data_mar[[j]]))
      mae_list[i, j] <- mae(original_data[missing_indices, j], imputed_data[missing_indices, j])
      rmse_list[i, j] <- rmse(original_data[missing_indices, j], imputed_data[missing_indices, j])
    }
    
    ## Calculate mean, standard deviation, and skewness for each variable
    desc_stats <- psych::describe(imputed_data, type = 2)
    for (j in 1:ncol(original_data)) {
      mean_list[i,] <- desc_stats$mean
      sd_list[i,] <- desc_stats$sd
      skewness_list[i,] <- desc_stats$skew
    }
    
  }
  
  # Calculate average metrics per feature
  average_mae <- colMeans(mae_list, na.rm = TRUE)
  average_rmse <- colMeans(rmse_list, na.rm = TRUE)
  average_mean <- colMeans(mean_list, na.rm = TRUE)
  average_sd <- colMeans(sd_list, na.rm = TRUE)
  average_skewness <- colMeans(skewness_list, na.rm = TRUE)
  
  ## Create histograms for each variable using facet_wrap
  imputed_data_long <- bind_rows(imputed_datasets_list) %>% 
    pivot_longer(cols = everything(), names_to = "name", values_to = "value")
  
  histogram <- ggplot(data = imputed_data_long, aes(x = value)) +
    geom_histogram(alpha = 0.7, fill = "black", color = "white") +
    facet_wrap(~ name, scales = "free") +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.25, "lines"),
      text = element_text(size = 14)
    ) +
    xlab("") +
    ylab("Frequency")
  
  return(list(imputed_datasets = imputed_datasets_list, 
              imputed_datasets_long = imputed_data_long,
              average_mae = average_mae, 
              average_rmse = average_rmse,
              average_mean = average_mean,
              average_sd = average_sd,
              average_skewness = average_skewness,
              histogram = histogram,
              average_total_variance = average_total_variance,
              average_ri_variance = average_ri_variance,
              average_fmi = average_fmi,
              coverage_rates = coverage_rates_named,
              coefficients = coefficients_list
  ))
}

## Function to calculate MAE, RMSE, mean, sd and skewness
metrics <- function(datasets_list, original_data) {
  mae_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  rmse_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  mean_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  sd_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  skewness_list <- matrix(NA, nrow = length(datasets_list), ncol = ncol(original_data))
  correlations_list <- list()
  
  for (i in 1:length(datasets_list)) {
    imputed_data <- datasets_list[[i]]
    
    ## Calculate MAE and RMSE for each variable
    for (j in 1:ncol(original_data)) {
      missing_indices <- which(is.na(data_mar[[j]]))
      mae_list[i, j] <- mae(original_data[missing_indices, j], imputed_data[missing_indices, j])
      rmse_list[i, j] <- rmse(original_data[missing_indices, j], imputed_data[missing_indices, j])
    }
    
    ## Calculate mean, standard deviation, and skewness for each variable
    desc_stats <- psych::describe(imputed_data, type = 2)
    for (j in 1:ncol(original_data)) {
      mean_list[i, ] <- desc_stats$mean
      sd_list[i, ] <- desc_stats$sd
      skewness_list[i, ] <- desc_stats$skew
    }
    
    cor_matrix <- cor(imputed_data, use = "pairwise.complete.obs")
    correlations_list[[i]] <- cor_matrix
    
  }
  
  # Calculate average metrics per feature
  average_mae <- colMeans(mae_list, na.rm = TRUE)
  average_rmse <- colMeans(rmse_list, na.rm = TRUE)
  average_mean <- colMeans(mean_list, na.rm = TRUE)
  average_sd <- colMeans(sd_list, na.rm = TRUE)
  average_skewness <- colMeans(skewness_list, na.rm = TRUE)
  
  ## Create histograms for each variable using facet_wrap
  imputed_data_long <- bind_rows(datasets_list) %>% 
    pivot_longer(cols = everything(), names_to = "name", values_to = "value")
  
  histogram <- ggplot(data = imputed_data_long, aes(x = value)) +
    geom_histogram(alpha = 0.7, fill = "black", color = "white") +
    facet_wrap(~ name, scales = "free") +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.spacing = unit(0.25, "lines"),
      text = element_text(size = 14)
    ) +
    xlab("") +
    ylab("Frequency")
  
  return(list(mea_list = mae_list,
              rmse_list = rmse_list,
              mean_list = mean_list,
              sd_list = sd_list,
              skewness_list = skewness_list,
              correlations_list = correlations_list,
              average_mae = average_mae,
              average_rmse = average_rmse,
              average_mean = average_mean,
              average_sd = average_sd,
              average_skewness = average_skewness,
              histogram = histogram,
              datasets_long = imputed_data_long))
  
}

# Imputations----

## Mean imputation----
mean_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "mean")
saveRDS(mean_imp, "mean_imp.rds")
mean_imp <- readRDS("mean_imp.rds")

### Reverse Yeo-Johnson transformartion 
mean_imp_rev <- reverse_yj(mean_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(mean_imp_rev, "mean_imp_rev.rds")
mean_imp_rev <- readRDS("mean_imp_rev.rds")

### Calculate metrics
mean_imp_rev_met <- metrics(mean_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
mean_imp_mae <- mean_imp_rev_met$average_mae
mean_imp_rmse <- mean_imp_rev_met$average_rmse
mean_imp_mean <- mean_imp_rev_met$average_mean
mean_imp_sd <- mean_imp_rev_met$average_sd
mean_imp_skewness <- mean_imp_rev_met$average_skewness
mean_imp_cor <- mean_imp_rev_met$correlations_list
mean_imp_histogram <- mean_imp_rev_met$histogram

### Calculate coverage
mean_cr <- coverage(mean_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(mean_cr, "mean_cr.rds")
mean_cr <- readRDS("mean_cr.rds")

mean_cr_rev <- coverage(mean_imp_rev$reverse_imputed_datasets, data)
saveRDS(mean_cr_rev, "mean_cr_rev.rds")
mean_cr_rev <- readRDS("mean_cr_rev.rds")

### Calculate efficiency
mean_mark <- microbenchmark(mean_imp = {mean_imp <- mice(data_mar_yj, m = 1, maxit = 1, printFlag = FALSE)
mean_imp <- complete(mean_imp)})
mean_mark <- summary(mean_mark)
saveRDS(mean_mark, "mean_mark.rds")
mean_mark <- readRDS("mean_mark.rds")

mean_mark2 <- mark(mean_imp = {mean_imp <- mice(data_mar_yj, m = 1, maxit = 1, printFlag = FALSE)
mean_imp <- complete(mean_imp)}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(mean_mark2, "mean_mark2.rds")
mean_mark2 <- readRDS("mean_mark2.rds")

## Median----
median_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "median")
saveRDS(median_imp, "median_imp.rds")
median_imp <- readRDS("median_imp.rds")

### Reverse Yeo-Johnson transformation
median_imp_rev <- reverse_yj(mean_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(median_imp_rev, "median_imp_rev.rds")
median_imp_rev <- readRDS("median_imp_rev.rds")

### Calculate metrics
median_imp_rev_met <- metrics(median_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
median_imp_mae <- median_imp_rev_met$average_mae
median_imp_rmse <- median_imp_rev_met$average_rmse
median_imp_mean <- median_imp_rev_met$average_mean
median_imp_sd <- median_imp_rev_met$average_sd
median_imp_skewness <- median_imp_rev_met$average_skewness
median_imp_cor <- median_imp_rev_met$correlations_list
median_imp_histogram <- median_imp_rev_met$histogram

### Calculate coverage
median_cr <- coverage(median_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(median_cr, "median_cr.rds")
median_cr <- readRDS("median_cr.rds")

median_cr_rev <- coverage(median_imp_rev$reverse_imputed_datasets, data)
saveRDS(median_cr_rev, "median_cr_rev.rds")
median_cr_rev <- readRDS("median_cr_rev.rds")

### Calculate efficiency
median_mark <- microbenchmark(median_imp = {median_imp <- data_mar_yj %>%
  mutate_all(~ impute(., fun = median))})
median_mark <- summary(median_mark)
saveRDS(median_mark, "median_mark.rds")
median_mark <- readRDS("median_mark.rds")

median_mark2 <- mark(median_imp = {median_imp <- data_mar_yj %>%
  mutate_all(~ impute(., fun = median))}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(median_mark2, "median_mark2.rds")
median_mark2 <- readRDS("median_mark2.rds")

## Zero----
zero_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "zero")
saveRDS(zero_imp, "zero_imp.rds")
zero_imp <- readRDS("zero_imp.rds")

### Reverse Yeo-Johnson transformation
zero_imp_rev <- reverse_yj(zero_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(zero_imp_rev, "zero_imp_rev.rds")
zero_imp_rev <- readRDS("zero_imp_rev.rds")

### Calculate metrics
zero_imp_rev_met <- metrics(zero_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
zero_imp_mae <- zero_imp_rev_met$average_mae
zero_imp_rmse <- zero_imp_rev_met$average_rmse
zero_imp_mean <- zero_imp_rev_met$average_mean
zero_imp_sd <- zero_imp_rev_met$average_sd
zero_imp_skewness <- zero_imp_rev_met$average_skewness
zero_imp_cor <- zero_imp_rev_met$correlations_list
zero_imp_histogram <- zero_imp_rev_met$histogram

### Calculate coverage
zero_cr <- coverage(zero_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(zero_cr, "zero_cr.rds")
zero_cr <- readRDS("zero_cr.rds")

zero_cr_rev <- coverage(zero_imp_rev$reverse_imputed_datasets, data)
saveRDS(zero_cr_rev, "zero_cr_rev.rds")
zero_cr_rev <- readRDS("zero_cr_rev.rds")

### Calculate efficiency
zero_mark <- microbenchmark(zero_imo = {zero_imp <- data_mar_yj %>%
  mutate_all(~ impute(., fun = 0))})
zero_mark <- summary(zero_mark)
saveRDS(zero_mark, "zero_mark.rds")
zero_mark <- readRDS("zero_mark.rds")

zero_mark2 <- mark(zero_imo = {zero_imp <- data_mar_yj %>%
  mutate_all(~ impute(., fun = 0))}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(zero_mark2, "zero_mark2.rds")
zero_mark2 <- readRDS("zero_mark2.rds")

## Regression imputation----
reg_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "regression")
saveRDS(reg_imp, "reg_imp.rds")
reg_imp <- readRDS("reg_imp.rds")

### Reverse Yeo-Johnson transformation
reg_imp_rev <- reverse_yj(reg_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(reg_imp_rev, "reg_imp_rev.rds")
reg_imp_rev <- readRDS("reg_imp_rev.rds")

##Calculate metrics
reg_imp_rev_met <- metrics(reg_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
reg_imp_mae <- reg_imp_rev_met$average_mae
reg_imp_rmse <- reg_imp_rev_met$average_rmse
reg_imp_mean <- reg_imp_rev_met$average_mean
reg_imp_sd <- reg_imp_rev_met$average_sd
reg_imp_skewness <- reg_imp_rev_met$average_skewness
reg_imp_cor <- reg_imp_rev_met$correlations_list
reg_imp_histogram <- reg_imp_rev_met$histogram

### Calculate coverage
reg_cr <- coverage(reg_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(reg_cr, "reg_cr.rds")
reg_cr <- readRDS("reg_cr.rds")

reg_cr_rev <- coverage(reg_imp_rev$reverse_imputed_datasets, data)
saveRDS(reg_cr_rev, "reg_cr_rev.rds")
reg_cr_rev <- readRDS("reg_cr_rev.rds")

### Calculate efficiency
methods_reg <- make.method(data_mar_yj)
methods_reg['quality'] <- 'polr'
methods_reg[names(data_mar_yj)[names(data_mar_yj) != 'quality']] <- 'norm.predict'
data_mar_yj$quality <- factor(data_mar_yj$quality, ordered = TRUE)

reg_mark <- microbenchmark(reg_imp = {reg_imp <- mice(data_mar_yj, method = methods_reg, m = 1, maxit = 1)
reg_imp <- complete(reg_imp)})
reg_mark <- summary(reg_mark)
saveRDS(reg_mark, "reg_mark.rds")
reg_mark <- readRDS("reg_mark.rds")

reg_mark2 <- mark(reg_imp = {reg_imp <- mice(data_mar_yj, method = methods_reg, m = 1, maxit = 1)
reg_imp <- complete(reg_imp)}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(reg_mark2, "reg_mark2.rds")
reg_mark2 <- readRDS("reg_mark2.rds")

## Stochastic regression imputation----
set.seed(850382)
sreg_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "stochastic regression")
saveRDS(sreg_imp, "sreg_imp.rds")
sreg_imp <- readRDS("sreg_imp.rds")

### Reverse Yeo-Johnson transformation
sreg_imp_rev <- reverse_yj(sreg_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(sreg_imp_rev, "sreg_imp_rev.rds")
sreg_imp_rev <- readRDS("sreg_imp_rev.rds")

### Calculate metrics
sreg_imp_rev_met <- metrics(sreg_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
sreg_imp_mae <- sreg_imp_rev_met$average_mae
sreg_imp_rmse <- sreg_imp_rev_met$average_rmse
sreg_imp_mean <- sreg_imp_rev_met$average_mean
sreg_imp_sd <- sreg_imp_rev_met$average_sd
sreg_imp_skewness <- sreg_imp_rev_met$average_skewness
sreg_imp_cor <- sreg_imp_rev_met$correlations_list
sreg_imp_histogram <- sreg_imp_rev_met$histogram

### Calculate coverage
sreg_cr <- coverage(sreg_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(sreg_cr, "sreg_cr.rds")
sreg_cr <- readRDS("sreg_cr.rds")

### Exclude dataset 16 as the polr() can't determine starting values
exclude <- c(7)
sreg_cr_rev <- coverage(sreg_imp_rev$reverse_imputed_datasets[-exclude], data)
saveRDS(sreg_cr_rev, "sreg_cr_rev.rds")
sreg_cr_rev <- readRDS("sreg_cr_rev.rds")

### Calculate efficiency
methods_sreg <- make.method(data_mar_yj)
methods_sreg['quality'] <- 'polr'
methods_sreg[names(data_mar_yj)[names(data_mar_yj) != 'quality']] <- 'norm.nob'
data_mar_yj$quality <- factor(data_mar_yj$quality, ordered = TRUE)

sreg_mark <- microbenchmark(sreg_imp = {sreg_imp <- mice(data_mar_yj, method = methods_sreg, m = 1, maxit = 1)
sreg_imp <- complete(sreg_imp)})
sreg_mark <- summary(sreg_mark)
saveRDS(sreg_mark, "sreg_mark.rds")
sreg_mark <- readRDS("sreg_mark.rds")

sreg_mark2 <- mark(sreg_imp = {sreg_imp <- mice(data_mar_yj, method = methods_sreg, m = 1, maxit = 1)
sreg_imp <- complete(sreg_imp)}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(sreg_mark2, "sreg_mark2.rds")
sreg_mark2 <- readRDS("sreg_mark2.rds")

## Hot deck based on a matching metric----
set.seed(87737)
knn_imp <- imputation(datasets_mar_yj, as.data.frame(data_yj), method = "knn")
saveRDS(knn_imp, "knn_imp.rds")
knn_imp<- readRDS("knn_imp.rds")

### Reverse Yeo-Johnson transformation
knn_imp_rev <- reverse_yj(knn_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(knn_imp_rev, "knn_imp_rev.rds")
knn_imp_rev <- readRDS("knn_imp_rev.rds")

### Calculate metrics
knn_imp_rev_met <- metrics(knn_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
knn_imp_mae <- knn_imp_rev_met$average_mae
knn_imp_rmse <- knn_imp_rev_met$average_rmse
knn_imp_mean <- knn_imp_rev_met$average_mean
knn_imp_sd <- knn_imp_rev_met$average_sd
knn_imp_skewness <- knn_imp_rev_met$average_skewness
knn_imp_cor <- knn_imp_rev_met$correlations_list
knn_imp_histogram <- knn_imp_rev_met$histogram

### Calculate coverage
knn_cr <- coverage(knn_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(knn_cr, "knn_cr.rds")
knn_cr <- readRDS("knn_cr.rds")

knn_cr_rev <- coverage(knn_imp_rev$reverse_imputed_datasets, data)
saveRDS(knn_cr_rev, "knn_cr_rev.rds")
knn_cr_rev <- readRDS("knn_cr_rev.rds")

### Calculate efficiency
knn_mark <- microbenchmark(knn_imp = {knn_imp <- kNN(data_mar_yj, 
               variable = colnames(data_mar_yj),
               k = 5,
               dist_var = colnames(data_mar_yj),
               trace = TRUE,
               imp_var = FALSE,
               addRF = FALSE,
               onlyRF = FALSE,
               addRandom = FALSE,
               useImputedDist = FALSE)})
## Unit seconds
saveRDS(knn_mark,"knn_mark.rds")
knn_mark <- readRDS("knn_mark.rds")
knn_mark_sum <- summary(knn_mark)


knn_mark2 <- mark(knn_imp = {knn_imp <- kNN(data_mar_yj, 
                                             variable = colnames(data_mar_yj),
                                             k = 5,
                                             dist_var = colnames(data_mar_yj),
                                             trace = TRUE,
                                             imp_var = FALSE,
                                             addRF = FALSE,
                                             onlyRF = FALSE,
                                             addRandom = FALSE,
                                             useImputedDist = FALSE)}, 
                  iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(knn_mark2, "knn_mark2.rds")
knn_mark2 <- readRDS("knn_mark2.rds")

## Multiple Imputation by Chained Equations----
### Specifiy methods
mice_data <- data_mar_yj
methods_mice <- make.method(mice_data)
methods_mice['quality'] <- 'polr'
methods_mice[names(mice_data)[names(mice_data) != 'quality']] <- 'norm.nob'
mice_data$quality <- factor(mice_data$quality, ordered = TRUE)

### Specification of the number of iterations
mice_imp_it <- mice(mice_data, method = methods_mice, m = 5 , maxit = 10, printFlag = TRUE, visitSequence = "monotone")
mice_imp_it_init <- mice(mice_data, method = methods_mice, m = 5, maxit = 10, printFlag = TRUE, visitSequence = "monotone", data.init = as.data.frame(knn_imp$imputed_datasets[[1]]))

mean_imp_plot <- convergence_plot(mice_imp_it$chainMean, method = "mean", n = 10)
var_imp_plot <- convergence_plot(mice_imp_it$chainVar, method = "var", n = 10)

mean_imp_init_plot <- convergence_plot(mice_imp_it_init$chainMean, method = "mean", n = 10)
var_imp_init_plot <- convergence_plot(mice_imp_it_init$chainVar, method = "var", n = 10)

### Perform imputation
set.seed(988398)
mice_imp <- multiple_imputation(datasets_mar_yj, as.data.frame(data_yj), knn_imp$imputed_datasets, method = "mice")
saveRDS(mice_imp, "mice_imp.rds")
mice_imp <- readRDS("mice_imp.rds")

### Extract imputation metrics
mice_imp_total_variance <- mice_imp$average_total_variance
mice_imp_ri_variance <- mice_imp$average_ri_variance
mice_imp_fmi <- mice_imp$average_fmi
mice_imp$coverage_rates

### Reverse Yeo-Johnson transformation
mice_imp_rev <- reverse_yj(mice_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(mice_imp_rev, "mice_imp_rev.rds")
mice_imp_rev <- readRDS("mice_imp_rev.rds")

### Calculate metrics
mice_imp_rev_met <- metrics(mice_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
mice_imp_mae <- mice_imp_rev_met$average_mae
mice_imp_rmse <- mice_imp_rev_met$average_rmse
mice_imp_mean <- mice_imp_rev_met$average_mean
mice_imp_sd <- mice_imp_rev_met$average_sd
mice_imp_skewness <- mice_imp_rev_met$average_skewness
mice_imp_cor <- mice_imp_rev_met$correlations_list
mice_imp_histogram <- mice_imp_rev_met$histogram

### Calculate coverage
mice_cr <- coverage(mice_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(mice_cr, "mice_cr.rds")
mice_cr <- readRDS("mice_cr.rds")

### Exclude datasets 39 and 43 because polr() can't determine starting values
exclude_mice <- c(13, 41, 45)
mice_cr_rev <- coverage(mice_imp_rev$reverse_imputed_datasets[-exclude_mice], data)
saveRDS(mice_cr_rev, "mice_cr_rev.rds")
mice_cr_rev <- readRDS("mice_cr_rev.rds")

### Calculate efficiency
mice_mark <- microbenchmark(mice_imp = {
  mice_imp <- mice(mice_data, method = methods_mice, m = 5, maxit = 10, printFlag = TRUE, visitSequence = "monotone", data.init = as.data.frame(knn_imp$imputed_datasets[[1]]))
  mice_imp <- complete(mice_imp)
})
### Unit is seconds
saveRDS(mice_mark, "mice_mark.rds")
mice_mark <- readRDS("mice_mark.rds")
mice_mark_sum <- summary(mice_mark)

mice_mark2 <- mark(mice_imp = {mice_imp <- mice(mice_data, method = methods_mice, m = 5, maxit = 10, printFlag = TRUE, visitSequence = "monotone", data.init = as.data.frame(knn_imp$imputed_datasets[[1]]))
mice_imp <- complete(mice_imp)}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(mice_mark2, "mice_mark2.rds")
mice_mark2 <- readRDS("mice_mark2.rds")

## Predictive mean matching----
### Specification of the number of iterations
pmm_imp_it <- mice(data_mar_yj, method = "pmm", m = 5 , maxit = 10, printFlag = TRUE, visitSequence = "monotone")
pmm_imp_it_init <- mice(data_mar_yj, method = "pmm", m = 5 , maxit = 10, printFlag = TRUE, visitSequence = "monotone", data_init = knn_imp$imputed_datasets[[1]])

mean_imp_plot2 <- convergence_plot(pmm_imp_it$chainMean, method = "mean", n = 10) #16
var_imp_plot2 <- convergence_plot(pmm_imp_it$chainVar, method = "var", n = 10) # 18

mean_imp_init_plot2 <- convergence_plot(pmm_imp_it_init$chainMean, method = "mean", n = 10) #17
var_imp_init_plot2 <- convergence_plot(pmm_imp_it_init$chainVar, method = "var", n = 10) #19

### Perform imputation
set.seed(736272)
pmm_imp <- multiple_imputation(datasets_mar_yj, as.data.frame(data_yj), knn_imp$imputed_datasets, method = "pmm")
saveRDS(pmm_imp, "pmm_imp.rds")
pmm_imp <- readRDS("pmm_imp.rds")

### Extract imputation metrics
pmm_imp_total_variance <- pmm_imp$average_total_variance
pmm_imp_ri_variance <- pmm_imp$average_ri_variance
pmm_imp_fmi <- pmm_imp$average_fmi

### Reverse Yeo-Johnson transformation
pmm_imp_rev <- reverse_yj(pmm_imp$imputed_datasets, datasets_amputed_yj$trans)
saveRDS(pmm_imp_rev, "pmm_imp_rev.rds")
pmm_imp_rev <- readRDS("pmm_imp_rev.rds")

### Calculate metrics
pmm_imp_rev_met <- metrics(pmm_imp_rev$reverse_imputed_datasets, data)
### Extract metrics
pmm_imp_mae <- pmm_imp_rev_met$average_mae
pmm_imp_rmse <- pmm_imp_rev_met$average_rmse
pmm_imp_mean <- pmm_imp_rev_met$average_mean
pmm_imp_sd <- pmm_imp_rev_met$average_sd
pmm_imp_skewness <- pmm_imp_rev_met$average_skewness
pmm_imp_cor <- pmm_imp_rev_met$correlations_list
pmm_imp_histogram <- pmm_imp_rev_met$histogram

### Calculcate coverage
pmm_cr <- coverage(pmm_imp$imputed_datasets, as.data.frame(data_yj))
saveRDS(pmm_cr, "pmm_cr.rds")
pmm_cr <- readRDS("pmm_cr.rds")

pmm_cr_rev <- coverage(pmm_imp_rev$reverse_imputed_datasets, data)
saveRDS(pmm_cr_rev, "pmm_cr_rev.rds")
pmm_cr_rev <- readRDS("pmm_cr_rev.rds")

### Calculate efficiency
pmm_mark <- microbenchmark(pmm_imp = {
  pmm_imp <- mice(data_mar_yj, method = "pmm", m = 5, maxit = 10, printFlag = TRUE, visitSequence = "monotone", data.init = as.data.frame(knn_imp$imputed_datasets[[1]]))
  pmm_imp <- complete(mice_imp)
})
### Unit it seconds
saveRDS(pmm_mark, "pmm_mark.rds")
pmm_mark <- readRDS("pmm_mark.rds")
pmm_mark_sum <- summary(pmm_mark)

pmm_mark2 <- mark(pmm_imp = {
  pmm_imp <- mice(data_mar_yj, method = "pmm", m = 5, maxit = 10, printFlag = TRUE, visitSequence = "monotone", data.init = as.data.frame(knn_imp$imputed_datasets[[1]]))
  pmm_imp <- complete(mice_imp)}, iterations = 100, check = TRUE, memory = TRUE, filter_gc = TRUE)
saveRDS(pmm_mark2, "pmm_mark2.rds")
pmm_mark2 <- readRDS("pmm_mark2.rds")

## Diagnostics----
### Kolmogorov-Smirnov test
ks_tests <- function(original_data, imputed_data_list) {
  # Initialize a list to store the results
  ks_test_results <- list()
  
  # Loop through each imputed dataset
  for (i in 1:length(imputed_data_list)) {
    imputed_data <- imputed_data_list[[i]]
    
    # Loop through each feature
    for (feature in names(original_data)) {
      # Check if the feature exists in the imputed dataset
      if (feature %in% names(imputed_data)) {
        ks_test <- ks.test(original_data[[feature]], imputed_data[[feature]])
        
        # Store the result if p-value < 0.05
        if (ks_test$p.value < 0.05) {
          ks_test_results[[paste0("imputed_", i, "_", feature)]] <- data.frame(
            imputed_dataset = i,
            feature = feature,
            p_value = ks_test$p.value
          )
        }
      }
    }
  }
  
  # Combine all results into a single data frame
  if (length(ks_test_results) > 0) {
    ks_test_summary <- do.call(rbind, ks_test_results)
    return(ks_test_summary)
  } else {
    return(data.frame(message = "No significant differences found"))
  }
}

### Perform Kolmogorov-Smirnov tests
ks_mean <- ks_tests(data, mean_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_median <- ks_tests(data, median_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_zero <- ks_tests(data, zero_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_reg <- ks_tests(data, reg_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_sreg <- ks_tests(data, sreg_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_knn <- ks_tests(data, knn_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_mice <- ks_tests(data, mice_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))

ks_pmm <- ks_tests(data, pmm_imp_rev$reverse_imputed_datasets) %>%
  group_by(imputed_dataset) %>%
  count() %>%
  filter(n < ncol(data))


# Performance----
## Preparation
imputation_methods <- c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM")
variables <- colnames(data)
variables_unit <- c(
  "Fixed Acidity (g/L)",
  "Volatile Acidity (g/L)",
  "Citric Acid (g/L)",
  "Residual Sugar (g/L)",
  "Chlorides (g/L)",
  "Free Sulfur Dioxide (mg/L)",
  "Total Sulfur Dioxide (mg/L)",
  "Density (g/ml)",
  "pH",
  "Sulphates (g/L)",
  "Alcohol (% vol)",
  "Quality")

### Determine and extend palette
palette <- wes_palette("Darjeeling2")
extended_palette <- colorRampPalette(palette, bias = 1.5)(8)
extended_palette2 <- colorRampPalette(palette, bias = 1.5)(9)

## Accuracy----
### Mean absolute error----
accuracy_mae <- tibble(
  variable = variables,
  variable_unit = variables_unit,
  Mean = mean_imp_mae,
  Median = median_imp_mae,
  Zero = zero_imp_mae,
  Regression = reg_imp_mae,
  Stochastic_regression = sreg_imp_mae,
  kNN = knn_imp_mae,
  MICE = mice_imp_mae,
  PMM = pmm_imp_mae
)
saveRDS(accuracy_mae, "accuracy_mae.rds")
accuracy_mae <- readRDS("accuracy_mae.rds")

### Transform the data
accuracy_mae_long <- accuracy_mae %>%
  dplyr::select(-variable) %>%
  pivot_longer(-variable_unit, names_to = "method", values_to = "mae") %>%
  mutate(method = ifelse(method == "Stochastic_regression", "Stochastic reg.", method))
accuracy_mae_long$method <- factor(accuracy_mae_long$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_mae <- accuracy_mae_long %>%
  dplyr::group_by(variable_unit) %>%
  dplyr::summarize(
    max_method = method[which.max(mae)],
    max = max(mae, na.rm = TRUE),
    min_method = method[which.min(mae)],
    min = min(mae, na.rm = TRUE))
    
### Plot
mae_plot <- ggplot(accuracy_mae_long, aes(x = method, y = mae, color = method, group = variable_unit)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable_unit, scales = "free_y") +
  scale_color_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Method",
       y = "Mean Absolute Erorr")


### Root mean square error----
accuracy_rmse <- tibble(
  variable = variables,
  variable_unit = variables_unit,
  Mean = mean_imp_rmse,
  Median = median_imp_rmse,
  Zero = zero_imp_rmse,
  Regression = reg_imp_rmse,
  Stochastic_regression = sreg_imp_rmse,
  kNN = knn_imp_rmse,
  MICE = mice_imp_rmse,
  PMM = pmm_imp_mae
)
saveRDS(accuracy_rmse, "accuracy_rmse.rds")
accuracy_rmse <- readRDS("accuracy_rmse.rds")

### Transform the data
accuracy_rmse_long <- accuracy_rmse %>%
  dplyr::select(-variable) %>%
  pivot_longer(-variable_unit, names_to = "method", values_to = "mae") %>%
  mutate(method = ifelse(method == "Stochastic_regression", "Stochastic reg.", method))
accuracy_rmse_long$method <- factor(accuracy_rmse_long$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_rmse <- accuracy_rmse_long %>%
  dplyr::group_by(variable_unit) %>%
  dplyr::summarize(
    max_method = method[which.max(mae)],
    max = max(mae, na.rm = TRUE),
    min_method = method[which.min(mae)],
    min = min(mae, na.rm = TRUE))

### Plot
rmse_plot <- ggplot(accuracy_rmse_long, aes(x = method, y = mae, color = method, group = variable_unit)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable_unit, scales = "free_y") +
  scale_color_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Method",
       y = "Root Mean Square Erorr")

### Coefficients----
#### Extract original data coefficients
original_coefficients_rev <- mean_cr_rev$coverage_rates[,1:2]
original_coefficients_yj <- mean_cr$coverage_rates[,1:2]
#### Extract coefficients per method
coefficients <- function(data, original_data, cols, method) {
  data <- as.data.frame(data)
  n <- ncol(original_data)
  colnames(data) <- colnames(original_data[,1:n])
  data_long <- data %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    mutate(method = method)
  return(data_long)
}

mean_coef <- coefficients(mean_cr_rev$coefficients[,1:11], data[,1:11], method = "Mean")
mean_coef_yj <- coefficients(mean_cr$coefficients[,1:11], data[,1:11], method = "Mean")
median_coef <- coefficients(median_cr_rev$coefficients[,1:11], data[,1:11], method = "Median")
median_coef_yj <- coefficients(median_cr$coefficients[,1:11], data[,1:11], method = "Median")
zero_coef <- coefficients(zero_cr_rev$coefficients[,1:11], data[,1:11], method = "Zero")
zero_coef_yj <- coefficients(zero_cr$coefficients[,1:11], data[,1:11], method = "Zero")
reg_coef <- coefficients(reg_cr_rev$coefficients[,1:11], data[,1:11], method = "Regression")
reg_coef_yj <- coefficients(reg_cr$coefficients[,1:11], data[,1:11], method = "Regression")
sreg_coef <- coefficients(sreg_cr_rev$coefficients[,1:11], data[,1:11], method = "Stochastic reg.")
sreg_coef_yj <- coefficients(sreg_cr$coefficients[,1:11], data[,1:11], method = "Stochastic reg.")
knn_coef <- coefficients(knn_cr_rev$coefficients[,1:11], data[,1:11], method = "kNN")
knn_coef_yj <- coefficients(knn_cr$coefficients[,1:11], data[,1:11], method = "kNN")
mice_coef <- coefficients(mice_cr_rev$coefficients[,1:11], data[,1:11], method = "MICE")
mice_coef_yj <- coefficients(mice_cr$coefficients[,1:11], data[,1:11], method = "MICE")
pmm_coef <- coefficients(pmm_cr_rev$coefficients[,1:11], data[,1:11], method = "PMM")
pmm_coef_yj <- coefficients(pmm_cr$coefficients[,1:11], data[,1:11], method = "PMM")

### Joint coefficients and transform data into long format for plotting
joint_coef <- rbind(mean_coef, median_coef, zero_coef, reg_coef, sreg_coef, knn_coef, mice_coef, pmm_coef)
joint_coef_rev <- left_join(joint_coef, original_coefficients_rev) 
joint_coef_rev <- joint_coef_rev %>%
  mutate(variable = recode(variable, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)"))
joint_coef_rev$method <- factor(joint_coef_rev$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

coef_plot <- ggplot(joint_coef_rev, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_coefficients, color = "Original coefficients")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original coefficients are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 16)) +
  labs(x = "Regression coefficients",
       y = "Method")

joint_coef_yj <- rbind(mean_coef_yj, median_coef_yj, zero_coef_yj, reg_coef_yj, sreg_coef_yj, knn_coef_yj, mice_coef_yj, pmm_coef_yj)
joint_coef_yj <- left_join(joint_coef_yj, original_coefficients_yj) 
joint_coef_yj <- joint_coef_yj %>%
  mutate(variable = recode(variable, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)"))
joint_coef_yj$method <- factor(joint_coef_yj$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

coef_plot_yj <- ggplot(joint_coef_yj, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_coefficients, color = "Original coefficients")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original coefficients are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 17)) +
  labs(x = "Regression coefficients",
       y = "Method")

### Coverage rate----
accuracy_crs <- tibble(
  coefficient = mean_cr_rev$coverage_rates$variable,
  Mean = mean_cr_rev$coverage_rates$coverage_rate,
  Median = median_cr_rev$coverage_rates$coverage_rate,
  Zero = zero_cr_rev$coverage_rates$coverage_rate,
  Regression = c(reg_cr_rev$coverage_rates$coverage_rate, rep(NA,1)),
  Stochastic_regression = c(sreg_cr_rev$coverage_rates$coverage_rate, rep(NA,1)),
  kNN = c(knn_cr_rev$coverage_rates$coverage_rate,rep(NA,1)),
  MICE = c(mice_cr_rev$coverage_rates$coverage_rate,rep(NA,1)),
  PMM = c(pmm_cr_rev$coverage_rates$coverage_rate,rep(NA,1))
)
saveRDS(accuracy_crs, "accuracy_crs.rds")
accuracy_crs <- readRDS("accuracy_crs.rds")


accuracy_crs_yj <- tibble(
  coefficient = mean_cr$coverage_rates$variable,
  Mean = mean_cr$coverage_rates$coverage_rate,
  Median = c(median_cr$coverage_rates$coverage_rate, rep(NA,1)),
  Zero = zero_cr$coverage_rates$coverage_rate,
  Regression = c(reg_cr$coverage_rates$coverage_rate, rep(NA, 1)),
  Stochastic_regression = c(sreg_cr$coverage_rates$coverage_rate, rep(NA,1)),
  kNN = c(knn_cr$coverage_rates$coverage_rate, rep(NA,1)),
  MICE = c(mice_cr$coverage_rates$coverage_rate, rep(NA,1)),
  PMM = c(pmm_cr$coverage_rates$coverage_rate, rep(NA,1))
)
saveRDS(accuracy_crs_yj, "accuracy_crs_yj.rds")
accuracy_crs_yj <- readRDS("accuracy_crs_yj.rds")

accuracy_crs_long <- accuracy_crs %>%
  select(-head(11))
  head(11) %>%
  pivot_longer(-coefficient, names_to = "method", values_to = "coverage_rate") %>%
  mutate(method = ifelse(method == "Stochastic_regression", "Stochastic reg.", method),
         variable = recode(coefficient, 
                                  "fixed.acidity" = "Fixed Acidity (g/L)",
                                  "volatile.acidity" = "Volatile Acidity (g/L)",
                                  "citric.acid" = "Citric Acid (g/L)",
                                  "residual.sugar" = "Residual Sugar (g/L)",
                                  "chlorides" = "Chlorides (g/L)",
                                  "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                                  "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                                  "density" = "Density (g/ml)",
                                  "pH" = "pH",
                                  "sulphates" = "Sulphates (g/L)",
                                  "alcohol" = "Alcohol (% vol)"))
accuracy_crs_long$method <- factor(accuracy_crs_long$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_crs <- accuracy_crs_long %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(
    max_method = method[which.max(coverage_rate)],
    max = max(coverage_rate, na.rm = TRUE),
    min_method = method[which.min(coverage_rate)],
    min = min(coverage_rate, na.rm = TRUE))

accuracy_crs_yj_long <- accuracy_crs_yj %>%
  head(11) %>%
  pivot_longer(-coefficient, names_to = "method", values_to = "coverage_rate") %>%
  mutate(method = ifelse(method == "Stochastic_regression", "Stochastic reg.", method),
         variable = recode(coefficient, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)"))
accuracy_crs_yj_long$method <- factor(accuracy_crs_yj_long$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_crs_yj <- accuracy_crs_yj_long %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(
    max_method = method[which.max(coverage_rate)],
    max = max(coverage_rate, na.rm = TRUE),
    min_method = method[which.min(coverage_rate)],
    min = min(coverage_rate, na.rm = TRUE))

### Plot
crs_plot <- ggplot(accuracy_crs_long, aes(x = method, y = coverage_rate, color = method, group = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable) +
  scale_color_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 16)) +
  labs(x = "Method",
       y = "Coverage rate")

crs_yj_plot <- ggplot(accuracy_crs_yj_long, aes(x = method, y = coverage_rate, color = method, group = variable)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ variable) +
  scale_color_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 17)) +
  labs(x = "Method",
       y = "Coverage rate")

### Intercepts
intercept_yj <- tibble(
  Intercept = mean_cr$coverage_rates$variable,
  Original = mean_cr$coverage_rates$original_coefficients,
  Mean = mean_cr$coverage_rates$mean_coefficients,
  Median = c(median_cr$coverage_rates$mean_coefficients, rep(NA,1)),
  Zero = zero_cr$coverage_rates$mean_coefficients,
  Regression = c(reg_cr$coverage_rates$mean_coefficients, rep(NA, 1)),
  Stochastic_reg = c(sreg_cr$coverage_rates$mean_coefficients, rep(NA,1)),
  kNN = c(knn_cr$coverage_rates$mean_coefficients, rep(NA,1)),
  MICE = c(mice_cr$coverage_rates$mean_coefficients, rep(NA,1)),
  PMM = c(pmm_cr$coverage_rates$mean_coefficients, rep(NA,1))) %>%
  tail(7)
saveRDS(intercept_yj, "intercept_yj.RDS")
intercept_yj <- readRDS("intercept_yj.RDS")

## Preservation of the data characteristics----
### Extract original data metrics----
original_mean <- tibble(variable = variables,
                       original_mean = psych::describe(data, type = 2)$mean)
original_median <- tibble(variable = variables,
                         original_median = psych::describe(data, type = 2)$median)
original_sd <- tibble(variable = variables,
                      original_sd = psych::describe(data, type = 2)$sd)
original_skewness <- tibble(variable = variables,
                            original_skewness = psych::describe(data, type = 2)$skew)

### Means----
### Extract average imputation means
preserv_mean <- tibble(variable = variables,
                       variable_unit = variables_unit,
                       Original = psych::describe(data, type = 2)$mean,
                       Mean = mean_imp_mean,
                       Median = median_imp_mean,
                       Zero = zero_imp_mean,
                       Regression = reg_imp_mean,
                       Stochastic_regression = sreg_imp_mean,
                       kNN = knn_imp_mean,
                       MICE = mice_imp_mean,
                       PMM = pmm_imp_mean)
saveRDS(preserv_mean, "preserv_mean.rds")
preserv_mean <- readRDS("preserv_mean.rds")

### Extract mean lists
mean_mean <- coefficients(mean_imp_rev_met$mean_list, data, method = "Mean")
median_mean <- coefficients(median_imp_rev_met$mean_list, data, method = "Median")
zero_mean <- coefficients(zero_imp_rev_met$mean_list, data, method = "Zero")
reg_mean <- coefficients(reg_imp_rev_met$mean_list, data, method = "Regression")
sreg_mean <- coefficients(sreg_imp_rev_met$mean_list, data, method = "Stochastic reg.")
knn_mean <- coefficients(knn_imp_rev_met$mean_list, data, method = "kNN")
mice_mean <- coefficients(mice_imp_rev_met$mean_list, data, method = "MICE")
pmm_mean <- coefficients(pmm_imp_rev_met$mean_list, data, method = "PMM")

joint_mean <- rbind(mean_mean, median_mean, zero_mean, reg_mean, sreg_mean, knn_mean, mice_mean, pmm_mean)
joint_mean <- left_join(joint_mean, original_mean)
joint_mean <- joint_mean %>%
  mutate(variable = recode(variable, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)",
                           "quality" = "Quality"))
joint_mean$method <- factor(joint_mean$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_mean <- joint_mean %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(
    max_method = method[which.max(value)],
    max = max(value, na.rm = TRUE),
    min_method = method[which.min(value)],
    min = min(value, na.rm = TRUE))

### Mean plot
mean_plot <- ggplot(joint_mean, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_mean, color = "Original mean")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original means are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Means",
       y = "Method")

method <- "MICE"
variable <- "Residual Sugar (g/L)"
max_value <- max(joint_mean$value[joint_mean$method == method & joint_mean$variable == variable])
joint_mean2 <- joint_mean[!(joint_mean$method == method & joint_mean$variable == variable & joint_mean$value == max_value),]

mean_plot2 <- ggplot(joint_mean2, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_mean, color = "Original mean")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original means are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Means",
       y = "Method")

### Standard deviations----
### Extract average imputation means
preserv_sd <- tibble(variable = variables,
                     Original = psych::describe(data, type = 2)$sd,
                     Mean = mean_imp_sd,
                     Median = median_imp_sd,
                     Zero = zero_imp_sd,
                     Regression = reg_imp_sd,
                     Stochastic_regression = sreg_imp_sd,
                     kNN = knn_imp_sd,
                     MICE = mice_imp_sd,
                     PMM = pmm_imp_sd)
saveRDS(preserv_sd, "preserv_sd.rds")
preserv_sd <- readRDS("preserv_sd.rds")

### Extract sd lists
mean_sd <- coefficients(mean_imp_rev_met$sd_list, data, method = "Mean")
median_sd <- coefficients(median_imp_rev_met$sd_list, data, method = "Median")
zero_sd <- coefficients(zero_imp_rev_met$sd_list, data, method = "Zero")
reg_sd <- coefficients(reg_imp_rev_met$sd_list, data, method = "Regression")
sreg_sd <- coefficients(sreg_imp_rev_met$sd_list, data, method = "Stochastic reg.")
knn_sd <- coefficients(knn_imp_rev_met$sd_list, data, method = "kNN")
mice_sd <- coefficients(mice_imp_rev_met$sd_list, data, method = "MICE")
pmm_sd <- coefficients(pmm_imp_rev_met$sd_list, data, method = "PMM")

joint_sd <- rbind(mean_sd, median_sd, zero_sd, reg_sd, sreg_sd, knn_sd, mice_sd, pmm_sd)
joint_sd <- left_join(joint_sd, original_sd)
joint_sd <- joint_sd %>%
  mutate(variable = recode(variable, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)",
                           "quality" = "Quality"))
joint_sd$method <- factor(joint_sd$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

method_sd <- "MICE"
variable_sd <- "Residual Sugar (g/L)"
max_value_sd <- max(joint_sd$value[joint_sd$method == method & joint_sd$variable == variable])
joint_sd2 <- joint_sd[!(joint_sd$method == method & joint_sd$variable == variable & joint_sd$value == max_value),]

### SD plot
sd_plot <- ggplot(joint_sd, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_sd, color = "Original standard deviation")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original standard deviations are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Standard Deviations",
       y = "Method")

sd_plot2 <- ggplot(joint_sd2, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_sd, color = "Original standard deviation")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original standard deviations are displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Standard Deviations",
       y = "Method")

### Skewness----
### Extract average skewness
preserv_skew <- tibble(variable = variables,
                       Original = psych::describe(data, type = 2)$skew,
                       Mean = mean_imp_skewness,
                       Median = median_imp_skewness,
                       Zero = zero_imp_skewness,
                       Regression = reg_imp_skewness,
                       Stochastic_regression = sreg_imp_skewness,
                       kNN = knn_imp_skewness,
                       MICE = mice_imp_skewness,
                       PMM = pmm_imp_skewness)
saveRDS(preserv_skew, "preserv_skew.rds")
preserv_skew <- readRDS("preserv_skew.rds")

### Extract skewness lists
mean_skew <- coefficients(mean_imp_rev_met$skewness_list, data, method = "Mean")
median_skew <- coefficients(median_imp_rev_met$skewness_list, data, method = "Median")
zero_skew <- coefficients(zero_imp_rev_met$skewness_list, data, method = "Zero")
reg_skew <- coefficients(reg_imp_rev_met$skewness_list, data, method = "Regression")
sreg_skew <- coefficients(sreg_imp_rev_met$skewness_list, data, method = "Stochastic reg.")
knn_skew <- coefficients(knn_imp_rev_met$skewness_list, data, method = "kNN")
mice_skew <- coefficients(mice_imp_rev_met$skewness_list, data, method = "MICE")
pmm_skew <- coefficients(pmm_imp_rev_met$skewness_list, data, method = "PMM")

joint_skew <- rbind(mean_skew, median_skew, zero_skew, reg_skew, sreg_skew, knn_skew, mice_skew, pmm_skew)
joint_skew <- left_join(joint_skew, original_skewness)
joint_skew <- joint_skew %>%
  mutate(variable = recode(variable, 
                           "fixed.acidity" = "Fixed Acidity (g/L)",
                           "volatile.acidity" = "Volatile Acidity (g/L)",
                           "citric.acid" = "Citric Acid (g/L)",
                           "residual.sugar" = "Residual Sugar (g/L)",
                           "chlorides" = "Chlorides (g/L)",
                           "free.sulfur.dioxide" = "Free Sulfur Dioxide (mg/L)",
                           "total.sulfur.dioxide" = "Total Sulfur Dioxide (mg/L)",
                           "density" = "Density (g/ml)",
                           "pH" = "pH",
                           "sulphates" = "Sulphates (g/L)",
                           "alcohol" = "Alcohol (% vol)",
                           "quality" = "Quality"))
joint_skew$method <- factor(joint_skew$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

### SD plot
skew_plot <- ggplot(joint_skew, aes(x=value, y = method, group = variable)) + 
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_skewness, color = "Original skewness")) +
  facet_wrap(~variable, scale = "free_x") +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original skewness is displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(x = "Skewness",
       y = "Method")

### Correlations----
### Function to flatten the correlations matrix
flatten_matrix <- function(matrix) {
  matrix_df <- data.frame(as.table(matrix))
  return(matrix_df)
}

### Original data
cor_original <- flatten_matrix(cor_matrix)
colnames(cor_original) <- c("variable1", "variable2", "original_correlation_coefficient")
col_original <- cor_original %>%
  filter(variable1 != variable2) %>% # Exclude diagonal values
  arrange(desc(abs(original_correlation_coefficient))) %>% # Sort by correlation in descending order
  head(5) %>%
  distinct(original_correlation_coefficient, .keep_all = TRUE) %>%
  unite("pair", variable1, variable2, sep = " - ") %>%
  mutate(pair = recode(pair,
                           "density - residual.sugar" = "Density (g/ml) - Residual Sugar (g/L)",
                           "alcohol - density" = "Alcohol (% vol) - Density (g/ml)",
                           "total.sulfur.dioxide - free.sulfur.dioxide" = "Total Sulfur Dioxde (mg/L) - Free Sulfur Dioxide (mg/L)")) %>%
  as_tibble()

### Data preparation
mean_cor_df <- do.call(rbind, lapply(seq_along(mean_imp_cor), function(i) flatten_matrix(mean_imp_cor[[i]]))) %>%
  mutate(method = "Mean")

median_cor_df <- do.call(rbind, lapply(seq_along(median_imp_cor), function(i) flatten_matrix(median_imp_cor[[i]]))) %>%
  mutate(method = "Median")

zero_cor_df <- do.call(rbind, lapply(seq_along(zero_imp_cor), function(i) flatten_matrix(zero_imp_cor[[i]]))) %>%
  mutate(method = "Zero")

reg_cor_df <- do.call(rbind, lapply(seq_along(reg_imp_cor), function(i) flatten_matrix(reg_imp_cor[[i]]))) %>%
  mutate(method = "Regression")

sreg_cor_df <- do.call(rbind, lapply(seq_along(sreg_imp_cor), function(i) flatten_matrix(sreg_imp_cor[[i]]))) %>%
  mutate(method = "Stochastic reg.")

knn_cor_df <- do.call(rbind, lapply(seq_along(knn_imp_cor), function(i) flatten_matrix(knn_imp_cor[[i]]))) %>%
  mutate(method = "kNN")

mice_cor_df <- do.call(rbind, lapply(seq_along(mice_imp_cor), function(i) flatten_matrix(mice_imp_cor[[i]]))) %>%
  mutate(method = "MICE")

pmm_cor_df <- do.call(rbind, lapply(seq_along(pmm_imp_cor), function(i) flatten_matrix(pmm_imp_cor[[i]]))) %>%
  mutate(method = "PMM")

### Combine all correlation coefficients across all methods
joint_cor_df <- rbind(mean_cor_df, median_cor_df, zero_cor_df, reg_cor_df, sreg_cor_df, knn_cor_df, mice_cor_df, pmm_cor_df)
joint_cor_df <- cbind(joint_cor_df, cor_original[,3])
colnames(joint_cor_df) <- c("variable1", "variable2", "correlation_coefficient", "method", "original_correlation_coefficient")
### Filter out the correlation coefficients of the three pairs
joint_cor <- joint_cor_df %>%
  filter((variable1 == "alcohol" & variable2 == "density") |
           (variable1 == "density" & variable2 == "residual.sugar") |
           (variable1 == "total.sulfur.dioxide" & variable2 == "free.sulfur.dioxide")) %>%
  unite("variable", variable1, variable2, sep = " - ") %>%
  mutate(variable = recode(variable,
                           "density - residual.sugar" = "Density (g/ml) - Residual Sugar (g/L)",
                           "alcohol - density" = "Alcohol (% vol) - Density (g/ml)",
                           "total.sulfur.dioxide - free.sulfur.dioxide" = "Total Sulfur Dioxde (mg/L) - Free Sulfur Dioxide (mg/L)")) %>%
  as_tibble()
joint_cor$method <- factor(joint_cor$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

summary_cor <- joint_cor %>%
  dplyr::group_by(variable) %>%
  dplyr::summarize(
    max_method = method[which.max(correlation_coefficient)],
    max = max(correlation_coefficient, na.rm = TRUE),
    min_method = method[which.min(correlation_coefficient)],
    min = min(correlation_coefficient, na.rm = TRUE))


joint_cor_sum <- joint_cor %>%
  group_by(variable, method) %>%
  summarise(
    min = min(correlation_coefficient),
    max = max(correlation_coefficient),
    average = mean(correlation_coefficient)) %>%
    as_tibble() %>%
  print(n = 8 * 3)

### Plots
cor1_plot <- joint_cor %>%
  filter(variable == "Alcohol (% vol) - Density (g/ml)") %>%
  ggplot(aes(x = correlation_coefficient, y = method, group = variable)) +
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_correlation_coefficient, color = "Original correlation coefficient")) +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original correlation coefficient is displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(title = "Alcohol (% vol) - Density (g/ml)",
       x = "Correlation coefficients",
       y = "Method")

cor2_plot <- joint_cor %>%
  filter(variable == "Density (g/ml) - Residual Sugar (g/L)") %>%
  ggplot(aes(x = correlation_coefficient, y = method, group = variable)) +
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_correlation_coefficient, color = "Original correlation coefficient")) +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original correlation coefficient is displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(title = "Density (g/ml) - Residual Sugar (g/L)",
       x = "Correlation coefficients",
       y = "Method")

cor3_plot <- joint_cor %>%
  filter(variable == "Total Sulfur Dioxde (mg/L) - Free Sulfur Dioxide (mg/L)") %>%
  ggplot(aes(x = correlation_coefficient, y = method, group = variable)) +
  geom_point(show.legend = FALSE, aes(color = method)) +
  geom_line(show.legend = FALSE, aes(x=original_correlation_coefficient, color = "Original correlation coefficient")) +
  scale_color_manual(" ", values = extended_palette2) +
  labs(caption = "Original correlation coefficient is displayed as the vertical line") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 15)) +
  labs(title = "Total Sulfur Dioxde (mg/L) - Free Sulfur Dioxide (mg/L)",
       x = "Correlation coefficients",
       y = "Method")

## Computational efficiency----
### Data preparation
mean_time <- c(mean_mark$mean, median_mark$mean,zero_mark$mean,reg_mark$mean,sreg_mark$mean,knn_mark$mean,mice_mark$mean,pmm_mark$mean)
median_time <- c(mean_mark$median, median_mark$median,zero_mark$median,reg_mark$median,sreg_mark$median,knn_mark_sum$median*1000,mice_mark_sum$median*1000,pmm_mark_sum$median*1000)
median_time_mark <- c(mean_mark2$median, median_mark2$median,zero_mark2$median,reg_mark2$median,sreg_mark2$median,knn_mark2$median,mice_mark2$median,pmm_mark2$median)
memory <- c(mean_mark2$mem_alloc, median_mark2$mem_alloc, zero_mark2$mem_alloc,reg_mark2$mem_alloc,sreg_mark2$mem_alloc,knn_mark2$mem_alloc,mice_mark2$mem_alloc,pmm_mark2$mem_alloc)
formatted_memory <- sapply(memory, scales::label_bytes())

efficiency <- tibble(method = imputation_methods,
                     median_time = median_time,
                     median_time_mark = median_time_mark*1000,
                     memory = memory*10^-6,
                     formatted_memory = formatted_memory)
efficiency$method <- factor(efficiency$method, levels = c("Mean", "Median", "Zero", "Regression", "Stochastic reg.", "kNN", "MICE", "PMM"))

### Plots
time_plot <- ggplot(efficiency, aes(x = method, y = median_time_mark, fill = method)) +
  geom_bar(show.legend = FALSE,stat = "identity") +
  scale_fill_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 20)) +
  labs(x = "Method",
       y = "Median time in miliseconds")

memory_plot <- ggplot(efficiency, aes(x = method, y = memory, fill = method)) +
  geom_bar(show.legend = FALSE, stat = "identity") +
  scale_fill_manual(" ", values = extended_palette) +
  theme_minimal() + 
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 20)) +
  labs(x = "Method",
       y = "Memory in Megabytes")
