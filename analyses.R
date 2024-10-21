library(tidyverse)
library(DescTools)

table_s1_data <- read_csv("data-raw/table_s1_healthy.csv")
table_s2_data <- read_csv("data-raw/table_s2_healthy_and_day0vaccinee.csv")

nrep <- 1000 # Number of bootstrap replicates


# Given the data, and for a single antigen, fits separate linear models predicting
# log titers based on age, year of birth, group1 imprinting probability, and H1N1
# imprinting probability
fit_linear_models <- function(data, antigen){
  
  data <- data %>% mutate(log2_titer = log(titer, base = 2)) %>%
    filter(antigen == !! antigen)
  
  age_model <- lm(log2_titer ~ age, data = data)
  yob_model <- lm(log2_titer ~ yob, data = data)
  group1_imprinting_model <- lm(log2_titer ~ imp_group1, data = data)
  h1n1_imprinting_model <- lm(log2_titer ~ imp_h1n1, data = data)
  
  AIC_table <- AIC(age_model, yob_model, group1_imprinting_model, h1n1_imprinting_model) %>%
    as_tibble(rownames = "model") %>%
    arrange(AIC) %>%
    mutate(loglik = sapply(model,
                             FUN = function(x){
                               logLik(get(x)) %>% as.numeric()
                             },
                           simplify = T, USE.NAMES = F)) %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    mutate(antigen = antigen) %>%
    select(antigen, everything()) %>%
    mutate(model = str_remove(model, "_model"))
  
  return(AIC_table)
}

# Compute Spearman correlations between titers to each antigen and various predictors
compute_spearman_correlations <- function(data, include_imp_h3n2 = F){
  
  candidate_correlates <- c("imp_h1n1", "imp_group1", "age", "yob")
  
  if(include_imp_h3n2){
    candidate_correlates <- c(candidate_correlates, "imp_h3n2")
  }
  
  
  data %>%
    group_by(antigen) %>%
    summarise(across(candidate_correlates, function(x){cor.test(x, titer, method = 'spearman')$estimate},
                     .names = "coef_{.col}"),
              across(candidate_correlates, function(x){cor.test(x, titer, method = 'spearman')$p.value},
                     .names = "pvalue_{.col}"),
              n_obs = n()
    ) %>%
    pivot_longer(cols = !any_of(c("antigen","n_obs") )) %>%
    mutate(value_type = str_extract(name, "[^_]+")) %>%
    mutate(name = str_remove(name, paste0(value_type,"_"))) %>%
    rename(predictor = name) %>%
    pivot_wider(names_from = value_type, values_from = value) %>%
    rowwise() %>%
    mutate(lower = CorCI(coef, n_obs)['lwr.ci'],
           upper = CorCI(coef, n_obs)['upr.ci']) %>%
    ungroup()
}

# Takes difference in Spearman correlation between predictors 1 and 2,#
# tests their significance using bootstrapping
run_bootstrap_correlation_test <- function(data, predictor1, predictor2, nrep){
  
  # Internal function: given table of Spearman correlation results, 
  # retrieves comparison of interest, puts in wide format
  get_pw_comparison <- function(spearman_corrs, predictor1, predictor2){
    spearman_corrs %>%
      filter(predictor %in% c(predictor1, predictor2)) %>%
      select(antigen, predictor, coef) %>%
      mutate(predictor = case_match(predictor,
                                    predictor1 ~ 1,
                                    predictor2 ~ 2 )) %>%
      pivot_wider(names_from = predictor, values_from = coef, names_prefix = "cor_predictor_") %>%
      # We care about the absolute value of the correlation coefficient
      mutate(across(!matches("antigen"), abs)) %>%
      select(antigen, cor_predictor_1, cor_predictor_2) %>%
      mutate(cor_diff = cor_predictor_1 - cor_predictor_2)
  }
  
  obs_corr_diff <- compute_spearman_correlations(data) %>%
   get_pw_comparison(predictor1 = predictor1, predictor2 = predictor2) %>%
    rename_with(.fn = function(x){paste0(x,"_obs")}, .cols = !matches("antigen"))
  
  bootstrap_corr_diff <- c()
  
  for(i in 1:nrep){
    for(antigen in unique(data$antigen)){
      
      row_indices <- (1:nrow(data))[data$antigen == antigen]
      
      resampled_rows <- sample(row_indices, replace = T)
      
      resampled_data <- data[resampled_rows,]
      
      
      resampled_corr <- compute_spearman_correlations(resampled_data) %>%
        get_pw_comparison(predictor1 = predictor1, predictor2 = predictor2)
      
      bootstrap_corr_diff <- bootstrap_corr_diff %>%
        bind_rows(resampled_corr)
      
    }

  }
  
  
  bootstrap_results <- bootstrap_corr_diff %>% 
    left_join(obs_corr_diff) %>%
    group_by(antigen, cor_predictor_1_obs, cor_predictor_2_obs, cor_diff_obs) %>%
    summarise(bootstrap_diff_mean = mean(cor_diff),
              bootstrap_diff_lower = quantile(cor_diff, 0.025),
              bootstrap_diff_upper = quantile(cor_diff, 0.975)) %>%
    mutate(predictor1 = predictor1, predictor2 = predictor2) %>%
    select(antigen, predictor1, predictor2, everything())
  
  return(bootstrap_results)
  
}



# For the blood bank data set, compute Spearman correlations between titers to
# each antigen and probability of imprinting with group1 and probability of
# imprinting with H1N1
spearman_cors_table_s1 <- compute_spearman_correlations(table_s1_data, include_imp_h3n2 = T) %>%
    filter(predictor != 'age', predictor != 'yob')

write_csv(spearman_cors_table_s1, "table_s1_spearman_corrs.csv")



# Run the models for the combined vaccination/blood bank dataset,
# export table with the results
model_comparison <- sapply(unique(table_s2_data$antigen), FUN = fit_linear_models,
       data = table_s2_data, simplify = F) %>% bind_rows() %>%
  rename(k = df)

# Annotate with Spearman correlations
spearman_cors_table_s2 <- compute_spearman_correlations(table_s2_data) %>%
  mutate(predictor = case_match(predictor,
                                "imp_h1n1" ~ "h1n1_imprinting",
                                "imp_group1" ~ "group1_imprinting",
                                "age" ~ "age",
                                "yob" ~ "yob")) %>%
  rename(spearman_coef = coef,
         spearman_p = pvalue,
         spearman_lower = lower,
         spearman_upper = upper)

model_comparison <- left_join(model_comparison,
                              spearman_cors_table_s2 %>% rename(model = predictor))

write_csv(model_comparison, "model_comparison.csv")

# As a complement to linear models, run bootstrap tests

# Is year of birth more strongly associated with titers than age is?
bootstrap_yob_vs_age <- run_bootstrap_correlation_test(table_s2_data,
                                                       predictor1 = "yob",
                                                       predictor2 = "age",
                                                       nrep = nrep)

# Is group1 imprinting probability more strongly associated with titers than age is?
bootstrap_group1_vs_age <- run_bootstrap_correlation_test(table_s2_data,
                                                            predictor1 = "imp_group1",
                                                            predictor2 = "age",
                                                            nrep = nrep)
# Does H1N1 imprinting beat age? 
bootstrap_h1N1_vs_age <- run_bootstrap_correlation_test(table_s2_data,
                                                        predictor1 = "imp_h1n1",
                                                        predictor2 = "age",
                                                        nrep = nrep)


# Does group 1 imprinting beat YOB?
bootstrap_group1_vs_yob <- run_bootstrap_correlation_test(table_s2_data,
                                                          predictor1 = "imp_group1",
                                                          predictor2 = "yob",
                                                          nrep = nrep)

bootstrap_results <- 
  bootstrap_yob_vs_age %>%
  bind_rows(bootstrap_group1_vs_age) %>%
  bind_rows(bootstrap_h1N1_vs_age) %>%
  bind_rows(bootstrap_group1_vs_yob)

write_csv(bootstrap_results,  "bootstrap_results.csv")




