library(tidyverse)
library(DescTools)

table_s1_data <- read_csv("data-raw/table_s1_healthy.csv")
table_s2_data <- read_csv("data-raw/table_s2_healthy_and_day0vaccinee.csv")
table_s4_data <- read_csv("data-raw/table_s3_vaccinated_cohort.csv")

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
compute_spearman_correlations <- function(data, include_imp_h3n2 = F, include_imp_h2n2 = F){
  
  candidate_correlates <- c("imp_h1n1", "imp_group1", "age", "yob")
  
  if(include_imp_h3n2){
    candidate_correlates <- c(candidate_correlates, "imp_h3n2")
  }
  if(include_imp_h2n2){
    candidate_correlates <- c(candidate_correlates, "imp_h2n2")
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

test_fold_change_difs <- function(data){
  
  d56_fold_changes <- data %>%
    select(-matches("_fold"), -specimen) %>%
    pivot_longer(cols = !any_of(c("protocol","specimen","subject","yob","age","day")),
                 names_to = "antigen", 
                 values_to = "titer") %>%
    pivot_wider(names_from = day, values_from = titer, names_prefix = "d") %>%
    mutate(d56_log_fold_change = log(d56, base = 2) - log(d0, base = 2)) %>%
    mutate(group = case_when(
      yob > 1990 ~ "children",
      # There's a gap from 1983 to 1996, so no need to handle the case yob == 1990
      yob < 1990 ~ "non-children"
    )) %>%
    group_by(antigen)
  
  tests <- lapply(as.list(unique(d56_fold_changes$antigen)),
         FUN = function(antigen){
           data_subset <- d56_fold_changes %>%
             filter(antigen == !!antigen)
           
           children_values <- data_subset %>% filter(group == 'children') %>%
             pull(d56_log_fold_change)
           non_children_values <- data_subset %>% filter(group == 'non-children') %>%
             pull(d56_log_fold_change)
           
           if(all(is.na(children_values)) & all(is.na(non_children_values))){
             return(c())
           }else{
             wilcox_test <- wilcox.test(children_values, non_children_values, alternative = "greater")
             
             return(tibble(
               antigen = antigen,
               children_GMFC = 2^(mean(children_values, na.rm = T)),
               adult_GMFC = 2^(mean(non_children_values, na.rm = T)),
               W = wilcox_test$statistic,
               p_value = wilcox_test$p.value
               
             ))
           }
         }) %>%
    bind_rows()
  
  # Adjust p-values using the Holm correction
  tests <- tests %>%
    mutate(p_value = p.adjust(p_value, method = "holm"))
  
  return(tests)
    
}


label_antigens <- function(data){
  data %>%
    mutate(antigen = case_match(
      antigen,
      "h1_headless" ~ "H1 stalk",
      "h1_michigan" ~ "h1_michigan",
      "h1_neut" ~ "h1_neut",
      "h2_headless" ~ "H2 stalk",
      "h2_japan" ~ "h2_japan",
      "h3_colorado" ~ "h3_colorado",
      "h3_headless" ~ "H3 stalk",
      "h3_neut" ~ "h3_neut",
      "h5_cattle" ~ "h5_cattle",
      "h5_headless" ~ "H5 stalk",
      "h5_pheasant" ~ "Clade 2.3.4.4b H5",
      "h5_pheasant_neut" ~ "h5_pheasant_neut",
      "h5_vietnam" ~ "Clade 1 H5",
      "h5_vietnam_neut" ~ "h5_vietnam_neut",
      "h5_vietnam_neut_tg" ~ "h5_vietnam_neut_tg",
      "h7_headless" ~ "H7 stalk",
      "h7_neut" ~ "h7_neut",
      "h7_shanghai" ~ "h7_shanghai"
    ))
  
}
  
label_predictor <- function(data){
  data %>%
    mutate(predictor = case_match(
      predictor,
      "imp_h1n1" ~ "H1N1 imprinting",
      "imp_group1" ~ "Group 1 imprinting",
      "imp_h3n2" ~ "H3N2 imprinting",
      "imp_h2n2" ~ "H2N2 imprinting",
      "yob" ~ "Birth year",
      "age" ~ "Age"
    ))
}
  
  

format_table_s1 <- function(spearman_cors_table_s1){
  spearman_cors_table_s1 %>%
    filter(antigen %in% c("h1_headless", "h2_headless",
                          "h3_headless", "h5_headless",
                          "h7_headless", "h5_pheasant",
                          "h5_vietnam")) %>%
    label_antigens() %>%
    label_predictor() %>%
    select(-n_obs) %>%
    mutate(Correlation = paste(antigen, predictor, sep = ' titers - '),
           `Coefficient (95% CI)` = paste0(
             round(coef, 2), " (", round(lower,2), ",", round(upper,2), ")"),
           `P value` = pvalue) %>%
    select(Correlation, `Coefficient (95% CI)`, `P value`)
    
}



# For the blood bank data set, compute Spearman correlations between titers to
# each antigen and probability of imprinting with group1 and probability of
# imprinting with H1N1
spearman_cors_table_s1 <- compute_spearman_correlations(
  table_s1_data %>%
    filter(age >= 18),
  include_imp_h3n2 = T,
  include_imp_h2n2 = T) %>%
    filter(predictor != 'age', predictor != 'yob') %>%
  format_table_s1()


write_csv(spearman_cors_table_s1, "table_s1_spearman_corrs.csv")



# Run the models for the combined vaccination/blood bank dataset,
# export table with the results
model_comparison <- sapply(unique(table_s2_data$antigen), FUN = fit_linear_models,
       data = table_s2_data, simplify = F) %>% bind_rows() %>%
  rename(k = df)

model_comparison <- model_comparison %>%
  label_antigens() %>%
  mutate(model = case_match(
    model,
    "yob" ~ "Birth year",
    "group1_imprinting" ~ "Group 1 imprinting",
    "age" ~ "Age",
    "h1n1_imprinting" ~ "H1N1 imprinting"
  )) %>%
  mutate(across(any_of(c('AIC','loglik','delta_AIC')), ~round(.x,2)))

write_csv(model_comparison, "table_s2_model_comparison.csv")

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

bootstrap_results <- bootstrap_results %>%
  label_antigens() %>%
  rename(predictor = predictor1) %>%
  label_predictor() %>%
  rename(predictor1 = predictor, predictor = predictor2) %>%
  label_predictor() %>%
  rename(predictor2 = predictor) %>%
  arrange(antigen) %>%
  mutate(across(where(is.numeric), ~round(.x,3))) %>%
  mutate(`Bootstrap difference (95% CI)` = 
           paste0(
             bootstrap_diff_mean, 
             " (",
             bootstrap_diff_lower,
             ",",
             bootstrap_diff_upper,
             ")"
           )) %>%
  rename(Antigen = antigen,
         `Predictor 1` = predictor1,
         `Predictor 2` = predictor2,
         `Titer correlation with predictor 1` = cor_predictor_1_obs,
         `Titer correlation with predictor 2` = cor_predictor_2_obs,
         `Observed difference` = cor_diff_obs) %>%
  select(Antigen, `Predictor 1`, `Predictor 2`, `Titer correlation with predictor 1`,
         `Titer correlation with predictor 2`, `Observed difference`,
         `Bootstrap difference (95% CI)`)

write_csv(bootstrap_results,  "table_s3_bootstrap_results.csv")

# Testing day 56 fold-change differences between kids and adults
fold_change_diffs_table_s4 <- test_fold_change_difs(table_s4_data) 
write_csv(fold_change_diffs_table_s4, "table_s4_fold_change_diffs.csv")
