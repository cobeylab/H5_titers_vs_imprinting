library(tidyverse)
library(DescTools)

table_s1_data <- read_csv("data-raw/table_s1_healthy.csv")
table_s2_data <- read_csv("data-raw/table_s2_healthy_and_day0vaccinee.csv")

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
compute_spearman_correlations <- function(data){
  data %>%
    group_by(antigen) %>%
    summarise(across(c("imp_h1n1", "imp_group1", "age", "yob"), function(x){cor.test(x, titer, method = 'spearman')$estimate},
                     .names = "coef_{.col}"),
              across(c("imp_h1n1", "imp_group1", "age", "yob"), function(x){cor.test(x, titer, method = 'spearman')$p.value},
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

# For the blood bank data set, compute Spearman correlations between titers to
# each antigen and probability of imprinting with group1 and probability of
# imprinting with H1N1
spearman_cors_table_s1 <- compute_spearman_correlations(table_s1_data) %>%
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





