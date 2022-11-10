# Regression models for article entitled
# Vaccination against COVID-19 in Geneva, Switzerland, and willingness
# of parents to get their children vaccinated: a cross-sectional population-based study
# Nick Pullen Nov 2022 Unité d'Epidémiologie Populationnelle, Hôpitaux Universitaires de Genève

library(tidyverse)
library(nnet)#for multinomial

dat = readRDS("PATH TO DATA.rds")

# Age and sex adjusted odds ratios for association of factors with COVID-19 vaccination status
calc_simpleORs = function(dataframe, ses_var, dropNA_in_count=FALSE) {
  # Exceptions for specific vars
  if (rlang::as_string(ensym(ses_var))=="hh_income_cat") {dataframe = dataframe %>% filter(hh_income_cat!="Ne sais pas, Ne souhaite pas répondre")}
  if (rlang::as_string(ensym(ses_var))=="i_education") {dataframe = dataframe %>% filter(i_education!="Autre")}
  tmp = dataframe %>%
    filter(i_sex != "Intersexe") %>%
    select(dose_yesno, age_at_Q_submission, i_sex, {{ses_var}}) %>%
    drop_na() %>%
    droplevels() %>%
    mutate(vac01 = if_else(dose_yesno == "Oui", 1L, 0L)) %>%
    glm(paste("vac01 ~ ", rlang::as_string(ensym(ses_var)), "+ age_at_Q_submission + i_sex"), data = ., family = "binomial")

  cbind(OR = exp(coef(tmp)), exp(confint(tmp)), pvalue = coef(summary(tmp))[,'Pr(>|z|)'], N_obs = nobs(tmp)) %>%
    as_tibble(rownames = "model_levels") %>%
    mutate(Variable = rlang::as_string(ensym(ses_var)), .before = 1) %>%
    mutate(across(3:5, round, digits=2), # ORs and CIs
           "Simple OR (95% CI)" = paste0(OR, " (", `2.5 %`, "-", `97.5 %`, ")")) %>%
    relocate(pvalue, .after = last_col()) %>%
    relocate(N_obs, .after = last_col())
}

simpleORs = map_dfr(quos(i_sex, age_at_Q_submission, i_education, hh_income_cat,
                         i_work_situation_recoded, i_hh_livewith,
                         i_chronic_disease, smoking_rec),
                    calc_simpleORs,
                    dataframe = dat)


# Multivariable (Age, sex, education, income) adjusted odds ratios for association of factors with COVID-19 vaccination status
calc_multivar_ORs = function(dataframe, ses_var, dropNA_in_count=FALSE) {
  tmp = dataframe %>%
    filter(i_sex != "Intersexe") %>%
    filter(hh_income_cat != "Ne sais pas, Ne souhaite pas répondre") %>%
    filter(i_education != "Autre") %>%
    select(dose_yesno, age_at_Q_submission, i_sex, i_education, hh_income_cat, {{ses_var}}) %>%
    drop_na() %>%
    droplevels() %>%
    mutate(vac01 = if_else(dose_yesno == "Oui", 1L, 0L)) %>%
    glm(paste("vac01 ~ ", rlang::as_string(ensym(ses_var)), "+ age_at_Q_submission + i_sex + i_education + hh_income_cat"), data = ., family = "binomial")

  cbind(OR = exp(coef(tmp)), exp(confint(tmp)), pvalue = coef(summary(tmp))[,'Pr(>|z|)'], N_obs = nobs(tmp)) %>%
    as_tibble(rownames = "model_levels") %>%
    mutate(Variable = rlang::as_string(ensym(ses_var)), .before = 1) %>%
    mutate(across(3:5, round, digits=2), # ORs and CIs
           "Multivariable OR (95% CI)" = paste0(OR, " (", `2.5 %`, "-", `97.5 %`, ")")) %>%
    relocate(pvalue, .after = last_col()) %>%
    relocate(N_obs, .after = last_col())
}

multivar_ORs = map_dfr(quos(i_sex, age_at_Q_submission, i_education, hh_income_cat,
                            i_work_situation_recoded, i_hh_livewith,
                            i_chronic_disease, smoking_rec),
                       calc_multivar_ORs,
                       dataframe = dat)


# Multinomial regression to assess relationship between variables (age, sex, education,
# income, child’s age group) and willingness of parents to get their children vaccinated
my_multinom = multinom(outcome_kids_3 ~ age_at_Q_submission + i_sex + i_education + hh_income_cat + bin_kids_factor,
                       data = dat %>% filter(i_sex != "Intersexe", hh_income_cat != "Ne sais pas, Ne souhaite pas répondre", i_education != "Autre", bin_kids_factor!="No kids") %>% droplevels())
# get number of observations (Complete case analysis) used by multinom # https://stackoverflow.com/a/47696331/3275826
nrow(residuals(my_multinom))
# or
model.frame(outcome_kids_3 ~ age_at_Q_submission + i_sex + i_education + hh_income_cat + bin_kids_factor,
            data = dat %>%
              filter(i_sex != "Intersexe",
                     hh_income_cat != "Ne sais pas, Ne souhaite pas répondre",
                     i_education != "Autre", bin_kids_factor!="No kids") %>%
              droplevels(), na.action = na.omit) %>%
  nrow
#summary(my_multinom)
z = summary(my_multinom)$coefficients/summary(my_multinom)$standard.errors

# 2-tailed z test
p = (1 - pnorm(abs(z), 0, 1)) * 2

## extract the coefficients from the model and exponentiate
rawORs = exp(coef(my_multinom)) %>% t %>%  as_tibble(rownames = "variable")
OR_CIs = exp(confint(my_multinom)) %>% as_tibble(rownames = "variable")
p_vals = as_tibble(t(p), rownames = "variable") %>%
  rename(Mixed_pval = Mixed, `Always no pval` = `Always no`)
multi_ORs_CIs = inner_join(rawORs, OR_CIs) %>%
  inner_join(p_vals) %>%
  mutate(Mixed_CI = paste0(round(`2.5 %.Mixed`, 2), "-", round(`97.5 %.Mixed`, 2)),
         `Always no CI` = paste0(round(`2.5 %.Always no`, 2), "-", round(`97.5 %.Always no`, 2))) %>%
  select(variable, Mixed, Mixed_CI, Mixed_pval, `Always no`, `Always no CI`, `Always no pval`)


# print session Info for reproducibility
sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
#
# Matrix products: default
#
# locale:
#   [1] LC_COLLATE=French_Switzerland.1252  LC_CTYPE=French_Switzerland.1252    LC_MONETARY=French_Switzerland.1252 LC_NUMERIC=C
# [5] LC_TIME=French_Switzerland.1252
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] nnet_7.3-17     forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.6    ggplot2_3.4.0   tidyverse_1.3.2
