# Setup ================
# Load libraries and raw data

library(readxl)
library(tidyverse)

el <- read_excel("Data/Banco de dados _ Idoso.xlsx")

meds <- read_excel("Data/Mini-73.xls")

# Data cleaning ================
# Rename and recode variables, create factors, handle missing values, etc.

colnames(el) <- gsub(" ", "_", colnames(el))
colnames(el) <- gsub("-", "_", colnames(el))
colnames(el) <- gsub("PANAS_POS", "PANASPOS", colnames(el))
colnames(el) <- gsub("PANAS_NEG", "PANASNEG", colnames(el))

# Exclude empty row at the end of the dataframe
el <- el[1:108, ]

el <- el %>%
  # Exclude variables which are part of MINI
  select(-starts_with("Risc_Sui"), -starts_with("Dados_epi")) %>% 
  rename(subject = "...1",
         GDS8_8 = GDS48_8) %>%
  mutate(
    condition = factor(Grupo...11, 
                       levels = c(0, 1),
                       labels = c("Sham TBS", "Active TBS")
    ),
    condition_rev = fct_rev(condition),
    across(starts_with("Avaliado_por"), as.factor),
    treatment_resistant_dep = as.factor(mini_69), 
    across(c(CGI2_0:G_CIRS12_0, BDNF_0:BDNF_6), as.numeric)) %>% 
  select(-c(Nome_completo, Iniciais, RGHC, Nascimento, E_mail, Telefone, Cód),
         subject, condition, everything())


el$treatment_resistant_dep[c(23, 61, 79)] <- 0

# Join dataframes
el <- el %>%
  full_join(meds %>%
              select(subject = ID, mini_73_total), 
            by = "subject")

# Transform to long format
outcomes <- el %>%
  select(
    c(starts_with("HDRS", ignore.case = F), 
      starts_with("GDS", ignore.case = F), 
      starts_with("CGI", ignore.case = F), 
      starts_with("PANAS", ignore.case = F), 
      starts_with("MADRS", ignore.case = F), 
      starts_with("YMRS", ignore.case = F),
      starts_with("TOTAL_TMA", ignore.case = F),
      starts_with("BDNF", ignore.case = F))
  ) %>%
  colnames

el_lg <- el %>%
  pivot_longer(cols = c(all_of(outcomes)),
               names_to = c(".value", "time"),
               names_pattern = "(.*)_([0-9]{1,2})",
               names_transform = list(time = as.integer),
  ) %>%
mutate(
  # Rename sum scores and recalculate
  hdrs_total_orig = HDRS,
  madrs_total_orig = MADRS,
  gds_total_orig = GDS,
  ymrs_total_orig = YMRS,
  panas_neg_total_orig = PANASPOS,
  panas_neg_total_orig = PANASNEG,
  time_log = log(time + 1), 
  time_log10 = log10(time + 1), 
  subject = as.factor(subject),
  # Harmonize evaluator names
  examiner = as.factor(
    recode(Avaliado_por,
           "Dr. Leonardo" = "Dr Leonardo Afonso dos Santos",
           "Leonardo Afonso dos Santos" = "Dr Leonardo Afonso dos Santos",
           "Dr Leonardo Afonso" = "Dr Leonardo Afonso dos Santos",
           "Dr Rafael Benatti" = "Dr Rafael Garcia Benatti",
           "Dr. Rafael" = "Dr Rafael Garcia Benatti",
           "Dr. Rafael Benatti" = "Dr Rafael Garcia Benatti",
           "Rafael Benatti" = "Dr Rafael Garcia Benatti",
           "rafael garcia benatti" = "Dr Rafael Garcia Benatti",
           "Rafael Garcia Benatti" = "Dr Rafael Garcia Benatti",
           "leandro" = "Dr Leandro Valiengo",
           "Leandro" = "Dr Leandro Valiengo",
           "LEANDRO" = "Dr Leandro Valiengo",
           "Dr. Leandro Valiengo" = "Dr Leandro Valiengo",
           "Dra. Jessica" = "Dra Jéssica Formiga e Silva",
           "Dra. Jéssica Formiga" = "Dra Jéssica Formiga e Silva",
           "Jéssica Formiga e Silva" = "Dra Jéssica Formiga e Silva",
           "Jéssica Formiga e Sulva" = "Dra Jéssica Formiga e Silva",
           "Dra Julia Loureiro" = "Dra Julia Cunha Loureiro",
           "Julia" = "Dra Julia Cunha Loureiro",
           "Julia Cunha Loureiro" = "Dra Julia Cunha Loureiro",
           "Julia Loureiro" = "Dra Julia Cunha Loureiro",
           "Dra. Lisiane" = "Dra Lisiane Martins",
           "Lisiane Martins" = "Dra Lisiane Martins",
           "Dra Maíra Lessa" = "Dra Maíra Lessa",
           "Dra. Maíra" = "Dra Maíra Lessa",
           "Maíra Pinheiro Maux Lessa" = "Dra Maíra Lessa",
           "Luara" = "Luara Tort",
           "Luara Tort" = "Luara Tort",
           "Renata" = "Renata Rocha Vaughan",
           "Renata rocha" = "Renata Rocha Vaughan",
           "Renata Rocha" = "Renata Rocha Vaughan",
           "Renata Rocha Vaughan" = "Renata Rocha Vaughan",
           "Roberta" = "Roberta Mattar",
           "Roberta de A M P F D Mattar" = "Roberta Mattar",
           "Roberta de Arruda Mendes Pereira Fiuze Dini Mattar" = "Roberta Mattar",
           "Roberta Mattar" = "Roberta Mattar",
           "Bianca e Bruna" = "Bianca",
           "Henriette" = "Henriette Baena"))
) %>%
  rowwise() %>%
  # create indicator variable whether wheather a scale is missing completely
  mutate(
    is_na_hdrs = sum(is.na(c_across(matches("^HDRS[0-9]+")))) == 17,
    is_na_madrs  = sum(is.na(c_across(matches("^MADRS[0-9]+")))) == 10,
    is_na_gds  = sum(is.na(c_across(matches("^GDS[0-9]+")))) == 15,
    is_na_panas  = sum(is.na(c_across(matches("^PANAS[0-9]+")))) == 20,
    is_na_ymrs  = sum(is.na(c_across(matches("^YMRS[0-9]+")))) == 11
  ) %>%
  ungroup() %>%
  select(-c(HDRS, MADRS, GDS, YMRS, PANASPOS, PANASNEG)) %>%
  # Re-arrange column order
  select(subject, condition, time, time_log, examiner,
         ends_with("_total"), starts_with("is_na_"),
         starts_with("HDRS", ignore.case = F), 
         starts_with("MADRS", ignore.case = F), 
         starts_with("GDS", ignore.case = F), 
         starts_with("CGI", ignore.case = F), 
         starts_with("PANAS", ignore.case = F), 
         starts_with("YMRS", ignore.case = F), 
         everything()) %>%
  mutate(
    gender = as.factor(case_when(
      Genero == 1 ~ "male",
      Genero == 2 ~ "female")),
    white_ethnicity = as.factor(ifelse(etnia == 1, "yes", "no")),
    unemployed = as.factor(ifelse(ocupacao == 2, "yes", "no")),
    wage = as.factor(ifelse(renda == 1, "low", "high")),
    relationship = as.factor(ifelse(estado_civil == 2, "yes", "no")),
    diabetes  = as.factor(ifelse(diabete == 1, "yes", "no")),
    physical_activity = as.factor(ifelse(atividade_fisica == 1, "low", "high")),
    in_psychotherapy = as.factor(mini_76 == 1),
    mini_74 = as.factor(mini_74),
    mini_79 = as.factor(mini_79),
    Benzodiazepines = if_else(
      is.na(mini_80) & is.na(medicacoes_atual), NA_real_,if_else(
        grepl("alprazolam|rivotril|lorazepam|clonazepam|bromazepam|dormonid", tolower(mini_80)) |
          grepl("lprazolam|rivotril|lorazepam|clonazepam|bromazepam|dormonid", tolower(medicacoes_atual)),
        1, 0), )
  ) %>%
  select(-c(Genero, Grupo...11, Grupo...900, diabete)) %>%
  relocate(gender:physical_activity, .after = Idade)

items_with_missings <- el_lg %>%
  # Count number of missings unless all scale items are missing
  summarize(
    across(matches("^HDRS[0-9]+"),
           ~sum(ifelse(is_na_hdrs, 0, is.na(.x)))),
    across(matches("^MADRS[0-9]+"),
           ~sum(ifelse(is_na_madrs, 0, is.na(.x)))),
    across(matches("^GDS[0-9]+"),
           ~sum(ifelse(is_na_gds, 0, is.na(.x)))),
    across(matches("^PANAS[0-9]+"),
           ~sum(ifelse(is_na_panas, 0, is.na(.x)))),
    across(matches("^YMRS[0-9]+"),
           ~sum(ifelse(is_na_ymrs, 0, is.na(.x))))
  ) %>%
  # Select items with at least one missing 
  select(where(~any(.x != 0))) %>%
  names()

# There is only one item with a missing value: subject 93, time 2, HDRS17.
# Mean imputation

el_lg <- el_lg %>%
  group_by(time) %>%
  mutate(
    across(all_of(items_with_missings) & starts_with("HDRS"),
           ~ifelse(is_na_hdrs, NA, replace_na(.x, round(mean(.x, na.rm=TRUE))))
           )
  ) %>%
  ungroup() %>%
  mutate(
    hdrs_total = HDRS1 + HDRS2 + HDRS3 + HDRS4 +
      HDRS5 + HDRS6 + HDRS7 + HDRS8 + HDRS9 +
      HDRS10 + HDRS11 + HDRS12 + HDRS13 + HDRS14 +
      HDRS15 + HDRS16 + HDRS17,
    madrs_total = MADRS1 + MADRS2 + MADRS3 + MADRS4 + MADRS5 +
      MADRS6 + MADRS7 + MADRS8 + MADRS9 + MADRS10,
    gds_total = 1-GDS1 + GDS2 + GDS3 + GDS4 + 1-GDS5 + GDS6 + 1-GDS7 + GDS8
    + GDS9 + GDS10 + 1-GDS11 + GDS12 + 1-GDS13 + GDS14 + GDS15,
    ymrs_total = YMRS1 + YMRS2 + YMRS3 + YMRS4 + YMRS5 + YMRS6 + YMRS7
    + YMRS8 + YMRS9 + YMRS10 + YMRS11,
    panas_pos_total = PANAS1 + PANAS3 + PANAS5 + PANAS9 + PANAS10
    + PANAS12 + PANAS14 + PANAS16 + PANAS17 + PANAS19,
    panas_neg_total = PANAS2 + PANAS4 + PANAS6 + PANAS7 + PANAS8
    + PANAS11 + PANAS13 + PANAS15 + PANAS18 + PANAS20
  ) %>%
  arrange(subject, time) %>%
  mutate(across(c(madrs_total, hdrs_total, examiner),
                ~ case_when(time == 0 ~ .x,
                            time == 1 ~ lag(.x, n = 1),
                            time == 2 ~ lag(.x, n = 2),
                            time == 4 ~ lag(.x, n = 3),
                            time == 6 ~ lag(.x, n = 4),
                            time == 8 ~ lag(.x, n = 5),
                            time == 12 ~ lag(.x, n = 6)),
                .names = "{col}_bl")) %>%
  ungroup() %>%
  mutate(hdrs_response = as.numeric(hdrs_total <= 0.5 * hdrs_total_bl),
         madrs_response = as.numeric(madrs_total <= 0.5 * madrs_total_bl),
         hdrs_remission = as.numeric(hdrs_total <= 7),
         madrs_remission = as.numeric(madrs_total <= 10))

dss <- el_lg

# Create wide format ===============
# Pivot long to wide format 

dss_wide <- dss %>%
  select(-time_log, -time_log10, -starts_with("is_na_")) %>%
  pivot_wider(names_from = time,
              values_from = c(starts_with("HDRS"),
                              starts_with("GDS"), 
                              starts_with("CGI"), 
                              starts_with("PANAS"), 
                              starts_with("MADRS"), 
                              starts_with("YMRS"),
                              starts_with("TOTAL_TMA"),
                              BDNF),
              names_sep = "_") %>%
  mutate(
    # diff is improvement = positive!
    diff_hdrs_2 = hdrs_total_0 - hdrs_total_2,
    diff_hdrs_6 = hdrs_total_0 - hdrs_total_6,
    diff_hdrs_8 = hdrs_total_0 - hdrs_total_8,
    diff_hdrs_12 = hdrs_total_0 - hdrs_total_12) %>%
  # count number of failed AD trials
  mutate(
    mini_73 = str_replace(str_replace(mini_73, ",$", ""), " e ", ","),
    failed_AD_trials = case_when(
      mini_73 == "-" ~ 0,
      mini_73 == "0" ~ 0,
      grepl("^não", tolower(mini_73)) ~ 1,
      grepl(",", mini_73) ~ str_count(mini_73, ",") + 1,
      mini_73 == "1" ~ 1,
      is.na(mini_73) ~ NA,
      mini_73 == "Teve melhora nos episódios de depressão prévia." ~ 0,
      TRUE ~ 1
    ),
    failed_AD_trials = ifelse(
      mini_73 == "todas",
      max(failed_AD_trials, na.rm=TRUE),
      failed_AD_trials
    )
  ) 

# Add new variables to long format as well
dss <- dss %>%
  full_join(
    dss_wide %>%
      select(subject, diff_hdrs_2:diff_hdrs_12, failed_AD_trials), 
    by = "subject")


# Save dataframes ------------

save(dss_wide, file = "Data/elderly_cleaned_wide.RData")

save(dss, file = "Data/elderly_cleaned_long.RData")
  
# Optionally create excel file for sharing: 
# library(writexl)
# write_xlsx(dss_wide, path = "Data/Idosos_wide.xlsx")


# Prepare machine learning dataset -----
# Dataset with only  baseline and outcome variables, or additionally data from
# weeks 1, 2

dss <- dss_wide %>%
  select(
    diff_hdrs_6, diff_hdrs_12, hdrs_total_6, hdrs_total_12, 
    hdrs_response_6, hdrs_response_12, hdrs_remission_6, hdrs_remission_12,
    everything(),
     -(ends_with("_1") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_2") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_4") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_6") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_6 & !hdrs_total_6 & !hdrs_response_6 & !hdrs_remission_6),
    -(ends_with("_8") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_12") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_12 & !hdrs_total_12 & !hdrs_response_12 & !hdrs_remission_12),
    -c(condition_rev, Nome_completo, Iniciais, RGHC, E_mail, Telefone)
  ) %>%
  filter(!is.na(hdrs_response_12))


save(dss, file = "elderly_ML.RData")

dss_w1 <- dss_wide %>%
  select(
    diff_hdrs_6, diff_hdrs_12, hdrs_total_6, hdrs_total_12, 
    hdrs_response_6, hdrs_response_12, hdrs_remission_6, hdrs_remission_12,
    everything(),
    # -(ends_with("_1") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_2") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_4") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_6") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_6 & !hdrs_total_6 & !hdrs_response_6 & !hdrs_remission_6),
    -(ends_with("_8") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_12") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_12 & !hdrs_total_12 & !hdrs_response_12 & !hdrs_remission_12),
    -c(condition_rev, Nome_completo, Iniciais, RGHC, E_mail, Telefone)
  ) %>%
  filter(!is.na(hdrs_response_12))

save(dss_w1, file = "Data/elderly_ML_inclWeek1.RData")

dss_w2 <- dss_wide %>%
  select(
    diff_hdrs_6, diff_hdrs_12, hdrs_total_6, hdrs_total_12, 
    hdrs_response_6, hdrs_response_12, hdrs_remission_6, hdrs_remission_12,
    everything(),
    # -(ends_with("_1") & !contains("EDINBURGH") & !contains("MINI")),
    # -(ends_with("_2") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_4") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_6") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_6 & !hdrs_total_6 & !hdrs_response_6 & !hdrs_remission_6),
    -(ends_with("_8") & !contains("EDINBURGH") & !contains("MINI")),
    -(ends_with("_12") & !contains("EDINBURGH") & !contains("MINI") & 
        !diff_hdrs_12 & !hdrs_total_12 & !hdrs_response_12 & !hdrs_remission_12),
    -c(condition_rev, Nome_completo, Iniciais, RGHC, E_mail, Telefone)
  ) %>%
  filter(!is.na(hdrs_response_12))

save(dss_w2, file = "Data/elderly_ML_inclWeek1And2.RData")
