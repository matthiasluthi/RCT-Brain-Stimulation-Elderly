# Setup ------------
# Load libraries and data
library(MASS)
library(tidyverse)
library(DescTools)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(emmeans)
library(flextable)
library(ftExtra)
library(ggpubr)
library(here)
library(infer)
library(ARTool)
library(gridExtra)
library(lemon)
library(effectsize)

load("Data/Idosos_Excitability_August2024.RData")

# Load custom functions
source("0-custom-functions.R")

select <- dplyr::select

## Variable definitions ------------
# Define a named vector mapping values with underscores to reader-friendly versions with spaces
rename_mapping <- c(
  "hdrs_total" = "HDRS-17",
  "Week_6" = "Week 6",
  "Week_12" = "Week 12",
  "Difference baseline - week 6" = "Difference baseline - week 6",
  "Difference baseline - week 12" = "Difference baseline - week 12",
  "SP_Limiar_direita" = "Right RMT",
  "SP_Limiar_esquerda" = "Left RMT",
  "PP_PEM_direita_0.5-1.5mV" = "Right MEP",
  "PP_PEM_esquerda_0.5-1.5mV" = "Left MEP",
  "PP_Intervalo_direita_2.4.6mseg" = "Right SICI",
  "PP_Intervalo_esquerda_2.4.6mseg" = "Left SICI",
  "PP_Intervalo_direita_10.15.20mseg" = "Right ICF",
  "PP_Intervalo_esquerda_10.15.20mseg" = "Left ICF"
)

mce_measures <- c(
  # Threshold for muscle response
  "SP_Limiar_esquerda", "SP_Limiar_direita", 
  # Motor evoked potential
  "PP_PEM_esquerda_0.5-1.5mV", "PP_PEM_direita_0.5-1.5mV",
  # Inhibition
  "PP_Intervalo_esquerda_2.4.6mseg", "PP_Intervalo_direita_2.4.6mseg",
  # Facilitation
  "PP_Intervalo_esquerda_10.15.20mseg", "PP_Intervalo_direita_10.15.20mseg"

  )

# _iqr: removal of outliers according to Tukey’s fences
mce_measures_iqr <- sapply(mce_measures, function(x) paste0(x, "_iqr"))
names(mce_measures_iqr) <- NULL

mce_measures_ext <- c(
  "SP_Limiar_esquerda", "SP_Limiar_direita", 
  "PP_PEM_esquerda_0.5-1.5mV", "PP_PEM_direita_0.5-1.5mV",
  "SP_Limiar_esquerda_iqr", "SP_Limiar_direita_iqr", 
  "PP_PEM_esquerda_0.5-1.5mV_iqr", "PP_PEM_direita_0.5-1.5mV_iqr",
  
  # Inhibition
  "PP_Intervalo_esquerda_2.4.6mseg_iqr", "PP_Intervalo_direita_2.4.6mseg_iqr",
  "PP_Intervalo_esquerda_2.4.6mseg_log", "PP_Intervalo_direita_2.4.6mseg_log",
  "PP_Intervalo_esquerda_2.4mseg_iqr", "PP_Intervalo_direita_2.4mseg_iqr",
  "PP_Intervalo_esquerda_2.4mseg_log", "PP_Intervalo_direita_2.4mseg_log",
  
  # Facilitation
  "PP_Intervalo_esquerda_10.15.20mseg_iqr", "PP_Intervalo_direita_10.15.20mseg_iqr",
  "PP_Intervalo_esquerda_10.15.20mseg_log", "PP_Intervalo_direita_10.15.20mseg_log",
  "PP_Intervalo_esquerda_10.15mseg_iqr", "PP_Intervalo_direita_10.15mseg_iqr",
  "PP_Intervalo_esquerda_10.15mseg_log", "PP_Intervalo_direita_10.15mseg_log"
)


all_mce_measures <- xct %>%
  select(ends_with("_iqr"), ends_with("_log"),
         -(PP_Intervalo_direita_1mseg_iqr:PP_Intervalo_direita_10mseg_iqr),
         -(PP_Intervalo_esquerda_1mseg_iqr:PP_Intervalo_esquerda_10mseg_iqr),
         -time_log
         ) %>%
  names()


## Additional DFs -----

xct0 <- filter(xct, time == 0)
xct_tms <- xct %>% filter(time == 0|time == 6|time == 12)

# Outliers & distributions ------------
## Raw/original data --------------
# Conclusion: They are not normally distributed. 
# Only thresholds look okay and pass the shapiro test (but barely)
# MEP: both directions
# Inhibition/facilitation: outliers to the right only

raw_histograms <- list()

# Histograms
for (var in mce_measures) {
  var_label <- rename_mapping[var]
  raw_histograms[[var]] <- 
    ggplot(xct, aes(.data[[var]])) +
      geom_histogram() + 
      xlab(var_label) +
      theme_pub2()
}


# Boxplots
raw_boxplots <- list()

for (var in mce_measures) {
  var_label <- rename_mapping[var]
  raw_boxplots[[var]] <- 
    ggplot(xct, aes(.data[[var]])) +
    geom_boxplot(coef = 1.5) + 
    xlab(paste(var_label, "- IQR coef: 1.5")) +
    theme_pub2()
  show(raw_boxplots[[var]])
}

raw_normality_tests <- vector()
# Normality test
for (var in mce_measures) {
  print(var)
  show(shapiro.test(xct[[var]]))
  raw_normality_tests <- c(raw_normality_tests, 
                           shapiro.test(xct[[var]])$p.value)
}

## Inter-quartile range (IQR) ---------------------

iqr_histograms <- list()
# Histograms
for (var in mce_measures_iqr) {
  var_label <- rename_mapping[sub("_iqr", "", var)]
  iqr_histograms[[var]] <- 
    ggplot(xct, aes(.data[[var]])) +
    geom_histogram() + 
    xlab(var_label) +
    # xlab(paste(var_label, "- IQR", IQR_SCALE)) +
    theme_pub2()
  # show(iqr_histograms[[var]])
}

ggarrange(
  iqr_histograms[[1]], iqr_histograms[[2]], 
  iqr_histograms[[3]], iqr_histograms[[4]],
  iqr_histograms[[5]], iqr_histograms[[6]],
  iqr_histograms[[7]], iqr_histograms[[8]],
  nrow=4, ncol = 2
)

iqr_normality_tests <- vector()
# Normality test
for (var in mce_measures_iqr) {
  print(var)
  show(shapiro.test(xct[[var]]))
  iqr_normality_tests <- c(iqr_normality_tests, 
                               shapiro.test(xct[[var]])$p.value)
}

## Raw vs. IQR comparisons ------------------------

# Table for normality tests
table_normality <- data.frame(
  MCE = rename_mapping[mce_measures],
  `P original data` = round(raw_normality_tests, 3),
  `P cleaned data` = round(iqr_normality_tests, 3),
  check.names = FALSE
) 

ft_normality <- table_normality %>%
  mutate(
    `P original data` = case_when(
      `P original data` == 1 ~ ">0.99",
      `P original data` == 0 ~ "<0.001",
      TRUE ~ as.character(round(`P original data`, 3))),
    `P cleaned data` = case_when(
      `P cleaned data` == 1 ~ ">0.99",
      `P cleaned data` == 0 ~ "<0.001",
      TRUE ~ as.character(round(`P cleaned data`, 3)))
  ) %>%
  flextable() %>%
  align(j = c(2, 3), align = "right", part = "all") %>%
  bold(i = table_normality$`P original data` < .05, j = 2) %>%
  bold(i = table_normality$`P cleaned data` < .05, j = 3) %>%
  fit_pagewidth()

ft_normality

save_as_docx(ft_normality, 
             path = here("Output", "Excitability",
                         "Table_Excitability_NormalityTests.docx"))

# Graphs for comparisons

ggarrange(
  raw_histograms[[1]] + xlim(c(33, 72)) +  scale_x_continuous(labels = ~paste0(., "%")) + ggtitle("Raw"), 
  iqr_histograms[[1]] + xlim(c(33, 72)) +  scale_x_continuous(labels = ~paste0(., "%")) + ylab(NULL) + ggtitle("Cleaned"),
  raw_histograms[[2]] + xlim(c(33, 72)) +  scale_x_continuous(labels = ~paste0(., "%")), 
  iqr_histograms[[2]] + xlim(c(33, 72)) +  scale_x_continuous(labels = ~paste0(., "%")) + ylab(NULL),
  
  raw_histograms[[3]] + xlim(c(0, 3000)), 
  iqr_histograms[[3]] + xlim(c(0, 3000)) + ylab(NULL),
  raw_histograms[[4]] + xlim(c(0, 1400)), 
  iqr_histograms[[4]] + xlim(c(0, 1400)) + ylab(NULL),
  
  ncol=2, nrow=4
  )
ggsave(here("Output", "Excitability", "Histograms_Group1.pdf"))

ggarrange(
  raw_histograms[[5]] + xlim(c(0, 8700)) + ggtitle("Raw"), 
  iqr_histograms[[5]] + xlim(c(0, 8700)) + ylab(NULL) + ggtitle("Cleaned"),
  raw_histograms[[6]] + xlim(c(0, 4300)), 
  iqr_histograms[[6]] + xlim(c(0, 4300)) + ylab(NULL),
  
  raw_histograms[[7]] + xlim(c(0, 3400)), 
  iqr_histograms[[7]] + xlim(c(0, 3400)) + ylab(NULL),
  raw_histograms[[8]] + xlim(c(0, 5600)), 
  iqr_histograms[[8]] + xlim(c(0, 5600)) + ylab(NULL),
  
  ncol=2, nrow=4
)

ggsave(here("Output", "Excitability", "Histograms_Group2.pdf"))


# Line-plots ------------
## HDRS-17 --------

xct %>%
  group_by(time, condition_rev) %>%  
  summarise(across(hdrs_total,
                   list(Mean = ~mean(.x, na.rm = TRUE), 
                        SD = ~sd(.x, na.rm = TRUE), 
                        SEM = ~sd(.x, na.rm = TRUE)/sqrt(sum(!is.na(.x))),
                        N = ~sum(!is.na(.x))
                   ),
                   .names = "{col}_{fn}")
  ) %>%
  ungroup() %>% 
  ggplot(aes(x = time, y = hdrs_total_Mean,
             color = condition_rev, shape = condition_rev)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3.5) +
    geom_errorbar(aes(
      ymin = hdrs_total_Mean - 1.96*hdrs_total_SEM,
      ymax = hdrs_total_Mean + 1.96*hdrs_total_SEM),
      width = 0.8,
      size = 0.75,
      alpha = 0.8
    ) +
    scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 8, 12),
                       labels = sapply(c("Base- line", 1, 2, 4, 6, 8, 12),
                                        function(x) str_wrap(x, width = 5))) +
                       # labels = c("Baseline", 1, 2, 4, 6, 8, 12)) +
    expand_limits(y = 0) +
    #scale_color_brewer(palette = "Dark2", direction = 1) +
    labs(y = "HDRS-17", x = "Week") +
    theme_pub2()

ggsave(here("Output", "Excitability", "Lineplot_HDRS.pdf"))

## TMS indices as medians & IQR --------
summaries <- xct_tms %>%
  # filter(time %in% c(0, 6, 12)) %>%
  group_by(time, condition) %>%  
  summarise(across(all_of(mce_measures_iqr),
                   list(Mean = ~mean(.x, na.rm = TRUE), 
                        SD = ~sd(.x, na.rm = TRUE), 
                        SEM = ~sd(.x, na.rm = TRUE)/sqrt(sum(!is.na(.x))),
                        N = ~sum(!is.na(.x)),
                        Median = ~median(.x, na.rm = TRUE), 
                        LowerQuartile = ~quantile(.x, 0.25, na.rm = TRUE),
                        UpperQuartile = ~quantile(.x, 0.75, na.rm = TRUE),
                        Med_LQ_UQ = ~paste0(
                          round(median(.x, na.rm = TRUE)), " (", 
                          round(quantile(.x, 0.25, na.rm = TRUE)), "-",
                          round(quantile(.x, 0.75, na.rm = TRUE)), ")")
                   ),
                   .names = "{col}_{fn}")
  ) %>%
  ungroup()

line_plots = list()

for (var in mce_measures_iqr) {
  current_median <- paste0(var, "_Median")
  current_LQ <- paste0(var, "_LowerQuartile")
  current_UQ <- paste0(var, "_UpperQuartile")
  var_label <- rename_mapping[sub("_iqr", "", var)]
  
  line_plots[[var]] <- summaries %>%
    ggplot(aes(x = time, y = .data[[current_median]],
               color = condition, shape = condition)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.3) +
    geom_errorbar(aes(
      ymin = .data[[current_LQ]],
      ymax = .data[[current_UQ]]),
      width = 0.8,
      size = 0.75,
      alpha = 0.8
    ) +
    scale_x_continuous(breaks = c(0, 6, 12),
                       labels = c("Baseline", 6, 12)) +
    expand_limits(y = 0) +
    # scale_color_discrete(direction = -1) +
    #scale_color_brewer(palette = "Dark2", direction = 1) +
    labs(y = str_wrap(var_label, width = 24), x = "Week") +
    theme_pub2()
  }

ggarrange(
  line_plots[[1]] + xlab(NULL) + ylim(0, 60) + theme(text = element_text(size = 12)), 
  line_plots[[2]] + xlab(NULL) + ylim(0, 60) + theme(text = element_text(size = 12)), 
  line_plots[[3]] + xlab(NULL) + ylim(0, 850) + theme(text = element_text(size = 12)), 
  line_plots[[4]] + xlab(NULL) + ylim(0, 850) + theme(text = element_text(size = 12)), 
  line_plots[[5]] + xlab(NULL) + ylim(0, 700) + theme(text = element_text(size = 12)), 
  line_plots[[6]] + xlab(NULL) + ylim(0, 700) + theme(text = element_text(size = 12)),
  line_plots[[7]] + ylim(0, 850) + theme(text = element_text(size = 12)), 
  line_plots[[8]] + ylim(0, 850) + theme(text = element_text(size = 12)),
  common.legend = T, nrow = 4, ncol = 2
)
ggsave(here("Output", "Excitability", "Lineplots_MCEMeasures.pdf"))
ggsave(here("Output", "Excitability", "Lineplots_MCEMeasures.png"))

# Tables ---------
## Table 1 ---------
# original data
table1 <- data.frame(matrix(nrow = 0, ncol = 4))
conditions_with_ns <- mapply(function(x, y) paste0(x, "---N = ", y), 
                             levels(xct0$condition_rev), 
                             xct0 %>%
                               count(condition_rev) %>%
                               pull(n))
colnames(table1) <- c("Characteristic",
                      conditions_with_ns, 
                      "P")
sample_ns <- xct0 %>% count(condition_rev) %>% pull(n) %>% unique()
vars <- c(
  list(c("Demographics", "", "header row"),
       c("gender", "Gender (% Female)", "cat", "female"),
       c("Idade", "Age, y - mean ± SD", "cont"),
       c("etnia", "White - no. (%)", "cat", 1),
       c("anos_estudo", "Years of education - mean ± SD", "cont"),
       c("renda", "At least 3 minimum wages - no. (%)", "cat", 2, 3, 4),
       c("ocupacao", "Retired - no. (%)", "cat", 4),
       c("Clinical characteristics", "", "header row"),
       c("hipertensao", "Hypertension - no. (%)", "cat", 1),
       c("hipotireoidismo", "Hypothyroidism - no. (%)", "cat", 1),
       c("mini_26", "Onset age of MDD, years - mean ± SD", "cont"),
       c("mini_27", "No. previous depressive episodes - median (IQR)", "count"),
       c("mini_28", "Duration of current episode, months - mean ± SD", "cont"),
       c("Scales at baseline - mean ± SD", "", "header row"),
       c("hdrs_total", "Hamilton Depression Rating Scale, 17-item version (HDRS-17)", "cont"),
       c("madrs_total", "Montgomery-Asberg Depression Rating Scale (MADRS)", "cont"),
       c("gds_total", "Geriatric Depression Scale (GDS)", "cont"),
       c("ymrs_total", "Young Mania Rating Scale (YMRS)", "cont"),
       c("panas_pos_total", "Positive affect (PANAS)", "cont"),
       c("panas_neg_total", "Negative affect (PANAS)", "cont"),
       c("CGI1", "Clinical Global Impression (CGI-1)", "cont"),
       c("G_CIRS_0", "Cumulative Illness Rating Scale-Geriatric (G-CIRS)", "cont"),
       c("Excitability measures at baseline - median (IQR)", "", "header row")),
       Map(function(value, name) c(paste0(value, "_iqr"), name, "non-normal"), names(rename_mapping)[6:13], rename_mapping[6:13])
)

for (var in vars) {
  var_name <- var[1]
  var_desc <- var[2]
  var_type <- var[3]
  print(var_name)
  if (grepl(var_type, "header row")) {
    table1[nrow(table1)+1, ] <- c(var_name, "", "", NA_real_)
    
  } else if (grepl(var_type, "cont")) {
    
    ns <- xct0 %>%
      group_by(condition_rev) %>%
      summarize(desc = paste0(round(mean(.data[[var_name]], na.rm=TRUE), 1),
                              " ± ",
                              round(sd(.data[[var_name]], na.rm=TRUE), 1))
      ) %>%
      pull(desc)
    
    current_formula <- as.formula(paste0("`", var_name, "`", " ~ condition"))
    
    table1[nrow(table1)+1, ] <- c(
      var_desc,
      ns,
      tidy(aov(current_formula, xct0))$p.value[1]
    )
    
  } else if (grepl(var_type, "non-normal")) {
    
    ns <- xct0 %>%
      group_by(condition_rev) %>%
      summarize(desc = paste0(
        round(median(.data[[var_name]], na.rm=TRUE), 0), " (", 
        round(quantile(.data[[var_name]], 0.25, na.rm=TRUE), 0), "-",
        round(quantile(.data[[var_name]], 0.75, na.rm=TRUE), 0), ")"
        )
      ) %>%
      pull(desc)

    current_formula <- as.formula(paste0("`", var_name, "`", " ~ condition"))
    
    table1[nrow(table1)+1, ] <- c(
      var_desc,
      ns,
      wilcox.test(current_formula, xct0)$p.value
    )  

  } else if (grepl(var_type, "cat")) {
    
    relevant_cats <- var[4:length(var)] 
    el_temp <- xct0 %>%
      mutate(current_var = .data[[var_name]] %in% relevant_cats) %>%
      select(condition_rev, current_var)
    
    ns <- el_temp %>%
      group_by(condition_rev) %>%
      summarize(desc = paste0(sum(current_var, na.rm=TRUE),
                              " (",
                              round(sum(current_var, na.rm=TRUE)/sum(!is.na(current_var))*100),
                              ")")
      ) %>%
      pull(desc)
    
    table1[nrow(table1)+1, ] <- c(
      var_desc,
      ns,
      test_independence(el_temp, "condition_rev", "current_var")
    )
    
  } else if (grepl(var_type, "count")) {
    ns <- xct0 %>%
      group_by(condition_rev) %>%
      summarize(desc = paste0(median(.data[[var_name]], na.rm=TRUE),
                              " (",
                              round(quantile(.data[[var_name]], na.rm=TRUE, probs = 0.25)),
                              "-",
                              round(quantile(.data[[var_name]], na.rm=TRUE, probs = 0.75)),
                              ")"
      )
      ) %>%
      pull(desc)
    
    table1[nrow(table1)+1, ] <- c(
      var_desc,
      ns,
      test_count(xct0, var_name, "condition_rev")
    )
  }
}

ft_table1  <- table1 %>%
  mutate(P = ifelse(P == "-", 100, round(as.numeric(P), 3)),
         P = case_when(
           P == 100 ~ "-",
           P == 1 ~ ">0.99",
           P == 0 ~ "<0.001",
           # P > 0.02 ~ as.character(round(P, 2)),
           TRUE ~ as.character(P))
  ) %>%
  flextable() %>%
  split_header(sep = "---") %>%
  italic(i = is.na(table1$P), j = 1) %>%
  merge_h(i = is.na(table1$P)) %>%
  autofit()

ft_table1

save_as_docx(ft_table1 %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Table1_Idosos_Excitability.docx"))

## Table MCE measures by group and time ------
tb_mce <- summaries %>%
  select(time, condition, ends_with("Med_LQ_UQ")) %>%
  pivot_longer(cols = SP_Limiar_esquerda_iqr_Med_LQ_UQ:PP_Intervalo_direita_10.15.20mseg_iqr_Med_LQ_UQ,
               names_to = "MCE measure",
               values_to = "Value"
  ) %>%
  pivot_wider(id_cols = `MCE measure`,
               names_from = c(time, condition),
               values_from = Value,
               names_sep = "_"         
  ) 

ft_mce <- tb_mce %>%
  mutate(`MCE measure` = rename_mapping[gsub("_iqr_Med_LQ_UQ", "", `MCE measure`)]) %>%
  flextable() %>%
  span_header() %>%
  autofit()

ft_mce

save_as_docx(ft_mce %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Table_MCEmeasures.docx"))

## Flowchart ---------

xct %>%
  group_by(condition, time) %>%
  summarize(
    N_HDRS = sum(!is.na(hdrs_total)),
    N_HDRS.MCE = sum(!is.na(hdrs_total) & (!is.na(PP_Intervalo_direita_1mseg_iqr) | !is.na(PP_Intervalo_direita_1mseg_iqr))),
    N_MCE1 = sum(!is.na(PP_Intervalo_direita_1mseg_iqr)),
    N_MCE2 = sum(!is.na(PP_Potência_direita_iqr)),
    N_MCE_alt = sum(!is.na(PP_Potência_direita_iqr) | !is.na(PP_Intervalo_direita_1mseg_iqr)),
  ) %>%
  rowwise() %>%
  mutate(
    N_MCE = max(c_across(N_MCE1:N_MCE2))
  ) %>%
  select(-N_MCE1, -N_MCE2)
  
xct %>% filter(time == 0|time==6|time==12) %>% select(subject, time, SP_Limiar_direita_iqr:PP_Intervalo_esquerda_10.15mseg_iqr) %>% View()

xct %>% filter(time == 0 & condition == "Active TBS") %>% select(subject, time, SP_Limiar_direita_iqr:PP_Intervalo_esquerda_10.15mseg_iqr) %>% View()



# Treatment effect on TMS indices -----------

## ART ANOVA ------
### Interaction time*group -----------
# Put in a table
tb_EffectsOnTMS <-  tibble(
  `MCE measure` = character(),
  `F` = numeric(),
  DF = character(),
  P = numeric()
)

for (week in c(6, 12)) {
  tb_EffectsOnTMS <- tb_EffectsOnTMS %>%
    add_row(`MCE measure` = paste0("Week_", week))
  
  for (tms in mce_measures_iqr) {
    current_formula <- as.formula(paste0("`", tms, "`", " ~ condition * time_fctr + (1|subject)"))
    current_model <-  tidy(anova(art(current_formula, data = xct_tms %>% 
                            filter(!is.na(.data[[tms]])) %>% 
                            filter(time <= week))))

    tb_EffectsOnTMS <- tb_EffectsOnTMS %>%
      add_row(`MCE measure` = tms,
              `F` = round(current_model$statistic[3], 2),
              DF = paste("1,", round(current_model$Df.res[3])),
              P = round(current_model$p.value[3], 3)
              )
    }
}

ft_EffectsOnTMS <- tb_EffectsOnTMS %>%
  mutate(
    `MCE measure` = rename_mapping[sub("_iqr", "", `MCE measure`)],
    "P(FDR)" = round(p.adjust(P, method = "fdr"), 3)
    ) %>%
  # rename(`MCE measure` = MCE_measure) %>%
  flextable() %>%
  italic(i = is.na(tb_EffectsOnTMS$P), j = 1) %>%
  align(j = 3, align = "right", part = "all") %>%
  autofit()

ft_EffectsOnTMS
save_as_docx(ft_EffectsOnTMS, 
             path = here("Output", "Excitability",
                         "Table_ARTAnova_EffectsOnMCE.docx"))

# alternative: Week on same row, FDR only for each endpoint
ft_EffectsOnTMS_alt <- cbind(
  # Predictor column
  tb_EffectsOnTMS[2:9, 1] %>%
    mutate(`MCE measure` = rename_mapping[gsub("_iqr", "", `MCE measure`)]),
  # Week 6
  tb_EffectsOnTMS[2:9, ] %>%
    mutate("Week6::P(FDR)" = round(p.adjust(P, method = "fdr"), 3)) %>%
    select("Week6::F"="F", "Week6::DF"=DF, "Week6::P"=P, "Week6::P(FDR)"),
  # Week 12
  tb_EffectsOnTMS[11:18, ] %>%
    mutate("Week12::P(FDR)" = round(p.adjust(P, method = "fdr"), 3)) %>%
    select("Week12::F"="F", "Week12::DF"=DF, "Week12::P"=P, "Week12::P(FDR)")
  )  %>%
  flextable() %>%
  # bold(c(1, 2), 4) %>%
  span_header(sep = "::") %>%
  align(part = "header", i = 1, align = "center") %>%
  autofit()

ft_EffectsOnTMS_alt

save_as_docx(ft_EffectsOnTMS_alt, 
             path = here("Output", "Excitability",
                         "Table_ARTAnova_EffectsOnMCE_Alt.docx"))

# Main Effect of time
# Put in a table
tb_ART_TimeEffect <- tibble(
  `MCE measure` = character(),
  `F` = numeric(),
  DF = character(),
  P = numeric()
)

for (week in c(6, 12)) {
  tb_ART_TimeEffect <- tb_ART_TimeEffect %>%
    add_row(`MCE measure` = paste0("Week_", week))
  
  for (tms in mce_measures_iqr) {
    current_formula <- as.formula(paste0("`", tms, "`", " ~ condition * time_fctr + (1|subject)"))
    current_model <-  tidy(anova(art(current_formula, data = xct_tms %>% 
                                       filter(!is.na(.data[[tms]])) %>% 
                                       filter(time <= week))))
    
    tb_ART_TimeEffect <- tb_ART_TimeEffect %>%
      add_row(`MCE measure` = tms,
              `F` = round(current_model$statistic[2], 2),
              DF = paste("1,", round(current_model$Df.res[2])),
              P = round(current_model$p.value[2], 3)
      )
  }
}

ft_ART_TimeEffect <- tb_ART_TimeEffect %>%
  mutate(
    `MCE measure` = rename_mapping[sub("_iqr", "", `MCE measure`)],
    "P(FDR)" = round(p.adjust(P, method = "fdr"), 3)
  ) %>%
  # rename("MCE measure" = MCE_measure) %>%
  flextable() %>%
  italic(i = is.na(tb_ART_TimeEffect$P), j = 1) %>%
  align(j = 3, align = "right", part = "all") %>%
  autofit()

ft_ART_TimeEffect
save_as_docx(ft_ART_TimeEffect, 
             path = here("Output", "Excitability",
                         "Table_ARTAnova_TimeEffectOnMCE.docx"))

# alternative: Week on same row, FDR only for each endpoint
ft_ART_TimeEffect_alt <- cbind(
  # Predictor column
  tb_ART_TimeEffect[2:9, 1] %>%
    mutate("MCE measure" = rename_mapping[gsub("_iqr", "", `MCE measure`)]),
  # Week 6
  tb_ART_TimeEffect[2:9, ] %>%
    mutate("Week6::P(FDR)" = round(p.adjust(P, method = "fdr"), 3)) %>%
    select("Week6::F"="F", "Week6::DF"=DF, "Week6::P"=P, "Week6::P(FDR)"),
  # Week 12
  tb_ART_TimeEffect[11:18, ] %>%
    mutate("Week12::P(FDR)" = round(p.adjust(P, method = "fdr"), 3)) %>%
    select("Week12::F"="F", "Week12::DF"=DF, "Week12::P"=P, "Week12::P(FDR)")
)  %>%
  flextable() %>%
  # bold(c(1, 2), 4) %>%
  span_header(sep = "::") %>%
  align(part = "header", i = 1, align = "center") %>%
  # hline(part = "header", i = 1) %>%
  autofit()

ft_ART_TimeEffect_alt

save_as_docx(ft_ART_TimeEffect_alt, 
             path = here("Output", "Excitability",
                         "Table_ARTAnova_TimeEffectOnMCE_Alt.docx"))

## Linear models -----
# Mainly done to compare results with logarithmic transformation to the Gordon 
# schizophrenia paper, but also to see how residuals behave.
# Also see what the addition of age gender covariates do. 

for (mce in mce_measures_ext) {
  print(mce)
  current_formula <- as.formula(paste0("`", mce, "`", " ~ condition * time + age + gender + (1|subject)"))
  current_model <-  lmer(current_formula, data = xct_tms %>% filter(!is.na(.data[[mce]])) %>% filter(time <8))
  if (tidy(current_model)[4,8] < 0.1){
    print(shapiro.test(residuals(current_model)))
    print(tidy(current_model)[4,3:8])
  }
}

for (mce in mce_measures_ext) {
  print(mce)
  current_formula <- as.formula(paste0("`", mce, "`", " ~ condition * time + age + gender + (1|subject)"))
  current_model <-  lmer(current_formula, data = xct_tms %>% filter(!is.na(.data[[mce]])))
  if (tidy(current_model)[4,8] < 0.1){
    print(shapiro.test(residuals(current_model)))
    print(tidy(current_model)[4,3:8])
  }
}

## Difference data & t-tests ----------
# Difference data are interesting because they can be used in t-tests, which 
# do not require strict normality if n is at least 30. As such, I did not try
# transformations. However, t-tests can still be disproportionally affected by
# outliers, so I used IQR to exclude outliers. 


iqr_diffs <- xct %>%
  select(contains("_iqr_diff")) %>%
  names()

# Create table of t-tests
for (var in iqr_diffs) {
  var_desc <- sub("_iqr_diff[0-9]{1,2}", "", var)
  current_formula <- as.formula(paste0("`", var, "`", " ~ condition"))
  current_t <- t_test(xct0, current_formula) %>%
    mutate(Measure = var_desc)
  
  if (var == iqr_diffs[1]) {
    table_diffs_iqr <- current_t %>%
      add_row(Measure = "Week_6", .before = 1) 
  } else if (var == iqr_diffs[9]) {
    table_diffs_iqr <- table_diffs_iqr %>%
      add_row(Measure = "Week_12") %>%
      add_row(current_t)
  } else {
    table_diffs_iqr <- table_diffs_iqr %>%
      add_row(current_t)
  }
}

# Change to flextable
ft_diffs_iqr <- table_diffs_iqr %>%
  mutate(
    Measure = rename_mapping[Measure], 
    t = round(statistic, 2), DF = round(t_df), P = round(p_value, 3), 
    .keep = "none") %>%
  flextable() %>%
  italic(i = grepl("Week", table_diffs_iqr$Measure), j = 1) %>%
  bold(i = table_diffs_iqr$p_value < 0.049, j = 4) %>%
  autofit()

ft_diffs_iqr

save_as_docx(ft_diffs_iqr %>% fit_pagewidth(), 
             path = "Output\\Excitability\\Excitability_TreatmentEffects_IQR.docx")


# Antidepressant effect in current sample -------

# No covariates: 
model_hdrs_6 <- lmer(hdrs_total ~ condition*time_log + (1|subject), 
                     data = filter(xct, time <= 6))
summary(model_hdrs_6)

model_hdrs_12 <- lmer(hdrs_total ~ condition*time_log + (1|subject), 
                      data = xct)
summary(model_hdrs_12)

# With covariates
model_hdrs_6 <- lmer(hdrs_total ~ condition*time_log + age + gender + handedness + years_study + (1|subject), 
                     data = filter(xct, time <= 6))
summary(model_hdrs_6)

model_hdrs_12 <- lmer(hdrs_total ~ condition*time_log + age + gender + handedness + years_study + (1|subject), 
                      data = xct)
summary(model_hdrs_12)



# Prediction of antidepressant effect --------

preds <- sapply(mce_measures_iqr, function(x) paste0(x, "_bl"))
# preds <- sapply(mce_measures_ext, function(x) paste0(x, "_bl"))
weeks <- c(6, 12)
pred_lmms <- list(Week_6 = list(), Week_12 = list())

# Create empty df that will contain model stats
tb_pred_lmm <- tibble(
  Predictor = character(),
  `Main effect::Estimate` = numeric(),
  `Main effect::SE` = numeric(),
  `Main effect::t` = numeric(),
  `Main effect::DF` = numeric(),
  `Main effect::P` = numeric(),
  `Interaction::Estimate` = numeric(),
  `Interaction::SE` = numeric(),
  `Interaction::t` = numeric(),
  `Interaction::DF` = numeric(),
  `Interaction::P` = numeric()
)

# Loop over weeks and predictors, fit models and add info to list and df 
for (week in weeks) {
  
  for (i in seq_along(preds)) {
    
    # Create formula
    current_formula <- as.formula(
      paste0("hdrs_total ~ condition * time_log * `", preds[i], "` + (1|subject)" ) #  + age + gender + years_study + handedness 
    )

    # Fit model
    current_model <- lmer(current_formula, filter(xct, time <= {{week}}))
    # Add model to list of models
    pred_lmms[[paste0("Week_", week)]][[preds[i]]] <- current_model
    # print(summary(current_model))
    # Create tidy version of model
    current_tidy <- tidy(current_model) %>% select(estimate:p.value)
    rownum_interaction <- grep("^conditionActive TBS:time_log:", 
                               tidy(current_model)$term) 

    if (current_tidy[rownum_interaction,]$p.value < 0.05) {
      print(paste("Significant three-way interaction for", preds[i], "until week", week))
      
      emm_formula <- as.formula(paste("~", preds[i], "* condition * time_log"))
      emm_list <- list(time_log = c(0, 1))
      emm_list[[preds[i]]] <- c(0, 1)
      current_contrasts <- emmeans(current_model, emm_formula,
                             at = emm_list) %>%
        contrast(interaction = "poly", by = "condition")
      print(current_contrasts)
      t_to_d(
        t = current_contrasts %>% tidy() %>% filter(condition == "Active TBS") %>% select(statistic) %>% pull(),
        df = current_contrasts %>% tidy() %>% filter(condition == "Active TBS") %>% select(df) %>% pull()
        ) %>%
        transmute(`Cohen's d` = paste0(
          round(d*log(week+1), 2), " [",
          round(CI_low*log(week+1), 2), ", ",
          round(CI_high*log(week+1), 2), "]"
        )) %>%
        print()
    }
   
    # Add week indicator if it changes
    if (i == 1) {
      tb_pred_lmm <- tb_pred_lmm %>%
        add_row(Predictor = paste0("Week_", week)) 
    }

    # Add current model stats to output dataframe
    tb_pred_lmm <- rbind(
      tb_pred_lmm, 
      setNames(
        c("Predictor" = preds[i], current_tidy[4,], current_tidy[rownum_interaction,]),
        names(tb_pred_lmm)
      )
    )
  }
}

# Transform model output df to flextable
ft_pred_lmm <- tb_pred_lmm %>%
  # filter(`Interaction::P` < 0.05 | is.na(`Interaction::P`)) %>%
  mutate(
    across(
      c(`Main effect::Estimate`:`Main effect::t`,
        `Interaction::Estimate`:`Interaction::t`),
      ~ round(., 3)),
    across(
      c(`Main effect::DF`, `Interaction::DF`),
      round),
    across(
      c(`Main effect::P`, `Interaction::P`),
      ~ round(., 3)),
    `Interaction::P` = ifelse(`Interaction::P` == 0.000, "<0.001", `Interaction::P`),
    Predictor = rename_mapping[gsub("_iqr_bl", "", Predictor)] #gsub("_", " ", Predictor)
  ) %>%
  flextable() %>%
  bold(~ `Main effect::P` < 0.05, 6) %>%
  bold(~ `Interaction::P` < 0.05, 11) %>%
  italic(~startsWith(Predictor, "Week"), 1) %>%
  span_header(sep = "::") %>%
  align(j = 11, align = "right", part = "all") %>%
  align(i = 1, align = 'center' ,part = "header")  %>%
  autofit() 

ft_pred_lmm

save_as_docx(ft_pred_lmm %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Prediction_HDRS_LMM_v2.docx"))

# Only interactions
ft_pred_lmm_onlyInteraction <- tb_pred_lmm %>%
  # filter(`Interaction::P` < 0.05 | is.na(`Interaction::P`)) %>%
  mutate(
    across(
      c(`Main effect::Estimate`:`Main effect::t`,
        `Interaction::Estimate`:`Interaction::t`),
      ~ round(., 3)),
    across(
      c(`Main effect::DF`, `Interaction::DF`),
      round),
    across(
      c(`Main effect::P`, `Interaction::P`),
      ~ round(., 3)),
    "P(FDR)" = round(p.adjust(`Interaction::P`, method = "fdr"), 3),
    `Interaction::P` = ifelse(`Interaction::P` == 0.000, "<0.001", `Interaction::P`),
    Predictor = rename_mapping[gsub("_iqr_bl", "", Predictor)] #gsub("_", " ", Predictor)
  ) %>%
  select(Predictor, "Interaction::Estimate":"P(FDR)") %>%
  rename(
    Estimate = `Interaction::Estimate`,
    SE = `Interaction::SE`,
    t = `Interaction::t`,
    DF = `Interaction::DF`,
    P = `Interaction::P`
  ) %>%
  flextable() %>%
  bold(~ P < 0.05, 6) %>%
  italic(~startsWith(Predictor, "Week"), 1) %>%
  align(j = 6, align = "right", part = "all") %>%
  align(i = 1, align = 'center' ,part = "header")  %>%
  autofit() 

ft_pred_lmm_onlyInteraction

save_as_docx(ft_pred_lmm_onlyInteraction %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Prediction_HDRS_LMM_OnlySignificantInteraction.docx"))

# alternative: Week on same row, FDR only for each endpoint
ft_pred_lmm_alt <- cbind(
  # Predictor column
  tb_pred_lmm[2:9, 1] %>%
    mutate( Predictor = rename_mapping[gsub("_iqr_bl", "", Predictor)]),
  # Week 6
  tb_pred_lmm[2:9, ] %>%
    mutate(
      "Week6::t" = round(`Interaction::t`, 2),
      "Week6::DF" = round(`Interaction::DF`),
      "Week6::P(FDR)" = round(p.adjust(`Interaction::P`, method = "fdr"), 3)) %>%
    select("Week6::t", "Week6::DF", "Week6::P(FDR)"),
  # Week 12
  tb_pred_lmm[11:18, ] %>%
    mutate(
      "Week12::t" = round(`Interaction::t`, 2),
      "Week12::DF" = round(`Interaction::DF`),
      "Week12::P(FDR)" = round(p.adjust(`Interaction::P`, method = "fdr"), 3)) %>%
    select("Week12::t", "Week12::DF", "Week12::P(FDR)")
  ) %>%
  flextable() %>%
  bold(c(1, 2), 4) %>%
  bold(1, 7) %>%
  span_header(sep = "::") %>%
  autofit()

ft_pred_lmm_alt
  
save_as_docx(ft_pred_lmm_alt, 
             path = here("Output", "Excitability",
                         "Table_Prediction_HDRS_LMM_Alt.docx"))


## Make graphs for significant predictions ----------
# Version Red-blue  
plot1 <- ggplot(xct0, aes(SP_Limiar_esquerda_iqr, -hdrs_total_diff6, 
                 color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  scale_color_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  scale_fill_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  scale_x_continuous(labels = ~paste0(., "%")) +
  theme_pub2() +
  labs(x = "Left RMT at baseline", y = "HDRS-17 improvement until week 6")

plot2 <- ggplot(xct0, aes(SP_Limiar_direita_iqr, -hdrs_total_diff6, 
                 color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  scale_x_continuous(labels = ~paste0(., "%")) +
  scale_color_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  scale_fill_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  theme_pub2() +
  labs(x = "Right RMT at baseline", y = "HDRS-17 improvement until week 6")

plot3 <- ggplot(xct0, aes(SP_Limiar_esquerda_iqr, -hdrs_total_diff12, 
                 color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  scale_x_continuous(labels = ~paste0(., "%")) +
  scale_color_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  scale_fill_manual(values = c("Sham TBS" = "red", "Active TBS" = "blue")) +
  theme_pub2() +
  labs(x = "Left RMT at baseline", y = "HDRS-17 improvement until week 12")

# Arrange 3 plots and a common legend bottom right with package gridExtra
legend <- g_legend(plot1 + theme(legend.position='right'))
g <- grid.arrange(plot1+theme(legend.position='hidden'), plot2+theme(legend.position='hidden'),
             plot3+theme(legend.position='hidden'), legend)
g


ggsave(here("Output", "Excitability", "Plots_SignificantPredictors_RedBlue.png"), g)
ggsave(here("Output", "Excitability", "Plots_SignificantPredictors_RedBlue.pdf"), g)



# alternative Version Ggplot default colors
plot1 <- ggplot(xct0, aes(SP_Limiar_esquerda_iqr, -hdrs_total_diff6, 
                            color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  # scale_color_discrete(direction = -1) +
  scale_x_continuous(labels = ~paste0(., "%")) +
  theme_pub2() +
  labs(x = "Left RMT at baseline", y = "HDRS-17 improvement until week 6")

plot2 <- ggplot(xct0, aes(SP_Limiar_direita_iqr, -hdrs_total_diff6, 
                          color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  scale_x_continuous(labels = ~paste0(., "%")) +
  # scale_color_discrete(direction = -1) +
  theme_pub2() +
  labs(x = "Right RMT at baseline", y = "HDRS-17 improvement until week 6")

plot3 <- ggplot(xct0, aes(SP_Limiar_esquerda_iqr, -hdrs_total_diff12, 
                          color = condition, fill = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  scale_x_continuous(labels = ~paste0(., "%")) +
  # scale_color_discrete(direction = -1) +
  theme_pub2() +
  labs(x = "Left RMT at baseline", y = "HDRS-17 improvement until week 12")

# Arrange 3 plots and a common legend bottom right with package gridExtra
legend <- g_legend(plot1 + theme(legend.position='right'))
g <- grid.arrange(plot1+theme(legend.position='hidden'), plot2+theme(legend.position='hidden'),
                  plot3+theme(legend.position='hidden'), legend)
g
  
ggsave(here("Output", "Excitability", "Plots_SignificantPredictors_ggplotColors.png"), g)
ggsave(here("Output", "Excitability", "Plots_SignificantPredictors_ggplotColors.pdf"), g)


# Old way (not possible to place legend in empty space)
ggarrange(plot1, plot2, plot3, common.legend = T, legend = "right")


# Correlation depression change - cortical excitability change -----
## Whole sample -------
weeks <- c(6, 12)

# Create empty df that will contain model stats
corr_table <- tibble(
  `Excitability measure` = character(),
  Coefficient = numeric(),
  P = numeric()
)

# Loop over weeks and excitability measures, and add correlation and df 
for (week in weeks) {
  
  current_hdrs <- paste0("hdrs_total_diff", week)
  
  for (i in seq_along(mce_measures_iqr)) {
    # browser()
    current_measure <- paste0(mce_measures_iqr[i], "_diff", week)
    current_cor <- cor.test(xct0[[current_hdrs]], xct0[[current_measure]])
    current_r <- round(current_cor$estimate, 2)
    current_P <- round(current_cor$p.value, 3)
    
    # Add week indicator if it changes
    if (i == 1) {
      corr_table <- corr_table %>%
        add_row(`Excitability measure` = paste0("Difference baseline - week ", week))
    }
    
    # Add current model stats to output dataframe
    corr_table <- rbind(
      corr_table,
      setNames(
        c(mce_measures_iqr[i], current_r, current_P),
        names(corr_table)
      )
    )
  }
}

corr_ft <- corr_table %>%
  mutate(
    `Excitability measure` = rename_mapping[gsub("_iqr", "", `Excitability measure`)] #gsub("_", " ", Predictor)
  ) %>%
  flextable() %>%
  bold(~ P < 0.05, 3) %>%
  italic(~startsWith(`Excitability measure`, "Week"), 1) %>%
  align(j = c(2, 3), align = "right", part = "all") %>%
  # align(i = 1, align = 'center' ,part = "header")  %>%
  italic(i = c(1, 10), j = 1) %>% 
  autofit() 

corr_ft

save_as_docx(corr_ft %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Table_CorrelationBetweenChanges_WholeSample.docx"))

## Condition Split -------
# Create empty df that will contain model stats
corr_split_table <- tibble(
  `Excitability measure` = character(),
  `Active TBS::Coefficient` = numeric(),
  `Active TBS::P` = numeric(),
  `Sham TBS::Coefficient` = numeric(),
  `Sham TBS::P` = numeric()
)

weeks <- c(6, 12)
xct0_active <- filter(xct0, condition == "Active TBS")
xct0_sham <- filter(xct0, condition == "Sham TBS")

# Loop over weeks and excitability measures, and add correlation and df 
for (week in weeks) {
  
  current_hdrs <- paste0("hdrs_total_diff", week)
  
  for (i in seq_along(mce_measures_iqr)) {
    # browser()
    current_measure <- paste0(mce_measures_iqr[i], "_diff", week)
    
    corr_active <- cor.test(xct0_active[[current_hdrs]], xct0_active[[current_measure]])
    r_active <- round(corr_active$estimate, 2)
    P_active <- round(corr_active$p.value, 3)
    
    corr_sham <- cor.test(xct0_sham[[current_hdrs]], xct0_sham[[current_measure]])
    r_sham <- round(corr_sham$estimate, 2)
    P_sham <- round(corr_sham$p.value, 3)
    
    # Add week indicator if it changes
    if (i == 1) {
      corr_split_table <- corr_split_table %>%
        add_row(`Excitability measure` = paste0("Difference baseline - week ", week))
    }
    
    # Add current model stats to output dataframe
    corr_split_table <- rbind(
      corr_split_table,
      setNames(
        c(mce_measures_iqr[i], r_active, P_active, r_sham, P_sham),
        names(corr_split_table)
      )
    )
  }
}

corr_split_ft <- corr_split_table %>%
  mutate(
    `Excitability measure` = rename_mapping[gsub("_iqr", "", `Excitability measure`)] #gsub("_", " ", Predictor)
  ) %>%
  flextable() %>%
  bold(~ `Sham TBS::P` < 0.05, 5) %>%
  bold(~ `Active TBS::P` < 0.05, 3) %>%
  italic(~startsWith(`Excitability measure`, "Week"), 1) %>%
  align(j = c(2, 3), align = "right", part = "all") %>%
  # align(i = 1, align = 'center' ,part = "header")  %>%
  italic(i = c(1, 10), j = 1) %>% 
  span_header(sep = "::") %>%
  autofit() 

corr_split_ft

save_as_docx(corr_split_ft %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Table_CorrelationBetweenChanges_ConditionSplit.docx"))

# TAble for active TBS only
corr_split_ActiveTBSonly_ft <- cbind(
    corr_split_table %>%
      select(`Excitability measure`, 
             `Change to week 6::Coefficient` = `Active TBS::Coefficient`,
             `Change to week 6::P` = `Active TBS::P`) %>%
      slice(2:9),
    corr_split_table %>%
      select(`Change to week 12::Coefficient` = `Active TBS::Coefficient`,
             `Change to week 12::P` = `Active TBS::P`) %>%
      slice(11:18) 
  ) %>%
  mutate(
    `Excitability measure` = rename_mapping[gsub("_iqr", "", `Excitability measure`)] #gsub("_", " ", Predictor)
  ) %>%
  flextable() %>%
  align(j = 2:5, align = "right", part = "all") %>%
  # align(i = 1, align = 'center' ,part = "header")  %>%
  span_header(sep = "::") %>%
  autofit() 

corr_split_ActiveTBSonly_ft

save_as_docx(corr_split_ActiveTBSonly_ft %>% fit_pagewidth(), 
             path = here("Output", "Excitability",
                         "Table_CorrelationBetweenChanges_ActiveTBSonly.docx"))
