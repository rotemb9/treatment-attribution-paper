# Analysis of the CVS attribution vignettes (processed data- after combination across weeks and participant exclusions)

## Required R package
library(tidyr)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(ggpubr)
library(ggeffects)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
set.seed(9)
data_path = "data/data_for_analysis"
fig_output_path = "figures"

### read data
data_filename = "vignettes_data_combined.csv"
data_wide = read.csv(file.path(data_path, data_filename))

### manipulation tests
data_wide$Q1_zscored = scale(data_wide$Q1, center = TRUE, scale = TRUE)
data_wide$Q3_zscored = scale(data_wide$Q3, center = TRUE, scale = TRUE)
data_wide$Q4_zscored = scale(data_wide$Q4, center = TRUE, scale = TRUE)
data_wide$Q6_zscored = scale(data_wide$Q6, center = TRUE, scale = TRUE)
data_wide$age_zscored = scale(data_wide$age, center = TRUE, scale = TRUE)
## the effect of symptom level on Q3
# model
data_wide$vignette_symptoms = factor(data_wide$vignette_symptoms, levels = c("Med", "Low", "High"))
manipulation_check_symptoms = lmer(Q3_zscored ~ vignette_symptoms + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = data_wide)
summary(manipulation_check_symptoms)

# summary statistics (mean and standard error)
data_wide$vignette_symptoms = factor(data_wide$vignette_symptoms, levels = c("Low", "Med", "High"))
summary_stats_symptoms = data_wide %>%
  group_by(vignette_symptoms) %>%
  summarise(
    mean_Q3 = mean(Q3, na.rm = TRUE),
    se_Q3 = sd(Q3, na.rm = TRUE) / sqrt(n())
  )

## The effect of alternative treatment level on Q6
# a = OTC drug; b = rest; c = vitamins
# model
data_wide$vignette_alternative_treatment = factor(data_wide$vignette_alternative_treatment, levels = c("c", "a", "b"))
manipulation_check_alt_treat = lmer(Q6_zscored ~ vignette_alternative_treatment + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = data_wide)
summary(manipulation_check_alt_treat)

# summary statistics (mean and standard error)
data_wide$vignette_alternative_treatment = factor(data_wide$vignette_alternative_treatment, levels = c("a", "b", "c"))
summary_stats_alt_treat <- data_wide %>%
  group_by(vignette_alternative_treatment) %>%
  summarise(
    mean_Q6 = mean(Q6, na.rm = TRUE),
    se_Q6 = sd(Q6, na.rm = TRUE) / sqrt(n())
  )

## the effect of outcome on Q1
# model
data_wide$vignette_outcome = factor(data_wide$vignette_outcome, levels = c("R", "NR"))
manipulation_check_outcome = lmer(Q1_zscored ~ vignette_outcome + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = data_wide)
summary(manipulation_check_outcome)

# summary statistics (mean and standard error)
summary_stats_outcome = data_wide %>%
  group_by(vignette_outcome) %>%
  summarise(
    mean_Q1 = mean(Q1, na.rm = TRUE),
    se_Q1 = sd(Q1, na.rm = TRUE) / sqrt(n())
  )


### Qualitative tests
# discounting model
recovery_data = data_wide[data_wide$vignette_outcome == "R",]
recovery_data$Q1_zscored = scale(recovery_data$Q1, center = TRUE, scale = TRUE)
recovery_data$Q3_zscored = scale(recovery_data$Q3, center = TRUE, scale = TRUE)
recovery_data$Q4_zscored = scale(recovery_data$Q4, center = TRUE, scale = TRUE)
recovery_data$Q6_zscored = scale(recovery_data$Q6, center = TRUE, scale = TRUE)
recovery_data$age_zscored = scale(recovery_data$age, center = TRUE, scale = TRUE)

discounting_model = lmer(Q1_zscored ~ Q3_zscored * Q4_zscored * Q6_zscored + vignette_gender + gender + age_zscored + (1 | src_subject_id),data = recovery_data)
summary(discounting_model)
discounting_model_no_int = lmer(Q1_zscored ~ Q3_zscored + Q4_zscored + Q6_zscored + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = recovery_data)
summary(discounting_model_no_int)
discounting_model_no_cov = lmer(Q1_zscored ~ Q3_zscored * Q4_zscored * Q6_zscored + (1 | src_subject_id),data = recovery_data)
summary(discounting_model_no_cov)

# Robustness check: add experimentally manipulated alternative-treatment type
# to the preregistered recovery-only reported-belief model.
# The qualitative pattern is unchanged.
recovery_data$vignette_alternative_treatment = factor(recovery_data$vignette_alternative_treatment, levels = c("c", "a", "b"))
discounting_model_no_int_with_alt_treat_level = lmer(Q1_zscored ~ Q3_zscored + Q4_zscored + Q6_zscored + vignette_alternative_treatment + vignette_gender + gender + age_zscored  + (1 | src_subject_id), data = recovery_data)
summary(discounting_model_no_int_with_alt_treat_level)
discounting_model_with_alt_treat_level = lmer(Q1_zscored ~ Q3_zscored * Q4_zscored * Q6_zscored + vignette_alternative_treatment + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = recovery_data)
summary(discounting_model_with_alt_treat_level)

## recovery vs. no recovery
int_outcome_q3 = lmer(Q1_zscored ~ Q3_zscored * vignette_outcome + vignette_gender + gender + age_zscored  + (1 | src_subject_id), data = data_wide)
summary(int_outcome_q3)

# plot
Q3_by_Q1 = ggplot(data = data_wide, aes(x = Q3, y = Q1)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) + # Points with transparency and custom size
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") + # Adds regression line with confidence interval
  labs(title = "Q1 as a Function of Q3 by Vignette Outcome",
       x = "Q3",
       y = "Q1") + # Adds labels for the title and axes
  # Adjust the y-axis breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  # Adjust the x-axis breaks
  scale_x_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  theme_minimal(base_size = 14) + # Applies a minimal theme with larger text
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black")
    ) +
      facet_wrap(~ vignette_outcome, labeller = labeller(vignette_outcome = c("R" = "Recovery", "NR" = "No Recovery")))

Q4_by_Q1 = ggplot(data = data_wide, aes(x = Q4, y = Q1)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) + # Points with transparency and custom size
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") + # Adds regression line with confidence interval
  labs(title = "Q1 as a Function of Q4 by Vignette Outcome",
       x = "Q4",
       y = "Q1") + # Adds labels for the title and axes
  # Adjust the y-axis breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  # Adjust the x-axis breaks
  scale_x_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  theme_minimal(base_size = 14) + # Applies a minimal theme with larger text
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  facet_wrap(~ vignette_outcome, labeller = labeller(vignette_outcome = c("R" = "Recovery", "NR" = "No Recovery")))


Q6_by_Q1 = ggplot(data = data_wide, aes(x = Q6, y = Q1)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) + # Points with transparency and custom size
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") + # Adds regression line with confidence interval
  labs(title = "Q1 as a Function of Q6 by Vignette Outcome",
       x = "Q6",
       y = "Q1") + # Adds labels for the title and axes
  # Adjust the y-axis breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  # Adjust the x-axis breaks
  scale_x_continuous(
    breaks = seq(0, 100, by = 10)
  ) +
  theme_minimal(base_size = 14) + # Applies a minimal theme with larger text
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  facet_wrap(~ vignette_outcome, labeller = labeller(vignette_outcome = c("R" = "Recovery", "NR" = "No Recovery")))



## Robustness analyses controlling for Q8
recovery_data$Q8_zscored = scale(recovery_data$Q8, center = TRUE, scale = TRUE)
discounting_model_with_q8 = lmer(Q1_zscored ~ Q3_zscored * Q4_zscored * Q6_zscored + Q8_zscored + vignette_gender + gender + age_zscored + (1|src_subject_id), data = recovery_data)
summary(discounting_model_with_q8)
# Singular fit: participant random-intercept variance is estimated at 0.
# The corresponding lm yields identical fixed-effect estimates and tests.
discounting_model_no_int_with_q8 = lmer(Q1_zscored ~ Q3_zscored + Q4_zscored + Q6_zscored + Q8_zscored + vignette_gender + gender + age_zscored + (1|src_subject_id), data = recovery_data)
summary(discounting_model_no_int_with_q8)
data_wide$Q8_zscored = scale(data_wide$Q8, center = TRUE, scale = TRUE)
int_outcome_q3_with_q8 = lmer(Q1_zscored ~ Q8_zscored + Q3_zscored * vignette_outcome + vignette_gender + gender + age_zscored  + (1 | src_subject_id), data = data_wide)
summary(int_outcome_q3_with_q8)

## Exploratory outcome-shift analyses
outcome_shift_data = data_wide
outcome_shift_data$main_treat_diff = outcome_shift_data$Q1 - outcome_shift_data$Q8

# Descriptive statistics for outcome-related belief shifts
outcome_shift_data %>%
  group_by(vignette_outcome) %>%
  summarise(
    mean_shift = mean(main_treat_diff, na.rm = TRUE),
    sd_shift = sd(main_treat_diff, na.rm = TRUE)
  )

# Center continuous covariates across the full dataset
outcome_shift_data$Q8_c = outcome_shift_data$Q8 - mean(outcome_shift_data$Q8, na.rm = TRUE)
outcome_shift_data$age_c = outcome_shift_data$age - mean(outcome_shift_data$age, na.rm = TRUE)

# Effect-code outcome
outcome_shift_data$outcome_effect =
  ifelse(outcome_shift_data$vignette_outcome == "R", 0.5, -0.5)

# Recode participant gender
outcome_shift_data$gender[outcome_shift_data$gender %in% c("Prefer not to say", "Other")] = "Other/Not reported"

outcome_shift_data$gender = factor(
  outcome_shift_data$gender,
  levels = c("Female", "Male", "Other/Not reported")
)

# Create and center gender dummy variables
outcome_shift_data$gender_male =
  as.numeric(outcome_shift_data$gender == "Male")
outcome_shift_data$gender_other =
  as.numeric(outcome_shift_data$gender == "Other/Not reported")

outcome_shift_data$gender_male_c =
  outcome_shift_data$gender_male - mean(outcome_shift_data$gender_male, na.rm = TRUE)
outcome_shift_data$gender_other_c =
  outcome_shift_data$gender_other - mean(outcome_shift_data$gender_other, na.rm = TRUE)

outcome_shift_data$vignette_male =
  as.numeric(outcome_shift_data$vignette_gender == "Male")
outcome_shift_data$vignette_male_c =
  outcome_shift_data$vignette_male - mean(outcome_shift_data$vignette_male, na.rm = TRUE)

# Combined model: recovery vs. non-recovery
main_treat_diff_outcome_model = lmer(
  main_treat_diff ~ outcome_effect + Q8_c +
    vignette_male_c + gender_male_c + gender_other_c + age_c +
    (1 | src_subject_id),
  data = outcome_shift_data
)
summary(main_treat_diff_outcome_model)

# Recovery: adjusted shift from zero
main_treat_diff_model_recovery = lmer(
  main_treat_diff ~ Q8_c +
    vignette_male_c + gender_male_c + gender_other_c + age_c +
    (1 | src_subject_id),
  data = outcome_shift_data[outcome_shift_data$vignette_outcome == "R", ]
)
summary(main_treat_diff_model_recovery)

# Non-recovery: adjusted shift from zero
main_treat_diff_model_no_recovery = lmer(
  main_treat_diff ~ Q8_c +
    vignette_male_c + gender_male_c + gender_other_c + age_c +
    (1 | src_subject_id),
  data = outcome_shift_data[outcome_shift_data$vignette_outcome == "NR", ]
)
summary(main_treat_diff_model_no_recovery)

## compare based on conditions (manipulated) rather than reported values
## since reported values are measured in a non-ideal way
data_wide$vignette_symptoms = factor(data_wide$vignette_symptoms, levels = c("Med", "Low", "High"))
data_wide$vignette_alternative_treatment = factor(data_wide$vignette_alternative_treatment, levels = c("c", "a", "b"))
recovery_data$vignette_symptoms = factor(recovery_data$vignette_symptoms, levels = c("Med", "Low", "High"))
recovery_data$vignette_alternative_treatment = factor(recovery_data$vignette_alternative_treatment, levels = c("c", "a", "b"))
discounting_model_conditions = lmer(Q1_zscored ~ vignette_symptoms + vignette_alternative_treatment + vignette_gender + gender + age_zscored + (1|src_subject_id), data = recovery_data)
summary(discounting_model_conditions)
discounting_model_conditions_int = lmer(Q1_zscored ~ vignette_symptoms * vignette_alternative_treatment + gender + age_zscored + (1|src_subject_id), data = recovery_data)
summary(discounting_model_conditions_int)
drop1(discounting_model_conditions_int, test = "F")
# Sensitivity analysis including vignette-patient gender.
# Because vignette gender is partially confounded with the factorial vignette
# design within the recovery subset, one interaction coefficient is not
# independently estimable and is dropped from the model matrix.
# The omnibus interaction remains non-significant.
discounting_model_conditions_int_sensitivity = lmer(Q1_zscored ~ vignette_symptoms * vignette_alternative_treatment + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = recovery_data)
summary(discounting_model_conditions_int_sensitivity)
drop1(discounting_model_conditions_int_sensitivity, test = "F")
discounting_model_conditions_no_alt = lmer(Q1_zscored ~ vignette_symptoms + vignette_gender + gender + age_zscored + (1|src_subject_id), data = recovery_data)
anova(discounting_model_conditions_no_alt, discounting_model_conditions)
# the same model controlling for Q8
discounting_model_conditions2 = lmer(Q1_zscored ~ vignette_symptoms + vignette_alternative_treatment + Q8_zscored + vignette_gender + gender + age_zscored + (1|src_subject_id), data = recovery_data)
summary(discounting_model_conditions2)
# experimental analogue of the outcome-relevance model
outcome_relevance_conditions = lmer(Q1_zscored ~ vignette_symptoms + vignette_outcome + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = data_wide)
outcome_relevance_conditions_int = lmer(Q1_zscored ~ vignette_symptoms * vignette_outcome + vignette_gender + gender + age_zscored + (1 | src_subject_id), data = data_wide)
summary(outcome_relevance_conditions_int)
drop1(outcome_relevance_conditions_int, test = "F")
#anova(outcome_relevance_conditions, outcome_relevance_conditions_int)

### manuscript figures
## figure 1:
# Panel A: model simulation bar plot c a int
p_model_int_c_a = readRDS("figures/fig_model_interaction_c_a.rds")

# Panel B: empirical data bar plot int symptom level by alternative treatment (manipulated variables)
int_symptoms_alternative = ggplot(
  recovery_data,
  aes(
    x = factor(vignette_symptoms, levels = c("Low", "Med", "High")),
    y = Q1,
    fill = factor(vignette_alternative_treatment, levels = c("a","c","b"))
  )
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "black",
    width = 0.6,
    alpha = 0.4,
    position = position_dodge(width = 0.7)
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    color = "black",
    position = position_dodge(width = 0.7)
  ) +
  geom_jitter(
    aes(color = factor(vignette_alternative_treatment, levels = c("a","c","b"))),
    size = 1,
    alpha = 0.5,
    show.legend = FALSE,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.7
    )
  ) +
  scale_fill_discrete(labels = c("a" = "OTC drug", "b" = "Rest", "c" = "Vitamins")) +
  scale_color_discrete(labels = c("a" = "OTC drug", "b" = "Rest", "c" = "Vitamins")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  labs(
    x = "Symptom level",
    y = "Main treatment\nbelief (Q1)",
    fill = "Alternative treatment",
    title = "Empirical data"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        legend.position = "right",
        panel.grid.minor = element_blank()
  )

# Combine into one figure
fig1 = ggarrange(p_model_int_c_a, int_symptoms_alternative,
                 labels = c("A", "B"),
                 ncol = 1, nrow = 2)

# Save
size_factor = 1
ggsave('fig1.pdf',
       plot = fig1,
       path = fig_output_path,
       width = 140*size_factor,
       height = 140*size_factor,
       units = "mm",
       dpi = 300)

## figure 2 - manipulation check results
## Panel A – Q3 by symptom severity
fig2_panelA = ggplot(data_wide, aes(x = factor(vignette_symptoms, levels = c("Low", "Med", "High")), y = Q3, fill = vignette_symptoms)) +
  # Add bar plot
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
  # Add individual data points
  geom_jitter(width = 0.2, size = 1, alpha = 0.2, color = "black") +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  # Customize axis labels and theme
  labs(
    x = "Vignette symptoms level ",
    y = "Belief patient has\nCOVID-19 (Q3)"
  ) +
  # Adjust the y-axis limits and breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)  # Specify tick marks every 10 units
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## Panel B – Q1 by outcome
fig2_panelB = ggplot(data_wide, aes(x = vignette_outcome, y = Q1, fill = vignette_outcome)) +
  # Add bar plot
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
  # Add individual data points
  geom_jitter(width = 0.2, size = 1, alpha = 0.2, color = "black") +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  # Customize axis labels and theme
  labs(
    x = "Vignette outcome",
    y = "Main treatment\nbelief (Q1)"
  ) +
  # Adjust the x-axis labels
  scale_x_discrete(
    labels = c("R" = "Recovery", "NR" = "No Recovery")
  ) +
  # Adjust the y-axis limits and breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)  # Specify tick marks every 10 units
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    #axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## Panel C – Q6 by alternative treatment
fig2_panelC = ggplot(data_wide, aes(x = factor(vignette_alternative_treatment, levels = c("a", "c", "b")), y = Q6, fill = vignette_alternative_treatment)) +
  # Add bar plot
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
  # Add individual data points
  geom_jitter(width = 0.2, size = 1, alpha = 0.2, color = "black") +
  # Add error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  # Customize axis labels and theme
  labs(
    x = "Alternative treatment",
    y = "Alternative treatment\nbelief (Q6)"
  ) +
  # Adjust the x-axis labels
  scale_x_discrete(
    labels = c("a" = "OTC drug", "b" = "Rest", "c" = "Vitamins")
  ) +
  # Adjust the y-axis limits and breaks
  scale_y_continuous(
    breaks = seq(0, 100, by = 10)  # Specify tick marks every 10 units
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    #axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

# Combine panels to figure 2
fig2 = ggarrange(fig2_panelA, fig2_panelB, fig2_panelC,
                 labels = c("A","B", "C"),
                 ncol = 3, nrow = 1,
                 widths = c(1, 0.8, 1))

size_factor=1
ggsave('fig2.pdf',
       plot = fig2,
       path = fig_output_path,
       width = 200*size_factor,
       height = 65*size_factor,
       units = "mm",
       dpi = 300)

bayes_benchmark_posteriors = readRDS("figures/bayes_benchmark_posteriors.rds")

# bottom raw - empirical results
discounting_model_no_int_no_z = lmer(
  Q1 ~ Q3 + Q4 + Q6 +
    vignette_gender + gender + age +
    (1 | src_subject_id),
  data = recovery_data
)
# Panel C – Q1 by Q3
pred_Q3 = ggpredict(
  discounting_model_no_int_no_z,
  terms = "Q3 [all]"
)

p_Q3 = ggplot() +
  # raw participant data
  geom_point(
    data = recovery_data,
    aes(x = Q3, y = Q1),
    alpha = 0.25
  ) +
  # model-based confidence band
  geom_ribbon(
    data = pred_Q3,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.25
  ) +
  # model-based regression line
  geom_line(
    data = pred_Q3,
    aes(x = x, y = predicted),
    linewidth = 1
  ) +
  labs(
    x = "Belief that the patient\nhas COVID-19 (Q3)",
    y = "Belief in main treatment effectiveness (Q1)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Panel  D – Q1 by Q4
pred_Q4 = ggpredict(
  discounting_model_no_int_no_z,
  terms = "Q4 [all]"
)

p_Q4 = ggplot() +
  # raw participant data
  geom_point(
    data = recovery_data,
    aes(x = Q4, y = Q1),
    alpha = 0.25
  ) +
  # model-based confidence band
  geom_ribbon(
    data = pred_Q4,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.25
  ) +
  # model-based regression line
  geom_line(
    data = pred_Q4,
    aes(x = x, y = predicted),
    linewidth = 1
  ) +
  labs(
    x = "Belief in natural recovery\n(Q4)",
    y = "Belief in main treatment effectiveness (Q1)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Panel E – Q1 by Q6
pred_Q6 = ggpredict(
  discounting_model_no_int_no_z,
  terms = "Q6 [all]"
)

p_Q6 = ggplot() +
  # raw participant data
  geom_point(
    data = recovery_data,
    aes(x = Q6, y = Q1),
    alpha = 0.25
  ) +
  # model-based confidence band
  geom_ribbon(
    data = pred_Q6,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.25
  ) +
  # model-based regression line
  geom_line(
    data = pred_Q6,
    aes(x = x, y = predicted),
    linewidth = 1
  ) +
  labs(
    x = "Belief in alternative treatment\neffectiveness (Q6)",
    y = "Belief in main treatment effectiveness (Q1)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Combine panels to figure 3
fig3_lower = ggarrange(p_Q3, p_Q4, p_Q6,
                       labels = c("C","D", "E"),
                       ncol = 3, nrow = 1)
fig3 = ggarrange(bayes_benchmark_posteriors, fig3_lower,
                 ncol = 1, nrow = 2,
                 heights = c(1.45, 1))

size_factor=1
ggsave('fig3.pdf',
       plot = fig3,
       path = fig_output_path,
       width = 200*size_factor,
       height = 200*size_factor,
       units = "mm",
       dpi = 300)

## figure S1
fig_S1 = readRDS("figures/fig_S1_c1_upper_bound.rds")
ggsave('figS1.pdf',
       plot = fig_S1,
       path = fig_output_path,
       width = 140,
       height = 100,
       units = "mm",
       dpi = 300)
