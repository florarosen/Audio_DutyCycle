

# Load necessary libraries
library(lme4)        # For linear mixed-effects models
library(car)         # For ANOVA
library(emmeans)     # For estimated marginal means (EMMs) and post-hoc analysis
library(ggplot2)     # For plotting
library(nortest)     # For normality tests
library(dplyr)
library(tidyr)
library(ggplot2)
library(extrafont)
library(MuMIn)
library(patchwork)
library(effectsize)
library(cowplot)
library(officer)
library(flextable)
library(stringr)

##------------------------------------ FFT Analysis on the data up to 5Hz----------------------------------------------------------------------------------

# Set working directory to the location of the CSV file (!!!!! IF DOESNT WORK, OPEN MANUALLY !!!)

setwd("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis")


# Colors for DutyCycle
matlab_colors <- c("25" = "#4DBEEE", "50" = "#0072BD", "75" = "#EDB120", "100" = "#D95319")

# Read the CSV files into data frames
#EEG
df <- read.csv("AggregatedEEGDataForR.csv")

#Stimuli
stimCochmod <- read.csv("StimZscoresCochModelAN.csv")
FreqAmps <- read.csv("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis/Frequency_Amplitude_LongFormat.csv")


# Fix factors concisely
df <- df %>%
  mutate(
    SubjectID = factor(SubjectID),
    DutyCycle = factor(DutyCycle, levels = c("25", "50", "75", "100")),
    Meter = factor(Meter, levels = c("Rel", "UnRel", "All")),
    FrequencySet = factor(FrequencySet, levels = c("AllAmps", "Set9", "Set13"))
  )

# Filter data as per requirement
data_sets <- list(
  data_all = df %>% filter(FrequencySet == "AllAmps"),
  data_upto5Hz = df %>% filter(FrequencySet == "Set9"),
  data_upto5Hz_MR = df %>% filter(FrequencySet == "Set9", Meter == "Rel"),
  data_upto5Hz_MUR = df %>% filter(FrequencySet == "Set9", Meter == "UnRel"),
  data_upto15Hz = df %>% filter(FrequencySet == "Set13"),
  data_upto15Hz_MR = df %>% filter(FrequencySet == "Set13", Meter == "Rel"),
  data_upto15Hz_MUR = df %>% filter(FrequencySet == "Set13", Meter == "UnRel")
)

# Process stimulus : Map TrackID -> DutyCycle (%)
FreqAmps <- FreqAmps %>%
  mutate(DutyCycle = case_when(
    str_detect(TrackID, "^50ms")  ~ "25",
    str_detect(TrackID, "^100ms") ~ "50",
    str_detect(TrackID, "^150ms") ~ "75",
    str_detect(TrackID, "^200ms") ~ "100",
    TRUE ~ TrackID
  )) %>%
  select(-TrackID) %>%
  mutate(
    DutyCycle = factor(DutyCycle, levels = c("25","50","75","100")),
    Frequency = as.numeric(Frequency),
    Amplitude = as.numeric(Amplitude)
  ) %>%
  filter(Frequency <= 45)

# ------------------------------------------------------------
# Overall response magnitude per subject x duty cycle
# (sum of all harmonic magnitudes up to 45 Hz)
# ------------------------------------------------------------
overall_mag <- data_sets$data_all %>%
  filter(Frequency <= 45) %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(
    amp_total = sum(Amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    DutyCycle = factor(as.character(DutyCycle), levels = c("25","50","75","100")),
    DutyNum   = as.numeric(as.character(DutyCycle))
  )

lowfrex_mag <- data_sets$data_all %>%
  filter(Frequency <= 5) %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(
    amp_total = sum(Amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    DutyCycle = factor(as.character(DutyCycle), levels = c("25","50","75","100")),
    DutyNum   = as.numeric(as.character(DutyCycle))
  )

# AMPLITUDE DIFFERENCES ACROSS CONDITIONS (overall response)

# Low frex amplitude 
model_amp<- lmer(amp_total ~ DutyCycle + (1 | SubjectID), data = lowfrex_mag)
summary(model_amp)
anova_amp <- Anova(model_amp, test = "F")
print(anova_amp)

res <- resid(model_amp)
shapiro.test(res)

qqnorm(res)
qqline(res)

emms_amp <- emmeans(model_amp, ~ DutyCycle)
pairs_amp <- contrast(emms_amp, method = "pairwise", adjust = "bonferroni")
pairs_amp

# Total amplitude 
model_amptot <- lmer(amp_total ~ DutyCycle + (1 | SubjectID), data = overall_mag)
summary(model_amptot)
anova_amptot <- Anova(model_amptot, test = "F")
print(anova_amptot)

restot <- resid(model_amptot)
shapiro.test(restot)

qqnorm(restot)
qqline(restot)

emms_amptot <- emmeans(model_amptot, ~ DutyCycle)
pairs_amptot <- contrast(emms_amptot, method = "pairwise", adjust = "bonferroni")
pairs_amptot


# Fit LMM
m1<- lmer(prop_low ~ DutyNum + (1|SubjectID), data = prop_0_5_df, REML = TRUE)
Anova(m1, type = 2, test.statistic = "F")  # KR F-test

# Trend analysis
summary(m1)
Anova(m1, type = 2) 


# -------------Calculate the proportion of amplitude below 5Hz (including first 2 harmonics and 5Hz)--------------------
# --------------------------
# PARAMETERS (align with methods)
# --------------------------
f0 <- 1/2.4             # pattern repetition rate (Hz) = 0.4166667
max_freq <- 45          # upper limit (Hz)
chunk_hz <- 5           # 5-Hz chunks

# ----------------------------
# 2) Define the same frequency bands as EEG (0–5 includes 5; next band is >5–10, etc.)
# ----------------------------
add_band <- function(df) {
  df %>%
    mutate(
      Band = case_when(
        Frequency <= 5 ~ "0–5",
        Frequency > 5  & Frequency <= 10 ~ "5–10",
        Frequency > 10 & Frequency <= 15 ~ "10–15",
        Frequency > 15 & Frequency <= 20 ~ "15–20",
        Frequency > 20 & Frequency <= 25 ~ "20–25",
        Frequency > 25 & Frequency <= 30 ~ "25–30",
        Frequency > 30 & Frequency <= 35 ~ "30–35",
        Frequency > 35 & Frequency <= 40 ~ "35–40",
        Frequency > 40 & Frequency <= 45 ~ "40–45",
        TRUE ~ NA_character_
      ),
      Band = factor(Band, levels = c("0–5","5–10","10–15","15–20","20–25","25–30","30–35","35–40","40–45"))
    ) %>%
    filter(!is.na(Band))
}




# ----------------------------
# 3) Stimulus band proportions (ONE distribution per DutyCycle)
# ----------------------------

FreqAmps_banded <- add_band(FreqAmps)

stim_band_props <- FreqAmps_banded %>%
  group_by(DutyCycle, Band) %>%
  summarise(band_amp = sum(Amplitude, na.rm = TRUE), .groups = "drop") %>%
  group_by(DutyCycle) %>%
  mutate(total_amp = sum(band_amp, na.rm = TRUE),
         prop = band_amp / total_amp,
         percent = 100 * prop) %>%
  ungroup()

# Wide table for "Table 1" (Stimulus)
stim_table1 <- stim_band_props %>%
  select(Band, DutyCycle, percent) %>%
  pivot_wider(names_from = DutyCycle, values_from = percent) %>%
  arrange(Band)

write.csv(stim_table1, "Table1_Stimulus_BandProportions_0to45.csv", row.names = FALSE)


# EEG GROUP LEVEL mean spectrum across subjects (AllAmps)
eeg_group_long <- data_sets$data_all %>%
  filter(Frequency <= 45) %>%
  group_by(DutyCycle, Frequency) %>%
  summarise(Amplitude = mean(Amplitude, na.rm = TRUE), .groups = "drop")

eeg_group_banded <- add_band(eeg_group_long)

eeg_band_props_group <- eeg_group_banded %>%
  group_by(DutyCycle, Band) %>%
  summarise(band_amp = sum(Amplitude, na.rm = TRUE), .groups = "drop") %>%
  group_by(DutyCycle) %>%
  mutate(total_amp = sum(band_amp, na.rm = TRUE),
         prop = band_amp / total_amp,
         percent = 100 * prop) %>%
  ungroup()

eeg_table1 <- eeg_band_props_group %>%
  select(Band, DutyCycle, percent) %>%
  pivot_wider(names_from = DutyCycle, values_from = percent) %>%
  arrange(Band)

write.csv(eeg_table1, "Table1_EEG_Group_BandProportions_0to45.csv", row.names = FALSE)

# SUBJECT LEVEL ANALYSIS 
prop_0_5_df <- data_sets$data_all %>%
  filter(Frequency <= 45) %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(
    amp_0_5   = sum(Amplitude[Frequency <= 5], na.rm = TRUE),
    amp_total = sum(Amplitude, na.rm = TRUE),
    prop_low  = amp_0_5 / amp_total,
    .groups = "drop"
  ) %>%
  mutate(DutyNum = as.numeric(as.character(DutyCycle)))

# Prop of all bands : 
prop_bands <- data_sets$data_all %>%
  filter(Frequency <= 45) %>%
  mutate(
    Band = case_when(
      Frequency <= 5 ~ "0–5",
      Frequency > 5  & Frequency <= 10 ~ "5–10",
      Frequency > 10 & Frequency <= 15 ~ "10–15",
      Frequency > 15 & Frequency <= 20 ~ "15–20",
      Frequency > 20 & Frequency <= 25 ~ "20–25",
      Frequency > 25 & Frequency <= 30 ~ "25–30",
      Frequency > 30 & Frequency <= 35 ~ "30–35",
      Frequency > 35 & Frequency <= 40 ~ "35–40",
      Frequency > 40 & Frequency <= 45 ~ "40–45",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Band)) %>%
  group_by(SubjectID, DutyCycle, Band) %>%
  summarise(band_amp = sum(Amplitude, na.rm = TRUE), .groups = "drop") %>%
  group_by(SubjectID, DutyCycle) %>%
  mutate(
    total_amp = sum(band_amp, na.rm = TRUE),
    prop = band_amp / total_amp
  ) %>%
  ungroup()

# Order bands nicely
prop_bands <- prop_bands %>%
  mutate(
    Band = factor(Band, levels = c("0–5","5–10","10–15","15–20","20–25","25–30","30–35","35–40","40–45"))
  )


# Fit LMM
m1<- lmer(prop_low ~ DutyNum + (1|SubjectID), data = prop_0_5_df, REML = TRUE)
Anova(m1, type = 2, test.statistic = "F")  # KR F-test

# Trend analysis
summary(m1)
Anova(m1, type = 2) 

# report the trend as a slope + 95% CI, not just a p-value
b  <- fixef(m1)["DutyNum"]
se <- sqrt(vcov(m1)["DutyNum","DutyNum"])
range_change <- 75 * b                          # predicted change in proportion, 25% -> 100%
ci <- 75 * (b + c(-1.96, 1.96) * se)            # 95% CI on that change
round(c(estimate = range_change, lower = ci[1], upper = ci[2]), 3)
effectsize::eta_squared(Anova(m1, type = 2, test.statistic = "F"), partial = TRUE)

# report the trend as a slope + 95% CI, not just a p-value
b  <- fixef(m1)["DutyNum"]
se <- sqrt(vcov(m1)["DutyNum","DutyNum"])
range_change <- 75 * b                          # predicted change in proportion, 25% -> 100%
ci <- 75 * (b + c(-1.96, 1.96) * se)            # 95% CI on that change
round(c(estimate = range_change, lower = ci[1], upper = ci[2]), 3)
effectsize::eta_squared(Anova(m1, type = 2, test.statistic = "F"), partial = TRUE)


# (Optional) As categorical fixed effect / not reported
prop_0_5_df <- prop_0_5_df %>%
  mutate(DutyCycle = factor(DutyCycle, levels = c("25","50","75","100")))

m_cat <- lmer(prop_low ~ DutyCycle + (1 | SubjectID), data = prop_0_5_df, REML = TRUE)
car::Anova(m_cat, type = 2, test.statistic = "F")  # KR df
emm <- emmeans(m_cat, ~ DutyCycle, lmer.df = "kenward-roger")
pairs(emm, adjust = "bonferroni")





# --- Build stimulus predictor at the same granularity: DutyCycle x Band ---
stim_prop_by_band <- stim_band_props %>%
  select(DutyCycle, Band, stim_prop = prop) %>%
  mutate(
    DutyCycle = factor(as.character(DutyCycle), levels = c("25","50","75","100")),
    Band = factor(as.character(Band), levels = c("0–5","5–10","10–15","15–20","20–25","25–30","30–35","35–40","40–45"))
  )

# --- Ensure prop_bands keys match exactly ---
prop_bands2 <- prop_bands %>%
  mutate(
    DutyCycle = factor(as.character(DutyCycle), levels = c("25","50","75","100")),
    Band = factor(as.character(Band), levels = c("0–5","5–10","10–15","15–20","20–25","25–30","30–35","35–40","40–45"))
  ) %>%
  left_join(stim_prop_by_band, by = c("DutyCycle","Band"))

prop_bands2_low <- prop_bands2 %>%
  filter(Band %in% c("0–5", "5–10", "10–15")) %>%
  droplevels()

eeg0_5 <- prop_bands2_low %>%
  filter(Band == "0–5") %>%
  droplevels()

# Sanity check: stim_prop must be present for all rows
stopifnot(!any(is.na(prop_bands2$stim_prop)))

# --- LMM: does stimulus band proportion predict subject-level EEG band proportion? ---
# Random intercept for subjects; Band as factor controls for baseline band differences.

m0_stim <- lmer(prop ~ stim_prop + (1 | SubjectID), data = eeg0_5, REML = FALSE)
summary(m0_stim)
car::Anova(m0_stim, type = 2)

## Print tables of results : 
# Helper functions for formating

format_table1 <- function(df) {
  
  # Rename first column
  names(df)[1] <- "Frequency band (Hz)"
  
  # Rename duty-cycle columns
  names(df) <- gsub("^25$",  "25% duty cycle",  names(df))
  names(df) <- gsub("^50$",  "50% duty cycle",  names(df))
  names(df) <- gsub("^75$",  "75% duty cycle",  names(df))
  names(df) <- gsub("^100$", "100% duty cycle", names(df))
  
  # Format numeric columns: round to 2 decimals + %
  df <- df %>%
    mutate(across(-`Frequency band (Hz)`,
                  ~ sprintf("%.2f%%", .)))
  
  df
}

format_ft_table1 <- function(ft, caption_text) {
  
  ft <- ft %>%
    set_caption(caption_text) %>%
    font(fontname = "Calibri", part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    italic(italic = FALSE, part = "all") %>%
    bold(bold = FALSE, part = "all") %>%
    align(align = "left", part = "all") %>%
    autofit() %>%
    theme_booktabs()
  
  # Increase space between caption and table
  ft <- line_spacing(ft, space = 1.4, part = "caption")
  
  ft
}

stim_table1 <- format_table1(stim_table1)
eeg_table1 <- format_table1(eeg_table1)


ft_stim <- flextable(stim_table1) %>%
  format_ft_table1(
    "Table 1a. Auditory nerve model outputs magnitude distribution across frequency bands (0–45 Hz), by duty cycle."
  )

ft_eeg <- flextable(eeg_table1) %>%
  format_ft_table1(
    "Table 1b. EEG group-level outputs magnitude distribution across frequency bands (0–45 Hz), by duty cycle."
  )


doc <- read_docx() %>%
  body_add_flextable(ft_stim) %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("", style = "Normal") %>%   # extra spacing
  body_add_par("", style = "Normal") %>%   # extra spacing
  body_add_flextable(ft_eeg)

print(doc, target = "Table1_Stimulus_EEG_GroupLevel.docx")

##----------------EXTRACT AUDITORY NERVE MODEL ZSCORES----------------------------------------------
# Pivot stimCochmod and retain only the last value for each DutyCycle
stimCochmod_long <- stimCochmod %>%
  pivot_longer(cols = starts_with("X"),  # Pivot columns that start with "X"
               names_to = "DutyCycle",    # Create a new column "DutyCycle" from the column names
               values_to = "stimcochZscore") %>%
  mutate(DutyCycle = as.numeric(gsub("X", "", DutyCycle)))  # Convert "DutyCycle" to numeric by removing "X"


#--------------------------------ANALYSIS UP TO 15 HZ----------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
mean_Amplitude_15Hz <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores_15Hz <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')

# ---------------------------
# Gain covariate: summed noise-subtracted magnitudes across ALL FOI ≤15 Hz
# ---------------------------
lowfreq15_mag <- data_sets$data_all %>%
  filter(Frequency <= 15) %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(
    amp_total = sum(Amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SubjectID = as.character(SubjectID),
    DutyCycle = factor(as.character(DutyCycle), levels = c("25","50","75","100"))
  )


stimCochmod_third <- stimCochmod_long %>%
  group_by(DutyCycle) %>%
  mutate(idx = row_number()) %>%
  filter(idx == 3) %>%
  ungroup()

# Merge data frames effectively
merged_data_15Hz <- merge(mean_Amplitude_15Hz, mean_zscores_15Hz, by = c("SubjectID", "DutyCycle"))
final_table_15Hz <- merge(merged_data_15Hz, stimCochmod_third, by = "DutyCycle")

df_model_15Hz <- final_table_15Hz %>%
  select(SubjectID, DutyCycle, meanZscore, stimcochZscore) %>%
  rename(
    participantID = SubjectID,
    condition = DutyCycle,
    participantZscore = meanZscore,
    stimulusZscore = stimcochZscore
  ) %>%
  mutate(condition = factor(condition, levels = c("25", "50", "75", "100")))

## ---------------!!!!!!!!!!!!! Measuring beat-prominence in EEG across conditions--------------------------------------
## ---- Data hygiene: condition as categorical factor ----
## Methods: duty cycle modeled as *categorical* fixed effect
df_model_15Hz <- df_model_15Hz %>%
  mutate(
    participantID = factor(participantID),
    condition = factor(condition, levels = c("25","50","75","100"))  # categorical; ordered is optional
  )

## ---- Fit LMM (random intercepts for participants) ----
model_condition_15Hz <- lmer(participantZscore ~ condition + (1 | participantID),
                            data = df_model_15Hz, REML = TRUE)

summary(model_condition_15Hz)

## ---- Assumption checks (as stated in Methods: QQ + residual plots) ----
resid_cond <- resid(model_condition_15Hz)

# Normality test (optional; you already do it)
shapiro.test(resid_cond)

par(mfrow = c(1, 2))
plot(resid_cond ~ fitted(model_condition_15Hz),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(resid_cond); qqline(resid_cond, col = "red")
par(mfrow = c(1, 1))

## ---- Omnibus effect of condition (Type II Wald F-tests, Kenward–Roger df) ----
anova_results15Hz <- car::Anova(model_condition_15Hz, type = 2, test.statistic = "F")
print(anova_results15Hz)

## ---- Effect size (ηp²) consistent with what you report ----
eta_sq15Hz <- effectsize::eta_squared(anova_results15Hz, partial = TRUE)
print(eta_sq15Hz)

## ---- Estimated marginal means (EMMs) per condition ----
## (these are what you report in Results as EMMs)
emm_15Hz <- emmeans(model_condition_15Hz, ~ condition,
                   lmer.df = "kenward-roger")  # aligns df method with your ANOVA table
emm_15Hz

## ---- Post-hoc pairwise comparisons (Bonferroni) ----
pairwise_15Hz <- pairs(emm_15Hz, adjust = "bonferroni")
pairwise_15Hz

# Convert pairwise output to a data frame
pw15 <- as.data.frame(pairwise_15Hz)

# Keep only significant comparisons (after Bonferroni)
significant_contrasts_15Hz <- pw15 %>%
  dplyr::filter(p.value < 0.05) %>%
  tidyr::separate(contrast, into = c("contrast1_raw", "contrast2_raw"), sep = " - ") %>%
  dplyr::mutate(
    # Remove the "condition" prefix to match your DutyCycle labels ("25","50","75","100")
    contrast1 = gsub("^condition", "", contrast1_raw),
    contrast2 = gsub("^condition", "", contrast2_raw)
  ) %>%
  dplyr::select(contrast1, contrast2, p.value)

significant_contrasts_15Hz

## ---- Post-hoc polynomial (linear trend) ----
## IMPORTANT (matches your Results statement only if you explicitly do it):
## emmeans’ "poly" contrast treats levels as equally spaced in the order given.
## If you want a true linear trend across 25/50/75/100, coding is OK because steps are equal.
poly_15Hz <- contrast(emm_15Hz, method = "poly")
poly_15Hz

## ---- Test each EEG EMM against the corresponding auditory-nerve model z-score (null) ----
## Methods 2.8.1: "each EMM tested against corresponding auditory nerve model z-score"
## Ensure the null vector is named with the same condition levels as emm_5Hz
stim_null <- stimCochmod_third$stimcochZscore
names(stim_null) <- as.character(stimCochmod_third$DutyCycle)

# Keep only the required levels and in the same order
stim_null <- stim_null[levels(df_model_15Hz$condition)]

test_results_15Hz <- test(emm_15Hz, null = stim_null)
test_results_15Hz

## Measuring stimulus beat-prominence impact on EEG beat-prominence in EEG across conditions
# Fit linear mixed-effects models
## ================== Methods 2.8.2: EEG ~ auditory-nerve model z-score (continuous) ==================
model_participantZscore_Lmer_intercept_15Hz <- lmer(participantZscore ~ stimulusZscore + (1 | participantID),
                                                   data = df_model_15Hz, REML = TRUE)

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_participantZscore_Lmer_intercept_15Hz) ~ fitted(model_participantZscore_Lmer_intercept_15Hz))
qqnorm(resid(model_participantZscore_Lmer_intercept_15Hz))
qqline(resid(model_participantZscore_Lmer_intercept_15Hz), col = "red")

# Summary and diagnostics
summary(model_participantZscore_Lmer_intercept_15Hz)
anova_results_stim15Hz <- car::Anova(model_participantZscore_Lmer_intercept_15Hz, type = 2, test.statistic = "F")
print(anova_results_stim15Hz)

# Compute partial eta-squared directly
eta_squared_results_stim15Hz <- effectsize::eta_squared(anova_results_stim15Hz, partial = TRUE)
print(eta_squared_results_stim15Hz)


# -------------------------------   PLOT Figure S2: Magnitude spectrum up to 15 Hz (with EMM-vs-model + pairwise stars) ---------------------------------------

dodge_width <- 0.15
jitter_dodge_position <- position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width)

# --- Build emm2_df_15Hz for plotting (EMMs + EEG-vs-model significance + model z-score) ---
emm2_df_15Hz <- as.data.frame(emm_15Hz)
emm2_df_15Hz$DutyCycle <- factor(as.character(emm2_df_15Hz$condition), levels = c("25","50","75","100"))
tr15 <- as.data.frame(test_results_15Hz)
emm2_df_15Hz$significance <- ifelse(tr15$p.value < 0.001, "***",
                             ifelse(tr15$p.value < 0.01, "**",
                                    ifelse(tr15$p.value < 0.05, "*", "")))
emm2_df_15Hz$stimcochZscore <- as.numeric(stim_null[as.character(emm2_df_15Hz$DutyCycle)])

plotFFT_15Hz <- ggplot(final_table_15Hz,
                         aes(x = DutyCycle, y = meanZscore, group = DutyCycle,
                             colour = DutyCycle, fill = DutyCycle)) +

  # Jittered points
  geom_point(position = jitter_dodge_position, size = 2,
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  # Points for stimcochZscore
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +

  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA,
               outlier.shape = NA, position = position_dodge(width = dodge_width)) +

  # Vertical lines for EMM vs stimcochZscore
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, y = emmean, yend = stimcochZscore),
               color = "black", size = 0.4, data = emm2_df_15Hz) +
  # Line ticks
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = emmean, yend = emmean),
               color = "black", size = 0.4, data = emm2_df_15Hz) +
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = stimcochZscore, yend = stimcochZscore),
               color = "black", size = 0.4, data = emm2_df_15Hz) +

  # Significance stars for EEG vs model (same as Figure 3A)
  geom_text(aes(x = as.numeric(DutyCycle) - 0.48, y = (emmean + stimcochZscore) / 2, label = significance),
            color = "black", size = 14, data = emm2_df_15Hz) +

  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +

  # Y-axis limits, breaks, and zero line
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.45)) +
  scale_x_discrete(expand = expansion(add = c(0.85, 0.6))) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +

  # Labels and title
  labs(title = "Magnitude Spectrum-Based Analysis up to 15Hz ",
       y = "EEG \nbeat-related z-score (a.u)",
       x = "Duty Cycle (%)") +

  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  )

# --- Pairwise significance bracket(s) for the 15 Hz analysis (same logic as Figure 3A) ---
x_positions_15 <- setNames(seq_along(levels(final_table_15Hz$DutyCycle)), levels(final_table_15Hz$DutyCycle))
contrast_height_15 <- -0.55 + 0.90 * (1.45 - (-0.55))
height_step_15 <- 0.03 * (1.45 - (-0.55))

for (i in seq_len(nrow(significant_contrasts_15Hz))) {
  c1_pos <- x_positions_15[as.character(significant_contrasts_15Hz$contrast1[i])]
  c2_pos <- x_positions_15[as.character(significant_contrasts_15Hz$contrast2[i])]
  label <- if (significant_contrasts_15Hz$p.value[i] < 0.001) {
    "***"
  } else if (significant_contrasts_15Hz$p.value[i] < 0.01) {
    "**"
  } else if (significant_contrasts_15Hz$p.value[i] < 0.05) {
    "*"
  } else { "" }
  plotFFT_15Hz <- plotFFT_15Hz +
    geom_segment(aes(x = c1_pos, xend = c2_pos, y = contrast_height_15, yend = contrast_height_15),
                 inherit.aes = FALSE, size = 0.4, color = "black") +
    geom_text(aes(x = (c1_pos + c2_pos) / 2, y = contrast_height_15 + 0.015 * (1.45 - (-0.55)), label = label),
              inherit.aes = FALSE, size = 18, color = "black")
  contrast_height_15 <- contrast_height_15 + height_step_15
}

print(plotFFT_15Hz)

ggsave("FFT_upto15Hz.svg", plot = plotFFT_15Hz, width = 25, height = 20, units = "cm", device = "svg")
ggsave("Figure_S2.pdf", plot = plotFFT_15Hz, width = 25, height = 20, units = "cm", device = cairo_pdf)
ggsave("Figure_S2.png", plot = plotFFT_15Hz, width = 25, height = 20, units = "cm", dpi = 300, device = "png")
ggsave("Figure_S2.tiff", plot = plotFFT_15Hz, width = 25, height = 20, units = "cm", dpi = 300, device = "tiff", compression = "lzw")


#--------------------------------ANALYSIS UP TO 5 HZ----------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
# Calculate mean amplitude and zscores
mean_Amplitude_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')


# Keep the first stimcochZscore for each DutyCycle for 5HZ ANALYSIS
stimCochmod_first <- stimCochmod_long %>%
  group_by(DutyCycle) %>%
  slice_head(n = 1)  # Keep the last value for each DutyCycle

# Merge data frames effectively
merged_data_5Hz <- merge(mean_Amplitude_5Hz, mean_zscores_5Hz, by = c("SubjectID", "DutyCycle"))
final_table_5Hz <- merge(merged_data_5Hz, stimCochmod_first, by = "DutyCycle")

df_model_5Hz <- final_table_5Hz %>%
  select(SubjectID, DutyCycle, meanZscore, stimcochZscore) %>%
  rename(
    participantID = SubjectID,
    condition = DutyCycle,
    participantZscore = meanZscore,
    stimulusZscore = stimcochZscore
  ) %>%
  mutate(condition = factor(condition, levels = c("25", "50", "75", "100")))

## ---------------!!!!!!!!!!!!! Measuring beat-prominence in EEG across conditions--------------------------------------
## ---- Data hygiene: condition as categorical factor ----
## Methods: duty cycle modeled as *categorical* fixed effect
df_model_5Hz <- df_model_5Hz %>%
  mutate(
    participantID = factor(participantID),
    condition = factor(condition, levels = c("25","50","75","100"))  # categorical; ordered is optional
  )

## ---- Fit LMM (random intercepts for participants) ----
model_condition_5Hz <- lmer(participantZscore ~ condition + (1 | participantID),
                            data = df_model_5Hz, REML = TRUE)

summary(model_condition_5Hz)

## ---- Assumption checks (as stated in Methods: QQ + residual plots) ----
resid_cond <- resid(model_condition_5Hz)

# Normality test (optional; you already do it)
shapiro.test(resid_cond)

par(mfrow = c(1, 2))
plot(resid_cond ~ fitted(model_condition_5Hz),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(resid_cond); qqline(resid_cond, col = "red")
par(mfrow = c(1, 1))

## ---- Omnibus effect of condition (Type II Wald F-tests, Kenward–Roger df) ----
anova_results <- car::Anova(model_condition_5Hz, type = 2, test.statistic = "F")
print(anova_results)

## ---- Effect size (ηp²) consistent with what you report ----
eta_sq <- effectsize::eta_squared(anova_results, partial = TRUE)
print(eta_sq)

## ---- Estimated marginal means (EMMs) per condition ----
## (these are what you report in Results as EMMs)
emm_5Hz <- emmeans(model_condition_5Hz, ~ condition,
                   lmer.df = "kenward-roger")  # aligns df method with your ANOVA table
emm_5Hz

## ---- Post-hoc pairwise comparisons (Bonferroni) ----
pairwise_5Hz <- pairs(emm_5Hz, adjust = "bonferroni")
pairwise_5Hz

## ---- Post-hoc polynomial (linear trend) ----
## IMPORTANT (matches your Results statement only if you explicitly do it):
## emmeans’ "poly" contrast treats levels as equally spaced in the order given.
## If you want a true linear trend across 25/50/75/100, coding is OK because steps are equal.
poly_5Hz <- contrast(emm_5Hz, method = "poly")
poly_5Hz

## ---- Test each EEG EMM against the corresponding auditory-nerve model z-score (null) ----
## Methods 2.8.1: "each EMM tested against corresponding auditory nerve model z-score"
## Ensure the null vector is named with the same condition levels as emm_5Hz
stim_null <- stimCochmod_first$stimcochZscore
names(stim_null) <- as.character(stimCochmod_first$DutyCycle)

# Keep only the required levels and in the same order
stim_null <- stim_null[levels(df_model_5Hz$condition)]

test_results_5Hz <- test(emm_5Hz, null = stim_null)
test_results_5Hz

# ----- Make a results table and export to word
# Create clean table
table_pub_5Hz <- as.data.frame(test_results_5Hz) %>%
  transmute(
    `Duty cycle (%)` = condition,
    `EEG beat-related z-score (EMM)` = sprintf("%.3f", emmean),
    `SE` = sprintf("%.3f", SE),
    `df` = sprintf("%.1f", df),
    `Auditory nerve model z-score` = sprintf("%.3f", null),
    `t` = sprintf("%.2f", t.ratio),
    `p-value` = "< 1 × 10⁻⁴"
  )

ft_5Hz <- flextable(table_pub_5Hz) %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  set_caption(
    caption = "Comparison of EEG beat-related z-scores against auditory nerve model predictions (magnitude spectrum-based analysis)."
  ) %>%
  add_footer_lines(
    values = "Table 1. Statistics are based on estimated marginal means from linear mixed-effects models with Kenward–Roger degrees of freedom. P-values correspond to tests of EEG beat-related z-scores against the corresponding auditory nerve model z-score."
  )

doc <- read_docx() %>%
  body_add_flextable(ft_5Hz)

print(doc, target = "Table_S1_EEG_vs_AuditoryNerve.docx")

## --- For Plotting : 
  # Prepare significance stars
    emm2_df_5Hz <- as.data.frame(emm_5Hz)
    emm2_df_5Hz$DutyCycle <- emm2_df_5Hz$condition
    emm2_df_5Hz$significance <- ifelse(test_results_5Hz$p.value < 0.001, "***",
                                   ifelse(test_results_5Hz$p.value < 0.01, "**",
                                          ifelse(test_results_5Hz$p.value < 0.05, "*", "")))


# Ensure that DutyCycle is of the same type (character) in both data frames
emm2_df_5Hz$DutyCycle <- as.character(emm2_df_5Hz$DutyCycle)
stimCochmod_first$DutyCycle <- as.character(stimCochmod_first$DutyCycle)


# Merge stimCochmod_last into emm2_df for consistent stimcochZscore usage
emm2_df_5Hz <- emm2_df_5Hz %>%
  left_join(stimCochmod_first %>% select(DutyCycle, stimcochZscore), by = "DutyCycle")
emm2_df_5Hz$DutyCycle <- as.factor(emm2_df_5Hz$DutyCycle)




# Calculate pairwise contrasts
contrast_results_5Hz <- contrast(emm_5Hz, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE)

# Filter out significant contrasts (p < 0.05)
significant_contrasts_5Hz <- contrast_results_5Hz %>%
  filter(p.value < 0.05)

# Convert contrast results to a data frame (if it's not already)
contrast_results_df_5Hz <- as.data.frame(contrast_results_5Hz)

# Manually split the 'contrast' column into two new columns: 'contrast1' and 'contrast2'
contrast_split_5Hz <- sapply(strsplit(contrast_results_df_5Hz$contrast, " - "), function(x) x)

# Create new columns 'contrast1' and 'contrast2'
contrast_results_df_5Hz$contrast1 <- as.numeric(gsub("condition", "", contrast_split_5Hz[1, ]))
contrast_results_df_5Hz$contrast2 <- as.numeric(gsub("condition", "", contrast_split_5Hz[2, ]))

# Filter out significant contrasts (p < 0.05)
significant_contrasts_5Hz <- contrast_results_df_5Hz %>%
  filter(p.value < 0.05)

# Ensure 'DutyCycle' is a factor in the relevant data frames
final_table_5Hz$DutyCycle <- as.factor(final_table_5Hz$DutyCycle)
emm2_df_5Hz$DutyCycle <- as.factor(emm2_df_5Hz$DutyCycle)  # Convert 'DutyCycle' in emm2_df to factor
# Reorder the levels of 'DutyCycle' in emm2_df
emm2_df_5Hz$DutyCycle <- factor(emm2_df_5Hz$DutyCycle, levels = c("25", "50", "75", "100"))


## Measuring stimulus beat-prominence impact on EEG beat-prominence in EEG across conditions
# Fit linear mixed-effects models
## ================== Methods 2.8.2: EEG ~ auditory-nerve model z-score (continuous) ==================
model_participantZscore_Lmer_intercept_5Hz <- lmer(participantZscore ~ stimulusZscore + (1 | participantID),
                          data = df_model_5Hz, REML = TRUE)

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_participantZscore_Lmer_intercept_5Hz) ~ fitted(model_participantZscore_Lmer_intercept_5Hz))
qqnorm(resid(model_participantZscore_Lmer_intercept_5Hz))
qqline(resid(model_participantZscore_Lmer_intercept_5Hz), col = "red")

# Summary and diagnostics
summary(model_participantZscore_Lmer_intercept_5Hz)
anova_results_stim <- car::Anova(model_participantZscore_Lmer_intercept_5Hz, type = 2, test.statistic = "F")
print(anova_results_stim)

# Compute partial eta-squared directly
eta_squared_results_stim <- effectsize::eta_squared(anova_results_stim, partial = TRUE)
print(eta_squared_results_stim)

## Test for EEG-AN gap directly : 
df_model_5Hz$deltaZ <- df_model_5Hz$participantZscore - df_model_5Hz$stimulusZscore
lmer(deltaZ ~ condition + (1|participantID), data = df_model_5Hz)


# -------------------------------   PLOTS ------------------------------------------------------------------------------------------
# Create the complete plot with boxplot, jitter, stimcochZscore points, vertical lines, and ticks

# Set plotting parameters
# Set dodge width for consistent spacing (optional if you want to adjust)
dodge_width <- 0.15
jitter_dodge_position <- position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width)

plotFFT <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +  # Now using meanZscore for y-axis
  
  
  # Jittered points
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = jitter_dodge_position, size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  # Points for stimcochZscore (black points)
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  
  # Vertical lines for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, y = emmean, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  
  # Line ticks for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = emmean, yend = emmean), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = stimcochZscore, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  
 
  # Add significance stars for EMM and stimcochZscore differences
geom_text(aes(x = as.numeric(DutyCycle) - 0.48, y = (emmean + stimcochZscore) / 2, label = significance),
          color = "black", size = 14, data = emm2_df_5Hz) +
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  ## Y-axis limits, breaks, and zero line
  #scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.25)) +
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.45)) +
  scale_x_discrete(expand = expansion(add = c(0.85, 0.6))) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
 
  
   # Labels and title
  labs(title = "Magnitude spectrum", y = "EEG \nbeat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  )
# Map DutyCycle levels to numeric positions on the x-axis
x_positions <- setNames(seq_along(levels(final_table_5Hz$DutyCycle)), levels(final_table_5Hz$DutyCycle))
max_y_fft <- max(final_table_5Hz$meanZscore, emm2_df_5Hz$stimcochZscore, na.rm = TRUE)
# Anchor pairwise brackets to a fixed fraction of the y-range (same fraction used
# in the ACF panel) so the bracket + star sit at the same visual level in both panels.
contrast_height <- -0.55 + 0.90 * (1.45 - (-0.55))   # base bracket at 90% of y-range
height_step <- 0.03 * (1.45 - (-0.55))               # 3% of range per stacked bracket

# Loop through significant contrasts
for (i in seq_len(nrow(significant_contrasts_5Hz))) {
  c1_label <- as.character(significant_contrasts_5Hz$contrast1[i])
  c2_label <- as.character(significant_contrasts_5Hz$contrast2[i])
  
  c1_pos <- x_positions[c1_label]
  c2_pos <- x_positions[c2_label]
  
  # **Corrected star logic**
  label <- if (significant_contrasts_5Hz$p.value[i] < 0.001) {
    "***"
  } else if (significant_contrasts_5Hz$p.value[i] < 0.01) {
    "**"
  } else if (significant_contrasts_5Hz$p.value[i] < 0.05) {
    "*"
  } else {
    ""
  }
  
  # Horizontal line
  plotFFT <- plotFFT +
    geom_segment(aes(x = c1_pos, xend = c2_pos, y = contrast_height, yend = contrast_height),
                 inherit.aes = FALSE, size = 0.4, color = "black") +
    geom_text(aes(x = (c1_pos + c2_pos) / 2, y = contrast_height + 0.015 * (1.45 - (-0.55)), label = label),
              inherit.aes = FALSE, size = 18, color = "black")
  
  contrast_height <- contrast_height + height_step
}
  

# Display the complete plot
print(plotFFT)

# -------------------------------   PLOTS (no stars, with stimulus squares) ---------------------------------------

dodge_width <- 0.15
jitter_dodge_position <- position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width)

plotFFT_nostars <- ggplot(final_table_5Hz,
                  aes(x = DutyCycle, y = meanZscore, group = DutyCycle,
                      colour = DutyCycle, fill = DutyCycle)) +
  
  # Jittered points
  geom_point(position = jitter_dodge_position, size = 2,
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA,
               outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  # Points for stimcochZscore (black points)
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  

  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits, breaks, and zero line
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.45)) +
  scale_x_discrete(expand = expansion(add = c(0.85, 0.6))) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(title = "Magnitude spectrum",
       y = "EEG \nbeat-related z-score (a.u)",
       x = "Duty Cycle (%)") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  )



# Display the complete plot
print(plotFFT_nostars)

# ------------------------STATS AND PLOTS FOR THE AUTOCORRELATION ANALYSIS ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------


# Read the CSV file into a data frame
df_acf <- read.csv("acf_z_snr_4stats.csv")


df_acf <- df_acf %>%
  rename(DutyCycle = condition,  # Renaming 'condition' to 'DutyCycle'
         SubjectID = subject    # Renaming 'subject' to 'SubjectID'
  )    

# Convert 'condition' to a factor with the new levels
df_acf$DutyCycle <- factor(df_acf$DutyCycle, levels = c("25%", "50%", "75%", "100%"), labels = c("25", "50", "75", "100"))
#Factor Subjects
df_acf$SubjectID <- as.factor(df_acf$SubjectID)

# Read the CSV file into a data frame
#stimCochmod <- read.csv('MR_z_stimCochmod.csv')
# Create the data frame
stimCochmod <- data.frame(
  `25%` = -0.3295543,
  `50%` = -0.3607533,
  `75%` = -0.3717696,
  `100%` = -0.3984415
)


# Pivot 'stimCochmod' from wide to long format effectively
pivot_long <- function(data, obs_row, prefix, value_name) {
  pivot_longer(data[obs_row, ], cols = starts_with(prefix), 
               names_to = "DutyCycle", values_to = value_name) %>%
    mutate(DutyCycle = as.numeric(gsub(prefix, "", DutyCycle)))
}

stimCochmod_long <- pivot_long(stimCochmod, 1, "X", "stimcochZscore")

# Merge data frames effectively
final_table_acf <- merge(df_acf, stimCochmod_long, by = "DutyCycle")


## ----------------------------LINEAR MIXED OF ACF DATA---------------------------------------------
## ---------------------------- Data prep (ACF) ---------------------------------------------
df_model_acf <- final_table_acf %>%
  select(SubjectID, DutyCycle, z_meter_rel, stimcochZscore) %>%
  rename(participantID = SubjectID,
         condition     = DutyCycle,
         participantZscore = z_meter_rel,
         stimulusZscore    = stimcochZscore) %>%
  mutate(
    participantID = factor(participantID),
    condition = factor(condition, levels = c("25","50","75","100"))
  )

## ---------------------------- 1) Duty cycle effect (categorical LMM) ----------------------
model_condition_acf <- lmer(participantZscore ~ condition + (1 | participantID),
                            data = df_model_acf, REML = TRUE)

summary(model_condition_acf)

## Assumption checks (same as FFT)
resid_acf <- resid(model_condition_acf)
shapiro.test(resid_acf)

par(mfrow = c(1,2))
plot(resid_acf ~ fitted(model_condition_acf),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
qqnorm(resid_acf); qqline(resid_acf, col = "red")
par(mfrow = c(1,1))

## Omnibus test (Type II Wald F with Kenward–Roger df)
anova_results_acf <- car::Anova(model_condition_acf, type = 2, test.statistic = "F")
anova_results_acf

## Effect size (ηp² from Wald F)
eta_sq_acf <- effectsize::eta_squared(anova_results_acf, partial = TRUE)
eta_sq_acf

## ---------------------------- 2) EMMs + post-hocs -----------------------------------------
emm_acf <- emmeans(model_condition_acf, ~ condition, lmer.df = "kenward-roger")
emm_acf

## Pairwise (Bonferroni)
pairs_acf <- pairs(emm_acf, adjust = "bonferroni")
pairs_acf

## Polynomial contrasts (linear / quadratic / cubic) on ordered levels
poly_acf <- contrast(emm_acf, method = "poly")
poly_acf

## ---------------------------- 3) EEG vs auditory-nerve model per condition ----------------
## IMPORTANT: use the ACF auditory nerve model z-scores that correspond to each duty cycle
stim_null_acf <- c("25" = -0.3295543,
                   "50" = -0.3607533,
                   "75" = -0.3717696,
                   "100"= -0.3984415)

## Ensure order matches factor levels
stim_null_acf <- stim_null_acf[levels(df_model_acf$condition)]

test_vs_stim_acf <- test(emm_acf, null = stim_null_acf)
test_vs_stim_acf

## ---------------------------- 4) Optional: test EMMs vs 0 (only if you report it) ----------
## Only include this if you actually discuss "above zero" in Results.
test_vs_zero_acf <- test(emm_acf, null = 0)
test_vs_zero_acf

## Measuring stimulus beat-prominence impact on EEG beat-prominence in EEG across conditions
# Fit linear mixed-effects models
## ================== Methods 2.8.2: EEG ~ auditory-nerve model z-score (continuous) ==================
model_participantZscore_Lmer_intercept_5Hz <- lmer(participantZscore ~ stimulusZscore + (1 | participantID),
                                                   data = df_model_5Hz, REML = TRUE)

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_participantZscore_Lmer_intercept_5Hz) ~ fitted(model_participantZscore_Lmer_intercept_5Hz))
qqnorm(resid(model_participantZscore_Lmer_intercept_5Hz))
qqline(resid(model_participantZscore_Lmer_intercept_5Hz), col = "red")

# Summary and diagnostics
summary(model_participantZscore_Lmer_intercept_5Hz)
anova_results_stim <- car::Anova(model_participantZscore_Lmer_intercept_5Hz, type = 2, test.statistic = "F")
print(anova_results_stim)

# Compute partial eta-squared directly
eta_squared_results_stim <- effectsize::eta_squared(anova_results_stim, partial = TRUE)
print(eta_squared_results_stim)

## Test for EEG-AN gap directly : 
df_model_5Hz$deltaZ <- df_model_5Hz$participantZscore - df_model_5Hz$stimulusZscore
lmer(deltaZ ~ condition + (1|participantID), data = df_model_5Hz)


#---------PLOTS WITH STARS AND STIM -----------------------------------------
# Create the first plot with only stim Z-score points (simplified plot)


# Create the complete plot
# --- Build emm2_df_acf for plotting (EMMs + significance + model z-score) ---
emm2_df_acf <- as.data.frame(emm_acf)
emm2_df_acf$DutyCycle <- factor(as.character(emm2_df_acf$condition), levels = c("25","50","75","100"))
tv_acf <- as.data.frame(test_vs_stim_acf)
emm2_df_acf$significance <- ifelse(tv_acf$p.value < 0.001, "***",
                             ifelse(tv_acf$p.value < 0.01, "**",
                                    ifelse(tv_acf$p.value < 0.05, "*", "")))
emm2_df_acf$stimcochZscore <- as.numeric(stim_null_acf[as.character(emm2_df_acf$DutyCycle)])

plotACF <- ggplot(final_table_acf, aes(x = DutyCycle, y = z_meter_rel, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +
  
  # Jittered points with fixed position
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = jitter_dodge_position, size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  # Points for stimcochZscore
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  
  # Vertical lines for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, y = emmean, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df_acf) +
  
  # Line ticks for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = emmean, yend = emmean), 
               color = "black", size = 0.4, data = emm2_df_acf) +
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = stimcochZscore, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df_acf) +
  geom_hline(yintercept = 0.0, linetype = "dashed", colour = "black", size = 0.5) +
  
  # Add significance stars for EMM and stimcochZscore differences
geom_text(aes(x = as.numeric(DutyCycle) - 0.48, y = (emmean + stimcochZscore) / 2 - 0.18, label = significance),
          color = "black", size = 14, data = emm2_df_acf) +
  
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-1.55, 1.85)) +
  scale_x_discrete(expand = expansion(add = c(0.85, 0.6))) +
  coord_cartesian(clip = "off") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(title = "Autocorrelation", y = "EEG \nbeat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  ) 

# **Add horizontal contrast bars for plotACF**
x_positions_acf <- setNames(seq_along(levels(final_table_acf$DutyCycle)), levels(final_table_acf$DutyCycle))
max_y_acf <- max(final_table_acf$z_meter_rel, emm2_df_acf$stimcochZscore, na.rm = TRUE)
# Same visual level and step as the magnitude panel (same fractions of the y-range).
contrast_height_acf <- -1.55 + 0.90 * (1.85 - (-1.55))   # base bracket at 90% of y-range
height_step_acf <- 0.03 * (1.85 - (-1.55))               # 3% of range per stacked bracket

# --- Build significant_contrasts_acf for the pairwise bracket (mirrors FFT) ---
contrast_results_acf <- as.data.frame(pairs_acf)
contrast_split_acf <- sapply(strsplit(contrast_results_acf$contrast, " - "), function(x) x)
contrast_results_acf$contrast1 <- as.numeric(gsub("condition", "", contrast_split_acf[1, ]))
contrast_results_acf$contrast2 <- as.numeric(gsub("condition", "", contrast_split_acf[2, ]))
significant_contrasts_acf <- contrast_results_acf %>% dplyr::filter(p.value < 0.05)

for (i in seq_len(nrow(significant_contrasts_acf))) {
  c1_label <- as.character(significant_contrasts_acf$contrast1[i])
  c2_label <- as.character(significant_contrasts_acf$contrast2[i])
  
  c1_pos <- x_positions_acf[c1_label]
  c2_pos <- x_positions_acf[c2_label]
  
  # **Corrected star logic**
  label <- if (significant_contrasts_acf$p.value[i] < 0.001) {
    "***"
  } else if (significant_contrasts_acf$p.value[i] < 0.01) {
    "**"
  } else if (significant_contrasts_acf$p.value[i] < 0.05) {
    "*"
  } else {
    ""
  }
  
  plotACF <- plotACF +
    geom_segment(aes(x = c1_pos, xend = c2_pos, y = contrast_height_acf, yend = contrast_height_acf),
                 inherit.aes = FALSE, size = 0.4, color = "black") +
    geom_text(aes(x = (c1_pos + c2_pos) / 2, y = contrast_height_acf + 0.015 * (1.85 - (-1.55)), label = label),
              inherit.aes = FALSE, size = 18, color = "black")
  
  contrast_height_acf <- contrast_height_acf + height_step_acf
}
# Display the plot
print(plotACF)

# -------------------------------   PLOT ACF (no EMM stars, with stimulus squares) ---------------------------------------

plotACF_nostars <- ggplot(final_table_acf,
                  aes(x = DutyCycle, y = z_meter_rel, group = DutyCycle,
                      colour = DutyCycle, fill = DutyCycle)) +
  
  # Jittered points
  geom_point(position = jitter_dodge_position, size = 2,
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA,
               outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  
  # Points for stimcochZscore
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  
  
  # Zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-1.55, 1.85)) +
  scale_x_discrete(expand = expansion(add = c(0.85, 0.6))) +
  coord_cartesian(clip = "off") +
  
  # Labels and title
  labs(title = "Autocorrelation",
       y = "EEG \nbeat-related z-score (a.u)",
       x = "Duty Cycle (%)") +
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  )



print(plotACF_nostars)


## -----------------------------------------------Overall plots----------------------------------------------
## ------------------------------------------------------------------------------------------------------------------------------------------

# Adjust theme
theme_adjusted <- theme(
  plot.margin = margin(t = 30, r = 10, b = 20, l = 40), # More space at the top
  plot.title = element_text(hjust = 0.5, size = 30, face = "bold", margin = margin(b = 20)), # Extra space below title
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  text = element_text(size = 28),
  axis.title = element_text(size = 28),
  axis.text = element_text(size = 28),
  legend.text = element_text(size = 28),
  axis.line = element_line(color = "black", size = 1.2)  # Add physical axis lines
)

# Apply theme adjustments
plotFFT <- plotFFT + theme_adjusted
plotACF <- plotACF + theme_adjusted
plotFFT_nostars <- plotFFT_nostars  + theme_adjusted
plotACF_nostars <- plotACF_nostars  + theme_adjusted

combined_plot <- plot_grid(
  plotFFT, 
  plotACF, 
  labels = c("A.", "B."), 
  label_size = 30, 
  ncol = 2, 
  align = "h",
  rel_widths = c(1, 1)  # Equal panel widths
)


combined_plot_nostars <- plot_grid(
  plotFFT_nostars, 
  plotACF_nostars, 
  labels = c("A.", "B."), 
  label_size = 30, 
  ncol = 2, 
  align = "h", 
  rel_widths = c(1, 1)  # Give ACF slightly more width to prevent overlap
)

# Display the combined plots
print(combined_plot)
print(combined_plot_nostars)


# Save the combined plot
ggsave("FFT_ACF_EEG_combined.svg", plot = combined_plot, width = 50, height = 20, units = "cm", device = "svg")
ggsave("FFT_ACF_EEG_combined.png", plot = combined_plot, width = 50, height = 20, units = "cm", dpi = 300, device = "png")
ggsave("FFT_ACF_EEG_combined_nostars.svg", plot = combined_plot_nostars, width = 50, height = 20, units = "cm", device = "svg")

# --- Production figure files for submission (this is Figure 3; no title baked in) ---
ggsave("Figure_3.pdf", plot = combined_plot, width = 50, height = 20, units = "cm", device = cairo_pdf)
ggsave("Figure_3.png", plot = combined_plot, width = 50, height = 20, units = "cm", dpi = 300, device = "png")
ggsave("Figure_3.tiff", plot = combined_plot, width = 50, height = 20, units = "cm", dpi = 300, device = "tiff", compression = "lzw")

