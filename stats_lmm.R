# ------------------------STATS AND PLOTS FOR THE FFT MAGNITUDE SPECTRUM ANALYSIS ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------
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
library(patchwork) 
library(MuMIn)
library(effectsize)

## -------------------------------Define all parameters
# Set working directory to the location of the CSV file (!!!!! IF DOESNT WORK, OPEN MANUALLY !!!)

setwd("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis")
source("R_rainclouds.R") # (!!!!! IF DOESNT FIND, OPEN THE DIRECTORY MANUALLY !!!)
source("summarySE.R")

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

## --------------------------------------------------------------------------------------------
### 1 - Analyze the distribution over the whole spectrum 
## Stimuli 
# Modify the TrackID column to DutyCycle with just the percentage as values
FreqAmps <- FreqAmps %>%
  mutate(DutyCycle = case_when(
    str_detect(TrackID, "^50ms") ~ "25",
    str_detect(TrackID, "^100ms") ~ "50",
    str_detect(TrackID, "^150ms") ~ "75",
    str_detect(TrackID, "^200ms") ~ "100",
    TRUE ~ TrackID  # Retain any TrackID values that don't match specified patterns
  )) %>%
  select(-TrackID)  # Remove the original TrackID column


# Define harmonic frequencies and categorize frequencies
harmonics <- seq(5, max(FreqAmps$Frequency, na.rm = TRUE), by = 5)  # Multiples of 5 Hz
harmonic_freqs <- seq(1.25, 45, by = 0.416)  # Detailed harmonic frequencies up to 45 Hz

FreqAmps <- FreqAmps %>%
  mutate(FrequencyGroup = case_when(
    Frequency < 5 ~ "Below 5Hz",
    Frequency > 5 & Frequency < 10 & !(Frequency %in% harmonics) ~ "Between 5 and 10 Hz",
    Frequency > 10 & Frequency < 15 & !(Frequency %in% harmonics) ~ "Between 10 and 15 Hz",
    Frequency > 15 & Frequency < 20 & !(Frequency %in% harmonics) ~ "Between 15 and 20 Hz",
    Frequency > 20 & Frequency < 25 & !(Frequency %in% harmonics) ~ "Between 20 and 25 Hz",
    Frequency > 25 & Frequency < 30 & !(Frequency %in% harmonics) ~ "Between 25 and 30 Hz",
    Frequency > 30 & Frequency < 35 & !(Frequency %in% harmonics) ~ "Between 30 and 35 Hz",
    Frequency > 35 & Frequency < 40 & !(Frequency %in% harmonics) ~ "Between 35 and 40 Hz",
    Frequency > 40 & Frequency < 45 & !(Frequency %in% harmonics) ~ "Between 40 and 45 Hz",
    Frequency %in% harmonics ~ "5Hz and Harmonics",
    TRUE ~ "Above 45Hz"  # For frequencies beyond the defined ranges
  ))

# Match with a larger tolerance to account for floating-point precision differences
FreqAmps <- FreqAmps %>%
  mutate(IsHarmonic = sapply(Frequency, function(x) any(abs(x - harmonic_freqs) < 1e-1)))  # Increased tolerance

# Calculate the total amplitude for all harmonic frequencies grouped by DutyCycle
total_harmonic_amplitude <- FreqAmps %>%
  filter(IsHarmonic) %>%
  group_by(DutyCycle) %>%
  reframe(TotalHarmonicAmplitude = sum(Amplitude, na.rm = TRUE))

# Summarize mean amplitudes by FrequencyGroup and DutyCycle
grouped_FreqAmps <- FreqAmps %>%
  filter(IsHarmonic) %>%
  group_by(DutyCycle, FrequencyGroup) %>%
  reframe(
    TotalAmplitude = sum(Amplitude, na.rm = TRUE),  # Sum of amplitudes for the group
    Proportion = TotalAmplitude / total_harmonic_amplitude$TotalHarmonicAmplitude[match(DutyCycle, total_harmonic_amplitude$DutyCycle)],  # Proportion
  ) %>%
  distinct()  # Ensure only unique rows

# Ensure proper ordering of FrequencyGroup and DutyCycle
grouped_FreqAmps <- grouped_FreqAmps %>%
  mutate(
    FrequencyGroup = factor(
      FrequencyGroup,
      levels = c(
        "Below 5Hz", 
        "Between 5 and 10 Hz", 
        "Between 10 and 15 Hz", 
        "Between 15 and 20 Hz", 
        "Between 20 and 25 Hz", 
        "Between 25 and 30 Hz", 
        "Between 30 and 35 Hz", 
        "Between 35 and 40 Hz", 
        "Between 40 and 45 Hz", 
        "5Hz and Harmonics"
      )
    ),
    DutyCycle = factor(DutyCycle, levels = c("25", "50", "75", "100"))  # Ensure proper order
  )

## EEG : for each Subject and Duty cycle condition
  # Filter data for AllAmps FrequencySet
all_amps_data <- df %>% filter(FrequencySet == "AllAmps")

# Define harmonic frequencies and categorize frequencies
harmonics <- seq(5, max(all_amps_data$Frequency, na.rm = TRUE), by = 5)  # Multiples of 5 Hz
harmonic_freqs <- seq(1.25, 45, by = 0.416)  # Detailed harmonic frequencies up to 45 Hz

all_amps_data <- all_amps_data %>%
  mutate(FrequencyGroup = case_when(
    Frequency < 5 ~ "Below 5Hz",
    Frequency > 5 & Frequency < 10 & !(Frequency %in% harmonics) ~ "Between 5 and 10 Hz",
    Frequency > 10 & Frequency < 15 & !(Frequency %in% harmonics) ~ "Between 10 and 15 Hz",
    Frequency > 15 & Frequency < 20 & !(Frequency %in% harmonics) ~ "Between 15 and 20 Hz",
    Frequency > 20 & Frequency < 25 & !(Frequency %in% harmonics) ~ "Between 20 and 25 Hz",
    Frequency > 25 & Frequency < 30 & !(Frequency %in% harmonics) ~ "Between 25 and 30 Hz",
    Frequency > 30 & Frequency < 35 & !(Frequency %in% harmonics) ~ "Between 30 and 35 Hz",
    Frequency > 35 & Frequency < 40 & !(Frequency %in% harmonics) ~ "Between 35 and 40 Hz",
    Frequency > 40 & Frequency < 45 & !(Frequency %in% harmonics) ~ "Between 40 and 45 Hz",
    Frequency %in% harmonics ~ "5Hz and Harmonics",
    TRUE ~ "Above 45Hz"  # For frequencies beyond the defined ranges
  ))

print(all_amps_data)


# Calculate the total amplitude per subject and duty cycle
subject_harmonic_amplitude <- all_amps_data %>%
  group_by(SubjectID, DutyCycle) %>%
  summarize(TotalHarmonicAmplitude = sum(Amplitude, na.rm = TRUE), .groups = "drop")

# Calculate the proportions for each subject, duty cycle, and frequency group
grouped_subject_amps <- all_amps_data %>%
  group_by(SubjectID, DutyCycle, FrequencyGroup) %>%
  summarize(
    TotalAmplitude = sum(Amplitude, na.rm = TRUE),  # Total amplitude for the group
    .groups = "drop"
  ) %>%
  left_join(subject_harmonic_amplitude, by = c("SubjectID", "DutyCycle")) %>%
  mutate(Proportion = TotalAmplitude / TotalHarmonicAmplitude)

# Ensure proper ordering of DutyCycle and FrequencyGroup
grouped_subject_amps <- grouped_subject_amps %>%
  mutate(
    DutyCycle = factor(DutyCycle, levels = c("25", "50", "75", "100")),
    FrequencyGroup = factor(
      FrequencyGroup,
      levels = c(
        "Below 5Hz", 
        "Between 5 and 10 Hz", 
        "Between 10 and 15 Hz", 
        "Between 15 and 20 Hz", 
        "Between 20 and 25 Hz", 
        "Between 25 and 30 Hz", 
        "Between 30 and 35 Hz", 
        "Between 35 and 40 Hz", 
        "Between 40 and 45 Hz", 
        "5Hz and Harmonics"
      )
    )
  )


## PLOT THE 2 TOGETHER 

# Combine EEG and stimulus data for plotting
# Scatterplot with jittered EEG points and black squares for Stimulus
scatter_plot <- ggplot(combined_data, aes(x = FrequencyGroup, y = Proportion)) +
  
  # EEG data: Jittered points
  geom_point(
    aes(color = DutyCycle),
    position = position_jitter(width = 0.2, height = 0),
    size = 2,
    alpha = 0.5,
    shape = 21,
    stroke = 1.5
  ) +
  
  # Stimulus data: Black squares
  geom_point(
    aes(y = StimulusProportion),
    shape = 22,
    color = "black",
    size = 3
  ) +
  
  # Facet by duty cycle
  facet_wrap(~ DutyCycle, ncol = 2) +
  
  # Customize plot appearance
  labs(
    title = "EEG vs. Stimulus Proportions Across Frequency Groups and Duty Cycles",
    x = "Frequency Band",
    y = "Proportion of Total Amplitude (up to 45Hz)",
    color = "EEG (Duty Cycle)",
    fill = "Stimulus (Duty Cycle)"
  ) +
  
  # Apply MATLAB colors for both fill and color
  scale_fill_manual(values = matlab_colors) +
  scale_color_manual(values = matlab_colors) +
  
  # Aesthetic adjustments
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Display the plot
print(scatter_plot)

# Save the plot
ggsave("Combined_EEG_vs_Stimulus_Proportions.svg",
       plot = scatter_plot,
       width = 16,
       height = 10,
       units = "cm",
       device = "svg")



## STATISTICAL ANALYSIS : model fixed effects : stimulus proportion, duty cycle, + interaction between the 2 

# Filter data for below 5Hz frequency group
below_5hz_data <- combined_data %>%
  filter(FrequencyGroup == "Below 5Hz")

stimulus_table <- below_5hz_data %>%
  select(DutyCycle, StimulusProportion) %>%
  distinct()

eeg_data <- below_5hz_data %>%
  select(-StimulusProportion)

combined_data_fixed <- eeg_data %>%
  left_join(stimulus_table, by = "DutyCycle")

# Aggregate StimulusProportion to one value per DutyCycle
stimulus_summary <- combined_data_fixed %>%
  group_by(DutyCycle) %>%
  summarise(StimulusProportion = unique(StimulusProportion)) # Unique values

# Add the unique StimulusProportion back to the original data
adjusted_data <- combined_data_fixed %>%
  select(-StimulusProportion) %>% # Remove the repeated StimulusProportion
  left_join(stimulus_summary, by = "DutyCycle")

# Fit the linear mixed-effects model
model <- lmer(
  Proportion ~ StimulusProportion + (1 | SubjectID),
  data = adjusted_data
)


residuals <- resid(model)
#test for normality of residuals
shapiro.test(residuals) 

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model) ~ fitted(model))
qqnorm(resid(model))
qqline(resid(model), col = "red")

# Check model summary
summary(model)

# Perform ANOVA using the car package
anova_results <- Anova(model, test = "F")
print(anova_results)

# Compute partial eta-squared directly
eta_squared_results <- eta_squared(anova_results, partial = TRUE)
print(eta_squared_results)


##---------------------------------------------
# Calculate mean amplitude and Z-scores concisely
mean_Amplitude <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')

# Do same until 5Hz only
mean_Amplitude_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')

# Pivot stimCochmod and retain only the last value for each DutyCycle
stimCochmod_long <- stimCochmod %>%
  pivot_longer(cols = starts_with("X"),  # Pivot columns that start with "X"
               names_to = "DutyCycle",    # Create a new column "DutyCycle" from the column names
               values_to = "stimcochZscore") %>%
  mutate(DutyCycle = as.numeric(gsub("X", "", DutyCycle)))  # Convert "DutyCycle" to numeric by removing "X"

# Keep the last stimcochZscore for each DutyCycle
stimCochmod_last <- stimCochmod_long %>%
  group_by(DutyCycle) %>%
  slice_tail(n = 1)  # Keep the last value for each DutyCycle

# Keep the first stimcochZscore for each DutyCycle for 5HZ ANALYSIS
stimCochmod_first <- stimCochmod_long %>%
  group_by(DutyCycle) %>%
  slice_head(n = 1)  # Keep the last value for each DutyCycle

# Merge data frames effectively
merged_data_5Hz <- merge(mean_Amplitude_5Hz, mean_zscores_5Hz, by = c("SubjectID", "DutyCycle"))
final_table_5Hz <- merge(merged_data_5Hz, stimCochmod_first, by = "DutyCycle")

merged_data <- merge(mean_Amplitude, mean_zscores, by = c("SubjectID", "DutyCycle"))
final_table <- merge(merged_data, stimCochmod_last, by = "DutyCycle")

##------------------------------------ Analysis on the data up to 5Hz----------------------------------------------------------------------------------
df_model_5Hz <- final_table_5Hz %>%
  select(SubjectID, DutyCycle, meanZscore, stimcochZscore) %>%
  rename(
    participantID = SubjectID,
    condition = DutyCycle,
    participantZscore = meanZscore,
    stimulusZscore = stimcochZscore
  ) %>%
  mutate(condition = factor(condition, levels = c("25", "50", "75", "100")))

## Measuring beat-prominence in EEG across conditions
# Fit linear mixed-effects models
model_condition_5Hz <- lmer(participantZscore ~ condition + (1|participantID), data = df_model_5Hz)
residuals <- resid(model_condition_5Hz )

#test for normality of residuals
shapiro.test(residuals)


# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_condition_5Hz) ~ fitted(model_condition_5Hz))
qqnorm(resid(model_condition_5Hz))
qqline(resid(model_condition_5Hz), col = "red")

# Summary and diagnostics
summary(model_condition_5Hz)
# Perform ANOVA using the car package
anova_results <- Anova(model_condition_5Hz, test = "F")
print(anova_results)

# Compute partial eta-squared directly
eta_squared_results <- eta_squared(anova_results, partial = TRUE)
print(eta_squared_results)

# Calculate Estimated Marginal Means and contrasts
emm2_5Hz <- emmeans(model_condition_5Hz, specs = "condition")
pairwise_results <- contrast(emm2_5Hz, method = "pairwise", adjust = "bonferroni")
pairwise_d <- effectsize::cohens_d(pairwise_results)

# Combine pairwise results with Cohen's d
pairwise_results_with_d <- cbind(as.data.frame(pairwise_results), d = pairwise_d$Cohens_d)
print(pairwise_results_with_d)

# Extract stimCochmod values for contrast testing from stimCochmod_first
stimCochmod_values_5Hz <- stimCochmod_first$stimcochZscore
names(stimCochmod_values_5Hz) <- stimCochmod_first$DutyCycle

# Perform contrast against dynamically extracted stimulus Z-scores
test_results_5Hz <- test(emm2_5Hz, null = stimCochmod_values_5Hz)

# Compute Cohen's d for comparisons against stimulus Z-scores
test_results_with_d <- test_results_5Hz %>%
  mutate(cohens_d = (emmean - null) / SE)  # Difference divided by SE

print(test_results_with_d)

print(test_results_5Hz)


# Test for a linear contrast
contrast(emm2_5Hz, method = "poly")

# Prepare significance stars
emm2_df_5Hz <- as.data.frame(emm2_5Hz)
emm2_df_5Hz$DutyCycle <- emm2_df_5Hz$condition
emm2_df_5Hz$significance <- ifelse(test_results_5Hz$p.value < 0.001, "***",
                                   ifelse(test_results_5Hz$p.value < 0.01, "**",
                                          ifelse(test_results_5Hz$p.value < 0.05, "*", "")))

# Add significance stars for EMMs against zero using test_zero_5Hz
emm2_df_5Hz$significance_zero <- ifelse(test_zero_5Hz$p.value < 0.001, "***",
                                        ifelse(test_zero_5Hz$p.value < 0.01, "**",
                                               ifelse(test_zero_5Hz$p.value < 0.05, "*", "")))


# Ensure that DutyCycle is of the same type (character) in both data frames
emm2_df_5Hz$DutyCycle <- as.character(emm2_df_5Hz$DutyCycle)
stimCochmod_first$DutyCycle <- as.character(stimCochmod_first$DutyCycle)


# Merge stimCochmod_last into emm2_df for consistent stimcochZscore usage
emm2_df_5Hz <- emm2_df_5Hz %>%
  left_join(stimCochmod_first %>% select(DutyCycle, stimcochZscore), by = "DutyCycle")
emm2_df_5Hz$DutyCycle <- as.factor(emm2_df_5Hz$DutyCycle)


# Set plotting parameters
# Set dodge width for consistent spacing (optional if you want to adjust)
dodge_width <- 0.15

# Colors for DutyCycle
matlab_colors <- c("25" = "#4DBEEE", "50" = "#0072BD", "75" = "#EDB120", "100" = "#D95319")

# Calculate pairwise contrasts
contrast_results_5Hz <- contrast(emm2_5Hz, method = "pairwise", adjust = "bonferroni") %>%
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
model_participantZscore_Lmer_intercept_5Hz <- lmer(participantZscore ~ stimulusZscore + (1|participantID), data = df_model_5Hz)

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_participantZscore_Lmer_intercept_5Hz) ~ fitted(model_participantZscore_Lmer_intercept_5Hz))
qqnorm(resid(model_participantZscore_Lmer_intercept_5Hz))
qqline(resid(model_participantZscore_Lmer_intercept_5Hz), col = "red")

# Summary and diagnostics
summary(model_participantZscore_Lmer_intercept_5Hz)
anova_results_stim <- Anova(model_participantZscore_Lmer_intercept_5Hz, test = "F")
print(anova_results_stim)

# Compute partial eta-squared directly
eta_squared_results_stim <- eta_squared(anova_results_stim, partial = TRUE)
print(eta_squared_results_stim)

# Calculate Marginal and Conditional R²
r2_values_5Hz <- r.squaredGLMM(model_participantZscore_Lmer_intercept_5Hz)
marginal_r2_5Hz <- r2_values_5Hz[1]
conditional_r2_5Hz <- r2_values_5Hz[2]





# -------------------------------   PLOTS ------------------------------------------------------------------------------------------
# Create the first plot with only stim Z-score points (simplified plot)
# Define the jittered position with a fixed seed for reproducibility
set.seed(123)  # Setting a seed for consistent jitter across plots
jitter_dodge_position <- position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width)

p5Hz_1 <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = stimcochZscore, group = DutyCycle, fill = DutyCycle)) +
  geom_point(size = 5, shape = 22, stroke = 2, color = "black", fill = "black") +
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-0.5,-0.25, 0, 0.25, 0.5, 0.75,1,1.25), limits = c(-0.55, 2))  +
  
  # Labels and title
  labs(title = "Magnitude spectrum (up to 5Hz)", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2),
    legend.position = "none",
    panel.grid = element_blank()  # Removes the grid
  )

# Create the complete plot with boxplot, jitter, stimcochZscore points, vertical lines, and ticks
p5Hz_2 <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +  # Now using meanZscore for y-axis
  
  
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
  # 
  # COULD NOT GET THIS TO WORK DESPITE BEST EFFORT = NEED TO ADD MANUALLY A CONTRAST BAR TO THE PLOT ON CANVAS !!! 
  # # Horizontal lines for significant contrasts using contrast1 and contrast2
  # geom_segment(data = significant_contrasts, aes(x = contrast1, xend = contrast2, y = 1.6, yend = 1.6), 
  #              color = "black", size = 0.5) +
  # 
  # # Stars for significant contrasts
  # geom_text(data = significant_contrasts, aes(x = (contrast1 + contrast2) / 2, y = 1.65, label = "***"), 
  #           color = "black", size = 8) +
  # 
  # Add significance stars for EMM and stimcochZscore differences
geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = (emmean + stimcochZscore) / 2, label = significance), 
          color = "black", size = 8, data = emm2_df_5Hz) +
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  ## Y-axis limits, breaks, and zero line
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Labels and title
  labs(title = "Magnitude spectrum", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
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
print(p5Hz_2)

# Plot with JUST EEG
p5Hz_3 <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +  # Now using meanZscore for y-axis
  
  
  # Jittered points
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = jitter_dodge_position, size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  ## Y-axis limits, breaks, and zero line
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Labels and title
  labs(title = "Magnitude spectrum", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
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
print(p5Hz_3)

# Stack the two plots vertically using patchwork
combined_plot_5Hz <- p5Hz_1 / p5Hz_2
print(combined_plot_5Hz)

# Create the new plot `p5Hz_4` using significance against zero
p5Hz_4 <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +  
  # Jittered points
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = jitter_dodge_position, size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  # Points for stimcochZscore (black points)
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  
  # Vertical lines for EMM and stimcochZscore
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, y = emmean, yend = 0), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  
  # Line ticks for EMM and stimcochZscore
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = emmean, yend = emmean), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = 0, yend = 0), 
               color = "black", size = 0.4, data = emm2_df_5Hz) +
  
  # Add significance stars for EMM against zero
  geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = emmean, label = significance_zero), 
            color = "black", size = 8, data = emm2_df_5Hz) +
  
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits, breaks, and zero line
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(title = "Magnitude spectrum", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
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

# Display the new plot
print(p5Hz_4)

# Save the plot
ggsave("FFT5Hz_YOUNG_Stim_Zscores_alone.svg", plot = p5Hz_1, width = 20, height = 15, units = "cm", device = 'svg')

# Save the plot
ggsave("FFT5Hz_YOUNG_EEG_Stim_Zscores.svg", plot = p5Hz_2, width = 20, height = 15, units = "cm", device = 'svg')

# Save the plot
ggsave("FFT5Hz_YOUNG_EEG_Zscores_alone.svg", plot = p5Hz_3, width = 20, height = 15, units = "cm", device = 'svg')

# Save the plot
ggsave("FFT5Hz_YOUNG_Stim_Zscores_2 plots.svg", plot = combined_plot_5Hz, width = 20, height = 30, units = "cm", device = 'svg')

# Save the new plot as "FFT5Hz_OLD_EEG_Stim_Zscores_againstZERO.svg"
ggsave("FFT5Hz_YOUNG_EEG_Stim_Zscores_againstZERO.svg", plot = p5Hz_4, width = 20, height = 15, units = "cm", device = 'svg')


# Add the 'group' column with the value 'YOUNG' for all rows
final_table_5Hz$group <- "YOUNG"

# Define the file path with the updated file name
file_path <- "/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_all/final_table_5Hz_YOUNG.csv"

# Export the updated data frame to a CSV file at the specified location
write.csv(final_table_5Hz, file_path, row.names = FALSE)

## ------------------------------------- Analysis on the data up to 15Hz --------------------------------------------------
# Prepare the model data with consistent variable names
df_model <- final_table %>%
  select(SubjectID, DutyCycle, meanZscore, stimcochZscore) %>%
  rename(
    participantID = SubjectID,
    condition = DutyCycle,
    participantZscore = meanZscore,
    stimulusZscore = stimcochZscore
  ) %>%
  mutate(condition = factor(condition, levels = c("25", "50", "75", "100")))

# Fit linear mixed-effects models
model_participantZscore_Lmer_intercept <- lmer(participantZscore ~ stimulusZscore + (1|participantID), data = df_model)
model_condition<- lmer(participantZscore ~ condition + (1|participantID), data = df_model)

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_participantZscore_Lmer_intercept) ~ fitted(model_participantZscore_Lmer_intercept))
qqnorm(resid(model_participantZscore_Lmer_intercept))
qqline(resid(model_participantZscore_Lmer_intercept), col = "red")

par(mfrow = c(1, 2))
plot(resid(model_condition) ~ fitted(model_condition))
qqnorm(resid(model_condition))
qqline(resid(model_condition), col = "red")

# Summary and diagnostics
summary(model_participantZscore_Lmer_intercept)
Anova(model_participantZscore_Lmer_intercept, test = "F")
summary(model_condition)
Anova(model_condition, test = "F")

# Calculate Marginal and Conditional R²
r2_values <- r.squaredGLMM(model_participantZscore_Lmer_intercept)
marginal_r2 <- r2_values[1]
conditional_r2 <- r2_values[2]

# Calculate Estimated Marginal Means and contrasts
emm2 <- emmeans(model_condition, specs = "condition")
contrast(emm2, method = "pairwise", adjust = "bonferroni")

# Summary and test
summary(emm2, infer = c(TRUE, TRUE))
test(emm2, null = 0)

# Extract stimCochmod values for contrast testing from stimCochmod_last
stimCochmod_values <- stimCochmod_last$stimcochZscore
names(stimCochmod_values) <- stimCochmod_last$DutyCycle

# Perform contrast against dynamically extracted stimulus Z-scores
test_results <- test(emm2, null = stimCochmod_values)

# Prepare significance stars
emm2_df <- as.data.frame(emm2)
emm2_df$DutyCycle <- emm2_df$condition
emm2_df$significance <- ifelse(test_results$p.value < 0.001, "***",
                               ifelse(test_results$p.value < 0.01, "**",
                                      ifelse(test_results$p.value < 0.05, "*", "")))

# Ensure that DutyCycle is of the same type (character) in both data frames
emm2_df$DutyCycle <- as.character(emm2_df$DutyCycle)
stimCochmod_last$DutyCycle <- as.character(stimCochmod_last$DutyCycle)


# Merge stimCochmod_last into emm2_df for consistent stimcochZscore usage
emm2_df <- emm2_df %>%
  left_join(stimCochmod_last %>% select(DutyCycle, stimcochZscore), by = "DutyCycle")
emm2_df$DutyCycle <- as.factor(emm2_df$DutyCycle)

# Set dodge width for consistent spacing (optional if you want to adjust)
dodge_width <- 0.15



# Calculate pairwise contrasts
contrast_results <- contrast(emm2, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE)

# Filter out significant contrasts (p < 0.05)
significant_contrasts <- contrast_results %>%
  filter(p.value < 0.05)

# Convert contrast results to a data frame (if it's not already)
contrast_results_df <- as.data.frame(contrast_results)

# Manually split the 'contrast' column into two new columns: 'contrast1' and 'contrast2'
contrast_split <- sapply(strsplit(contrast_results_df$contrast, " - "), function(x) x)

# Create new columns 'contrast1' and 'contrast2'
contrast_results_df$contrast1 <- as.numeric(gsub("condition", "", contrast_split[1, ]))
contrast_results_df$contrast2 <- as.numeric(gsub("condition", "", contrast_split[2, ]))

# Filter out significant contrasts (p < 0.05)
significant_contrasts <- contrast_results_df %>%
  filter(p.value < 0.05)

# Ensure 'DutyCycle' is a factor in the relevant data frames
final_table$DutyCycle <- as.factor(final_table$DutyCycle)
emm2_df$DutyCycle <- as.factor(emm2_df$DutyCycle)  # Convert 'DutyCycle' in emm2_df to factor
# Reorder the levels of 'DutyCycle' in emm2_df
emm2_df$DutyCycle <- factor(emm2_df$DutyCycle, levels = c("25", "50", "75", "100"))


# -------------------------------   PLOTS ------------------------------------------------------------------------------------------
# Create the first plot with only stim Z-score points (simplified plot)
p1 <- ggplot(final_table, aes(x = DutyCycle, y = stimcochZscore, group = DutyCycle, fill = DutyCycle)) +
  geom_point(size = 5, shape = 22, stroke = 2, color = "black", fill = "black") +
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5), limits = c(-0.3, 0.7)) +
  
  # Labels and title
  labs(title = "Magnitude spectrum", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2),
    legend.position = "none",
    panel.grid = element_blank()  # Removes the grid
  )

# Create the complete plot with boxplot, jitter, stimcochZscore points, vertical lines, and ticks
p2 <- ggplot(final_table, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +  # Now using meanZscore for y-axis
  
  # Violin plot
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), adjust = 0.4, trim = FALSE, alpha = 0.3) +
  
  # Jittered points
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = position_jitterdodge(dodge.width = dodge_width), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot
  geom_boxplot(alpha = 0.3, width = 0.2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  
  # Points for stimcochZscore (black points)
  geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2) +
  
  # Vertical lines for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, y = emmean, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df) +
  
  # Line ticks for EMM and stimcochZscore (Use 'DutyCycle')
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = emmean, yend = emmean), 
               color = "black", size = 0.4, data = emm2_df) +
  geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, y = stimcochZscore, yend = stimcochZscore), 
               color = "black", size = 0.4, data = emm2_df) +
  # 
  # COULD NOT GET THIS TO WORK DESPITE BEST EFFORT = NEED TO ADD MANUALLY A CONTRAST BAR TO THE PLOT ON CANVAS !!! 
  # # Horizontal lines for significant contrasts using contrast1 and contrast2
  # geom_segment(data = significant_contrasts, aes(x = contrast1, xend = contrast2, y = 1.6, yend = 1.6), 
  #              color = "black", size = 0.5) +
  # 
  # # Stars for significant contrasts
  # geom_text(data = significant_contrasts, aes(x = (contrast1 + contrast2) / 2, y = 1.65, label = "***"), 
  #           color = "black", size = 8) +
  # 
  # Add significance stars for EMM and stimcochZscore differences
  geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = (emmean + stimcochZscore) / 2, label = significance), 
            color = "black", size = 8, data = emm2_df) +
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
    # Y-axis limits and breaks
    scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5), limits = c(-0.3, 0.7)) +
    
     # Labels and title
     labs(title = "Magnitude spectrum", y = "EEG beat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
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
print(p2)

# Stack the two plots vertically using patchwork
combined_plot <- p1 / p2
print(combined_plot)

# Save the plot
ggsave("FFT15Hz_Stim_Zscores_alone.svg", plot = p1, width = 20, height = 15, units = "cm", device = 'svg')

# Save the plot
ggsave("FFT15Hz_EEG_Stim_Zscores.svg", plot = p2, width = 20, height = 15, units = "cm", device = 'svg')

# Save the plot
ggsave("FFT15Hz_Stim_Zscores_2 plots.svg", plot = combined_plot, width = 20, height = 30, units = "cm", device = 'svg')