

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
#pairwise_d <- effectsize::cohens_d(pairwise_results)

# Combine pairwise results with Cohen's d
#pairwise_results_with_d <- cbind(as.data.frame(pairwise_results), d = pairwise_d$Cohens_d)
#print(pairwise_results_with_d)

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


# Ensure that DutyCycle is of the same type (character) in both data frames
emm2_df_5Hz$DutyCycle <- as.character(emm2_df_5Hz$DutyCycle)
stimCochmod_first$DutyCycle <- as.character(stimCochmod_first$DutyCycle)


# Merge stimCochmod_last into emm2_df for consistent stimcochZscore usage
emm2_df_5Hz <- emm2_df_5Hz %>%
  left_join(stimCochmod_first %>% select(DutyCycle, stimcochZscore), by = "DutyCycle")
emm2_df_5Hz$DutyCycle <- as.factor(emm2_df_5Hz$DutyCycle)




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
geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = (emmean + stimcochZscore) / 2, label = significance), 
          color = "black", size = 13, data = emm2_df_5Hz) +
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  ## Y-axis limits, breaks, and zero line
  #scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.25)) +
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.55, 1.45)) +
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
contrast_height <- max_y_fft + 0.1
height_step <- 0.05

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
    geom_text(aes(x = (c1_pos + c2_pos) / 2, y = contrast_height + 0.02, label = label),
              inherit.aes = FALSE, size = 15, color = "black")
  
  contrast_height <- contrast_height + height_step
}
  

# Display the complete plot
print(plotFFT)







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
# Prepare df_model for clarity
df_model_acf <- final_table_acf %>%
  select(SubjectID, DutyCycle, z_meter_rel, stimcochZscore) %>%
  rename(participantID = SubjectID, condition = DutyCycle,
         participantZscore = z_meter_rel, stimulusZscore = stimcochZscore) %>%
  mutate(condition = factor(condition, levels = c("25", "50", "75", "100")))


## Measuring beat-prominence in EEG across conditions
# Fit linear mixed-effects models
model_condition_acf <- lmer(participantZscore ~ condition + (1|participantID), data = df_model_acf)
residuals_acf <- resid(model_condition_acf)

#test for normality of residuals
shapiro.test(residuals_acf)

#Check assumptions visually for model
par(mfrow=c(1,2))
plot(resid(model_condition_acf ) ~ fitted(model_condition_acf ))
qqnorm(resid(model_condition_acf )); qqline(resid(model_condition_acf ), col = "red")

# Summary and diagnostics
summary(model_condition_acf )
anova_results_acf <- Anova(model_condition_acf , test="F")
print(anova_results_acf)

# Compute partial eta-squared directly
eta_squared_results_acf <- eta_squared(anova_results_acf, partial = TRUE)
print(eta_squared_results_acf)





# Posthoc testing
# Calculate Estimated Marginal Means and contrasts
emm2_acf <- emmeans(model_condition_acf, specs = "condition")
pairwise_results_acf <- contrast(emm2_acf, method = "pairwise", adjust = "bonferroni")
#pairwise_d <- effectsize::cohens_d(pairwise_results)

# Combine pairwise results with Cohen's d
#pairwise_results_with_d <- cbind(as.data.frame(pairwise_results), d = pairwise_d$Cohens_d)
#print(pairwise_results_with_d)

# Calculate pairwise contrasts and convert to data frame
contrast_results_acf <- contrast(emm2_acf, method = "pairwise", adjust = "bonferroni") %>%
  as.data.frame()  # Convert contrast results to a data frame


# Perform contrast against stimulus Z-scores
stimCochmod <- c(-0.3295543, -0.3607533, -0.3717696, -0.3984415)
names(stimCochmod) <- c("25", "50", "75", "100")
test_results_acf <- test(emm2_acf, null = stimCochmod)

# Compute Cohen's d for comparisons against stimulus Z-scores
test_results_with_d_acf <- test_results_acf %>%
  mutate(cohens_d = (emmean - null) / SE)  # Difference divided by SE

print(test_results_with_d_acf)

# Perform contrast against zero
summary(emm2_acf, infer = c(TRUE, TRUE))
test_zero_acf <- test(emm2_acf, null = 0)

# Add the null value explicitly to the output
test_zero_acf <- test_zero_acf %>%
  mutate(null = 0)  # Add a column with the null value

# Compute Cohen's d for comparisons against stimulus Z-scores
test_zero_with_d_acf <- test_zero_acf  %>%
  mutate(cohens_d = (emmean - null) / SE)  # Difference divided by SE

print(test_zero_with_d_acf )


# Test for a linear contrast
contrast(emm2_acf, method = "poly")



# Extract significance levels based on p-values
emm2_df_acf <- as.data.frame(emm2_acf)
emm2_df_acf$DutyCycle <- emm2_df_acf$condition  # Ensure DutyCycle is in emm2_df
emm2_df_acf$significance <- ifelse(test_results_acf$p.value < 0.001, "***",
                               ifelse(test_results_acf$p.value < 0.01, "**",
                                      ifelse(test_results_acf$p.value < 0.05, "*", "")))


# Add significance stars for EMMs against zero using results_zero
emm2_df_acf$significance_zero <- ifelse(test_zero_acf$p.value < 0.001, "***",
                                    ifelse(test_zero_acf$p.value < 0.01, "**",
                                           ifelse(test_zero_acf$p.value < 0.05, "*", "")))


# Merge stimcochZscore from final_table into emm2_df based on matching DutyCycle/condition
stimcochZscore_map <- final_table_acf %>% select(DutyCycle, stimcochZscore) %>% distinct()
emm2_df_acf <- merge(emm2_df_acf, stimcochZscore_map, by = "DutyCycle")


# Calculate pairwise contrasts
contrast_results_acf <- contrast(emm2_acf, method = "pairwise", adjust = "bonferroni") %>%
  summary(infer = TRUE)

# Filter out significant contrasts (p < 0.05)
significant_contrasts_acf <- contrast_results_acf %>%
  filter(p.value < 0.05)

# Convert contrast results to a data frame (if it's not already)
contrast_results_df_acf <- as.data.frame(contrast_results_acf)

# Make sure the 'contrast' column is a character
contrast_results_df_acf$contrast <- as.character(contrast_results_df_acf$contrast)

# Manually split the 'contrast' column into two new columns: 'contrast1' and 'contrast2'
# Using strsplit() and sapply to extract the two parts
contrast_split <- sapply(strsplit(contrast_results_df_acf$contrast, " - "), function(x) x)

# Create new columns 'contrast1' and 'contrast2'
contrast_results_df_acf$contrast1 <- contrast_split[1, ]
contrast_results_df_acf$contrast2 <- contrast_split[2, ]

# Now we need to remove the 'condition' prefix to just get the numeric values
contrast_results_df_acf <- contrast_results_df_acf %>%
  mutate(contrast1 = gsub("condition", "", contrast1),
         contrast2 = gsub("condition", "", contrast2))

# Filter out significant contrasts (p < 0.05)
significant_contrasts_acf <- contrast_results_df_acf %>%
  filter(p.value < 0.05)

# Calculate pairwise contrasts and convert to a data frame
contrast_results_acf <- as.data.frame(contrast(emm2_acf, method = "pairwise", adjust = "bonferroni"))

# Manually split the 'contrast' column into 'contrast1' and 'contrast2'
contrast_split <- sapply(strsplit(contrast_results_acf$contrast, " - "), function(x) x)

# Create new columns 'contrast1' and 'contrast2'
contrast_results_acf$contrast1 <- as.numeric(gsub("condition", "", contrast_split[1, ]))
contrast_results_acf$contrast2 <- as.numeric(gsub("condition", "", contrast_split[2, ]))

# Filter out significant contrasts (p < 0.05)
significant_contrasts_acf <- contrast_results_acf %>%
  filter(p.value < 0.05)

# Convert contrast1 and contrast2 to numeric if not already
significant_contrasts_acf$contrast1 <- as.numeric(significant_contrasts_acf$contrast1)
significant_contrasts_acf$contrast2 <- as.numeric(significant_contrasts_acf$contrast2)

#---------PLOTS -----------------------------------------
# Create the first plot with only stim Z-score points (simplified plot)


# Create the complete plot
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
geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = (emmean + stimcochZscore) / 2, label = significance), 
          color = "black", size = 13, data = emm2_df_acf) +
  
  # Custom color scale
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-1.55, 1.85)) +
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
contrast_height_acf <- max_y_acf + 0.1
height_step_acf <- 0.05

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
    geom_text(aes(x = (c1_pos + c2_pos) / 2, y = contrast_height_acf + 0.02, label = label),
              inherit.aes = FALSE, size = 15, color = "black")
  
  contrast_height_acf <- contrast_height_acf + height_step_acf
}
# Display the plot
print(plotACF)


## -----------------------------------------------Overall plot----------------------------------------------
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

combined_plot <- plot_grid(
  plotFFT, 
  plotACF, 
  labels = c("A.", "B."), 
  label_size = 30, 
  ncol = 2, 
  align = "h", 
  rel_widths = c(1, 1.1)  # Give ACF slightly more width to prevent overlap
)

# Display the combined plot
print(combined_plot)


# Save the combined plot
ggsave("FFT_ACF_EEG_combined.svg", plot = combined_plot, width = 50, height = 20, units = "cm", device = "svg")
