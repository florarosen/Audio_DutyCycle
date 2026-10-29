library(ggplot2)
library(dplyr)
library(lme4)        # For linear mixed-effects models
library(car)         # For ANOVA
library(emmeans)     # For estimated marginal means (EMMs) and post-hoc analysis
library(ggplot2)     # For plotting
library(nortest)     # For normality tests
library(dplyr)
library(tidyr)
library(extrafont)
library(knitr)
library(PMCMRplus)
library(nlme)
library(cowplot)
library(readr)
library(forcats)
library(lmerTest)
library(effectsize)
library(MuMIn)
library(cowplot)
library(grid)
library(ggplot2)
library(rsvg)      # For handling SVG files
library(gridExtra) # For rasterizing the image

source("R_rainclouds.R")
source("summarySE.R")

matlab_colors <- c("25" = "#4DBEEE", "50" = "#0072BD", "75" = "#EDB120", "100" = "#D95319")

# Set working directory to the location of the CSV files
setwd("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis")

# Load the ITI data into a dataframe
itiData <- read.csv("itiDataLongFormat.csv")


df <- read.csv("AggregatedTAPDataForR_Individual.csv")  # All Individual Data 
stimCochmod <- read.csv("StimZscoresCochModelAN.csv") # Read stimulus data for cochlear model


# Load acf data into dataframe
acf_data <-read.csv("acf_tapping_indtrials.csv")
#Hard coded acf cochlear model values !!
stimCochACF <- data.frame(
  `25` = -0.3295543,
  `50` = -0.3607533,
  `75` = -0.3717696,
  `100` = -0.3984415
)


## ---------------------------------------------- FFT ANALYSIS -------------------------------------------------------------------------------
# Fix factors concisely
df <- df %>%
  mutate(
    SubjectID = factor(SubjectID),
    DutyCycle = factor(DutyCycle, levels = c("25", "50", "75", "100")),
    Meter = factor(Meter, levels = c("Rel", "UnRel")),
    FrequencySet = factor(FrequencySet, levels = c("Set9", "Set13","AllAmps"))
  )

# Filter data as per requirement
data_sets <- list(
  data_upto5Hz = df %>% filter(FrequencySet == "Set9"),
  data_upto5Hz_MR = df %>% filter(FrequencySet == "Set9", Meter == "Rel"),
  data_upto5Hz_MUR = df %>% filter(FrequencySet == "Set9", Meter == "UnRel"),
  data_upto15Hz = df %>% filter(FrequencySet == "Set13"),
  data_upto15Hz_MR = df %>% filter(FrequencySet == "Set13", Meter == "Rel"),
  data_upto15Hz_MUR = df %>% filter(FrequencySet == "Set13", Meter == "UnRel")
)

# Calculate mean amplitude and z-scores concisely for Set9 (up to 5Hz)
mean_Amplitude_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores_5Hz <- data_sets$data_upto5Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')

# Calculate mean amplitude and z-scores concisely for Set13 (up to 15Hz)
mean_Amplitude_15Hz <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanAmplitude = mean(Amplitude, na.rm = TRUE), .groups = 'drop')

mean_zscores_15Hz <- data_sets$data_upto15Hz_MR %>%
  group_by(SubjectID, DutyCycle) %>%
  summarise(meanZscore = mean(ZScore, na.rm = TRUE), .groups = 'drop')

# Get stimulus values for Set9 and Set13
stim_values_5Hz <- as.numeric(stimCochmod[1, ])
stim_values_15Hz <- as.numeric(stimCochmod[3, ])


# Merge data frames effectively for Set9
merged_data_5Hz <- merge(mean_Amplitude_5Hz, mean_zscores_5Hz, by = c("SubjectID", "DutyCycle"))

final_table_5Hz <- merged_data_5Hz %>%
  mutate(stimcochZscore = case_when(
    DutyCycle == "25" ~ stim_values_5Hz[1],
    DutyCycle == "50" ~ stim_values_5Hz[2],
    DutyCycle == "75" ~ stim_values_5Hz[3],
    DutyCycle == "100" ~ stim_values_5Hz[4]
  ))%>%
  arrange(SubjectID)


# Merge data frames effectively for Set13
merged_data_15Hz <- merge(mean_Amplitude_15Hz, mean_zscores_15Hz, by = c("SubjectID", "DutyCycle"))

final_table_15Hz <- merged_data_15Hz %>%
  mutate(stimcochZscore = case_when(
    DutyCycle == "25" ~ stim_values_15Hz[1],
    DutyCycle == "50" ~ stim_values_15Hz[2],
    DutyCycle == "75" ~ stim_values_15Hz[3],
    DutyCycle == "100" ~ stim_values_15Hz[4]
  ))%>%
  arrange(SubjectID)


##------------------------------------ Analysis on the data up to 5Hz----------------------------------------------------------------------------------


## Measuring beat-prominence in tapping across conditions
# Fit linear mixed-effects models
model_condition_5Hz <- lmer(meanZscore ~ DutyCycle + (1|SubjectID), data =  final_table_5Hz)
residuals <- resid(model_condition_5Hz )

#test for normality of residuals
shapiro.test(residuals)


# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_condition_5Hz) ~ fitted(model_condition_5Hz))
qqnorm(resid(model_condition_5Hz))
qqline(resid(model_condition_5Hz), col = "red")

# not normal !!!!!! do Friedman test instead !!!
#If normality of residuals is not respected : 
friedman.test(meanZscore ~ DutyCycle | SubjectID, data = final_table_5Hz)


## Measuring stimulus beat-prominence impact on EEG beat-prominence in EEG across conditions
# Fit linear mixed-effects models
model_stim_5Hz <- lmer(meanZscore ~ stimcochZscore + (1|SubjectID), data =  final_table_5Hz)
residuals_modelstim <- resid(model_stim_5Hz )

#test for normality of residuals
shapiro.test(residuals_modelstim )

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_stim_5Hz) ~ fitted(model_stim_5Hz))
qqnorm(resid(model_stim_5Hz))
qqline(resid(model_stim_5Hz), col = "red")

# not normal !!!!!! do Friedman test instead !!!
#If normality of residuals is not respected : 
friedman_result<-friedman.test(meanZscore ~ stimcochZscore | SubjectID, data = final_table_5Hz)

# Extract Kendall’s W (W = chi-squared / (n * (k - 1)))
n <- length(unique(final_table_5Hz$SubjectID))  # Number of subjects
k <- length(unique(final_table_5Hz$DutyCycle))  # Number of conditions
kendalls_W <- friedman_result$statistic / (n * (k - 1))

# Print results
print(friedman_result)
cat("Kendall’s W:", kendalls_W, "\n")

# Perform Wilcoxon signed-rank test for each DutyCycle level
wilcoxon_results <- final_table_5Hz %>%
  group_by(DutyCycle) %>%
  summarise(
    test_result = list(wilcox.test(meanZscore - stimcochZscore, mu = 0))
  )

# Print the results
wilcoxon_results$test_result

# Plotting 
# Convert Wilcoxon results into a tidy data frame
wilcoxon_results_df <- wilcoxon_results %>%
  mutate(
    W = sapply(test_result, function(x) x$statistic),
    p_value = sapply(test_result, function(x) x$p.value)
  ) %>%
  select(DutyCycle, W, p_value)  # Keep relevant columns

# Define significance thresholds for stars
wilcoxon_results_df$significance <- ifelse(wilcoxon_results_df$p_value < 0.001, "***",
                                               ifelse(wilcoxon_results_df$p_value < 0.01, "**",
                                                      ifelse(wilcoxon_results_df$p_value < 0.05, "*", "")))

# Print the cleaned results table
print(wilcoxon_results_df)

# Calculate median z_meter_rel for each DutyCycle to position the lines
medians <- final_table_5Hz %>%
  group_by(DutyCycle) %>%
  summarize(median_z = median(meanZscore))

# Merge median values into final_table_5Hz 
final_table_5Hz <- merge(final_table_5Hz , medians, by = "DutyCycle")

# Merge significance data into final_table_5Hz 
final_table_5Hz <- merge(final_table_5Hz , wilcoxon_results_df[, c("DutyCycle", "significance")], by = "DutyCycle", all.x = TRUE)


# Plotting 
# Set a dodge width for consistent spacing of points and boxplots
dodge_width <- 0.15

# Ensure EMM plot looks well-spaced with minimal nudge and consistent spacing
plotFFT <- ggplot(final_table_5Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +
 
  # Jittered points for meanZscore with slight spacing using position_jitterdodge
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = position_jitterdodge(dodge.width = dodge_width), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  # Boxplot with no fill, using position_dodge for separation
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  # # Points for stimcochZscore with black color fill and darker outline
   geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2, show.legend = TRUE) +
  # 
  # # Vertical lines positioned with subtle left offset for minimal separation
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, 
  #                  y = median_z, yend = stimcochZscore), color = "black", size = 0.4) +
  # # Line ticks at the EMM and stimcochZscore
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = median_z, yend = median_z), color = "black", size = 0.4) +  # Tick at EMM
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = stimcochZscore, yend = stimcochZscore), color = "black", size = 0.4) +  # Tick at stimcochZscore
  geom_hline(yintercept = 0.0, linetype = "dashed", colour = "black", size = 0.5) +

  
  # Custom color scale for different DutyCycle conditions
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.6, 1.4)) +
  
  # Update the labels and title
  labs(title = "Magnitude spectrum", y = "Tapping \nbeat-related z-score", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  # Custom theme settings
  # Custom theme settings (matching Magnitude Spectrum plot)
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
    axis.line = element_line(color = "black", size = 1.2)
  ) 
# +
  # guides(fill = "none", color = "none", 
  #        shape = guide_legend(override.aes = list(size = 5, shape = 22, fill = "black", color = "black"),
  #                             title = "Stim Z-score")) +
  # 
  # # Add stars for significance levels with small left offset to avoid overlap
  # geom_text(aes(x = as.numeric(DutyCycle) - 0.35, y = (median_z + stimcochZscore) / 2, label = significance), 
  #           color = "black", size = 13, fontface = "plain", family = "sans")

# Display the plot
print(plotFFT)

# Save the plot as requested
ggsave("FFT_tapping.svg", plot = plotFFT, width = 20, height = 15, units = "cm", device = "svg")

##------------------------------------ Analysis on the data up to 15Hz----------------------------------------------------------------------------------


## Measuring beat-prominence in tapping across conditions
# Fit linear mixed-effects models
model_condition_15Hz <- lmer(meanZscore ~ DutyCycle + (1|SubjectID), data =  final_table_15Hz)
residuals <- resid(model_condition_15Hz )

#test for normality of residuals
shapiro.test(residuals)


# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_condition_15Hz) ~ fitted(model_condition_15Hz))
qqnorm(resid(model_condition_15Hz))
qqline(resid(model_condition_15Hz), col = "red")

# not normal !!!!!! do Friedman test instead !!!
#If normality of residuals is not respected : 
friedman.test(meanZscore ~ DutyCycle | SubjectID, data = final_table_15Hz)


## Measuring stimulus beat-prominence impact on EEG beat-prominence in EEG across conditions
# Fit linear mixed-effects models
model_stim_15Hz <- lmer(meanZscore ~ stimcochZscore + (1|SubjectID), data =  final_table_15Hz)
residuals_modelstim <- resid(model_stim_15Hz )

#test for normality of residuals
shapiro.test(residuals_modelstim )

# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_stim_15Hz) ~ fitted(model_stim_15Hz))
qqnorm(resid(model_stim_15Hz))
qqline(resid(model_stim_15Hz), col = "red")

# not normal !!!!!! do Friedman test instead !!!
#If normality of residuals is not respected : 
friedman_result<-friedman.test(meanZscore ~ stimcochZscore | SubjectID, data = final_table_15Hz)

# Extract Kendall’s W (W = chi-squared / (n * (k - 1)))
n <- length(unique(final_table_5Hz$SubjectID))  # Number of subjects
k <- length(unique(final_table_5Hz$DutyCycle))  # Number of conditions
kendalls_W <- friedman_result$statistic / (n * (k - 1))

# Print results
print(friedman_result)
cat("Kendall’s W:", kendalls_W, "\n")

# Perform Wilcoxon signed-rank test for each DutyCycle level
wilcoxon_results <- final_table_15Hz %>%
  group_by(DutyCycle) %>%
  summarise(
    test_result = list(wilcox.test(meanZscore - stimcochZscore, mu = 0))
  )

# Print the results
wilcoxon_results$test_result

# Plotting 
# Convert Wilcoxon results into a tidy data frame
wilcoxon_results_df <- wilcoxon_results %>%
  mutate(
    W = sapply(test_result, function(x) x$statistic),
    p_value = sapply(test_result, function(x) x$p.value)
  ) %>%
  select(DutyCycle, W, p_value)  # Keep relevant columns

# Define significance thresholds for stars
wilcoxon_results_df$significance <- ifelse(wilcoxon_results_df$p_value < 0.001, "***",
                                           ifelse(wilcoxon_results_df$p_value < 0.01, "**",
                                                  ifelse(wilcoxon_results_df$p_value < 0.05, "*", "")))

# Print the cleaned results table
print(wilcoxon_results_df)

# Calculate median z_meter_rel for each DutyCycle to position the lines
medians <- final_table_15Hz %>%
  group_by(DutyCycle) %>%
  summarize(median_z = median(meanZscore))

# Merge median values into final_table_5Hz 
final_table_15Hz <- merge(final_table_15Hz , medians, by = "DutyCycle")

# Merge significance data into final_table_5Hz 
final_table_15Hz <- merge(final_table_15Hz , wilcoxon_results_df[, c("DutyCycle", "significance")], by = "DutyCycle", all.x = TRUE)


# Plotting 
# Set a dodge width for consistent spacing of points and boxplots
dodge_width <- 0.15

# Ensure EMM plot looks well-spaced with minimal nudge and consistent spacing
plotFFT15Hz <- ggplot(final_table_15Hz, aes(x = DutyCycle, y = meanZscore, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +
  
  # Jittered points for meanZscore with slight spacing using position_jitterdodge
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = position_jitterdodge(dodge.width = dodge_width), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  # Boxplot with no fill, using position_dodge for separation
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  # # Points for stimcochZscore with black color fill and darker outline
   geom_point(aes(y = stimcochZscore), color = "black", fill = "black", size = 5, shape = 22, stroke = 2, show.legend = TRUE) +
  # 
  # # Vertical lines positioned with subtle left offset for minimal separation
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, 
  #                  y = median_z, yend = stimcochZscore), color = "black", size = 0.4) +
  # # Line ticks at the EMM and stimcochZscore
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = median_z, yend = median_z), color = "black", size = 0.4) +  # Tick at EMM
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = stimcochZscore, yend = stimcochZscore), color = "black", size = 0.4) +  # Tick at stimcochZscore
  geom_hline(yintercept = 0.0, linetype = "dashed", colour = "black", size = 0.5) +
  
  
  # Custom color scale for different DutyCycle conditions
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks
  scale_y_continuous(breaks = c(-0.5, 0, 1.2), limits = c(-0.6, 1.4)) +
  
  # Update the labels and title
  labs(title = "Magnitude spectrum", y = "Tapping \nbeat-related z-score (a.u)", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  # Custom theme settings
  # Custom theme settings (matching Magnitude Spectrum plot)
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
    axis.line = element_line(color = "black", size = 1.2)
  ) 

# Display the plot
print(plotFFT15Hz)

# Save the plot as requested
ggsave("FFT_tapping15Hz.svg", plot = plotFFT15Hz, width = 20, height = 15, units = "cm", device = "svg")

## ---------ACF data analysis -------------------------------------------------------
# Fix factors concisely
# Apply factor conversion
acf_data$SubjectID <- factor(acf_data$SubjectID)
acf_data$DutyCycle <- factor(acf_data$DutyCycle, levels = c("25", "50", "75", "100"))

# Apply factor conversion
# Add the long format stimCochACF to your data frame
stimCochACF_long <- rep(stimCochACF, each = nrow(acf_data) / length(unique(acf_data$DutyCycle)))
acf_data <- acf_data %>%
  mutate(stimCochACF = stimCochACF_long)
# Ensure stimCochACF is numeric
acf_data$stimCochACF <- as.numeric(acf_data$stimCochACF)

## Measuring beat-prominence in tapping across conditions
# Fit linear mixed-effects models
model_acf <- lmer(z_meter_rel ~ DutyCycle + (1|SubjectID), data =  acf_data)
residuals_acf <- resid(model_acf )

#test for normality of residuals
shapiro.test(residuals_acf)


# Check assumptions visually for model with continuous stim Z-score
par(mfrow = c(1, 2))
plot(resid(model_acf ) ~ fitted(model_acf ))
qqnorm(resid(model_acf ))
qqline(resid(model_acf ), col = "red")

# not normal !!!!!! do Friedman test instead !!!
#If normality of residuals is not respected : 
friedman_result_acf<-friedman.test(z_meter_rel ~ DutyCycle  | SubjectID, data = acf_data)

# Extract Kendall’s W (W = chi-squared / (n * (k - 1)))
n <- length(unique(acf_data$SubjectID))  # Number of subjects
k <- length(unique(acf_data$DutyCycle))  # Number of conditions
kendalls_W <- friedman_result_acf$statistic / (n * (k - 1))

# Print results
print(friedman_result_acf)
cat("Kendall’s W:", kendalls_W, "\n")

# Perform Wilcoxon signed-rank test for each DutyCycle level against stimulus
wilcoxon_results_acf <- acf_data%>%
  group_by(DutyCycle) %>%
  summarise(
    test_result_acf = list(wilcox.test(z_meter_rel - stimCochACF, mu = 0))
  )

# Print the results
wilcoxon_results_acf$test_result_acf

# Perform Wilcoxon signed-rank test for each DutyCycle level against ZERO
wilcoxon_results_acf_zero <- acf_data%>%
  group_by(DutyCycle) %>%
  summarise(
    test_result_acf_zero = list(wilcox.test(z_meter_rel, mu = 0))
  )

# Print the results
wilcoxon_results_acf_zero$test_result_acf_zero



# Plotting 
# Convert Wilcoxon results into a tidy data frame
wilcoxon_results_acf_df <- wilcoxon_results_acf %>%
  mutate(
    W = sapply(test_result_acf, function(x) x$statistic),
    p_value = sapply(test_result_acf, function(x) x$p.value)
  ) %>%
  select(DutyCycle, W, p_value)  # Keep relevant columns

# Define significance thresholds for stars
wilcoxon_results_acf_df$significance <- ifelse(wilcoxon_results_acf_df$p_value < 0.001, "***",
                                                ifelse(wilcoxon_results_acf_df$p_value < 0.01, "**",
                                                       ifelse(wilcoxon_results_acf_df$p_value < 0.05, "*", "")))

# Print the cleaned results table
print(wilcoxon_results_acf_df)

# Calculate median z_meter_rel for each DutyCycle to position the lines
medians <- acf_data %>%
  group_by(DutyCycle) %>%
  summarize(median_z = median(z_meter_rel))

# Merge median values into acf_data
acf_data <- merge(acf_data, medians, by = "DutyCycle")

# Merge significance data from wilcoxon_results_acf_df into acf_data
acf_data <- merge(acf_data, wilcoxon_results_acf_df[, c("DutyCycle", "significance")], by = "DutyCycle", all.x = TRUE)


# Set a dodge width for consistent spacing of points and boxplots
dodge_width <- 0.15

# Updated plotACF to match the layout of the Magnitude Spectrum plot
plotACF <- ggplot(acf_data, aes(x = DutyCycle, y = z_meter_rel, group = DutyCycle, colour = DutyCycle, fill = DutyCycle)) +
  
  # Jittered points for meanZscore with consistent spacing (matching Magnitude Spectrum plot)
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = position_jitterdodge(dodge.width = dodge_width), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  # Boxplot with consistent position and appearance (matching Magnitude Spectrum plot)
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = dodge_width)) +
  # # Points for stimCochACF with black color fill and darker outline (matching Magnitude Spectrum plot)
   geom_point(aes(y = stimCochACF), color = "black", fill = "black", size = 5, shape = 22, stroke = 2, show.legend = TRUE) +
  # 
  # # Vertical lines positioned consistently with Magnitude Spectrum plot
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.2, 
  #                  y = median_z, yend = stimCochACF), color = "black", size = 0.4) +
  # # # Line ticks at the median and stimCochACF (matching Magnitude Spectrum plot)
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = median_z, yend = median_z), color = "black", size = 0.4) +  # Tick at median
  # geom_segment(aes(x = as.numeric(DutyCycle) - 0.2, xend = as.numeric(DutyCycle) - 0.15, 
  #                  y = stimCochACF, yend = stimCochACF), color = "black", size = 0.4) +  # Tick at stimCochACF
  geom_hline(yintercept = 0.0, linetype = "dashed", colour = "black", size = 0.5) +


  # Custom color scale matching Magnitude Spectrum plot
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  
  # Y-axis limits and breaks (consistent with the original Autocorrelation plot)
  scale_y_continuous(breaks = c(-1.5, 0, 1.5), limits = c(-1.6, 1.9)) +
  
  # Add space between x-axis categories (matching Magnitude Spectrum plot)
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  
  # Keep the original labels and title
  labs(title = "Autocorrelation", y = "Tapping \nbeat-related z-score", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  # Custom theme settings (matching Magnitude Spectrum plot)
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    axis.title.y = element_text(hjust = 0.5, vjust = 0.5), # Centers vertically and horizontally
    legend.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)
  ) 
# +
  #guides(fill = "none", color = "none", 
  #shape = guide_legend(override.aes = list(size = 5, shape = 22, fill = "black", color = "black"),
  #                          title = "Stim Z-score")) +
  
  # Add stars for significance levels with matching position and size
  #geom_text(data = acf_data, aes(x = as.numeric(DutyCycle) - 0.35, y = (median_z + stimCochACF) / 2, label = significance), 
  #         color = "black", size = 13, fontface = "plain", family = "sans")

# Display the plot
print(plotACF)

# Save the plot
ggsave("ACF_tapping.svg", plot = plotACF, width = 20, height = 15, units = "cm", device = "svg")

## Overall plot
# Create a third plot that displays both plots side by side
combined_plot <- plot_grid(plotFFT, plotACF, labels = c("B.", "C."), label_size = 24, ncol = 2, align = "h")

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("FFT_ACF_tapping_combined.svg", plot = combined_plot, width = 40, height = 20, units = "cm", device = "svg")



#------------------------------- ITI ANALYSIS -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#Factor all variables 
itiData$Condition <- factor(itiData$Condition)
itiData$Subject <- factor(itiData$Subject)
itiData$Trial <- factor(itiData$Trial)

head(itiData)

# Calculate the median ITIValue within each trial for each participant and condition
itiData <- itiData %>%
  group_by(Subject, Condition, Trial) %>%
  mutate(medITI = median(ITIValue)) %>%
  ungroup()

# Calculate the median ITIValue within each trial for each participant and condition
itiData <- itiData %>%
  group_by(Subject, Condition, Trial) %>%
  mutate(medError = median(ITIError)) %>%
  ungroup()



# Group by 'Condition' and calculate the IQR, range, and MAD
dispersion_by_condition <- itiData %>%
  group_by(Condition) %>%
  summarise(
    IQR = IQR(ITIValue),
    Range = max(ITIValue) - min(ITIValue),
    MAD = mad(ITIValue)
  )

# View the dispersion measures by condition
print(dispersion_by_condition)

# Test for homegeneity of variance using Levene test
levene_resultALL <- leveneTest(ITIValue ~ Condition, itiData)
print(levene_resultALL)

# Fligner-Killeen test of homogeneity of variances
fligner_test_result <- fligner.test(ITIValue ~ Condition, data = itiData)
print(fligner_test_result)

#
#model <- lmer(ITIValue ~ Condition + (1|Subject), data=itiData)
#summary(model)

#model <- lmer(ITIValue ~ Condition + (1 + Condition|Subject), data=itiData)
#summary(model)

#modelError <- lmer(ITIError ~ Condition + (1 + Condition|Subject), data=itiData)
#summary(modelError)
#Anova(modelError, test="F")

#group by condition and participant to get only one median ITI per trial+ other metrics
MedITI <- itiData  %>%
  group_by(Subject, Condition,Trial) %>%
  summarise(medITI  = median(ITIValue, na.rm = TRUE),
            medError  = median(ITIError, na.rm = TRUE),
            meanITI = mean(ITIValue, na.rm = TRUE),
            sdITI  =  sd(ITIValue, na.rm = TRUE)
  ) %>%
  mutate(CV = sdITI / meanITI) %>%
  ungroup()


#group by condition and participant to get only one ITI per participant + other metrics
allMedITI <- MedITI %>%
  group_by(Subject, Condition) %>%
  summarise(avgMedITI   = median(medITI),
            avgMedError = median(medError),
            avgMeanITI  = mean(meanITI),
            avgsdITI    = mean(sdITI),
            avgCV       = mean(CV),
            .groups = 'drop')

#--------------------------------------------------------------------------------------------------------------------------------
# Statistical Analysis
# ITI error
#Compute descriptive stats
boxplot_stats_ERROR <- allMedITI%>%
  group_by(Condition) %>%
  summarise(
    Min = min(avgMedError, na.rm = TRUE),
    Q1 = quantile(avgMedError, 0.25, na.rm = TRUE),
    Median = median(avgMedError, na.rm = TRUE),
    Q3 = quantile(avgMedError, 0.75, na.rm = TRUE),
    Max = max(avgMedError, na.rm = TRUE)
  ) %>%
  ungroup()


#MEDIAN ITIT error 
results_normality_Error <- allMedITI %>%
  group_by(Condition) %>%
  summarise(
    shapiro_p_value_medError = shapiro.test(avgMedError)$p.value
  )

print(results_normality_Error)

# Test for homegeneity of variance using Levene test
levene_result_error <- leveneTest(avgMedError~ Condition, allMedITI)
print(levene_result_error)

# Perform the Friedman test for median ITIs
friedman_test_result_error <- with(allMedITI, friedman.test(avgMedError, Condition, Subject))

# Number of subjects
num_subjects <- length(unique(allMedITI$Subject))

# Number of conditions - make sure this is correct for your dataset
num_conditions <- length(unique(allMedITI$Condition))

# Calculate Kendall's W
kendalls_w_error <- friedman_test_result_error$statistic / (num_subjects * (num_conditions - 1))

# Print Kendall's W
print(kendalls_w_error)
print(friedman_test_result_error)




## ITIs

#Extract boxplot stats 
boxplot_stats_ITI <- allMedITI%>%
  group_by(Condition) %>%
  summarise(
    Min = min(avgMedITI, na.rm = TRUE),
    Q1 = quantile(avgMedITI, 0.25, na.rm = TRUE),
    Median = median(avgMedITI, na.rm = TRUE),
    Q3 = quantile(avgMedITI, 0.75, na.rm = TRUE),
    Max = max(avgMedITI, na.rm = TRUE)
  ) %>%
  ungroup()

print( boxplot_stats_ITI)

#Check for normality (most probably not normally distributed)
results_normality <- allMedITI %>%
  group_by(Condition) %>%
  summarise(
    shapiro_p_value_median = shapiro.test(avgMedITI)$p.value,
    shapiro_p_value_mean = shapiro.test(avgMeanITI)$p.value,
    shapiro_p_value_sd = shapiro.test(avgsdITI)$p.value,
    shapiro_p_value_cv = shapiro.test(avgCV)$p.value
  )

print(results_normality)

# Ensure Condition is treated as a factor and its levels are in the correct order if specific ordering is desired
MedITI$Condition <- as.factor(MedITI$Condition) # Adjust levels as needed
allMedITI$Condition <- as.factor(allMedITI$Condition) # Adjust levels as needed

# Test for homegeneity of variance using Levene test
levene_result <- leveneTest(medITI ~ Condition, MedITI)
print(levene_result)

# Perform the Friedman test for median ITIs
friedman_test_result_median <- friedman.test(avgMedITI ~ Condition | Subject, data = allMedITI)


# Number of subjects
num_subjects <- length(unique(allMedITI$Subject))

# Number of conditions - make sure this is correct for your dataset
num_conditions <- length(unique(allMedITI$Condition))

# Calculate Kendall's W
kendalls_w <- friedman_test_result_median$statistic / (num_subjects * (num_conditions - 1))





# Obtain the number of subjects and conditions
num_subjects <- length(unique(allMedITI$Subject))
num_conditions <- length(unique(allMedITI$Condition))

# Calculate Kendall's W
kendalls_w <- friedman_test_result_median$statistic / (num_subjects * (num_conditions - 1))

# Print the results
print(friedman_test_result_median)
print(paste("Kendall's W for Median ITIs:", kendalls_w))

# Print Kendall's W
print(kendalls_w)
print(friedman_test_result_median)

# Test for homegeneity of variance using Levene test
levene_result <- leveneTest(avgMedITI ~ Condition, allMedITI)
print(levene_result)

# Perform the Friedman test for mean
friedman_test_result_mean <- with(allMedITI, friedman.test(avgMeanITI, Condition, Subject))
print(friedman_test_result_mean)

# Perform the Friedman test for CV
friedman_test_result_cv <- with(allMedITI, friedman.test(avgCV, Condition, Subject))
print(friedman_test_result_cv)


if (FALSE) {
#Raincloud plot with ALL ITI Values
plotALLITIs <- ggplot(itiData, aes(x = Condition, y = ITIValue, fill = Condition, colour = Condition)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE) +
  geom_point(position = position_jitter(width = .15), size = .25) +
  geom_boxplot(aes(x = as.numeric(Condition) + 0.25, y = ITIValue), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('ITI') + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
  ggtitle("ALL ITIs") +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0, 1.6)) + # Adjusting y-axis limits and breaks
  scale_colour_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors)+
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", size = 0.5) +  # Dashed line at y = 0.8
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black", size = 0.5)   # Dashed line at y = 0.4


# Display the plot
print(plotALLITIs)
ggsave("ALLITIs.svg", plot = plotALLITIs, width = 8, height = 6, device = 'svg')

#Raincloud plot with ALL Error Values
plotALLITIErrors <- ggplot(itiData, aes(x = Condition, y = ITIError, fill = Condition, colour = Condition)) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE) +
  geom_point(position = position_jitter(width = .20), size = .25) +
  geom_boxplot(aes(x = as.numeric(Condition) + 0.15, y = ITIError), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('ITI Error') + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  ggtitle("ALL ITI Errors")



# Display the plot
print(plotALLITIErrors)
ggsave("ALLITIs.svg", plot = plotALLITIs, width = 8, height = 6, device = 'svg')

# Plot for each subject 
# Get unique subjects
unique_subjects <- unique(itiData$Subject)

# Loop over each subject
for(subject in unique_subjects) {
  # Subset the data for the current subject
  subject_data <- itiData[itiData$Subject == subject, ]
  
  # Create the plot for the current subject
  plotALLITIs <- ggplot(subject_data, aes(x = Condition, y = ITIValue, fill = Condition, colour = Condition)) +
    geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 2, trim = FALSE, alpha = 0.3) +
    geom_point(position = position_jitter(width = .15), size = .25) +
    geom_boxplot(aes(x = as.numeric(as.factor(Condition)) + 0.25, y = ITIValue), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab('ITI') + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
    ggtitle(paste("ALL ITIs - Subject", subject)) +
    theme(plot.title = element_text(hjust = 0.5)) +  # Ensure title is centered 
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0, 1.6)) +
    scale_colour_manual(values = matlab_colors) +
    scale_fill_manual(values = matlab_colors) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", size = 0.5) +
    geom_hline(yintercept = 0.6, linetype = "dashed", colour = "black", size = 0.5) +
    geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black", size = 0.5)
  
  # Display the plot (optional)
  print(plotALLITIs)
  
  # Save the plot
  ggsave(paste0("ALLITIs_Subject_", subject, ".svg"), plot = plotALLITIs, width = 6, height = 6, device = 'svg')
}

# PLOT MEDIAN ITIs of each trial
plotMedITIsTrials <- ggplot(MedITI, aes(x = Trial, y = medITI, fill = Condition, colour = Condition)) +
  geom_flat_violin(adjust = 2, trim = FALSE, alpha = 0.3) +
  geom_point(position = position_jitter(width = .15), size = 1, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = .1, position = position_dodge(0.9), colour = "BLACK") +
  ylab('Median ITI') +
  xlab('Trial') +
  theme_cowplot() +
  guides(fill = FALSE, colour = FALSE) +
  ggtitle("Median ITIs Across Conditions and Trials") +
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1.2), limits = c(0, 1.6)) +
  facet_wrap(~ Condition, ncol = 1, scales = "free_x") +
  facet_wrap(~ Condition, ncol = 1, scales = "free_x") +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", size = 0.5) +
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black", size = 0.5)+
  scale_fill_manual(values = matlab_colors) +  # Assigning colors based on Condition
  scale_colour_manual(values = matlab_colors)  # Consistent color for lines and points

print(plotMedITIsTrials)
ggsave("MedITIsTrials.svg", plot = plotMedITIsTrials, width = 20, height = 15, units = "cm", device = 'svg')

# PLOT MEDIAN ITIs of each trial per condition (4 Median ITIs value/condition/ participant)
plotMedITIs <- ggplot(MedITI, aes(x = Condition, y = medITI, fill = Condition, colour = Condition)) +
  geom_flat_violin(position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE,, alpha = 0.3) +
  geom_point(position = position_jitter(width = .08), size = 1) +
  geom_boxplot(aes(x = as.numeric(Condition) + 0.25, y = medITI), outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('Median ITI of each trial') + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
  ggtitle("Median ITIs Across Conditions and Trials") +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), limits = c(0, 1.6)) + # Adjusting y-axis limits and breaks
  scale_colour_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors)+
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", size = 0.5) +  # Dashed line at y = 0.8
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black", size = 0.5)   # Dashed line at y = 0.4

print(plotMedITIs)
ggsave("MedITIs.svg", plot = plotMedITIs, width = 8, height = 6, device = 'svg')
}

# -------------------GOOD PLOT -----------------------------------
#-------------------------------------------------------------------------------------------------------------------------


# PLOT and save Average MEDIAN ITIs for each condition (4 Median ITIs value/condition/ participant)
plotAvgITIs<- 
  
  
  # Plot setup
  ggplot(allMedITI, aes(x = Condition, y = avgMedITI, fill = Condition, colour = Condition, )) +
  # Violin plot
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), adjust = 0.4, trim = FALSE, alpha = 0.3) +
  # Jittered points
  geom_point(aes(color = Condition, fill = Condition), position = position_jitterdodge(dodge.width = 0.15), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  # Boxplot
  geom_boxplot(alpha = 0.3, width = .2, color = "black", fill = NA, outlier.shape = NA, position = position_dodge(width = 0.15)) +
  
  ylab("Inter-tap Interval (s)") + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
  ggtitle("Average Median ITI") +
  scale_y_continuous(breaks = c(0, 0.4, 0.6,0.8, 1.2), limits = c(0, 1.25)) +
  scale_colour_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", size = 0.5) +  # Dashed line at y = 0.8
  #geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black", size = 0.5) +
  theme_minimal(base_family = "Arial") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 24),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 24),
    legend.text = element_text(size = 24),
    axis.line = element_line(color = "black", size = 1.2)  # Add physical axis lines
  ) +
  guides(color = guide_legend(title = "Stim Coch Zscore", override.aes = list(size = 4)))

# Print the plot
print(plotAvgITIs)
ggsave("AvgMedianITIs.svg", plot = plotAvgITIs,  width = 20, height = 15, units = "cm", device = 'svg')

if (FALSE) {
# PLOT and save Average MEDIAN ITI ERRORs for each condition (4 Median ITIs value/condition/ participant)
plotAvgITIError<- ggplot(allMedITI, aes(x = Condition, y =avgMedError, fill = Condition, colour = Condition)) +
  geom_flat_violin(position = position_nudge(x = .3, y = 0), adjust = 0.5, trim = TRUE,, alpha = 0.3) +
  geom_point(aes(color = Condition, fill = Condition), position = position_jitter(width = 0.1), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1, alpha = 0.5) +
  geom_boxplot(aes(x = as.numeric(Condition) , y = avgMedError), alpha = 0.1,width = .18, position = position_nudge(x = 0.3, y = 0), colour = "BLACK") +
  ylab("Inter-tap Interval Error (s)") + xlab('Duty Cycle (%)') + theme_cowplot() + guides(fill = FALSE, colour = FALSE) +
  ggtitle('Median ITI Error') +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_y_continuous(breaks = c(0, 5,10,15, 20), limits = c(0, 20)) + # Adjusting y-axis limits and breaks
  scale_colour_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors)+
  theme_minimal(base_family = "Arial") +  # Using Arial as a fallback to Calibri
  # Custom theme settings
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center and format title
    legend.position = "none",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    text = element_text(size = 20),  # Global text size
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 20),  # Axis text size
    legend.text = element_text(size = 16),  # Adjust legend text size
    axis.line = element_line(color = "black", size = 1.2)  # Add physical axis lines
    
  ) +
  
  guides(color = guide_legend(title = "Stim Coch Zscore", override.aes = list(size = 4)))


print(plotAvgITIError)
ggsave("AvgITIError.svg", plot = plotAvgITIError,  width = 20, height = 15, units = "cm", device = 'svg')


#--------------------------------------------------------------------------------------------------------------------------------
# Statistical Analysis

#Extract boxplot info
boxplot_stats_stability <- stability_data %>%
  group_by(DutyCycle) %>%
  summarise(
    Min = min(Rvalue, na.rm = TRUE),
    Q1 = quantile(Rvalue, 0.25, na.rm = TRUE),
    Median = median(Rvalue, na.rm = TRUE),
    Mean = mean(Rvalue),
    Q3 = quantile(Rvalue, 0.75, na.rm = TRUE),
    Max = max(Rvalue, na.rm = TRUE)
  ) %>%
  ungroup()

#Check for normality (most probably not normally distributed)
results_normality_stab <- stability_data %>%
  group_by(DutyCycle) %>%
  summarise(
    shapiro_p_value_stability = shapiro.test(Rvalue)$p.value)

print(results_normality_stab)

# Perform the Friedman test for median
friedman_test_result_stability <- with(stability_data, friedman.test(Rvalue, DutyCycle, SubjectID))

# Calculate Kendall's W
kendalls_w_stability <- friedman_test_result_stability$statistic / (num_subjects * (num_conditions - 1))

print(friedman_test_result_stability)
print(kendalls_w_stability)

# Test for homegeneity of variance using Levene test
levene_result_stab <- leveneTest(Rvalue ~ DutyCycle, stability_data)
print(levene_result_stab)


# Perform the Conover posthoc test
all_pval <- frdAllPairsConoverTest(stability_data$Rvalue,stability_data$DutyCycle,stability_data$SubjectID,p.adjust = "holm")
print(all_pval)















}

## ----------------------------------------------------  FINAL PUBLICATION FIGURE WITH 3 PANELS --------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Adjust margins to prevent label overlap
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
plotAvgITIs <- plotAvgITIs + theme_adjusted



# Create bottom row (B & C side by side with equal width)
bottom_row <- plot_grid(
  plotFFT, 
  plotACF, 
  labels = c("B.", "C."), 
  label_size = 30, 
  ncol = 2, 
  align = "h", 
  rel_widths = c(1, 1.1)  # Give ACF slightly more width to prevent overlap
)

# Center A (plotAvgITIs) with respect to B & C
top_row <- plot_grid(
  NULL, plotAvgITIs, NULL,  
  ncol = 3,  
  labels = c("", "A.", ""),  
  label_size = 30,  
  label_x = c(0, 0.05, 0),  # Adjust label placement
  rel_widths = c(0.5, 1, 0.5)  # Center A with empty space on both sides
)

# Arrange the final layout
final_plot <- plot_grid(
  top_row, bottom_row,  
  ncol = 1,               
  rel_heights = c(1, 1)  # Ensure equal height proportions
)


# Display the final plot
print(final_plot)




# Save the final figure
ggsave("FFT_ACF_combined_with_ITIs.svg", plot = final_plot, width = 50, height = 35, units = "cm", device = "svg")
