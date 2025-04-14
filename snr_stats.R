# Clear and load packages
rm(list = ls())
library(tidyverse) 
library(ggpubr)
library(rstatix)
library(ggplot2)
library(ggstatsplot)
library(lme4)
library(emmeans)
library(car)
library(MuMIn)
library(effectsize)
library(nortest)


# Set working directory to the location of the CSV file
setwd("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis")

# Read the CSV file into a data frame
data <- read.csv("z_snr_young.csv")
final_table_15Hz_EEG <- read.csv("final_table_15Hz_EEG_ordered.csv") # Averaged Individual Data 


# Rename variables with "X" prefix
names(data) <- gsub("^X", "", names(data))

# Replace '.' at the end with "%"
names(data) <- gsub("\\.$", "%", names(data))

# Add a row identifier based on row numbers
#data$Participant <- seq_len(nrow(data))

# Convert SNR datafrom wide to long format
long_data <- gather(data, key = "Duty_Cycle", value = "zscore_snr", -Subject)

# Convert Participant and Condition to factors
long_data$SubjectID <- as.factor(long_data$Subject)
long_data$DutyCycle <- factor(long_data$Duty_Cycle, levels = c("25", "50", "75", "100"))

final_table_15Hz_EEG$SubjectID <- as.factor(final_table_15Hz_EEG$SubjectID )
final_table_15Hz_EEG$DutyCycle <- factor(final_table_15Hz_EEG$DutyCycle, levels = c("25", "50", "75", "100"))

# Define a z-SNR threshold (e.g., z > 2.0 or z > 3.0)
z_threshold <- 3.0

# Add a column indicating whether the z-score exceeds the threshold
long_data$significant <- ifelse(long_data$zscore_snr > z_threshold, TRUE, FALSE)

# Optionally filter the data to include only significant z-scores
significant_data <- long_data %>% filter(significant == TRUE)

# Check normality for each group using Shapiro-Wilk test
normality_tests<-lapply(split(long_data$zscore_snr, long_data$DutyCycle), shapiro.test)

# Extracting the results into a data frame
normality_results <- data.frame(
  DutyCycle = names(normality_tests),
  W = sapply(normality_tests, function(x) x$statistic),
  p_value = sapply(normality_tests, function(x) x$p.value)
)


# Perform Levene's Test for homogeneity of variance
levene_test<- leveneTest(meanZscore ~ DutyCycle, data = final_table_15Hz)


# Extracting the results into a data frame
levene_results <- data.frame(
  Df = levene_test$`Df`,
  `F value` = levene_test$`F value`,
  `Pr(>F)` = levene_testz$`Pr(>F)`
)

# Fit the linear mixed-effects model
model <- lmer(zscore_snr ~ DutyCycle + (1|SubjectID), data = long_data)
residuals <- resid(model)

#test for normality of residuals
shapiro.test(residuals)

#Check assumptions visually for model
par(mfrow=c(1,2))
plot(resid(model) ~ fitted(model))
qqnorm(resid(model)); qqline(resid(model), col = "red")


#If normality of residuals is not respected : 
friedman.test(zscore_snr ~ DutyCycle | SubjectID, data = long_data)

#If normality of residuals IS respected 
# ANova
# Summary and diagnostics
summary(model)
anova_results <- Anova(model, test="F")
summary(anova_results)
print(anova_results)

# Compute partial eta-squared directly
eta_squared_results <- eta_squared(anova_results, partial = TRUE)
print(eta_squared_results)




# Calculate Marginal and Conditional R²
r2_values <- r.squaredGLMM(model)
marginal_r2 <- r2_values[1]  # Marginal R²
conditional_r2 <- r2_values[2]  # Conditional R²

# Calculate the pairwise comparisons for post-hoc tests
emms <- emmeans(model, ~ DutyCycle)
contrast_results <- contrast(emms, method = "pairwise", adjust = "bonferroni")

# Compare each DutyCycle's mean z-score to the z-threshold
contrast_to_threshold <- contrast(emms, method = list(
  `25% vs Threshold` = c(1, 0, 0, 0) - z_threshold,
  `50% vs Threshold` = c(0, 1, 0, 0) - z_threshold,
  `75% vs Threshold` = c(0, 0, 1, 0) - z_threshold,
  `100% vs Threshold` = c(0, 0, 0, 1) - z_threshold
), adjust = "bonferroni")

# Display results
print(contrast_to_threshold)

contrast_results_df <- data.frame(
  group1 = c("Duty_Cycle25", "Duty_Cycle25", "Duty_Cycle25", "Duty_Cycle50", "Duty_Cycle50", "Duty_Cycle75"),
  group2 = c("Duty_Cycle50", "Duty_Cycle75", "Duty_Cycle100", "Duty_Cycle75", "Duty_Cycle100", "Duty_Cycle100"),
  estimate = c(-1.0230848, -1.9400788, -0.5726399, -0.9169941, 0.4504448, 1.3674389),
  SE = c(0.3848844, 0.3848844, 0.3848844, 0.3848844, 0.3848844, 0.3848844),
  df = c(69, 69, 69, 69, 69, 69),
  t.ratio = c(-2.658161, -5.040679, -1.487823, -2.382518, 1.170338, 3.552856),
  p.value = c(0.05853941, 2.144046e-05, 0.8481215, 0.1197289, 1.0000000, 0.004154226),
  y.position = c(8, 10, 8, 8, 8, 11),
  stars = c("ns", "***", "ns", "ns", "ns", "**")
)

# Return the data frame for kable
list(
  anova_results = as.data.frame(anova_results),
  contrast_results_df = contrast_results_df
)

# Filter only significant results (excluding "ns")
significant_contrasts <- contrast_results_df %>%
  filter(stars != "ns")

# Define colors for each condition (without '%')
matlab_colors <- c("25" = "#4DBEEE", "50" = "#0072BD", "75" ="#EDB120" , "100" = "#D95319")

# Create the base plot
p2 <- ggplot(long_data, aes(x = DutyCycle, y = zscore_snr, group = DutyCycle)) +
  
  # Jittered points for z_meter_rel with a darker circle outline and lighter fill
  geom_point(aes(color = DutyCycle, fill = DutyCycle), position = position_jitter(width = 0.1), size = 2, 
             show.legend = FALSE, shape = 21, stroke = 1.5, alpha = 0.5) +
  
  # Boxplot with no fill
  geom_boxplot(alpha = 0.3, width = .3, color = "black", outlier.shape = NA) +
  
  # Custom color scale for different DutyCycle conditions
  scale_color_manual(values = matlab_colors) +
  scale_fill_manual(values = matlab_colors) +
  scale_y_continuous(breaks = c(0,2,4, 6,8,10), limits = c(0, 12)) +  # Set Y-axis limits
  
  # Labels and title
  labs(title = "SNR", y = "Z-snr (a.u)", x = "Duty Cycle (%)") +
  
  theme_minimal() +
  # Custom theme settings
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),  # Center and format title
    legend.position = "top",
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    text = element_text(size = 20),  # Global text size
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 20),  # Axis text size
    legend.text = element_text(size = 16),  # Adjust legend text size
    axis.line = element_line(color = "black", size = 1.2)  # Add physical axis lines
  )

# Add significance bars for significant results only
p2 <- p2 +
  geom_signif(
    comparisons = list(
      c("25", "75"),
      c("75", "100")
    ),
    map_signif_level = TRUE,
    annotations = significant_contrasts$stars,
    y_position = significant_contrasts$y.position,
    textsize = 6,
    tip_length = 0.01
  )

# Display the plot
print(p2)


# Correlation between meter related Zscore and z-snr--------------------------------
# Function to calculate Spearman's correlation
get_spearman_correlation <- function(df) {
  correlation <- cor(df$meanZscore, df$zscore_snr, method = "spearman")
  p_value <- cor.test(df$meanZscore, df$zscore_snr, method = "spearman")$p.value
  
  return(data.frame(correlation = correlation, p_value = p_value))
}

# Combine the data into one data frame for analysis
# Assuming final_table_15Hz_EEG and long_data have a common identifier like 'SubjectID' and 'DutyCycle'
combined_data <- final_table_15Hz_EEG %>%
  inner_join(long_data, by = c("SubjectID", "DutyCycle"))

# Calculate Spearman's correlation for each Duty Cycle
correlation_results <- combined_data %>%
  group_by(DutyCycle) %>%
  do(get_spearman_correlation(.))

print(correlation_results)
# Define your plot for EEG vs SNR correlation
plot_eeg_snr <- ggplot(combined_data, aes(x = meanZscore, y = zscore_snr, color = DutyCycle, group = DutyCycle)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid") +
  labs(title = "Correlation between EEG FFT Meter-related Z-score and SNR Z-scores by Duty Cycle",
       x = "Mean Meter-related Z-score EEG",
       y = "SNR Z-score",
       color = "Duty Cycle (%)") +
  scale_color_manual(values = matlab_colors) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

# Display the plot
print(plot_eeg_snr)