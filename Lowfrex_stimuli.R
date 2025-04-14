# Load necessary libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(rstatix)

# Set working directory to the location of the CSV file (!!!!! IF DOESNT WORK, OPEN MANUALLY !!!)


# Read the CSV files into data frames
FreqAmps <- read.csv("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis/Frequency_Amplitude_LongFormat.csv")

#---------------------------------------------------------------------------------------------------------------------------------------------
## AMPLITUDE ANALYSIS
# Modify the TrackID column to DutyCycle with just the percentage as values
FreqAmps <- FreqAmps %>%
  mutate(DutyCycle = case_when(
    str_detect(TrackID, "^50ms") ~ "25%",
    str_detect(TrackID, "^100ms") ~ "50%",
    str_detect(TrackID, "^150ms") ~ "75%",
    str_detect(TrackID, "^200ms") ~ "100%",
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
    DutyCycle = factor(DutyCycle, levels = c("25%", "50%", "75%", "100%"))  # Ensure proper order
  )

# Create faceted bar plot
distrib_proportions <- ggplot(grouped_FreqAmps, aes(x = FrequencyGroup, y = Proportion, fill = FrequencyGroup)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ DutyCycle, ncol = 1, strip.position = "top") +  # Facet by DutyCycle, one per row
  labs(
    title = "Proportion of Total Amplitude by Frequency Group and Duty Cycle",
    x = "Frequency Group",
    y = "Proportion of Total Amplitude",
    fill = "Frequency Group"
  ) +
  scale_fill_brewer(palette = "Spectral") +  # Use a distinct color palette
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),  # Customize facet titles
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )

# Save the plot
ggsave("Stim_amplitude_distrib_proportion.svg", plot = distrib_proportions, width = 20, height = 30, units = "cm", device = 'svg')

below_5hz_data <- grouped_FreqAmps %>%
  filter(FrequencyGroup == "Below 5Hz")

print(below_5hz_data)

harm5hz_data <- grouped_FreqAmps %>%
  filter(FrequencyGroup == "5Hz and Harmonics")

print(harm5hz_data)

# Perform Kruskal-Wallis test for "Below 5Hz"
kruskal_test_result <- kruskal.test(Proportion ~ DutyCycle, data = below_5hz_data)
print(kruskal_test_result)

# Perform pairwise Wilcoxon test with Bonferroni correction
pairwise_result <- pairwise_wilcox_test(below_5hz_data, Proportion ~ DutyCycle, p.adjust.method = "bonferroni")
print(pairwise_result)
print(kruskal_test_result)tage_Lowfrex.svg", plot = p4, width = 20, height = 15, units = "cm", device = 'svg')