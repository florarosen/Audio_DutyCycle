# --------------------------------------------
# Stimulus: energy ratios (harmonics only, up to 45 Hz)
# What you asked:
# 1) Band ratios per duty cycle (<=5 Hz included)
# 2) Ratio of energy contained only in f3–f11 per duty cycle
# --------------------------------------------

library(dplyr)
library(stringr)

# Load
FreqAmps <- read.csv("/Users/florarosenzweig/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/Postprocess_R_analysis/Frequency_Amplitude_LongFormat.csv")

# Map TrackID -> DutyCycle (%)
FreqAmps <- FreqAmps %>%
  mutate(DutyCycle = case_when(
    str_detect(TrackID, "^50ms")  ~ "25%",
    str_detect(TrackID, "^100ms") ~ "50%",
    str_detect(TrackID, "^150ms") ~ "75%",
    str_detect(TrackID, "^200ms") ~ "100%",
    TRUE ~ TrackID
  )) %>%
  select(-TrackID) %>%
  mutate(DutyCycle = factor(DutyCycle, levels = c("25%", "50%", "75%", "100%")))

# ---- Harmonic grid (base = 1/2.4 Hz), up to 45 Hz ----
f0 <- 1/2.4
harmonic_freqs <- f0 * (1:floor(45 / f0))

# Robust matching of each Frequency to its harmonic index k (1,2,3,...)
# (Your CSV contains only harmonics, but this makes it safe.)
tol_hz <- 1e-3
FreqAmps <- FreqAmps %>%
  mutate(
    HarmonicIndex = vapply(Frequency, function(x) {
      diffs <- abs(x - harmonic_freqs)
      k <- which.min(diffs)
      if (length(k) == 0 || diffs[k] > tol_hz) NA_integer_ else k
    }, integer(1))
  ) %>%
  filter(!is.na(HarmonicIndex)) %>%
  filter(Frequency <= 45)

# ======================================================
# 1) Ratios of energy (amplitude) bands per duty cycle
#    with 5 Hz included in the <=5 Hz band.
#    Since you only have harmonics, bands are defined
#    by harmonic frequencies falling in each range.
# ======================================================

FreqAmps <- FreqAmps %>%
  mutate(Band = case_when(
    Frequency <= 5  ~ "0–5 Hz (incl. 5)",
    Frequency > 5  & Frequency <= 10 ~ "5–10 Hz",
    Frequency > 10 & Frequency <= 15 ~ "10–15 Hz",
    Frequency > 15 & Frequency <= 20 ~ "15–20 Hz",
    Frequency > 20 & Frequency <= 25 ~ "20–25 Hz",
    Frequency > 25 & Frequency <= 30 ~ "25–30 Hz",
    Frequency > 30 & Frequency <= 35 ~ "30–35 Hz",
    Frequency > 35 & Frequency <= 40 ~ "35–40 Hz",
    Frequency > 40 & Frequency <= 45 ~ "40–45 Hz",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Band))

# Total energy (sum amplitude) per DutyCycle over 0–45 Hz harmonics
# Harmonic index corresponding to ~5 Hz
harmonic_5hz_index <- which.min(abs(harmonic_freqs - 5))

# Sanity check (optional but recommended)
harmonic_freqs[harmonic_5hz_index]

totals_dc <- FreqAmps %>%
  group_by(DutyCycle) %>%
  summarise(
    TotalAmp_0_45 = sum(Amplitude, na.rm = TRUE),
    TotalAmp_0_45_no5Hz = sum(
      Amplitude[!(HarmonicIndex %in% k5_harmonics)],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Band energy + ratio per DutyCycle
band_ratios <- FreqAmps %>%
  group_by(DutyCycle, Band) %>%
  summarise(BandAmp = sum(Amplitude, na.rm = TRUE), .groups = "drop") %>%
  left_join(totals_dc, by = "DutyCycle") %>%
  mutate(BandRatio = BandAmp / TotalAmp_0_45) %>%
  arrange(DutyCycle, Band)

print(band_ratios)

# Optional: wide format for a clean table
band_ratios_wide <- band_ratios %>%
  select(DutyCycle, Band, BandRatio) %>%
  tidyr::pivot_wider(names_from = Band, values_from = BandRatio)

print(band_ratios_wide)

# ======================================================
# 2) Ratio of energy contained only in f3–f11 per duty cycle
#    Here "f3–f11" is interpreted as harmonic indices 3..11.
#    (i.e., 3*f0 through 11*f0)
#    Ratio is computed relative to TOTAL energy (0–45 Hz harmonics).
# ======================================================

f3_f11_ratios <- FreqAmps %>%
  group_by(DutyCycle) %>%
  summarise(
    TotalAmp_0_45 = sum(Amplitude, na.rm = TRUE),
    Amp_f3_f11 = sum(Amplitude[HarmonicIndex >= 3 & HarmonicIndex <= 11], na.rm = TRUE),
    Ratio_f3_f11 = Amp_f3_f11 / TotalAmp_0_45,
    .groups = "drop"
  ) %>%
  arrange(DutyCycle)

print(f3_f11_ratios)


# index of the harmonic closest to 5 Hz (should be 12)
k5 <- which.min(abs(harmonic_freqs - 5))

# all indices corresponding to 5 Hz and its harmonics: k5, 2*k5, 3*k5, ...
k5_harmonics <- which((1:length(harmonic_freqs)) %% k5 == 0)

# sanity check (optional)
print(harmonic_freqs[k5_harmonics])

f3_f11_ratios_no5Hz <- FreqAmps %>%
  group_by(DutyCycle) %>%
  summarise(
    TotalAmp_0_45_no5Hz = sum(
      Amplitude[HarmonicIndex != harmonic_5hz_index],na.rm = TRUE),
    Amp_f3_f11 = sum(Amplitude[HarmonicIndex >= 3 & HarmonicIndex <= 11], na.rm = TRUE),
    Ratio_f3_f11 = Amp_f3_f11 / TotalAmp_0_45_no5Hz,
    .groups = "drop"
  ) %>%
  arrange(DutyCycle)

print(f3_f11_ratios_no5Hz)

f3_f11_ratios_no5HzHarm <- FreqAmps %>%
  group_by(DutyCycle) %>%
  summarise(
    TotalAmp_0_45_no5HzHarm = sum(Amplitude[!(HarmonicIndex %in% k5_harmonics)], na.rm = TRUE),
    Amp_f3_f11 = sum(Amplitude[HarmonicIndex >= 3 & HarmonicIndex <= 11], na.rm = TRUE),
    Ratio_f3_f11_no5HzHarmDen = Amp_f3_f11 / TotalAmp_0_45_no5HzHarm,
    .groups = "drop"
  ) %>%
  arrange(DutyCycle)

print(f3_f11_ratios_no5HzHarm)



print(f3_f11_within_0_15)