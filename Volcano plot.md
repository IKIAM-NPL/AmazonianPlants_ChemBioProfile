Characterization of four Amazonian plant extracts by liquid chromatography-mass spectrometry (LC-MS) data analysis:
================
Jhesel Almeida
2025-22-02

# Install necessary packages
install.packages("readxl")    # To read Excel files
install.packages("ggplot2")   # To create plots
install.packages("dplyr")     # For data manipulation
install.packages("ggrepel")   # To avoid overlapping labels in plots

# Load the installed packages
library(readxl)    # To read Excel files
library(ggplot2)   # To create plots
library(dplyr)     # For data manipulation
library(ggrepel)   # To avoid overlapping labels in plots

# Load data from the Excel file
file <- "Data/Volcano_POS.xlsx"  # Ensure the file is in your working directory
data <- read_excel(file)  # Read the Excel file into a dataframe

# Convert commas to dots in numeric values and ensure they are numeric
data <- data %>%
  mutate(across(Kp_D1:Pa_F3, ~ as.numeric(gsub(",", ".", .))))
# This step replaces commas with dots in numeric columns and converts them to numeric type.

# Calculate means for each condition
data <- data %>%
  rowwise() %>%
  mutate(
    Kp_D_mean = mean(c(Kp_D1, Kp_D2, Kp_D3), na.rm = TRUE),
    Kp_F_mean = mean(c(Kp_F1, Kp_F2, Kp_F3), na.rm = TRUE),
    Ws_D_mean = mean(c(Ws_D1, Ws_D2, Ws_D3), na.rm = TRUE),
    Ws_F_mean = mean(c(Ws_F1, Ws_F2, Ws_F3), na.rm = TRUE),
    Gn_D_mean = mean(c(Gn_D1, Gn_D2, Gn_D3), na.rm = TRUE),
    Gn_F_mean = mean(c(Gn_F1, Gn_F2, Gn_F3), na.rm = TRUE),
    Pa_D_mean = mean(c(Pa_D1, Pa_D2, Pa_D3), na.rm = TRUE),
    Pa_F_mean = mean(c(Pa_F1, Pa_F2, Pa_F3), na.rm = TRUE)
  ) %>%
  ungroup()
# This calculates the mean of triplicate measurements for each condition, ignoring NA values.

# Calculate log2FoldChange (log2FC) and p-values using t-tests
data <- data %>%
  rowwise() %>%
  mutate(
    log2FC_Kp = log2(Kp_F_mean / Kp_D_mean),
    log2FC_Ws = log2(Ws_F_mean / Ws_D_mean),
    log2FC_Gn = log2(Gn_F_mean / Gn_D_mean),
    log2FC_Pa = log2(Pa_F_mean / Pa_D_mean),
    
    # Check explicitly if there is enough data to perform t-tests
    p_value_Kp = {
      kp_f_values <- c(Kp_F1, Kp_F2, Kp_F3)
      kp_d_values <- c(Kp_D1, Kp_D2, Kp_D3)
      if (sum(!is.na(kp_f_values)) > 1 && sum(!is.na(kp_d_values)) > 1) {
        t.test(kp_f_values, kp_d_values)$p.value
      } else {
        NA_real_
      }
    },
    
    p_value_Ws = {
      ws_f_values <- c(Ws_F1, Ws_F2, Ws_F3)
      ws_d_values <- c(Ws_D1, Ws_D2, Ws_D3)
      if (sum(!is.na(ws_f_values)) > 1 && sum(!is.na(ws_d_values)) > 1) {
        t.test(ws_f_values, ws_d_values)$p.value
      } else {
        NA_real_
      }
    },
    
    p_value_Gn = {
      gn_f_values <- c(Gn_F1, Gn_F2, Gn_F3)
      gn_d_values <- c(Gn_D1, Gn_D2, Gn_D3)
      if (sum(!is.na(gn_f_values)) > 1 && sum(!is.na(gn_d_values)) > 1) {
        t.test(gn_f_values, gn_d_values)$p.value
      } else {
        NA_real_
      }
    },
    
    p_value_Pa = {
      pa_f_values <- c(Pa_F1, Pa_F2, Pa_F3)
      pa_d_values <- c(Pa_D1, Pa_D2, Pa_D3)
      if (sum(!is.na(pa_f_values)) > 1 && sum(!is.na(pa_d_values)) > 1) {
        t.test(pa_f_values, pa_d_values)$p.value
      } else {
        NA_real_
      }
    }
  ) %>%
  ungroup()
# This calculates log2 fold change and p-values for each condition using t-tests, ensuring there are enough non-NA values.

# Calculate -log10(p-value)
data <- data %>%
  mutate(
    neg_log10_p_Kp = -log10(p_value_Kp),
    neg_log10_p_Ws = -log10(p_value_Ws),
    neg_log10_p_Gn = -log10(p_value_Gn),
    neg_log10_p_Pa = -log10(p_value_Pa)
  )
# This converts p-values to -log10 scale for better visualization in volcano plots.

# Define colors based on significance
data <- data %>%
  mutate(
    Significance_Kp = case_when(
      p_value_Kp < 0.05 & abs(log2FC_Kp) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    Significance_Ws = case_when(
      p_value_Ws < 0.05 & abs(log2FC_Ws) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    Significance_Gn = case_when(
      p_value_Gn < 0.05 & abs(log2FC_Gn) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    Significance_Pa = case_when(
      p_value_Pa < 0.05 & abs(log2FC_Pa) > 1 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )
# This assigns significance labels based on p-value and log2 fold change thresholds.

# Combine data for Kp, Ws, Gn, and Pa into a single dataframe
combined_data <- data %>%
  select(Metabolite_name, log2FC_Kp, neg_log10_p_Kp, Significance_Kp) %>%
  rename(log2FC = log2FC_Kp, neg_log10_p = neg_log10_p_Kp, Significance = Significance_Kp) %>%
  mutate(Condition = "Kp") %>%
  bind_rows(
    data %>%
      select(Metabolite_name, log2FC_Ws, neg_log10_p_Ws, Significance_Ws) %>%
      rename(log2FC = log2FC_Ws, neg_log10_p = neg_log10_p_Ws, Significance = Significance_Ws) %>%
      mutate(Condition = "Ws")
  ) %>%
  bind_rows(
    data %>%
      select(Metabolite_name, log2FC_Gn, neg_log10_p_Gn, Significance_Gn) %>%
      rename(log2FC = log2FC_Gn, neg_log10_p = neg_log10_p_Gn, Significance = Significance_Gn) %>%
      mutate(Condition = "Gn")
  ) %>%
  bind_rows(
    data %>%
      select(Metabolite_name, log2FC_Pa, neg_log10_p_Pa, Significance_Pa) %>%
      rename(log2FC = log2FC_Pa, neg_log10_p = neg_log10_p_Pa, Significance = Significance_Pa) %>%
      mutate(Condition = "Pa")
  )
# This combines all conditions into one dataframe for easier plotting.

# Verify the combined data
print(head(combined_data))  # Check the combined data

# Create the combined Volcano Plot with separate labels
p <- ggplot(combined_data, aes(x = log2FC, y = neg_log10_p, color = Significance, shape = Condition)) +
  geom_point(alpha = 0.8, size = 5) +
  geom_label_repel(
    aes(label = ifelse(Significance == "Significant", Metabolite_name, NA)), 
    hjust = 0.5, vjust = 0.5,  # Center labels on points
    size = 3, 
    box.padding = 0.5,  # Spacing around text
    point.padding = 0.5,  # Spacing around points
    max.overlaps = 30,   # Maximum overlaps before relocating
    min.segment.length = 0.5,  # Minimize the length of connecting segments
    label.padding = 0.2,  # Spacing between text and frame
    label.size = 0.5,  # Thickness of frame border
    fill = "white",    # Background color of frame
    color = "black"    # Border color of frame
  ) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
  scale_shape_manual(values = c("Kp" = 16, "Ws" = 17, "Gn" = 18, "Pa" = 15)) +  # Different shapes for each condition
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Combined Volcano",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  ggtitle("Volcano Plot") +
  theme(legend.title = element_blank())
  
print(p)

# Save the plot as a JPG file in the Results folder
ggsave("Results/Volcano_POS.jpg", plot = p, width = 8, height = 6, dpi = 300, units = "in")
# This saves the volcano plot as a high-resolution JPG file.