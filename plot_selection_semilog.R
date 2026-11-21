library(tidyverse)
library(ggrepel)
library(scales)

# 1. Load and verify columns
df_dnds <- read_csv("all_hotspots_dnds.csv")

# 2. Semi-log Selection Plot
ggplot(df_dnds, aes(x = dS, y = dN)) +
  geom_point(aes(color = Selection, size = Ratio), alpha = 0.7) +
  scale_color_manual(values = c("Positive" = "#e74c3c", "Purifying" = "#3498db")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.5) +
  # Apply the log transformation safely
  scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0, 1, 10, 50, 100, 300)) +
  geom_text_repel(data = subset(df_dnds, Ratio > 3),
                  aes(label = gsub(".*_", "", GeneID)), 
                  size = 3, fontface = "bold") +
  theme_bw() +
  labs(title = "MPXV Selection Pressure (Semi-Log)",
       x = "Synonymous (dS)", y = "Non-Synonymous (dN) - Log Scale")

ggsave("R_mpxv_selection_log.png", width = 10, height = 7, dpi = 300)

