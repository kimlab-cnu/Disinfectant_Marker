getwd()

setwd("./input/")  # Change the directory to your data file storage

Sys.glob("*")

traj <- as.data.frame(read.csv("HD_clinical_traj.CSV"))

library(dplyr); library(tidyr); library(ggplot2)

df_long <- traj %>%
  pivot_longer(cols = starts_with("z"), names_to   = c("measure_id", "time"),
    names_pattern = "(.+)([0-9]{2})$", values_to  = "zscore"
  ) %>%
  mutate(
    time = factor(
      time,
      levels = c("01", "02", "03", "04"),
      labels = c("First", "time2", "time3", "Last")
    ),
    measure = case_when(
      measure_id == "zfev1fvc" ~ "FEV1/FVC",
      measure_id == "zfvc" ~ "FVC",
      measure_id == "zfev" ~ "FEV1",
      measure_id == "zfef" ~ "FEF25–75",
      TRUE ~ measure_id
    )
  )

# summary
summary_df <- df_long %>%
  group_by(group, measure, time) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE),
            .groups = "drop")

# visualization
plot_traj <- function(summary_df, group_value,
                                  y_limits = c(-4, 1.5),
                                  y_breaks = seq(-3, 1, by = 1),
                                  base_size = 14,
                                  show_legend = FALSE) {
  
  p <- summary_df %>%
    filter(group == group_value) %>%
    ggplot(aes(x = time, y = mean_z,
               group = measure, color = measure)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    labs(x = NULL, y = "Z-score") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = 13),
      axis.title.y = element_text(
        margin = margin(r = 5),
        face = "plain"
      ),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      plot.title = element_blank()
    )
  
  return(p)
}

nor_p <- plot_traj(summary_df, "Normal")
obs_p <- plot_traj(summary_df, "Obstructive")
res_p <- plot_traj(summary_df, "Restrictive")
sev_p <- plot_traj(summary_df, "Severe")

nor_p
obs_p
res_p
sev_p

getwd()
setwd("../output/")

ggsave("z-score_exposed_traj_normal.png", plot=nor_p, dpi=300)
ggsave("z-score_exposed_traj_obstructive.png", plot=obs_p, dpi=300)
ggsave("z-score_exposed_traj_restrictive.png", plot=res_p, dpi=300)
ggsave("z-score_exposed_traj_severe.png", plot=sev_p, dpi=300)
