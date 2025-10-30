library(ggplot2)
library(latex2exp)
library(ggpubr)

# Load the lifespan sensitivity estimation results
set_1 <- read.csv("Lifespan_sensitivity_estimation_gammahalf_omega5.csv", header=TRUE)
set_2 <- read.csv("Lifespan_sensitivity_estimation_gamma2_omega5.csv", header=TRUE)
set_3 <- read.csv("Lifespan_sensitivity_estimation_gammahalf_omega15.csv", header=TRUE)
set_4 <- read.csv("Lifespan_sensitivity_estimation_gamma2_omega15.csv", header=TRUE)

# combine sets and remove any row that did not converge during optimization
real_MLE_free <- rbind(set_1,set_2,set_3,set_4)
real_MLE_free <- real_MLE_free[,-1]
real_MLE_free <- subset(real_MLE_free, code==0)

# Plot the Log-likelihood
ggplot(real_MLE_free, aes(x = nloglik)) +
  geom_histogram(
    aes(fill = set, color = set),
    binwidth = 2,
    alpha = 0.6,
    position = "identity",
    size = 0.5
  ) + 
  labs(
    x = "Negative Log-Likelihood Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  xlim(100, 350) + 
  scale_fill_manual(
    values = c(
      "combination 1" = "steelblue4", 
      "combination 2" = "cadetblue3",
      "combination 3" = "khaki3",
      "combination 4" = "palevioletred2"
    ),
    labels = list(
      expression(gamma == 1/0.5 ~ "," ~ omega == 1/5),
      expression(gamma == 1/2 ~ "," ~ omega == 1/5),
      expression(gamma == 1/0.5 ~ "," ~ omega == 1/15),
      expression(gamma == 1/2 ~ "," ~ omega == 1/15)
    )
  ) +
  scale_color_manual(
    values = c(
      "combination 1" = "steelblue4", 
      "combination 2" = "cadetblue3",
      "combination 3" = "khaki3",
      "combination 4" = "palevioletred2"
    ),
    labels = list(
      expression(gamma == 1/0.5 ~ "," ~ omega == 1/5),
      expression(gamma == 1/2 ~ "," ~ omega == 1/5),
      expression(gamma == 1/0.5 ~ "," ~ omega == 1/15),
      expression(gamma == 1/2 ~ "," ~ omega == 1/15)
    )
  ) +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))+
  theme(
    legend.position = c(0.99, 0.99),       # Top-right inside the plot (normalized coordinates)
    legend.justification = c(1, 1),  # Anchor the legend's top-right corner
    legend.background = element_rect(fill = "white", color = "black"),
    legend.box.background = element_rect(color = "black")
  )




