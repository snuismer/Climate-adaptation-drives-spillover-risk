library(tidyverse)
library(ggtext)
library(optimx)
library(scales)
library(tinytable)

#==============================================================================


# Import village-level summary data and individual-level human serostatus data

v <- read_csv("Data/Raw/BruMERS_Human_Animals_Village_level_Summary_09_22_2025.csv")
h <- read_csv("Data/Raw/BruMers_Human_Data_09_23_2025.csv")

# Modify column names in both data frames

colnames(v) <- c(
  "county", "village", 
  "size_livestock", "size_bovine", "size_camel", "size_goat", "size_sheep",
  "n_tested_livestock", "n_tested_bovine", "n_tested_camel", "n_tested_goat", "n_tested_sheep",
  "n_positive_bovine", "n_positive_camel", "n_positive_goat", "n_positive_sheep",
  "n_negative_bovine", "n_negative_camel", "n_negative_goat", "n_negative_sheep",
  "n_tested_human", "n_positive_human", "n_negative_human"
)

colnames(h) <- c(
  "participant_ID", "household", "village", "manyatta", "county",
  "human_age", "human_age_group", "blood_sample", "serostatus"
)

# Modify data frames

v <- v %>%
  # Add proportion variables for each livestock species
  mutate(
    prop_bovine = size_bovine/size_livestock,
    prop_camel = size_camel/size_livestock,
    prop_goat = size_goat/size_livestock,
    prop_sheep = size_sheep/size_livestock
  ) %>%
  # Arrange by county and village
  arrange(county, village)

h <- h %>%
  # Add a numeric serostatus variable
  mutate(
    serostatus_numeric = case_when(
      serostatus == "Negative" ~ 0,
      serostatus == "Positive" ~ 1
    )
  ) %>%
  # Reorder columns
  select(
    county, village, manyatta, household, participant_ID, blood_sample,
    serostatus, serostatus_numeric, human_age, human_age_group
  ) %>%
  arrange(county, village, manyatta, household, participant_ID)

#==============================================================================


# Data validation

assertthat::assert_that(
  sum(
    v$size_livestock ==
      v$size_bovine + v$size_camel + v$size_goat + v$size_sheep
  ) == nrow(v)
)

assertthat::assert_that(
  sum(
    v$n_tested_livestock ==
      v$n_tested_bovine + v$n_tested_camel + v$n_tested_goat + v$n_tested_sheep
  ) == nrow(v)
)

assertthat::assert_that(
  sum(
    v$n_tested_human == v$n_positive_human + v$n_negative_human
  ) == nrow(v)
)

assertthat::assert_that(
  sum(v$n_tested_human) == nrow(h)
)

assertthat::assert_that(
  sum(v$n_positive_human) == sum(h$serostatus == "Positive")
)

assertthat::assert_that(
  sum(v$n_negative_human) == sum(h$serostatus == "Negative")
)

all.equal(
  current = unique(v$village), target = unique(h$village)
)

#==============================================================================


# Add expected herd lifespan values to the village data frame

bovine.lifespan <- 6
camel.lifespan <- 16
goat.lifespan <- 2
sheep.lifespan <- 3

v <- v %>%
  mutate(
    expected_herd_lifespan = 
      (size_bovine * bovine.lifespan + 
         size_camel * camel.lifespan + 
         size_goat * goat.lifespan +
         size_sheep * sheep.lifespan) / 
      size_livestock,
    expected_herd_lifespan_s = rethinking::standardize(expected_herd_lifespan),
  )

#==============================================================================


# Generate village-level Brucella FOI estimates using human serostatus data

# Define the Log Likelihood function of Parms where Parms[1] is lambda

LogLike <- function(Parms) {
  
  P <- 0
  
  for(i in 1:count) {
    
    P <- P + 
      (sub$serostatus_numeric[i]) * log(1 - exp(-sub$human_age[i]*Parms[1])) + 
      (1 - sub$serostatus_numeric[i]) * log(exp(-sub$human_age[i]*Parms[1]))
  }
  
  return(-P)
}

# Define a data frame in which to store the ML solutions

ml <- data.frame(matrix(ncol = 2, nrow = 0))

# Conduct ML analysis

for(site in unique(h$village)) {
  
  # Select the site and size it
  sub <- filter(h, village == site)
  county <- unique(sub$county)
  count <- length(sub$serostatus_numeric)
  seroprevalence <- mean(sub$serostatus_numeric)
  
  # Now maximize the Log Likelihood for lambda
  ml.sol <- optimx(
    par = c(0.01), fn = LogLike, method = "L-BFGS-B",
    lower = c(0.0001), upper = c(0.2), control = list(maxit = 10000)
  )
  row <- c(county, site, seroprevalence, ml.sol$p1)
  ml <- rbind(ml, row)
}

colnames(ml) <- c("county", "village", "seroprevalence", "lambda")

ml <- ml %>%
  mutate(
    seroprevalence = as.numeric(seroprevalence),
    lambda = as.numeric(lambda)
  ) %>%
  arrange(county, village)

ml

# Plot FOI against seroprevalence

ml %>%
  ggplot(aes(x = seroprevalence, y = lambda, color = county)) + 
  geom_point() +
  xlab("Brucella seroprevalence in humans") +
  ylab("Brucella force of infection into humans") +
  scale_x_continuous(limits = c(0, 0.3), labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = c(0, 0.015), labels = label_number(accuracy = 0.001)) +
  theme_minimal()

# Add these lambda estimates to the village data frame

v <- v %>%
  left_join(
    ., ml,
    by = c("county", "village")
  )

# Save village-level data with FOI info attached

write_csv(v, file = "Data/Clean/village_level_data_w_FOI.csv")

#==============================================================================


# Make a table of primary data

# Import village-level data with FOI information

v <- read_csv("Data/Clean/village_level_data_w_FOI.csv")

# Generate table

v %>%
  mutate(
    livestock_seroprevalence = 
      (n_positive_bovine + n_positive_camel + n_positive_goat + n_positive_sheep) / n_tested_livestock
  ) %>%
  mutate(
    size_camel = paste0(size_camel, " (", round(prop_camel, digits = 2), ")"),
    size_bovine = paste0(size_bovine, " (", round(prop_bovine, digits = 2), ")"),
    size_goat = paste0(size_goat, " (", round(prop_goat, digits = 2), ")"),
    size_sheep = paste0(size_sheep, " (", round(prop_sheep, digits = 2), ")"),
    livestock_seroprevalence = paste0(
      round(livestock_seroprevalence, digits = 3), " (", n_tested_livestock, ")"
    ),
    human_seroprevalence = paste0(
      round(seroprevalence, digits = 3), " (", n_tested_human, ")"
    ),
    lambda = round(lambda, digits = 3)
  ) %>%
  rename(
    Village = village,
    Camel = size_camel,
    Cattle = size_bovine,
    Goat = size_goat,
    Sheep = size_sheep,
    "*Brucella* seroprevalence in livestock (no. tested)" = livestock_seroprevalence,
    "*Brucella* seroprevalence in humans (no. tested)" = human_seroprevalence,
    "Estimated *Brucella* FOI into humans" = lambda
  ) %>%
  select(
    Village, 
    Camel, Cattle, Goat, Sheep,
    `*Brucella* seroprevalence in livestock (no. tested)`,
    `*Brucella* seroprevalence in humans (no. tested)`,
  ) %>%
  tt() %>%
  group_tt(
    i = list(
      "*Kajiado County*" = 1,
      "*Marsabit County*" = 9
    ),
    j = list(
      "Livestock count (prop. of village herd)" = 2:5
    )
  ) %>%
  style_tt(i = c(1, 10), background = "tan") %>%
  style_tt(j = 2:7, align = "c") %>%
  format_tt(markdown = TRUE) %>%
  save_tt("Figures/table.tex", overwrite = TRUE)

#==============================================================================


# Linear modeling

# Fit simpler linear models

m.l <- glm(lambda ~ expected_herd_lifespan_s, data = v)
m.c <- glm(lambda ~ county, data = v)
m.lc <- glm(lambda ~ expected_herd_lifespan_s + county, data = v)

AICcmodavg::aictab(list(m.l, m.c, m.lc), modnames = c("lifespan", "county", "lifespan + county"))

# Fit the model with an interaction between expected herd lifespan and county,
# setting Marsabit county to the base case

m.interaction <- glm(
  lambda ~ expected_herd_lifespan_s * county, 
  data = v %>% mutate(county = factor(county, levels = c("Marsabit", "Kajiado")))
)

summary(m.interaction)

# Sample the quadratic posterior to visualize the expected herd lifespan effects
# in both counties

# Sample the posterior
s <- rethinking::sample.qa.posterior(m.interaction)

# Generated the implied slope values in both Marsabit and Kajiado counties
post.samples <- data.frame(
  county = rep(c("Marsabit", "Kajiado"), each = nrow(s)),
  # Slopes in Marsabit, then Kajiado
  slope = c(
    s$expected_herd_lifespan_s,
    s$expected_herd_lifespan_s + s$`expected_herd_lifespan_s:countyKajiado`
  )
)

# Plot
post.samples %>%
  ggplot(aes(x = slope, color = county)) +
  geom_density() +
  xlab("Slope coefficient for expected herd lifespan") +
  theme_minimal()

#==============================================================================


# Create plots based on the interaction model

# Generate predictions for Kajiado and Marsabit, only considering observed
# values of expected herd lifespan

n.preds.per.county <- 100
min.Kajiado <- v %>%
  filter(county == "Kajiado") %>%
  pull(expected_herd_lifespan_s) %>%
  min()
max.Kajiado <- v %>%
  filter(county == "Kajiado") %>%
  pull(expected_herd_lifespan_s) %>%
  max()
min.Marsabit <- v %>%
  filter(county == "Marsabit") %>%
  pull(expected_herd_lifespan_s) %>%
  min()
max.Marsabit <- v %>%
  filter(county == "Marsabit") %>%
  pull(expected_herd_lifespan_s) %>%
  max()

main.text.size <- 18
linewidth <- 1.5


# Generate panel a, which will show the proportion of livestock types in
# each village herd as well as the expected herd lifespan in each village

a <- v %>%
  select(county, village, prop_bovine, prop_goat, prop_sheep, prop_camel) %>%
  tidyr::pivot_longer(
    cols = c("prop_bovine", "prop_goat", "prop_sheep", "prop_camel"),
    names_to = "species"
  ) %>%
  left_join(
    ., select(v, county, village, expected_herd_lifespan), by = c("county", "village")
  ) %>%
  mutate(
    species = case_when(
      species == "prop_bovine" ~ "cattle",
      species == "prop_goat" ~ "goats",
      species == "prop_sheep" ~ "sheep",
      species == "prop_camel" ~ "camels"
    ),
    village = factor(village, levels = v %>% arrange(expected_herd_lifespan) %>% pull(village)),
    expected_herd_lifespan = format(round(expected_herd_lifespan, 1), nsmall = 1)
  ) %>%
  ggplot(aes(x = village, y = value)) +
  geom_col(aes(fill = species)) +
  geom_label(aes(y = 0.08, label = expected_herd_lifespan), size = 5) +
  ylab("Proportion of herd") +
  facet_wrap(~county, scales = "free_x") +
  scale_fill_manual(values = c("burlywood", "black", "gray", "antiquewhite")) +
  theme_minimal() +
  theme(
    text = element_text(size = main.text.size),
    strip.text.x = element_text(size = main.text.size + 2),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = main.text.size - 2),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )


# Generate panel b, which will show the relationship between FOI and expected
# herd lifespan along with county-level confidence intervals

# Generate a data frame with expected lifespan values for each county for 
# which we want model predictions of FOI
data.for.preds <- data.frame(
  county = rep(c("Kajiado", "Marsabit"), each = n.preds.per.county),
  expected_herd_lifespan_s = c(
    seq(from = min.Kajiado, to = max.Kajiado, length.out = n.preds.per.county),
    seq(from = min.Marsabit, to = max.Marsabit, length.out = n.preds.per.county)
  )
) %>%
  mutate(
    expected_herd_lifespan = 
      (expected_herd_lifespan_s * sd(v$expected_herd_lifespan)) + mean(v$expected_herd_lifespan)
  )

# Generate model-based predictions for these expected lifespan values
predictions1 <- predict(m.interaction, newdata = data.for.preds, se.fit = TRUE)

predictions1 <- data.frame(
  county = data.for.preds$county,
  expected_herd_lifespan_s = data.for.preds$expected_herd_lifespan_s,
  expected_herd_lifespan = data.for.preds$expected_herd_lifespan,
  lambda = predictions1$fit,
  lower_ci = predictions1$fit - (1.96 * predictions1$se.fit),
  upper_ci = predictions1$fit + (1.96 * predictions1$se.fit)
)

# Plot
b <- ggplot(
  data = predictions1,
  aes(x = expected_herd_lifespan, y = lambda, group = county, color = county)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci, group = county),
    color = NA,
    fill = alpha("gainsboro", 0.5)
  ) +
  geom_line(
    linewidth = linewidth
  ) +
  geom_point(
    data = v,
    size = 4
  ) +
  scale_color_manual(
    values = c("darkseagreen3", "firebrick")
  ) +
  xlab("Expected herd lifespan") +
  ylab("*Brucella* FOI into humans") +
  xlim(3, 5.3) +
  ylim(0, 0.015) +
  theme_minimal() +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(
    text = element_text(size = main.text.size),
    axis.title.y = element_markdown(size = main.text.size -2),
    legend.position = "inside",
    legend.position.inside = c(0.83, 0.25),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = main.text.size + 2),
    legend.key.size = unit(0.5, "inch"),
    panel.grid.minor = element_blank()
  )

ggsave("Outputs/human_FOI_vs_lifespan.jpg", plot = b, width = 8, height = 6)


# Generate panel c, which will show the impact of increasing herd lifespan 
# in Kajiado beyond that which has been currently observed

# Generate a data frame with expected lifespan values for 
# which we want model predictions of FOI
data.for.preds <- data.frame(
  county = rep(c("Kajiado"), each = n.preds.per.county),
  expected_herd_lifespan_s = 
    seq(from = min.Kajiado, to = max.Marsabit, length.out = n.preds.per.county)
) %>%
  mutate(
    expected_herd_lifespan = 
      (expected_herd_lifespan_s * sd(v$expected_herd_lifespan)) + mean(v$expected_herd_lifespan)
  )

# Generate model-based predictions for these expected lifespan values
predictions2 <- predict(m.interaction, newdata = data.for.preds, se.fit = TRUE)

predictions2 <- data.frame(
  county = data.for.preds$county,
  expected_herd_lifespan_s = data.for.preds$expected_herd_lifespan_s,
  expected_herd_lifespan = data.for.preds$expected_herd_lifespan,
  lambda = predictions2$fit,
  lower_ci = predictions2$fit - (1.96 * predictions2$se.fit),
  upper_ci = predictions2$fit + (1.96 * predictions2$se.fit)
)

# What is the average livestock herd composition in each county?
v %>%
  group_by(county) %>%
  summarize(
    mean_prop_bovine = mean(prop_bovine),
    mean_prop_camel = mean(prop_camel),
    mean_prop_goat = mean(prop_goat),
    mean_prop_sheep = mean(prop_sheep)
  )

# Generate a table used for labeling panel c, which will be based on
# adding camels to the average livestock herd composition currently observed
# in Kajiado
herd.composition.table <- data.frame(
  n_bovine = c(32, 27, 22), # reduce by 5% each time
  n_goats = c(36, 36, 36), # stays the same, at average observed value
  n_sheep = c(32, 32, 32), # stays the same, at average observed value
  n_camel = c(0, 5, 10), # increase by 5% each time
  text = c(
    "0% camels\n32% cattle\n36% goats\n32% sheep",
    "5% camels\n27% cattle\n36% goats\n32% sheep",
    "10% camels\n22% cattle\n36% goats\n32% sheep"
  )
) %>%
  mutate(
    n_tot = n_bovine + n_camel + n_goats + n_sheep,
    expected_herd_lifespan = 
      (n_bovine * bovine.lifespan + 
         n_camel * camel.lifespan + 
         n_goats * goat.lifespan +
         n_sheep * sheep.lifespan) / n_tot,
    expected_herd_lifespan_s = (expected_herd_lifespan - mean(v$expected_herd_lifespan)) / sd(v$expected_herd_lifespan)
  )

# Add model-predicted Brucella FOI values for each scenario in the herd
# composition table
df <- herd.composition.table %>%
  select(expected_herd_lifespan_s) %>%
  mutate(county = "Kajiado") %>%
  mutate(
    lambda = predict(m.interaction, newdata = .)
  )

herd.composition.table$lambda <- df$lambda

# Plot
c <- ggplot(
  data = predictions2,
  aes(x = expected_herd_lifespan, y = lambda)
) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci, group = county),
    color = NA,
    fill = alpha("gainsboro", 0.5)
  ) +
  geom_line(
    linewidth = linewidth,
    lty = 2,
    color = "darkseagreen3"
  ) +
  geom_line(
    data = predictions1 %>% filter(county == "Kajiado"),
    linewidth = linewidth,
    color = "darkseagreen3"
  ) +
  geom_segment(
    data = herd.composition.table,
    aes(x = expected_herd_lifespan, y = 0.011, yend = lambda + 0.0003)
  ) +
  geom_label(
    data = herd.composition.table,
    aes(x = expected_herd_lifespan, y = 0.011, label = text),
    size = main.text.size/4
  ) +
  xlab("Expected herd lifespan") +
  ylab("*Brucella* FOI into humans") +
  xlim(3, 5.3) +
  ylim(0, 0.015) +
  theme_minimal() +
  theme(
    text = element_text(size = main.text.size),
    axis.title.y = element_markdown(size = main.text.size -2),
    panel.grid.minor = element_blank()
  )

cowplot::plot_grid(
  a, b, c,
  labels = "auto",
  label_size = 20,
  ncol = 1
)

ggsave(
  filename = "Outputs/human_FOI_vs_lifespan_multipanel.jpg",
  height = 12,
  width = 10,
  units = "in"
)

#==============================================================================


# Sensitivity analysis for expected lifespan

# Import randomly generated lifespans for each animal group
sensitivity.lifespans <- read_csv("Data/Sensitivity/Random_Lifespan.csv")

# Set up a data frame that will hold the parameter values and p-values for
# the four parameters in the interaction model, given each unique set of
# lifespan values for the livestock species
d.params <- data.frame(
  replicate = rep(
    1:nrow(sensitivity.lifespans), 
    each = 4
  ),
  parameter = rep(
    c("intercept", "slope", "county_effect", "interaction_effect"),
    times = nrow(sensitivity.lifespans)
  ),
  value = rep(NA, each = nrow(sensitivity.lifespans)),
  p_value = rep(NA, each = nrow(sensitivity.lifespans))
) %>%
  mutate(
    parameter = factor(parameter, levels = c("intercept", "slope", "county_effect", "interaction_effect"))
  )

# Set up a data frame that will hold the model predictions (i.e., Brucella FOI
# vs. expected herd lifespan relationship in each county), give each unique
# set of lifespan values for the livestock species
d.preds <- data.frame()

# Loop through all sets of lifespan values
for(i in 1:nrow(sensitivity.lifespans)) {
  
  # Generate a data frame that calculates the expected herd lifespan in each
  # village, given this unique set of lifespan values
  v.replicate <- v %>%
    mutate(
      expected_herd_lifespan_replicate = 
        (size_bovine * sensitivity.lifespans$Cow[i] +
           size_camel * sensitivity.lifespans$Camel[i] +
           size_goat * sensitivity.lifespans$Goat[i] +
           size_sheep * sensitivity.lifespans$Sheep[i]) /
        size_livestock,
      expected_herd_lifespan_replicate_s = rethinking::standardize(expected_herd_lifespan_replicate)
    )
  
  # Fit the interaction model, given these expected herd lifespan values
  m.replicate <- glm(
    lambda ~ expected_herd_lifespan_replicate_s * county, 
    data = v.replicate %>% 
      mutate(county = factor(county, levels = c("Marsabit", "Kajiado")))
  )
  
  # Populate the d.params data frame with the parameter values and p-values
  # from the interaction model fit on this iteration of lifespan values
  d.params[d.params$replicate == i, "value"] <- coef(m.replicate)
  d.params[d.params$replicate == i, "p_value"] <- summary(m.replicate)$coefficients[,4]
  
  # Generate model-based predictions, given this iteration of data and model fit
  min.Kajiado <- v.replicate %>%
    filter(county == "Kajiado") %>%
    pull(expected_herd_lifespan_replicate_s) %>%
    min()
  max.Kajiado <- v.replicate %>%
    filter(county == "Kajiado") %>%
    pull(expected_herd_lifespan_replicate_s) %>%
    max()
  min.Marsabit <- v.replicate %>%
    filter(county == "Marsabit") %>%
    pull(expected_herd_lifespan_replicate_s) %>%
    min()
  max.Marsabit <- v.replicate %>%
    filter(county == "Marsabit") %>%
    pull(expected_herd_lifespan_replicate_s) %>%
    max()
  
  data.for.preds <- data.frame(
    county = rep(c("Kajiado", "Marsabit"), each = n.preds.per.county),
    expected_herd_lifespan_replicate_s = c(
      seq(from = min.Kajiado, to = max.Kajiado, length.out = n.preds.per.county),
      seq(from = min.Marsabit, to = max.Marsabit, length.out = n.preds.per.county)
    )
  ) %>%
    mutate(
      expected_herd_lifespan_replicate = 
        (expected_herd_lifespan_replicate_s * sd(v.replicate$expected_herd_lifespan_replicate)) + mean(v.replicate$expected_herd_lifespan_replicate)
    )
  
  preds <- predict(m.replicate, newdata = data.for.preds, se.fit = TRUE)
  
  # Bind relevant outputs into the d.preds data frame
  temp <- data.for.preds %>%
    mutate(
      lambda = preds$fit,
      replicate = rep(i, nrow(.))
    )
  
  d.preds <- bind_rows(d.preds, temp)
}

# Organize the d.params data frame
d.params <- d.params %>%
  mutate(
    parameter = case_when(
      parameter == "intercept" ~ "Intercept",
      parameter == "slope" ~ "Expected herd lifespan effect",
      parameter == "county_effect" ~ "Kajiado effect",
      parameter == "interaction_effect" ~ "Expected herd lifespan x Kajiado interaction"
    ),
    parameter = factor(
      parameter, 
      levels = c(
        "Intercept", 
        "Expected herd lifespan effect", 
        "Kajiado effect", 
        "Expected herd lifespan x Kajiado interaction")
    )
  )

# Generate a data frame with parameter values and p-values for the actual
# observed dataset, for comparison with all the sensitivity iteraions
observed.params <- data.frame(
  parameter = factor(
    c(
      "Intercept", 
      "Expected herd lifespan effect", 
      "Kajiado effect", 
      "Expected herd lifespan x Kajiado interaction"
    )
  ),
  value = coef(m.interaction),
  p_value = summary(m.interaction)$coefficients[,4]
)

# Table summarizing range of parameter values and p-values observed using the
# sensitivity analysis data
d.params %>%
  group_by(parameter) %>%
  summarize(
    min_param = min(value),
    max_param = max(value),
    min_pvalue = min(p_value),
    max_pvalue = max(p_value)
  )

# Plot histogram of parameter values for each of the four parameters
d.params %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "gray") +
  xlab("parameter value") +
  geom_vline(data = observed.params, (aes(xintercept = value)), color = "red") +
  facet_wrap(~parameter, scales = "free_x") +
  theme_minimal() +
  theme(
    text = element_text(size = 16)
  )

ggsave(
  filename = "Outputs/sensitivity_parameter_value_histogram.jpg",
  height = 6,
  width = 8,
  units = "in"
)

# Plot histogram of p-values values for each of the four parameters
d.params %>%
  ggplot(aes(x = p_value)) +
  geom_histogram(fill = "gray") +
  xlab("p-value") +
  geom_vline(data = observed.params, (aes(xintercept = p_value)), color = "red") +
  facet_wrap(~parameter, scales = "free_x") +
  theme_minimal() +
  theme(
    text = element_text(size = 16)
  )

ggsave(
  filename = "Outputs/sensitivity_p-value_value_histogram.jpg",
  height = 6,
  width = 8,
  units = "in"
)
