#packages
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dagitty)
library(ggtext)


#####Define DAG#####
dag <- dagitty("dag {human -> food;
                     human -> house;
                     house -> green;
                     house -> diversity;
                     food -> diversity
                     green -> diversity
                     green -> protected
                     diversity -> protected;
}")

#####Identify adjustment sets#####

#direct effect of housing density on bird diversity
adjustmentSets(dag, exposure = "house", outcome = "diversity", effect="direct")

adjustmentSets(dag, exposure = "house", outcome = "diversity", effect="total")


#####Simualtion setttings##### 

n    <- 100 #sample size
sims <- 1000 #simulations

#direction and magnitude of direct effects
b.fp <- 0.9   # human -> food
b.hp <- 0.8   # human -> house
b.gh <- -1.0  # house -> green
b.sh <- 0.3   # house -> diversity
b.sf <- 0.4   # food -> diversity
b.sg <- 0.8   # green -> diversity
b.lg <- 1.0   # green -> protected
b.ls <- 0.8   # diversity -> protected

#derive true total effect of housing density on species diversity
b.tot.sh <- b.sh + b.gh * b.sg

#####Regression models to be fit#####
models <- c(
  "s ~ h + p",             # valid total
  "s ~ h + f",             # valid total
  "s ~ h + p + g",         # valid direct
  "s ~ h + f + g",         # valid direct
  "s ~ h",                 # biased (confounder)
  "s ~ h + p + f + g + l"  # biased (collider)
)
models <- lapply(models, as.formula)

# Display labels
heading_models <- toupper(models)
heading_accuracy <- c(
  "Valid total effect", "Valid total effect",
  "Valid direct effect", "Valid direct effect",
  "Biased by confounder", "Biased by collider"
)

#####Simulate#####
out <- matrix(NA, nrow = sims, ncol = length(models))
colnames(out) <- paste0("model", 1:length(models))

set.seed(999)
for (j in 1:sims) {
  p <- rnorm(n, 0, 1)
  f <- b.fp * p + rnorm(n, 0, 0.5)
  h <- b.hp * p + rnorm(n, 0, 0.5)
  g <- b.gh * h + rnorm(n, 0, 0.5)
  s <- b.sf * f + b.sh * h + b.sg * g + rnorm(n, 0, 0.5)
  l <- b.lg * g + b.ls * s + rnorm(n, 0, 0.5)
  df <- data.frame(p, f, h, g, s, l)
  
  out[j, ] <- vapply(models, function(form) {
    fit <- lm(form, data = df)
    co  <- coef(fit)
    if ("h" %in% names(co)) unname(co["h"]) else NA
  }, numeric(1))
}

#####Organize data for plot#####

#colors
pal <- brewer.pal(n = 8, "Set2")
col_valid  <- pal[5]  # purple-blue for valid
col_biased <- pal[2]  # orange for biased

#convert output to long format
out_long <-
  as.data.frame(out) |>
  mutate(sim = row_number()) |>
  pivot_longer(
    cols = starts_with("model"),
    names_to = "model",
    values_to = "beta_h"
  )

#table with panel headings
labels_tbl <- tibble(
  model = paste0("model", seq_along(models)),
  formula_disp = heading_models,
  intent_disp  = heading_accuracy,
  validity     = ifelse(grepl("^Valid", heading_accuracy), "Valid", "Biased")
)

labels_tbl$formula_disp <- c(
  "Diversity ~ Housing + Human",
  "Diversity ~ Housing + Food",
  "Diversity ~ Housing + Human + Green",
  "Diversity ~ Housing + Food + Green",
  "Diversity ~ Housing",
  "Diversity ~ Housing + Human + Food + Green + Protected"
)

levels_vec <- with(labels_tbl, paste0("Model: ", formula_disp, "\n", intent_disp))

tag_map <- setNames(paste0(LETTERS[seq_along(levels_vec)], ") "), levels_vec)

out_long <- out_long |>
  left_join(labels_tbl, by = "model") |>
  mutate(
    facet_lab = paste0("Model: ", formula_disp, "\n", intent_disp),
    facet_lab = factor(facet_lab, levels = levels_vec),
    validity  = factor(validity, levels = c("Valid", "Biased"))
  )

#axis limits
x_limits <- range(out_long$beta_h, na.rm = TRUE)

#colors for histograms
fill_vals <- c("Valid" = col_valid, "Biased" = col_biased)

# ---------------- ggplot figure ----------------
p <- ggplot(out_long, aes(x = beta_h, fill = validity)) +
  geom_histogram(bins = 40, color = NA) +
  geom_segment(aes(x = b.sh, xend = b.sh, y = 0, yend = 250),
               linetype = "dashed", linewidth = 0.7) +
  geom_segment(aes(x = b.tot.sh, xend = b.tot.sh, y = 0, yend = 250),
               linetype = "dotted", linewidth = 0.7) +
  annotate("text", x = b.sh,     y = 302, label = "True\n direct effect",
           vjust = 1.2, size = 3.3, lineheight = 0.9) +
  annotate("text", x = b.tot.sh, y = 302, label = "True\n total effect",
           vjust = 1.2, size = 3.3, lineheight = 0.9) +
  facet_wrap(
    ~ facet_lab, ncol = 2, scales = "fixed",
    labeller = labeller(facet_lab = function(l) {
      sapply(l, function(lbl) {
        tagged  <- paste0(tag_map[[lbl]], lbl)  # e.g., "A) Model: ...\nValid ..."
        indent  <- strrep(" ", max(nchar(tag_map[[lbl]]) + nchar("Model: ") - 5, 0))
        sub("\n", paste0("\n", indent), tagged, fixed = TRUE)
      })
    })
  ) +
  scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0.05, 0.02))) +
  scale_fill_manual(values = fill_vals, name = NULL) +
  coord_cartesian(ylim = c(0, 300)) +   
  labs(x = "Estimated housing density effect", y = "Simulations") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust=0), 
        legend.position = "none",
        legend.box = "horizontal")

print(p)
