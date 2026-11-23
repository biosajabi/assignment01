###1 - OVERVIEW----

#####Question----
#How does the transposable element differ between Trichinella spiralis and Trichinella murrelli in terms of overall repeat abundance, major element class composition, and element length distribution?

####Data Access----
# This code file should run on everyones system assuming they have downloaded the required packages and are in the latest software version. 

# Please make sure to download datasets from this NCBI link for T. Spirallis and T. Murrelli and save with the appropriate name and data/ directory:
#T. Spirallis: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001447595.1/ 
#T. Murrelli: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001447425.1/ 

####Installation---- 
#Package installation & load libraries + Data exploration for both species (T. spiralis & T. murrelli)

library(tidyverse)
library(patchwork)


###2 - Edit 1 - Data Exploration ----
#Changes:
#1. Replaced repetitive data exploration and cleaning code with a for loop to improve efficiency and make it easier to add more species.
#2. Added column validation using setdiff() and stop() to catch missing columns early and prevent errors in downstream analysis (ie. if file is read into script and does not contain the required columns, the loop will terminate)

files <- c(
  "T. spiralis" = "../data/trichinella_spiralis_te.csv",
  "T. murrelli" = "../data/trichinella murrelli_te.csv"
)

# Columns required for analysis
required_cols <- c("qry_id", "qry_start", "qry_end",
                   "qry_width", "qry_left",
                   "matching_class", "perc_div", "sw_score")

# Empty list to store cleaned datasets
ls_data <- list()

for (i in names(files)) {
  
  # Read dataset
  df <- read_csv(files[[i]])
  
  # Check for missing columns so loop terminates if a dataset lacks the info needed for analysis
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in", i, ":",
               paste(missing_cols, collapse = ", ")))
  }
  # Calculate contig length by standardizing qry_left to a positive numeric value before calculating contig length for both species.
  df <- df |>
    mutate(
      qry_left = abs(as.numeric(gsub("[()]", "", qry_left))),
      contig_length = qry_end + qry_left
    )
  # Exploration
  print(class(df))
  print(names(df))
  print(summary(df))
  print(dim(df))
  
  # Select columns and add species label
  ls_data[[i]] <- df |>
    select(all_of(required_cols),(contig_length)) |>
    mutate(species = i)
}


#####Combination----

df_tstm <- bind_rows(ls_data)
unique(df_tstm$species)
dim(df_tstm)

###3 (OG)- Broad Classification----
# Group detailed repeat annotations into broad, interpretable classes using str_detect() from stringr.

unique(df_tstm$matching_class)
#Question - Would grouping repeats into broad classes hide potentially important differences between TE families or subfamilies?
df_tstm <- df_tstm |>
  mutate(
    broad_class = case_when(
      str_detect(matching_class, "DNA") ~ "DNA Transposon",
      str_detect(matching_class, "LINE") ~ "LINE",
      str_detect(matching_class, "SINE") ~ "SINE",
      str_detect(matching_class, "LTR") ~ "LTR Retrotransposon",
      str_detect(matching_class, "Simple_repeat") ~ "Simple Repeat",
      str_detect(matching_class, "simple_repeat") ~ "Simple Repeat",
      str_detect(matching_class, "Satellite") ~ "Satellite",
      str_detect(matching_class, "PLE") ~ "PLE",
      str_detect(matching_class, "Helitron") ~ "RC/Helitron",
      str_detect(matching_class, "Low_complexity") ~ "Low Complexity",
      str_detect(matching_class, "Retroposon") ~ "Other Retrotransposon",
      TRUE ~ "Other/Unknown"))

#check unique groups in broad_class
unique(df_tstm$broad_class)

####Non-TE Filtering----

# Filter out non-TEs (decided to keep satellites as a good non transposable element comparison cohort for discussion)

combined_te <- df_tstm |>
  filter(!broad_class %in% c("Simple Repeat", "Low Complexity", "Other/Unknown")) |>
  mutate(broad_class = fct_drop(broad_class)) |>
  mutate(broad_class = fct_infreq(broad_class))

df_tstmclean <- combined_te |>
  mutate(df_telength = qry_end - qry_start +1)

df_tstmclean |>
  select(qry_start, qry_end, df_telength) |>
  head()

###4 - SUMMARY STATS----

# summary statistics for results slides

####Metric 1----

#Overall repeat statistics per species

overall_te_summary <- combined_te |>
  group_by(species) |>
  summarise(
    total_te_count = n(),
    median_length_bp = median(qry_width),
    mean_length_bp = mean(qry_width),
    total_te_bp = sum(qry_width)
  )

cat("\n--- Transposable Element Summary per Species ---\n")
print(overall_te_summary)

####Metric 2----

#Abundance and percentage per repeat class

te_class_summary <- combined_te |>
  count(species, broad_class) |>
  group_by(species) |>
  mutate(percentage = (n / sum(n)) * 100) |>
  ungroup()

cat("\n--- Transposable Element Class Proportions ---\n")
print(te_class_summary)

###5 - DATA VISUALIZATION----

####QrylengthXspecies----
#figure 1: Element Length Distribution (EDA / Quality Check)

fig1 <- ggplot(combined_te, aes(x = species, y = qry_width, fill = species)) + geom_boxplot(outlier.alpha = 0.15) +
  scale_y_log10() + 
  scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) + labs(
    title = "Comparison of Transposable Element Lengths",
    x = "Species",
    y = "TE Length in bp (log10 scale)"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

#print to Plots viewer and save directly to output/
print(fig1)
ggsave("../output/figure1_length_distribution.png", plot = fig1, width = 6, height = 4.5)

####Class composition----
#figure 2: Major Repeat Class Composition (Main Comparative Result)

fig2 <- ggplot(te_class_summary, aes(x = broad_class, y = percentage, fill = species)) + geom_col(position = position_dodge(width = 0.8), width = 0.7) + scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) + labs(
  title = "Proportion of Transposable Element Classes",
  x = "Transposable Element Class",
  y = "Percentage of Total TEs (%)",
  fill = "Species"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

# Original Figure 2
print(fig2)
ggsave("../output/figure2_class_composition.png",
       plot = fig2, width = 7, height = 5)

####Contig repeat density----

#step 8: figure 3: Contig-Level Repeat Density & Wilcoxon Test (Novelty Component)
#count number of repeats per contig using qry_id

contig_te_counts <- combined_te |>
  count(species, qry_id) #longer contigs will naturally have more TE annotations, may affect raw repeat counts
#perform Wilcoxon rank-sum test

wilcox_te_res <- wilcox.test(n ~ species, data = contig_te_counts)
p_val_label <- paste0("Wilcoxon rank-sum p = ", format.pval(wilcox_te_res$p.value, digits = 3))

fig3 <- ggplot(contig_te_counts, aes(x = species, y = n, fill = species)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
  scale_y_log10() +
  scale_fill_manual(values = c("T. spiralis" = "#2b8cbe", "T. murrelli" = "#e6550d")) +
  annotate("text", x = 1.5, y = max(contig_te_counts$n) * 0.8, label = p_val_label, fontface = "italic") +
  labs(
    title = "Repeat Annotation Counts per Contig",
    subtitle = "Non-normalized counts across assembly contigs",
    x = "Species",
    y = "Repeat Annotations per Contig (log10 scale)"
  ) + 
  theme_classic() +
  theme(legend.position = "none")

print(fig3)
ggsave("../output/figure3_contig_density.png", plot = fig3, width = 6, height = 4.5)

cat("\n--- Statistical Test on TE Contig Density ---\n")
print(wilcox_te_res)

###

###6 - Edit 2 - Length-Normalized TE Density----
#Normalized repeat counts by contig length to make fairer comparisons between different sized contigs.
# Get one length per contig
contig_lengths <- df_tstm |>
  distinct(species, qry_id, contig_length)

stopifnot(
  nrow(contig_lengths) ==
    nrow(distinct(contig_lengths, species, qry_id))
) #checks that each contig has one length, makes life easier when combining data.

# Calculate annotations per kilobase
contig_density <- contig_te_counts |>
  left_join(contig_lengths, by = c("species", "qry_id")) |>
  mutate(te_per_kb = n / (contig_length / 1000))

# Validate calculated densities
stopifnot(
  all(is.finite(contig_density$te_per_kb)),
  all(contig_density$contig_length > 0)
)

# Compare species
density_summary <- contig_density |>
  group_by(species) |>
  summarise(
    n_contigs = n(),
    median_density = median(te_per_kb),
    mean_density = mean(te_per_kb),
    .groups = "drop"
  )

print(density_summary)

# Updated Figure 3
fig3_updated <- ggplot(
  contig_density,
  aes(x = species, y = te_per_kb, fill = species)
) +
  geom_boxplot(outlier.alpha = 0.15) +
  scale_fill_manual(values = c(
    "T. spiralis" = "#2b8cbe",
    "T. murrelli" = "#e6550d"
  )) +
  labs(
    title = "Length-Normalized TE Density by Species",
    x = "Species",
    y = "TE Annotations per kb"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(fig3_updated)

ggsave("../output/figure3_updated.png",
       plot = fig3_updated, width = 7, height = 5)

# Statistical comparison
density_test <- wilcox.test(
  te_per_kb ~ species,
  data = contig_density,
  exact = FALSE
)

print(density_test)


###7 - Edit 3.1 -  Divergence Summary----
#I wanted to see what the full distribution and class-specific comparison would look like if the 5% divergence split was removed.

# Summarize all divergence values without a 5% cutoff
div_summary <- combined_te |>
  filter(!is.na(matching_class), !is.na(perc_div)) |>
  group_by(species, matching_class) |>
  summarise(
    n = n(),
    mean_div = mean(perc_div),
    median_div = median(perc_div),
    sd_div = sd(perc_div),
    .groups = "drop"
  )

print(div_summary)


# Visualize the full distribution of divergence values
fig_div_updated <- combined_te |>
  filter(!is.na(matching_class), !is.na(perc_div)) |>
  ggplot(aes(x = matching_class,
             y = perc_div,
             fill = species)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(
    "T. spiralis" = "#2b8cbe",
    "T. murrelli" = "#e6550d")) +
  labs(
    title = "TE Divergence by Repeat Class",
    x = "Repeat Class",
    y = "Percent Divergence (%)",
    fill = "Species"
  ) +
  theme_classic() +
  coord_flip()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(fig_div_updated)

ggsave("../output/figure_divergence_updated.png",
       plot = fig_div_updated,
       width = 11, height = 10)

##

###8 - Edit 3.2 - Statistical Comparison----

# Keep repeat classes with at least 10 observations from each species (practical cut off, not a guarantee of statistical power)
div_counts <- combined_te |>
  filter(!is.na(matching_class), !is.na(perc_div)) |>
  count(matching_class, species) |>
  group_by(matching_class) |>
  filter(n_distinct(species) == 2, all(n >= 10)) |>
  ungroup()

# Keep only eligible repeat classes
div_test_data <- combined_te |>
  filter(!is.na(perc_div)) |>
  semi_join(div_counts, by = c("matching_class", "species"))

# Run a Wilcoxon test for each repeat class and then correct to account for multiple comparisons.
div_tests <- div_test_data |>
  group_by(matching_class) |>
  summarise(
    n_spiralis = sum(species == "T. spiralis"),
    n_murrelli = sum(species == "T. murrelli"),
    median_spiralis = median(perc_div[species == "T. spiralis"]),
    median_murrelli = median(perc_div[species == "T. murrelli"]),
    p_value = wilcox.test(
      perc_div ~ species,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) |>
  mutate(
    median_difference = median_spiralis - median_murrelli,
    p_adjusted = p.adjust(p_value, method = "BH")
  ) |>
  arrange(p_adjusted)

div_tests |>
  select(matching_class, median_difference, p_value, p_adjusted) |>
  print(n = Inf)
#NOTE: Significant results suggest differences in repeat divergence distributions, but do not establish differences in TE age or activity. Non-significant adjusted p-values indicate insufficient evidence of a difference, not that the distributions are identical.

###9 - Multi-Variable Bubble Plots----
# what if we look at sequence length x percent divergence? 

mv_combined <- ggplot(data = df_tstmclean, aes(x = perc_div, y  = df_telength, size = sw_score, fill = broad_class, shape = species)) +
  geom_point(alpha=0.5, color = "black") +
  scale_shape_manual(values = c(
    "T. spiralis" = 21,
    "T. murrelli" = 24)) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21, 
        color = "black", 
        alpha = 1
      )
    )
  ) +
  labs(
    fill = "Repeat Class",
    shape = "Species",
    size = "SW_score") +
  labs(
    title = "Multivariable Comparison between\n T. spiralis vs. T.murelli",
    x = "Percent Divergence (%)",
    y = "Sequence Length (bp)",
  ) +
  theme_minimal()

mv_edited <- ggplot(data = subset(df_tstmclean, broad_class != "LTR Retrotransposon"), aes(x = perc_div, y  = df_telength, size = sw_score, fill = broad_class, shape = species)) +
  geom_point(alpha=0.5, color = "black") +
  scale_shape_manual(values = c(
    "T. spiralis" = 21,
    "T. murrelli" = 24)) +
  
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21, 
        color = "black", 
        alpha = 1
      )
    )
  ) +
  labs(
    fill = "Repeat Class",
    shape = "Species",
    size = "SW_score") +
  labs(
    title = "Multivariable Comparison between\n T. spiralis vs. T.murelli (without LTR)",
    x = "Percent Divergence (%)",
    y = "Sequence Length (bp)",
  ) +
  theme_minimal()

Bubble_plot <- mv_combined + mv_edited +
  plot_layout(ncol = 2)
Bubble_plot

ggsave(
  filename = "../output/fig_45_Comp.png",
  plot = Bubble_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)


