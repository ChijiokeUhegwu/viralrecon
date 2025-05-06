#Name: Chijioke Uhegwu
#Project: Viral Recon analysis
#Date: 24-04-2025

# Loading packages ----
if(!require(pacman)) install.packages("pacman")
  pacman::p_load(tidyverse, here, RcolorBrewer, plotly)

# importing the data  ----
var <- read_csv(here("data/variants_long_table.csv"))
head (var)  #to display the first 6 rows of the dataset

#visualization of the data ----
#The distribution of DP values per sample
DP_per_sample <- ggplot(data = var_tb) + 
  aes(x=SAMPLE, y=DP, fill=SAMPLE) + 
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_brewer(palette="RdYlBu") +
  theme(legend.position="right") +
  ggtitle("DP per Sample")
DP_per_sample
ggsave(DP_per_sample, filename = here("outputs/DP per sample.jpg"), width=10, height=8)

#The distribution of DP values per chromosome and per sample
DP_per_chromosome <- ggplot(data = var_tb) + 
  aes(x=CHROM, y=DP, fill= SAMPLE) + 
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_brewer(palette="RdYlBu") + 
  labs(title="DP_per_Chromosome") + facet_grid(. ~ SAMPLE)
DP_per_chromosome
ggsave(DP_per_chromosome, filename=here("outputs/DP per chromosome.jpg"), width=10, height=8)

#Variants effects per sample
p_effect <- ggplot(data = var_tb) +
  aes(y=EFFECT, fill= SAMPLE) +
  scale_fill_brewer(palette="Spectral") + 
  labs(title="Effect_per_Sample") + 
  theme(legend.position="bottom") +
  geom_bar()
p_effect
ggsave(p_effect, filename=here("outputs/variants effects per sample.jpg"), width=10, height=8)

#Variants per Gene (Tells you which genes are most frequently mutated)
var_per_gene <- var_tb %>%
  count(GENE, sort = TRUE) %>%
  ggplot(aes(x = reorder(GENE, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Variants per Gene",
       x = "Gene",
       y = "Count of Variants") +
  theme_minimal()
var_per_gene
ggsave(var_per_gene, filename=here("outputs/variants per gene.jpg"), width=10, height=8)

# Read depth per position
p_DP_POS <- ggplot(data = var_tb) +
  aes(x=POS, y=DP, fill= SAMPLE) + 
  scale_fill_brewer(palette="RdBu") + 
  labs(title="DP_per_Position") + 
  theme(legend.position="bottom") + 
  geom_point(shape = 21, size = 5, alpha = 0.7)
p_DP_POS
ggsave(p_DP_POS, filename=here("outputs/DP per position.jpg"), width=10, height=8)

# Mutation effect types (Overview of functional consequences)
mut_effects <- var_tb %>%
  count(EFFECT) %>%
  ggplot(aes(x = "", y = n, fill = EFFECT)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Mutation Effects") +
  theme_void() +
  theme(legend.position = "right")
mut_effects
ggsave(mut_effects, filename=here("outputs/proportion of mutation effects.png"), width=10, height=8)

# Heatmap of Variant Positions by Sample
#Visualizes where mutations occur in the genome across all samples.
pos_matrix <- var_tb %>%
  select(SAMPLE, POS) %>%
  mutate(POS = as.numeric(POS)) %>%
  mutate(POS_bin = cut(POS, breaks = seq(0, 30000, by = 1000))) %>%
  count(SAMPLE, POS_bin)

pos_matrix_heatmap <- ggplot(pos_matrix, aes(x = POS_bin, y = SAMPLE, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  labs(title = "Heatmap of Mutation Density Across Genome by Sample",
       x = "Genome Position (binned)",
       y = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pos_matrix_heatmap
ggsave(pos_matrix_heatmap, filename=here("outputs/heatmap of variant positions by sample.jpg"), width=10, height=8)

#Amino Acid Changes – Frequency of Notable Mutations
amino_changes <- var_tb %>%
  filter(!is.na(HGVS_P_1LETTER)) %>%
  count(HGVS_P_1LETTER, sort = TRUE) %>%
  top_n(10) %>%
  ggplot(aes(x = reorder(HGVS_P_1LETTER, n), y = n)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(title = "Top Amino Acid Changes",
       x = "Amino Acid Change",
       y = "Count") +
  theme_minimal()
amino_changes
ggsave(amino_changes, filename=here("outputs/amino acid changes.png"), width=10, height=8)

#Allele Frequency Distribution
#Determine how confident your variant calls are (AF close to 1 are high confidence).
allele_freq <- var_tb %>%
  mutate(AF = as.numeric(AF)) %>%
  ggplot(aes(x = AF)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
  labs(title = "Allele Frequency Distribution of Variants",
       x = "Allele Frequency (AF)",
       y = "Number of Variants") +
  theme_minimal()
allele_freq
ggsave(allele_freq, filename=here("outputs/allele frequency distribution.jpg"), width=10, height=8)

#Variants per Sample (Compare mutation loads across different samples)
var_per_sample <- var_tb %>%
  count(SAMPLE) %>%
  ggplot(aes(x = reorder(SAMPLE, n), y = n)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  coord_flip() +
  labs(title = "Variant Count per Sample",
       x = "Sample",
       y = "Number of Variants") +
  theme_minimal()
var_per_sample
ggsave(var_per_sample, filename=here("outputs/variant count per sample.jpg"), width=10, height=8)

#plotly for interactivity 
var_plotly <- ggplot(var_tb %>% count(EFFECT), aes(x = EFFECT, y = n, fill = EFFECT)) +
  geom_col()
ggplotly(var_plotly)
