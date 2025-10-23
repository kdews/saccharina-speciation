# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(ggpubr)
if (require(showtext)) {
  showtext_auto()
  if (interactive())
    showtext_opts(dpi = 100)
  else
    showtext_opts(dpi = 300)
}
# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/scratch1/kdeweese/latissima_scratch/saccharina_speciation/"
  setwd(wd)
  # Input VCF
  vcf_file <- "master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.filt_QUAL30_bial.vcf.gz"
  # Filter name
  filter_name <- "MEAN_DP"
  # Output directory
  outdir <- "saccharina-speciation/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  vcf_file <- line_args[1]
  filter_name <- line_args[2]
  outdir <- line_args[3]
}
# Set number of cores for multi-threading
my_cores <- system("echo $SLURM_CPUS_PER_TASK", intern = T)
if (my_cores != "") {
  my_cores <- as.integer(my_cores)
} else {
  my_cores <- 1
}
# Parse input
vcf_base <- gsub("\\.vcf.*", "", (basename(vcf_file)))
vcf_base_prefix <- gsub("\\..*", "", vcf_base)
vcf_base_suffix <- gsub(paste0(vcf_base_prefix, "."), "", vcf_base)
depth_dist_file <- paste0(vcf_base, ".ldepth.mean")
bcf_depth_file <- paste0(vcf_base, ".per_sample.stats")
# Output
filt_depth_base <- paste(
  vcf_base_prefix,
  filter_name,
  vcf_base_suffix,
  sep = "."
)
filt_depth_tab_file <- paste0(filt_depth_base, ".tsv")
filt_depth_vcf_file <- paste0(filt_depth_base, ".vcf.gz")
filt_depth_plot_file <- paste0("depth_dist.", vcf_base_suffix, ".png")

# Create list of sites to filter out with extreme low/high coverage (DP)
# Read VCFtools --site-mean-depth file
depth_dist <- read_tsv(depth_dist_file, num_threads = my_cores)
# Remove outliers (>100x)
rm_outliers_depth_dist <- depth_dist %>%
  filter(MEAN_DEPTH <= 100)
# Filter mean depth distribution with bias toward higher coverage
min_cutoff <- round(quantile(rm_outliers_depth_dist$MEAN_DEPTH, probs = 0.1))
max_cutoff <- round(quantile(rm_outliers_depth_dist$MEAN_DEPTH, probs = 0.99985))
# Plot distribution of read depth (log10-scaled) to select DP cutoffs
p1 <- ggplot(rm_outliers_depth_dist, aes(x = MEAN_DEPTH)) +
  geom_histogram() +
  geom_vline(xintercept = min_cutoff, color = "red", linetype = "dashed") +
  geom_vline(xintercept = max_cutoff, color = "red", linetype = "dashed") +
  annotate(
    geom = "text",
    label = min_cutoff,
    x = min_cutoff,
    y = (dim(rm_outliers_depth_dist)[1])/4,
    hjust = 1.5,
    color = "red",
    size = rel(7)
  ) +
  annotate(
    geom = "text",
    label = max_cutoff,
    x = max_cutoff,
    y = (dim(rm_outliers_depth_dist)[1])/4,
    hjust = 1.5,
    color = "red",
    size = rel(7)
  ) +
  labs(
    x = "Mean site read depth",
    y = "Sites"
  ) +
  scale_x_log10()

# Read BCFtools stats depth table
# Parse for only lines containing depth information
all_lines <- readLines(bcf_depth_file)
start_line_index <- grep("# DP, Depth distribution", all_lines)[1]
# BCFtools standard depth table is 500 bins + >500 bin = 501 lines
bcf_depth <- read_tsv(bcf_depth_file, skip = start_line_index, n_max = 501)
colnames(bcf_depth) <- gsub("# ", "", colnames(bcf_depth))
colnames(bcf_depth) <- gsub("\\[.*\\]", "", colnames(bcf_depth))
# Convert DP bin to numeric, handling ">500"
bcf_depth <- bcf_depth %>%
  mutate(
    bin_numeric = as.numeric(ifelse(grepl("^>", bin), sub(">", "", bin), bin)),
    bin_numeric = ifelse(is.na(bin_numeric), max(bin_numeric, na.rm = TRUE), bin_numeric)
  )
# Plot number of genotypes vs. depth
p2 <- ggplot(
  bcf_depth, 
  aes(x = bin_numeric, y = `number of genotypes` / `number of sites`)
) +
  geom_line(color = "steelblue") +
  geom_point(size = 0.5, color = "steelblue") +
  scale_y_log10() +
  labs(
    x = "Binned read depth",
    y = "Genotypes"
  ) +
  geom_vline(xintercept = max_cutoff, color = "red", linetype = "dashed") +
  annotate(
    geom = "text",
    label = max_cutoff,
    x = max_cutoff,
    y = median(bcf_depth$`number of genotypes`),
    hjust = -1,
    color = "red",
    size = rel(7)
  )
p <- ggarrange(p1, p2)
# Write plots to output directory
ggsave(filt_depth_plot_file, p, path = outdir, height = 7, width = 10)

# Create table of sites to keep that pass filters for site mean depth
filt_depth_tab <- depth_dist %>%
  filter(MEAN_DEPTH >= min_cutoff) %>%
  filter(MEAN_DEPTH <= max_cutoff) %>%
  select(CHROM, POS)
write_tsv(
  x = filt_depth_tab,
  file = filt_depth_tab_file,
  col_names = F,
  num_threads = my_cores
)

# Report variant stats and filter parameters
msg1 <- paste("Filter:", min_cutoff, "< SITE_MEAN_DP <", max_cutoff)
msg2 <- paste("Before filter:", dim(depth_dist)[1], "sites")
msg3 <- paste("After filter:", dim(filt_depth_tab)[1], "sites")
cat(msg1, msg2, msg3, sep = "\n")
cat("\n")

