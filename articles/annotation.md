# Running POTATO on Genomes

## Introduction

This vignette demonstrates the complete POTATO workflow:

1.  Loading genomes and potatoes
2.  Running annotation tools (kofam, BLAST, HMM)
3.  Scoring pathway presence/absence
4.  Analyzing and visualizing results

## Prerequisites

Assumes you have: - A configured sack (see “Setting Up a Potato
Environment” vignette) - Protein FASTA files (.faa) - Annotation
databases configured in `potato_config.yaml` - Conda environment with
bioinformatics tools

## Load Sack and Add Genomes

``` r

library(potato)

# Load existing sack
sack <- create_sack("my_potato_analysis")

# Add genomes
sack <- add_genomes(
  sack,
  genome_paths = "genomes/*.faa"
)

# Load potatoes
sack <- load_potatoes(
  sack,
  potato_dir = system.file("potatoes", package = "potato")
)

# Check setup
print(sack)
```

## Running Annotation Tools

POTATO uses three annotation approaches. You can run all three or select
specific tools.

### Kofam Annotation

Annotate proteins against KEGG Orthology (KO) database:

``` r

# Run kofam
sack <- run_kofam(
  sack,
  conda_env = "potato",
  workers = 8
)
```

**Parameters:** - `conda_env` - Name of conda environment with
kofamscan - `workers` - Number of parallel processes (default: 4)

**Output:** Results stored in `sack@results$kofam` with columns: -
`genome` - Genome name - `gene_name` - Protein ID - `ko` - KO
identifier - `score` - Kofam score - `threshold` - Per-gene threshold
from KEGG

### BLAST Annotation

Search against custom protein database:

``` r

# Run BLAST
sack <- run_blast(
  sack,
  conda_env = "potato",
  workers = 8
)
```

**Output:** Results in `sack@results$blast` with: - `genome`,
`gene_name` - Query info - `subject_id` - Hit identifier - `evalue`,
`bitscore` - Alignment metrics - `pident` - Percent identity

### HMM Annotation

Profile HMM searches (PFAM, TIGRfam, custom profiles):

``` r

# Run HMM
sack <- run_hmm(
  sack,
  conda_env = "potato",
  workers = 8
)
```

**Output:** Results in `sack@results$hmm` with: - `genome`,
`gene_name` - Query info - `target_name` - HMM profile name - `evalue`,
`score` - Search metrics - `tc_threshold` - Trusted cutoff (if
available)

### Run All Tools

``` r

# Sequential execution
sack <- run_kofam(sack, conda_env = "potato", workers = 8)
sack <- run_blast(sack, conda_env = "potato", workers = 8)
sack <- run_hmm(sack, conda_env = "potato", workers = 8)

# Save after annotation
saveRDS(sack, "my_potato_analysis/sack_annotated.rds")
```

## Scoring Pathways

Apply quality thresholds and calculate pathway completeness:

``` r

# Score pathways with default thresholds
sack <- score_pathways(sack)

# Or customize thresholds
sack <- score_pathways(
  sack,
  kofam_threshold = NULL,      # Use per-gene thresholds
  blast_evalue = 1e-10,
  blast_bitscore = 50,
  hmm_evalue = 1e-10
)
```

**Scoring logic:** - **Kofam:** Uses per-gene threshold from KEGG
(stored in results) - **HMM:** Uses trusted cutoff (TC) if available,
else e-value - **BLAST:** Global e-value and bitscore filters

**Output:** `sack@scores` tibble with columns: - `genome`, `potato`,
`pathway` - Identifiers - `total_steps_detected`, `total_steps` - Gene
counts (all genes) - `fraction` - Completeness fraction
(detected/total) - `present` - Pathway present/absent (passes
threshold) - `present_fraction`, `present_gaps` - Dual scoring results -
`gap_count`, `max_gaps` - Gap-based scoring metrics - `essential_*` -
Metrics for required genes only

### Dual Scoring

Pathways can use **fraction-based** and/or **gap-based** scoring:

**Fraction-based:** Requires ≥X% of genes detected

    fraction ≥ min_fraction → present

**Gap-based:** Can path from input→output with ≤N missing genes

    gap_count ≤ max_gaps → present

Pathway passes if **EITHER** method succeeds.

## Analyzing Results

### Export Results

``` r

# Gene-level results (all hits with thresholds)
gene_results <- get_gene_results(sack)
write.csv(gene_results, "gene_results.csv", row.names = FALSE)

# Pathway-level scores
pathway_scores <- get_pathway_scores(sack)
write.csv(pathway_scores, "pathway_scores.csv", row.names = FALSE)
```

### Summary Statistics

``` r

# Pathway prevalence across genomes
prevalence <- pathway_scores %>%
  group_by(potato, pathway) %>%
  summarize(
    n_genomes = n(),
    n_present = sum(present),
    prevalence = n_present / n_genomes
  )

# Find near-miss pathways (just below threshold)
near_misses <- find_near_miss_pathways(
  sack,
  max_distance = 0.15  # Within 15% of threshold
)
```

### Missing Gene Analysis

``` r

# Identify systematically missing genes
missing_genes <- summarize_missing_genes(sack)

# Shows:
# - Genes that never pass thresholds
# - Genes with hits below threshold
# - Genes with no hits at all
```

### Inspect Thresholds

``` r

# View annotation quality for specific potato
inspect_gene_thresholds(
  sack,
  potato_id = "glyoxylate_cycle",
  genome_name = "my_genome_001"
)
```

## Visualizing Results

### Pathway Heatmap

``` r

# Heatmap of pathway presence across genomes
plot_all_pathways_heatmap(sack)
```

### Pathway Prevalence

``` r

# Bar plot of pathway prevalence
plot_pathway_prevalence(sack)
```

### Annotation Coverage

``` r

# How many genes detected per genome
plot_annotation_coverage(sack)
```

### Interactive Pathway View

``` r

# Interactive network plot
potato <- sack@potatoes[[1]]
plot_potato(potato, interactive = TRUE)

# HTML detail view in browser
view_pathway_detail(potato, pathway = "main")
```

## Saving and Exporting

``` r

# Save complete sack (with all results)
saveRDS(sack, "my_potato_analysis/sack_complete.rds")

# Export provenance (reproducibility metadata)
print_provenance(sack)

# Get KEGG URL for pathway
kegg_url <- get_kegg_url(
  potato = sack@potatoes[[1]],
  pathway_id = "main",
  genes = TRUE,
  compounds = TRUE
)
browseURL(kegg_url)
```

## Troubleshooting

### Annotation failed for some genomes

Check the error logs in `results/annotations/*/`:

``` r

# Find failed genomes
failed <- sack@results %>%
  filter(is.na(kofam) | length(kofam) == 0)
```

### No pathways detected

Check thresholds are appropriate:

``` r

# View distribution of scores
gene_results %>%
  filter(tool == "kofam") %>%
  ggplot(aes(x = score, color = passed)) +
  geom_density()

# Try relaxing thresholds
sack <- score_pathways(sack, kofam_threshold = 50)
```

### Out of memory

Reduce parallel workers:

``` r

sack <- run_kofam(sack, workers = 2, conda_env = "potato")
```

Or process genomes in batches:

``` r

# Split genomes into batches
genomes1 <- sack@genomes[1:10]
genomes2 <- sack@genomes[11:20]

sack1 <- sack
sack1@genomes <- genomes1
sack1 <- run_kofam(sack1, workers = 4, conda_env = "potato")

sack2 <- sack
sack2@genomes <- genomes2
sack2 <- run_kofam(sack2, workers = 4, conda_env = "potato")

# Combine results
# (merge sack1@results and sack2@results)
```

## Next Steps

- **Create custom potatoes** - See vignette: “Creating Custom Potatoes”
- **Compare across studies** - Load multiple sacks and compare
- **Publication figures** - Use ggplot2 to customize visualizations
