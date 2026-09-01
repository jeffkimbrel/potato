# Setting Up a Potato Environment

## Introduction

This vignette walks through setting up a POTATO analysis environment,
including:

- Installing dependencies
- Configuring annotation tools
- Creating a new sack
- Loading pathway definitions (potatoes)

## Prerequisites

### R Package Installation

Install POTATO from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("jeffkimbrel/potato")
library(potato)
```

### Conda Environment

POTATO uses bioinformatics tools that are best installed via Conda.
Create the environment from the included `environment.yaml`:

``` bash
# From the POTATO package directory
conda env create -f environment.yaml

# Activate the environment
conda activate potato
```

The environment includes:

- **kofamscan** - KEGG Orthology annotation
- **HMMER3** - Profile HMM searches
- **BLAST+** - Sequence similarity searches

## Creating a New Sack

A “sack” is POTATO’s project container that holds:

- Configuration (database paths, tool settings)
- Loaded pathway definitions (potatoes)
- Genome files
- Annotation results
- Pathway scores

### Initialize Directory Structure

``` r

# Create a new sack in the current directory
initialize_potato_sack(
  sack_root = "my_potato_analysis",
  sack_id = "marine_mags_2024"
)
```

This creates:

    my_potato_analysis/
    ├── config/
    │   └── potato_config.yaml
    ├── potatoes/
    ├── genomes/
    └── results/

### Configure Annotation Databases

Edit `config/potato_config.yaml` to point to your local databases:

``` yaml
databases:
  kofam:
    type: kofam
    profiles: /path/to/kofam/profiles
    ko_list: /path/to/kofam/ko_list

  blast:
    type: blast
    database: /path/to/blast/db

  hmm:
    type: hmm
    profiles: /path/to/hmm/profiles
```

**Database locations:**

- **KOFAM**: Download from [KEGG
  FTP](https://www.genome.jp/ftp/db/kofam/)
- **BLAST**: Custom databases or use NCBI nr/refseq
- **HMM**: PFAM, TIGRfam, or custom profiles

### Create Sack Object

Load the initialized sack into R:

``` r

# From within the sack directory
sack <- create_sack()

# Or specify path
sack <- create_sack("path/to/my_potato_analysis")
```

## Loading Potatoes

Potatoes are pathway definitions in JSON format. POTATO includes several
built-in pathways.

### Built-in Pathways

``` r

# List available built-in potatoes
list.files(system.file("potatoes", package = "potato"),
           pattern = "\\.json$")
```

### Load Specific Potatoes

``` r

# Load from built-in potatoes
sack <- load_potatoes(
  sack,
  potato_dir = system.file("potatoes", package = "potato"),
  pattern = "glyoxylate|nitrogen"  # Filter by name
)

# Or load all built-in potatoes
sack <- load_potatoes(
  sack,
  potato_dir = system.file("potatoes", package = "potato")
)
```

### Load Custom Potatoes

``` r

# Copy custom potatoes to your sack's potatoes/ directory
# Then load them
sack <- load_potatoes(
  sack,
  potato_dir = "my_potato_analysis/potatoes"
)
```

### Inspect Loaded Potatoes

``` r

# View sack contents
print(sack)

# List loaded potato IDs
names(sack@potatoes)

# View a specific potato
potato <- sack@potatoes[[1]]
print(potato)

# Get verification status
get_pathway_verification(sack)
```

## Adding Genomes

Add protein FASTA files (.faa) to your sack:

``` r

# Add genomes from directory
sack <- add_genomes(
  sack,
  genome_paths = "path/to/genomes/*.faa"
)

# Check added genomes
sack@genomes
```

**IMPORTANT**: POTATO requires **protein** FASTA files (.faa), not
nucleotide sequences.

## Saving Your Sack

The sack object can be saved for later:

``` r

# Save to RDS
saveRDS(sack, "my_potato_analysis/sack.rds")

# Reload later
sack <- readRDS("my_potato_analysis/sack.rds")
```

## Next Steps

Now that your environment is set up, you can:

1.  **Run annotation** - See vignette: “Running POTATO on Genomes”
2.  **Create custom potatoes** - See vignette: “Creating Custom
    Potatoes”
3.  **Analyze results** - Explore pathway presence/absence patterns

## Troubleshooting

### Conda environment not found

``` r

# Specify conda environment explicitly
sack <- run_kofam(sack, conda_env = "potato")
```

### Database paths not recognized

Check that paths in `potato_config.yaml` are absolute and the files
exist:

``` r

# Validate config
validate_database_configs(sack@config)
```

### Memory issues with large genome sets

Process genomes in batches:

``` r

# Run annotation with fewer parallel workers
sack <- run_kofam(sack, workers = 4, conda_env = "potato")
```
