# Setting Up a Potato Environment

## Introduction

This vignette walks through setting up a POTATO analysis environment,
including:

- Installing dependencies
- Configuring annotation tools
- Creating a new project
- Loading the project into a sack (PotatoSack S7 object)
- Adding genomes for annotation

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

## Creating a New Project

A POTATO project is an on-disk folder structure that holds configuration
files and input data. The in-memory “sack” (PotatoSack S7 object) is
created from this project structure.

### Initialize Directory Structure

``` r

# Create a new project
initialize_potato_project(
  path = "my_potato_analysis",
  copy_potatoes = TRUE  # Copies bundled potatoes from package
)
```

This creates:

    my_potato_analysis/
    ├── potato_config.yaml      # Configuration file
    ├── my_potato_analysis.Rproj # RStudio project file
    ├── potatoes/                # Pathway definitions (JSON files)
    ├── genomes/                 # Input genome files
    ├── results/                 # Output files
    │   ├── annotations/        # Gene-level annotation results
    │   └── scores/             # Pathway-level scoring results
    ├── .claude/                 # Claude Code agents for building potatoes
    ├── .gitignore
    └── README.md

### Configure Annotation Databases

Edit `potato_config.yaml` to point to your local databases:

``` yaml
databases:
  kofam:
    type: kofam
    profiles_dir: /path/to/kofam/profiles/
    ko_list: /path/to/kofam/ko_list
    executable: exec_annotation
    thresholds:
      score_ratio: 1.0

  blast:
    type: blast
    path: /path/to/blast_db.faa
    executable: blastp
    thresholds:
      evalue: 1e-7
      bitscore: 50

  hmm:
    type: hmm
    path: /path/to/profiles.hmm
    executable: hmmsearch
    thresholds:
      evalue: 1e-10
```

**Database locations:**

- **KOFAM**: Download from [KEGG
  FTP](https://www.genome.jp/ftp/db/kofam/)
- **BLAST**: Custom databases or use NCBI nr/refseq
- **HMM**: PFAM, TIGRfam, or custom profiles

### Create Sack Object

Load the initialized project into an R PotatoSack object:

``` r

# From within the project directory
sack <- create_sack()

# Or specify path
sack <- create_sack("my_potato_analysis")
```

**What
[`create_sack()`](https://jeffkimbrel.github.io/potato/reference/create_sack.md)
does:** - Loads `potato_config.yaml` and validates database
configuration - Loads all potatoes from the `potatoes/` directory -
Creates an empty PotatoSack S7 object ready for genome annotation

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

### Using a Template Configuration

If you already have a working project with configured databases, you can
reuse that configuration:

``` r

# Create new project using existing config as template
initialize_potato_project(
  path = "new_project",
  template_config = "existing_project/potato_config.yaml"
)
```

This copies the `databases:`, `annotation:`, and `scoring:` sections
from your template, so you don’t need to reconfigure database paths.

## Adding Genomes

Add protein FASTA files (.faa or .fasta) to your sack:

``` r

# Add genomes from directory (glob pattern)
sack <- add_genomes(sack, "path/to/genomes/*.faa")

# Add specific files
sack <- add_genomes(sack, c("genome1.faa", "genome2.faa"))

# Check added genomes
sack@genomes
```

**IMPORTANT**: - POTATO requires **protein** FASTA files (.faa or
.fasta), not nucleotide sequences - Genomes are registered by path (not
copied), so keep them in place - Use
[`add_genomes()`](https://jeffkimbrel.github.io/potato/reference/add_genomes.md)
multiple times to add more genomes incrementally

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
exist. The config is validated automatically when you call
[`create_sack()`](https://jeffkimbrel.github.io/potato/reference/create_sack.md),
so you’ll see errors immediately if paths are wrong.

### Memory issues with large genome sets

Process genomes in batches:

``` r

# Run annotation with fewer parallel workers
sack <- run_kofam(sack, workers = 4, conda_env = "potato")
```
