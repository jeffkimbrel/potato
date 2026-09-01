# Potato Sack S7 Class

An S7 class representing a complete potato annotation project, including
configuration, potatoes, genomes, annotation results, and provenance
tracking.

## Usage

``` r
PotatoSack(
  sack_id = character(0),
  sack_root = character(0),
  config = NULL,
  potatoes = list(),
  genomes = list(),
  results = NULL,
  scores = NULL,
  completed_stages = character(0),
  provenance = list(),
  metadata = list()
)
```

## Fields

- `sack_id`:

  Unique identifier for this sack (hash or user-provided)

- `sack_root`:

  Path to the sack directory

- `config`:

  Loaded potato configuration

- `potatoes`:

  List of Potato objects

- `genomes`:

  List of GenomeFile objects

- `results`:

  Annotation results table (tibble)

- `scores`:

  Pathway scoring results (tibble)

- `completed_stages`:

  Vector of completed annotation stages

- `provenance`:

  Provenance tracking for each stage

- `metadata`:

  Additional metadata
