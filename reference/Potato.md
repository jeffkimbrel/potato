# Potato S7 Class

Represents a pathway definition with genes and pathway structure

## Usage

``` r
Potato(
  id = character(0),
  name = character(0),
  genes = list(),
  edges = list(),
  tags = character(0),
  source = character(0),
  notes = character(0),
  scoring = list(),
  input = list(),
  output = list(),
  json_path = character(0),
  graph = NULL,
  compound_coordinates = list()
)
```

## Arguments

- id:

  Character. Unique identifier for the potato

- name:

  Character. Human-readable name

- genes:

  List. Genes (enzymes) in the pathway. Each gene can have:

  - `marker`: Logical. If TRUE, this gene is a diagnostic marker for the
    pathway

  - `required`: Logical. If TRUE, pathway cannot be complete without
    this gene

  - `type`: Character. "enzyme", "compound", or "transporter"

  - `databases`: List. Detection methods by database type (e.g., kofam,
    blast, hmm)

- edges:

  List. Edges connecting genes or pathways object for multi-pathway
  networks

- tags:

  Character vector. Tags for organizing potatoes

- source:

  Character. Source/origin of the pathway definition

- notes:

  Character. Additional notes about the pathway

- scoring:

  List. Scoring parameters including:

  - `min_fraction`: Minimum fraction of genes required (default 0.75)

  - `marker_mode`: "any" or "all" - how many marker genes needed for
    positive call

- json_path:

  Character. Path to the source JSON file

- graph:

  igraph object or NULL. Cached graph representation
