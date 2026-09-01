# Creating Custom Potatoes

## Introduction

This vignette covers:

- Creating custom pathway definitions (potatoes)
- V2 schema structure and design principles
- Building potatoes from KEGG modules
- Validating and testing potatoes
- Best practices for pathway definition

## Why Create Custom Potatoes?

Built-in potatoes cover common pathways, but you may need custom
definitions for:

- Novel or recently characterized pathways
- Organism-specific metabolic variants
- Custom gene sets or biomarkers
- Pathways not in KEGG
- Alternative gene annotations

## V2 Schema Overview

POTATO uses a **network schema** where related pathways share genes and
metabolic context.

### Key Components

**1. Genes** - Detection methods and enzyme information

``` json
{
  "id": "aceA",
  "name": "isocitrate lyase",
  "databases": {
    "kofam": ["K01637"],
    "hmm": ["PF00463"]
  },
  "ec": ["4.1.3.1"],
  "type": "enzyme"
}
```

**2. Compounds** - Metabolites that flow through pathway

``` json
{
  "id": "C00026",
  "name": "2-Oxoglutarate",
  "kegg_compound": "C00026"
}
```

**3. Pathways** - Network topology and scoring

``` json
"pathways": {
  "main": {
    "name": "Glyoxylate Cycle",
    "edges": [
      {"from": "C00026", "to": "aceA", "required": true, "marker": true}
    ],
    "scoring": {
      "min_fraction": 0.75,
      "max_gaps": 1
    }
  }
}
```

## Creating a Potato from KEGG

The easiest way to create a potato is from a KEGG module.

### Using the build-potato Skill

``` r

# Use the interactive agent
# In Claude Code, type:
# /build-potato
```

This launches an interactive agent that: 1. Fetches KEGG module data 2.
Converts to V2 schema 3. Adds detection methods 4. Generates validated
JSON

### Manual KEGG Conversion

``` r

# Example: Glyoxylate cycle from KEGG M00012
# https://www.genome.jp/entry/M00012
```

**KEGG module format:**

    M00012  Glyoxylate cycle
      K01637  isocitrate lyase [EC:4.1.3.1]
      K01638  malate synthase [EC:2.3.3.9]

**Convert to V2 JSON:**

``` json
{
  "schema_version": "v2",
  "id": "glyoxylate_cycle",
  "name": "Glyoxylate Cycle",
  "source": "KEGG M00012",
  "tags": ["metabolism", "carbon"],

  "genes": [
    {
      "id": "aceA",
      "name": "isocitrate lyase",
      "databases": {"kofam": ["K01637"]},
      "ec": ["4.1.3.1"],
      "type": "enzyme"
    },
    {
      "id": "aceB",
      "name": "malate synthase",
      "databases": {"kofam": ["K01638"]},
      "ec": ["2.3.3.9"],
      "type": "enzyme"
    }
  ],

  "compounds": [
    {"id": "C00026", "name": "2-Oxoglutarate", "kegg_compound": "C00026"},
    {"id": "C00149", "name": "Malate", "kegg_compound": "C00149"}
  ],

  "pathways": {
    "main": {
      "name": "Glyoxylate Cycle",
      "type": "independent",
      "verified": false,
      "edges": [
        {"from": "C00026", "to": "aceA", "required": true, "marker": true},
        {"from": "aceA", "to": "C00149", "required": true, "marker": true},
        {"from": "C00149", "to": "aceB", "required": true, "marker": false}
      ],
      "scoring": {
        "min_fraction": 1.0,
        "marker_mode": "all"
      }
    }
  }
}
```

## Schema Details

### Required Fields

**Root level:** - `schema_version`: “v2” - `id`: Unique identifier
(letters, numbers, underscores) - `name`: Human-readable name - `genes`:
Array of gene definitions - `compounds`: Array of compounds -
`pathways`: Object with pathway definitions

**Gene:** - `id`: Gene identifier - `name`: Gene name - `databases`:
Detection methods (kofam, blast, hmm) - `type`: “enzyme”, “transport”,
“chaperone”, or “complex”

**Pathway:** - `name`: Pathway name - `type`: “variant” or
“independent” - `verified`: Always `false` for new pathways - `edges`:
Array of connections - `scoring`: Scoring parameters

### Database Types

**kofam** - KEGG Orthology

``` json
"databases": {
  "kofam": ["K01637", "K01638"]
}
```

**blast** - Protein sequences

``` json
"databases": {
  "blast": ["WP_000123456", "NP_987654"]
}
```

**hmm** - Profile HMMs (includes PFAM)

``` json
"databases": {
  "hmm": ["PF00463", "TIGR01234"]
}
```

### Protein Complexes

Multi-subunit complexes where **all** components are required:

``` json
{
  "genes": [
    {"id": "cutA", "databases": {"kofam": ["K18020"]}, "type": "enzyme"},
    {"id": "cutB", "databases": {"kofam": ["K18021"]}, "type": "enzyme"},
    {"id": "cutC", "databases": {"kofam": ["K18022"]}, "type": "enzyme"},
    {
      "id": "cutABC",
      "type": "complex",
      "components": ["cutA", "cutB", "cutC"],
      "notes": "All three subunits required"
    }
  ],
  "pathways": {
    "main": {
      "edges": [
        {"from": "C00577", "to": "cutABC", "required": true}
      ]
    }
  }
}
```

### Alternative Genes

Use separate edges for alternatives (OR logic):

``` json
"edges": [
  {"from": "C00026", "to": "geneA", "required": true},
  {"from": "C00026", "to": "geneB", "required": true}
]
```

Either geneA **OR** geneB satisfies the requirement.

### Scoring Parameters

**Fraction-based:**

``` json
"scoring": {
  "min_fraction": 0.75  // 75% of genes required
}
```

**Gap-based:**

``` json
"scoring": {
  "max_gaps": 1,              // Allow 1 missing gene
  "input": ["C00026"],        // Start compound
  "output": ["C00149"]        // End compound
}
```

**Dual scoring** (pathway passes if EITHER succeeds):

``` json
"scoring": {
  "min_fraction": 0.67,
  "max_gaps": 1,
  "input": ["C00026"],
  "output": ["C00149"]
}
```

## Validating Potatoes

Always validate before using:

``` r

# Load and validate
potato <- load_potato("my_custom_potato.json")

# Returns list with:
# $valid - TRUE/FALSE
# $errors - Critical issues (must fix)
# $warnings - Non-critical suggestions

result <- validate_potato(potato)

if (!result$valid) {
  cat("Errors:\n")
  print(result$errors)
}

if (length(result$warnings) > 0) {
  cat("Warnings:\n")
  print(result$warnings)
}
```

### Common Validation Errors

**Missing required fields:**

    Error: Missing required field: 'id'

**Invalid references:**

    Error: Edge references gene 'aceX' not found in genes array

**Cycle detection:**

    Warning: contains cycles (metabolic cycles are OK, but check carefully)

**Unreachable compounds:**

    Error: Input 'C00026' cannot reach any output

## Testing Potatoes

Test on known genomes before production use:

``` r

# Create test sack
sack <- create_sack()

# Add test potato
sack@potatoes <- list(my_potato = potato)

# Add test genome (known to have/lack pathway)
sack <- add_genomes(sack, "test_genomes/ecoli.faa")

# Run annotation
sack <- run_kofam(sack, conda_env = "potato")

# Score
sack <- score_pathways(sack)

# Check results match expectations
results <- get_pathway_scores(sack)
print(results)
```

## Visualizing Potatoes

Check pathway structure:

``` r

# Print text view
print_potato(potato)

# Interactive network plot
plot_potato(potato, interactive = TRUE)

# HTML detail view
view_pathway_detail(potato, pathway = "main")
```

## Multi-Pathway Networks

Group related pathways that share genes:

``` json
{
  "id": "nitrogen_fixation",
  "name": "Nitrogen Fixation Network",

  "genes": [
    {"id": "nifH", "databases": {"kofam": ["K02588"]}},
    {"id": "nifD", "databases": {"kofam": ["K02586"]}},
    {"id": "vnfH", "databases": {"kofam": ["K22896"]}}
  ],

  "pathways": {
    "mo_nitrogenase": {
      "name": "Mo-Nitrogenase",
      "type": "variant",
      "edges": [
        {"from": "nifH", "to": "nifD", "required": true}
      ]
    },
    "v_nitrogenase": {
      "name": "V-Nitrogenase",
      "type": "variant",
      "edges": [
        {"from": "vnfH", "to": "vnfG", "required": true}
      ]
    }
  }
}
```

## Best Practices

### 1. Start Simple

Begin with core genes only. Add alternatives and variants later.

### 2. Use Existing Potatoes as Templates

``` r

# Copy a similar potato as starting point
file.copy(
  system.file("potatoes/glyoxylate_cycle.json", package = "potato"),
  "my_potato.json"
)
```

### 3. Document Thoroughly

- `notes` at pathway level for biological context
- `notes` on genes for annotation sources
- Clear pathway `name` and `source`

### 4. Test Incrementally

- Validate after each change
- Test on positive/negative control genomes
- Compare scores to expectations

### 5. Never Set verified: true

Only humans verify. Start with `"verified": false` and update manually
after validation.

### 6. Version Control

Track potato JSON files in git to maintain history of changes.

## Example: Custom Biomarker Pathway

``` json
{
  "schema_version": "v2",
  "id": "my_biomarker",
  "name": "Custom Biomarker Set",
  "source": "Lab study 2024",
  "tags": ["biomarker", "custom"],
  "notes": "Genes identified in study XYZ as diagnostic markers",

  "genes": [
    {
      "id": "markerA",
      "name": "Biomarker A",
      "databases": {
        "kofam": ["K12345"],
        "blast": ["WP_123456789"]
      },
      "type": "enzyme",
      "notes": "Identified via comparative genomics"
    },
    {
      "id": "markerB",
      "name": "Biomarker B",
      "databases": {
        "hmm": ["custom_profile_B"]
      },
      "type": "enzyme",
      "notes": "Custom HMM built from alignment"
    }
  ],

  "compounds": [],

  "pathways": {
    "main": {
      "name": "Biomarker Presence",
      "type": "independent",
      "verified": false,
      "notes": "Presence of both markers indicates trait X",
      "edges": [],
      "scoring": {
        "min_fraction": 1.0,
        "marker_mode": "all",
        "notes": "Both markers required"
      }
    }
  }
}
```

## Getting Help

- **Schema reference:** See ROADMAP.md in package repo
- **Examples:** Browse `inst/potatoes/` directory
- **Interactive builder:** Use `/build-potato` skill
- **Validation:** Run
  [`validate_potato()`](https://jeffkimbrel.github.io/potato/reference/validate_potato.md)
  frequently

## Contributing Potatoes

To contribute a potato to the POTATO package:

1.  Create and validate your potato JSON
2.  Test on diverse genomes
3.  Document sources and verification
4.  Submit via GitHub pull request
