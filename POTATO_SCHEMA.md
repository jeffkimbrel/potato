# POTATO Multi-Pathway Network Schema

**Version:** 0.9.3  
**Date:** 2026-08-10  
**Last Updated:** 2026-08-10

## Overview

POTATO uses multi-pathway network JSON files to represent related metabolic pathways that share genes or metabolic context. Each potato defines:
1. **Global genes** - Detection methods (KO IDs, BLAST refs, HMM profiles) and metadata
2. **Pathways** - Pathway-specific topology, gene roles, and scoring parameters

## Top-Level Structure

```json
{
  "id": "pathway_network_id",
  "name": "Human Readable Network Name",
  "source": "KEGG M00XXX, M00YYY",
  "tags": ["metabolism", "carbohydrate"],
  "notes": "Brief description of the network",
  
  "nodes": [
    /* Global gene definitions */
  ],
  
  "pathways": {
    "pathway_id_1": { /* pathway 1 definition */ },
    "pathway_id_2": { /* pathway 2 definition */ }
  }
}
```

## Global Nodes Array

Defines **detection methods** for all genes in the network. Each gene appears **once** in this array.

```json
{
  "id": "geneSymbol",
  "name": "enzyme name",
  "databases": {
    "kofam": ["K00001"],
    "blast": ["ref_seq_1"],
    "hmm": ["PF00001"]
  },
  "ec": ["1.1.1.1"],
  "notes": "Biological context",
  "type": "enzyme",
  "x": 100,
  "y": 200,
  "x_compounds": 150,
  "y_compounds": 250
}
```

**Key points:**
- `id` must be unique across all genes
- `databases` lists detection methods by tool type
- Coordinates (`x`, `y`, `x_compounds`, `y_compounds`) are optional for visualization
- **NO step, required, or marker fields** - these are pathway-specific

## Pathways Object

Each pathway key maps to a pathway definition with:

```json
"pathway_id": {
  "name": "Pathway Display Name",
  "type": "variant",  // or "independent"
  "kegg_module": "M00123",
  "notes": "Pathway-specific notes",
  "verified": false,  // ALWAYS false unless manually validated
  
  "nodes": {
    "geneSymbol": {
      "step": 1,
      "required": true,
      "marker": false
    }
  },
  
  "edges": [
    {
      "from": "geneA",
      "to": "geneB",
      "compound": "metabolite name",
      "kegg_compound": "C00001"
    }
  ],
  
  "input": {
    "compound": "substrate",
    "kegg_compound": "C00031",
    "targets": ["geneA"]
  },
  
  "output": {
    "compound": "product",
    "kegg_compound": "C00118",
    "sources": ["geneZ"]
  },
  
  "scoring": {
    "min_fraction": 0.75,
    "marker_mode": "any"
  }
}
```

## Pathway-Specific Attributes

These fields are defined **per pathway** in the `pathways.pathway_id.nodes` object:

| Field | Type | Description |
|-------|------|-------------|
| `step` | integer | Sequential reaction step (1, 2, 3...) |
| `required` | boolean | Must be present for pathway completion? |
| `marker` | boolean | Diagnostic gene for this pathway? |

**Important:** The same gene can have different `step`, `required`, and `marker` values in different pathways.

**Example:**
```json
// In pathway "classic"
"gnaD": {"step": 3, "required": true, "marker": true}

// In pathway "non_phosphorylative"  
"gnaD": {"step": 1, "required": true, "marker": true}
```

## Edge Format

Edges define pathway topology using **gene IDs** (not step-based node IDs):

```json
{
  "from": "geneA",      // Source gene ID
  "to": "geneB",        // Target gene ID
  "compound": "glucose",
  "kegg_compound": "C00031"
}
```

- `from` and `to` must reference gene IDs in the global `nodes` array
- `compound` and `kegg_compound` are optional but recommended
- For OR branches (alternatives), create separate edges from previous step to each alternative
- **Compound name normalization:** Multi-part compounds (e.g., "pyruvate + G3P") have parts sorted alphabetically during graph construction to merge equivalent representations ("A + B" == "B + A")
- **Empty edges arrays:** Transporter pathways may have empty `edges: []` arrays since their topology is defined solely by input → genes → output connections

## Pathway Types

### variant
Alternative routes to the **same biological outcome**. Examples:
- Mo-nitrogenase vs V-nitrogenase (both fix N₂)
- ED classic vs ED semi-phosphorylative (both degrade glucose)

**Interpretation:** Finding one variant explains the absence of another. "Uses pathway A, not pathway B."

### independent
Different purposes but share metabolic space. Examples:
- TCA cycle vs glyoxylate shunt (respiration vs acetate assimilation)
- Forward TCA vs reverse TCA (oxidation vs carbon fixation)

**Interpretation:** Can have both, either, or neither. Each serves a distinct function.

## Shared Genes

A gene used by multiple pathways is:
1. Defined **once** in global `nodes`
2. Referenced in multiple pathway `nodes` with pathway-specific attributes

**Detection happens once** - if geneA is detected, it satisfies all pathways that use it.

## Scoring

Each pathway is scored independently:
- Count genes detected at each step
- Calculate `fraction = detected_steps / total_steps`
- Compare to `min_fraction` threshold
- Pathway is `present` if `fraction >= min_fraction`

## Input/Output Fields

Define substrate entry and product exit for each pathway:

```json
"input": {
  "compound": "D-glucose",
  "kegg_compound": "C00031",
  "targets": ["glk"]  // First step gene(s)
}
```

**For transporters:** Use location qualifiers:
- `"NH4_external"` → `"NH4_internal"`
- `"phosphate_periplasm"` → `"phosphate_cytoplasm"`

## Verification Status

```json
"verified": false  // Per-pathway field
```

- **ALWAYS `false`** for new/unverified pathways
- **NEVER set to `true`** programmatically
- Only humans manually verify after validation
- Verified pathways should not be edited without explicit approval

## Coordinate Fields

For visualization:

| Field | Usage |
|-------|-------|
| `x`, `y` | Enzyme-only layout coordinates |
| `x_compounds`, `y_compounds` | Bipartite layout with compound nodes |

Compound nodes can be positioned explicitly via the top-level `compound_coordinates` object:

```json
{
  "compound_coordinates": {
    "COMPOUND_C00031": {
      "id": "COMPOUND_C00031",
      "x": 100,
      "y": 200
    }
  }
}
```

Compound nodes are auto-positioned if not specified in `compound_coordinates`.

## Example: Bifunctional Enzymes

Gene that catalyzes reactions at multiple steps:

```json
// Global nodes
{
  "id": "sskdgK",
  "name": "bifunctional 2-keto-3-deoxygluconokinase",
  "databases": {"kofam": ["K18126"]},
  "notes": "Bifunctional kinase; phosphorylates KDG to KDPG (C00204→C04442) and phosphorylates D-glycerate to 2-phospho-D-glycerate (C01216→C01286)"
}

// In pathway
"nodes": {
  "sskdgK": {"step": 2, "required": true, "marker": true}
}

// Edges - gene appears in multiple reactions
{"from": "gnaD", "to": "sskdgK", "compound": "KDG"},
{"from": "sskdgK", "to": "kdpgA", "compound": "KDPG"}
```

## Common Patterns

### OR Branches (Alternative Enzymes)
Multiple genes at the same step:
```json
"nodes": {
  "geneA": {"step": 2, "required": true},
  "geneB": {"step": 2, "required": true}
}
"edges": [
  {"from": "gene1", "to": "geneA"},
  {"from": "gene1", "to": "geneB"},
  {"from": "geneA", "to": "gene3"},
  {"from": "geneB", "to": "gene3"}
]
```

### Optional Steps
Gene not required for pathway completion:
```json
"pgl": {"step": 3, "required": false, "marker": false}
```

### Marker Genes
Diagnostic genes specific to this pathway:
```json
"mlrA": {"step": 1, "required": true, "marker": true}
```

## Validation Rules

1. All `edges.from` and `edges.to` must reference genes in global `nodes`
2. All pathway `nodes` keys must reference genes in global `nodes`
3. Step numbers should be sequential (1, 2, 3...)
4. At least one gene should have `marker: true` per pathway
5. No cycles allowed in edge topology (must be DAG)
6. `input.targets` and `output.sources` must reference genes in pathway `nodes`

## Design Principles

1. **Genes defined once** - Detection methods in global `nodes`, context in pathway `nodes`
2. **Pathway-specific roles** - Same gene can be step 1 in one pathway, step 3 in another
3. **Independent scoring** - Each pathway evaluated separately
4. **Biological context** - Pathway types explain relationships (variant vs independent)
5. **Shared detection** - Gene detected once, counts for all pathways

## Further Reading

- Full implementation details: `CLAUDE.md`
- Development roadmap: `ROADMAP.md`
- Workflow diagrams: `WORKFLOW.md`
