# POTATO V2 Schema Reference

**Schema Version:** `v2` (frozen/stable)  
**Package Version:** v0.14.5+  
**Last Updated:** 2026-11-02

This document defines the V2 schema for POTATO pathway definition files (potatoes). The V2 schema is frozen and stable - breaking changes would require V3.

---

## Overview

A potato is a JSON file defining:
- **Genes** - Detection methods (KEGG, PFAM, BLAST, HMM) and protein complexes
- **Compounds** - Metabolic intermediates (KEGG compound IDs)
- **Pathways** - One or more pathway variants with topology (edges) and scoring parameters
- **Metadata** - Source, tags, notes for organization

Potatoes are **self-contained** - all information needed to detect and score a pathway is in one file.

---

## Top-Level Structure

```json
{
  "schema_version": "v2",
  "id": "pathway_name",
  "name": "Human Readable Name",
  "source": "KEGG M00123 / custom",
  "tags": ["metabolism", "energy"],
  "notes": "Brief biological description",
  
  "genes": [ /* gene definitions */ ],
  "compounds": [ /* compound definitions */ ],
  
  "pathways": {
    "pathway_id": { /* pathway definition */ }
  }
}
```

### Required Fields

- `schema_version` - Must be `"v2"`
- `id` - Unique identifier (snake_case, alphanumeric + underscore/hyphen)
- `name` - Human-readable display name
- `genes` - Array of gene/complex definitions (at least one)
- `pathways` - Object with at least one pathway

### Optional Fields

- `source` - Origin (KEGG module, literature citation, custom)
- `tags` - Array of strings for organization
- `notes` - Biological context, caveats, references
- `compounds` - Array of compound definitions (required if edges reference compounds)

---

## Gene Definitions

Genes define **detection methods only** - no pathway-specific attributes (step, required, marker) at this level.

### Regular Gene

```json
{
  "id": "nifH",
  "name": "nitrogenase iron protein",
  "type": "enzyme",
  "databases": {
    "kofam": ["K02588"],
    "blast": ["ref_seq_id"],
    "hmm": ["PF00142", "custom_hmm"]
  },
  "ec": ["1.18.6.1"],
  "reactions": ["R05185"],
  "notes": "Highly conserved marker gene",
  "x": 100,
  "y": 200
}
```

**Fields:**
- `id` - **Required**. Unique gene identifier (referenced by edges, complexes)
- `name` - **Required**. Human-readable enzyme name
- `type` - **Recommended**. One of: `"enzyme"`, `"transport"`, `"chaperone"`, `"complex"`
- `databases` - **Required for non-complex genes**. Object with database types as keys:
  - `kofam` - Array of KEGG Orthology IDs (K##### format)
  - `blast` - Array of reference sequence IDs
  - `hmm` - Array of HMM profile NAMEs (includes PFAM PF##### and custom HMMs)
- `ec` - Optional. Array of EC numbers
- `reactions` - Optional. Array of KEGG reaction IDs
- `notes` - Optional. Biological context, caveats
- `x`, `y` - Optional. Curated visualization coordinates

**Important:**
- **PFAM profiles go in `hmm` field** (PFAM is a type of HMM database)
- Standard database types only: `kofam`, `blast`, `hmm` (no custom names)
- NO `step`, `required`, `marker` at gene level (these go on edges or pathway-level)

### Protein Complex

```json
{
  "id": "nifHDK",
  "name": "molybdenum nitrogenase complex",
  "type": "complex",
  "components": ["nifH", "nifD", "nifK"],
  "notes": "All three components required for N2 fixation"
}
```

**Fields:**
- `id` - **Required**. Complex identifier (referenced by edges)
- `name` - **Required**. Complex name
- `type` - **Required**. Must be `"complex"`
- `components` - **Required**. Array of gene IDs (all must exist as separate genes)
- `notes` - Optional. Biological context

**Important:**
- NO `databases` field on complexes (detection via components)
- All components must exist as separate gene entries with their own detection methods
- Edges reference complex ID (appears as single node in visualization)
- Pathway's `genes` array lists components (not complex ID) for scoring

**Naming convention:** Use concatenated component names (e.g., `ureABC`, `nifHDK`, `napAB`)

### Gene Types

- **`enzyme`** - Catalytic gene (most common)
- **`transport`** - Transporter/permease
- **`chaperone`** - Non-catalytic accessory gene (electron carriers, regulators, maturation factors)
- **`complex`** - Multi-component protein complex

---

## Compound Definitions

```json
{
  "id": "C00697",
  "name": "Nitrogen",
  "kegg_compound": "C00697"
}
```

**Fields:**
- `id` - **Required**. Compound identifier (referenced by edges, input/output)
- `name` - **Required**. Human-readable compound name
- `kegg_compound` - Optional. KEGG compound ID

**Location qualifiers:** For transport pathways, use suffixes like `"NH4_external"`, `"NH4_internal"`

---

## Pathway Definitions

Each potato can have multiple pathways (alternative variants, independent functions).

```json
"pathways": {
  "pathway_id": {
    "name": "Pathway Display Name",
    "type": "variant",
    "verified": false,
    "notes": "Pathway-specific biological context",
    "tags": ["nitrogen", "anaerobic", "diagnostic"],
    
    "genes": ["geneA", "geneB", "geneC"],
    "markers": ["geneA", "geneB"],
    "input": ["C00697"],
    "output": ["C00014"],
    
    "edges": [ /* bipartite graph */ ],
    
    "scoring": {
      "min_fraction": 0.75,
      "max_gaps": 0
    }
  }
}
```

### Required Fields

- `name` - Human-readable pathway name
- `type` - Pathway type: `"variant"` or `"independent"`
  - `variant` - Alternative routes to same outcome (Mo-nif vs V-nif)
  - `independent` - Different functions, shared metabolic space (TCA vs glyoxylate)
- `verified` - **ALWAYS `false` for new pathways**. Only humans set to `true` after validation.
- `edges` - Array of edges connecting genes and compounds (can be empty for transport pathways)
- `scoring` - Scoring parameters (see Scoring Methods below)

### Optional Fields

- `notes` - Pathway-specific biological context
- `kegg_module` - KEGG module ID if applicable
- `tags` - Array of pathway-specific classification tags (e.g., ["nitrogen", "anaerobic", "diagnostic"])
- `genes` - Explicit gene list (includes non-catalytic genes like chaperones). If omitted, genes inferred from edges.
- `markers` - Array of gene IDs that are diagnostic markers (individual genes, NOT complexes)
- `input` - Array of compound IDs defining pathway start (required for gap-based scoring)
- `output` - Array of compound IDs defining pathway end (required for gap-based scoring)

---

## Edges (Bipartite Graph)

Edges connect genes AND compounds in a bipartite graph. The `required` attribute lives on edges (not genes) because same gene can be required in one pathway but optional in another.

```json
{
  "from": "C00697",
  "to": "nifH",
  "required": true,
  "reaction": "R05185"
}
```

**Fields:**
- `from` - **Required**. Gene ID or compound ID
- `to` - **Required**. Gene ID or compound ID
- `required` - **Required**. Boolean. Is this edge required for pathway functionality?
- `reaction` - Optional. KEGG reaction ID

**DEPRECATED:** `marker` field on edges (v0.14.5+). Use pathway-level `markers` array instead.

**Important:**
- Edges can connect: compound→gene, gene→compound, or gene→gene
- Same gene can have different `required` values in different pathways
- Edges can reference complex IDs (e.g., `"to": "nifHDK"`)
- Empty edges array is valid for transport pathways with no internal reactions

**Why edge-level `required`?** Same gene can act on different substrates in different pathways. Edge-level attributes capture this pathway-specific context.

---

## Marker Genes

**Marker genes** are diagnostic - their presence strongly indicates pathway presence.

### Pathway-Level Markers (v0.14.5+)

**RECOMMENDED METHOD:** Define markers at pathway level:

```json
"pathways": {
  "nitrogen_fixation": {
    "markers": ["nifH", "vnfH"]
  }
}
```

Gene IDs must exist in pathway's `genes` array (explicit or inferred from edges).

**Guidelines:**
- List individual genes (NOT complex IDs) - e.g., `["urtA"]` not `["urtABCDE"]`
- For complexes: list most diagnostic component(s) - e.g., `["nirB"]` from nirBD complex
- Markers = diagnostic genes only, not comprehensive lists
- Use biological judgment: which genes are pathway-specific?

### Edge-Level Marker Attribute (DEPRECATED)

```json
"edges": [
  {"from": "substrate", "to": "nifH", "marker": true}  // DEPRECATED - use pathway-level markers
]
```

**Status:** Edge-level `marker` attribute deprecated in v0.14.5. Validation warns if present. Use pathway-level `markers` array instead.

### Marker Mode

`marker_mode` is a **runtime parameter** for `run_kofam()`, NOT a potato field:
- `run_kofam(marker_mode = "any")` - Screen for pathway presence using markers only
- `run_kofam(marker_mode = "all")` - Require all markers detected

**Removed from schema (v0.14.5):** `marker_mode` field in potato scoring blocks. Validation warns if present.

---

## Scoring Methods

Two scoring methods available (can use both simultaneously):

### Fraction-Based Scoring

```json
"scoring": {
  "min_fraction": 0.75
}
```

**How it works:** Counts detected genes / total genes. Pathway passes if fraction ≥ threshold.

**Limitations:**
- Float precision issues (0.67 threshold fails 2/3 genes at 0.6666...)
- Ignores topology (80% with broken middle step still scores 0.8)
- Doesn't verify connectivity

**Use when:** Gene count matters more than connectivity, or working with incomplete MAGs where any gene presence is informative.

### Gap-Based Scoring (v0.11.0+)

```json
"input": ["C00697"],
"output": ["C00014"],
"scoring": {
  "max_gaps": 0
}
```

**How it works:** DAG traversal from `input` → `output` compounds. Finds best path (fewest missing genes) considering OR branches. Pathway passes if gap_count ≤ max_gaps.

**Benefits:**
- Biologically meaningful ("allow 1 missing gene")
- Verifies connectivity (ensures complete route exists)
- Handles alternative routes correctly (finds best path)
- No float precision issues

**Requirements:**
- `input` and `output` must be defined (arrays of compound IDs)
- Compounds must exist in `compounds` array
- Path must be reachable (no disconnected subgraphs)

**Use when:** Pathway connectivity matters, or pathway is complex with branching.

### Dual Scoring

```json
"input": ["C00697"],
"output": ["C00014"],
"scoring": {
  "min_fraction": 0.75,
  "max_gaps": 1
}
```

Pathway passes if **EITHER** method succeeds. Both metrics reported in results.

---

## Complete Example

```json
{
  "schema_version": "v2",
  "id": "nitrogen_fixation_network",
  "name": "Nitrogen Fixation Network",
  "source": "KEGG M00175",
  "tags": ["nitrogen", "energy"],
  "notes": "Two nitrogenase systems: Mo-dependent (common) and V-dependent (Mo-scarce).",
  
  "genes": [
    {
      "id": "nifH",
      "name": "nitrogenase iron protein",
      "type": "enzyme",
      "databases": {"kofam": ["K02588"]},
      "notes": "Fe protein, highly conserved marker"
    },
    {
      "id": "nifD",
      "name": "nitrogenase molybdenum-iron protein alpha chain",
      "type": "enzyme",
      "databases": {"kofam": ["K02586"]},
      "ec": ["1.18.6.1"]
    },
    {
      "id": "nifK",
      "name": "nitrogenase molybdenum-iron protein beta chain",
      "type": "enzyme",
      "databases": {"kofam": ["K02591"]},
      "ec": ["1.18.6.1"]
    },
    {
      "id": "nifHDK",
      "name": "molybdenum nitrogenase complex",
      "type": "complex",
      "components": ["nifH", "nifD", "nifK"],
      "notes": "Complete Mo-dependent nitrogenase"
    },
    {
      "id": "vnfH",
      "name": "vanadium nitrogenase iron protein",
      "type": "enzyme",
      "databases": {"kofam": ["K22899"]}
    },
    {
      "id": "vnfD",
      "name": "vanadium-dependent nitrogenase alpha chain",
      "type": "enzyme",
      "databases": {"kofam": ["K22896"]}
    },
    {
      "id": "vnfK",
      "name": "vanadium-dependent nitrogenase beta chain",
      "type": "enzyme",
      "databases": {"kofam": ["K22897"]}
    },
    {
      "id": "vnfG",
      "name": "vanadium nitrogenase delta subunit",
      "type": "enzyme",
      "databases": {"kofam": ["K22898"]}
    },
    {
      "id": "vnfHDKG",
      "name": "vanadium nitrogenase complex",
      "type": "complex",
      "components": ["vnfH", "vnfD", "vnfK", "vnfG"],
      "notes": "Complete V-dependent nitrogenase"
    },
    {
      "id": "nifF",
      "name": "flavodoxin",
      "type": "chaperone",
      "databases": {"kofam": ["K00531"]},
      "notes": "Electron carrier, optional/accessory"
    }
  ],
  
  "compounds": [
    {"id": "C00697", "name": "Nitrogen", "kegg_compound": "C00697"},
    {"id": "C00014", "name": "Ammonia", "kegg_compound": "C00014"}
  ],
  
  "pathways": {
    "nitrogen_fixation": {
      "name": "Nitrogen Fixation (Mo or V dependent)",
      "type": "variant",
      "verified": false,
      "notes": "Either complete system indicates N2 fixation capability",
      
      "genes": ["nifH", "nifD", "nifK", "vnfH", "vnfD", "vnfK", "vnfG", "nifF"],
      "markers": ["nifH", "vnfH"],
      "input": ["C00697"],
      "output": ["C00014"],
      
      "edges": [
        {"from": "C00697", "to": "nifHDK", "required": false, "reaction": "R05185"},
        {"from": "nifHDK", "to": "C00014", "required": false, "reaction": "R05185"},
        {"from": "C00697", "to": "vnfHDKG", "required": false},
        {"from": "vnfHDKG", "to": "C00014", "required": false}
      ],
      
      "scoring": {
        "max_gaps": 0,
        "notes": "Gap-based finds complete Mo-nif OR V-nif route"
      }
    }
  }
}
```

---

## Validation Rules

The `validate_potato()` function checks:

### Structural
- `schema_version` is `"v2"`
- Required top-level fields present
- Gene/compound IDs unique
- No duplicate pathway IDs

### Gene Definitions
- Regular genes have `databases` field
- Complexes have `components` array (not `databases`)
- Complex components exist as genes
- NO `step`, `required`, `marker` on genes (these go on edges/pathways)

### Pathways
- Each pathway has `name`, `type`, `verified`, `edges`, `scoring`
- `verified` should be `false` for new pathways
- Edge `from`/`to` reference valid gene/compound IDs
- No cycles in DAG (warning only)
- If `input`/`output` defined, compounds exist
- If `markers` field exists, gene IDs exist in pathway

### Scoring
- If `max_gaps` present, `input` and `output` required
- `min_fraction` not between valid gene fractions (float precision trap)
- `marker_mode` field deprecated (v0.14.5) - validation warns if present

### Markers
- Edge-level `marker` attribute deprecated (v0.14.5) - validation warns if present
- Pathway-level `markers` array: gene IDs must exist in pathway's genes
- Markers must be individual genes, NOT complex IDs

### Orphaned Entities
- All defined genes used in at least one pathway (edges OR explicit `genes` array)
- Complex components counted as "used" when complex is used
- All defined compounds used in at least one pathway

---

## Package Version Features

Features tied to specific package versions:

- **v0.9.0+** - Multi-pathway networks (multiple pathways per potato)
- **v0.11.0+** - Gap-based scoring (`max_gaps`, `input`, `output`)
- **v0.11.0+** - Simplified input/output (arrays of IDs, not verbose objects)
- **v0.14.4+** - Pathway-level `markers` field
- **v0.14.4+** - Complex components automatically marked as "used" in validation
- **v0.14.5+** - Edge-level `marker` deprecated (use pathway-level `markers`)
- **v0.14.5+** - `marker_mode` removed from potato schema (runtime parameter only)
- **v0.14.5+** - Pathway-level and potato-level `tags` for filtering
- **v0.14.5+** - Dark mode theme support for visualization (`view_pathway_detail`, interactive plots)

Schema version stays `v2` across all these releases - only implementation features evolve.

---

## See Also

- **ROADMAP.md** - Planned schema enhancements (V3 would require breaking changes)
- **CLAUDE.md** - Project context and design decisions
- **inst/.claude/commands/build-potato.md** - Agent instructions for building potatoes
- **inst/potatoes/** - Example potato files

---

**Questions or Issues?** https://github.com/jeffkimbrel/potato/issues
