# POTATO - Project Context for Claude

## Project Overview

**POTATO** (Pathway annOTATOr) is an R package for annotating MAGs (metagenome-assembled genomes) against curated metabolic pathways. It's the successor to GATOR (Genome annotATOR), redesigned around self-contained pathway definitions (potatoes) as DAG structures in JSON.

**Current Status:** v0.9.3 (2026-08-05) - Multi-pathway networks with full scoring support, essential-only scoring metrics, result export functions, 27% test coverage (91 passing tests).

**Key Innovation:** Each "potato" (pathway) is a self-contained JSON file defining:
- Genes with multi-tool detection methods (KEGG, PFAM, BLAST, HMM)
- Pathway structure as a DAG (handles complex branching)
- Confidence scoring via gene specificity weighting
- LLM-assisted building and interpretation

---

## Important Files

- **[ROADMAP.md](ROADMAP.md)** - Complete implementation plan, phases, file formats, future enhancements
- **[WORKFLOW.md](WORKFLOW.md)** - Mermaid diagrams of current workflow and architecture
- **[environment.yaml](environment.yaml)** - Conda environment with bioinformatics tools
- **R/zzz.R** - Reticulate setup, loads Python backend on package load
- **inst/python/** - Python backend code (will be rewritten for v1)
- **inst/potatoes/** - Example potato JSON files (7 examples)
- **inst/.claude/commands/build-potato.md** - Agent instructions for building potato JSONs (updated v0.7.1 - expert guide, pushes back on vague requests, suggests PFAM/BLAST proactively)

---

## Key Design Decisions

### 1. Self-Contained Potatoes
Each potato JSON defines its own genes + pathway logic. No central gene registry (unless user wants optional canonical genes for consistency).

**Why:** Portability, flexibility, version control, independent evolution

### 2. DAG Structure
Pathway logic is a directed acyclic graph in JSON, not a string in Excel.

**Why:** Handles KEGG module complexity, visually graphable, LLM can reason about structure

### 3. Gene Specificity Weighting
Pre-computed score: `specificity = 1 / (num potatoes containing gene)`

**Why:** Pathway-specific genes (mlrA, nifH) are more informative than ubiquitous ones (gapA). Confidence = weighted sum of found genes.

### 4. Tool Configuration
Separate `tools.json` for user-specific database paths. Continue with available tools if some are disabled.

**Why:** Users have different setups (KEGG is expensive). Graceful degradation.

**Future enhancement (Phase 5):** Optionally embed HMM/BLAST sequences directly in potato JSON for ultimate portability. Makes tools.json optional if all data is embedded.

### 5. Compounds as Edge Metadata
Metabolites decorate edges, not structural nodes.

**Why:** Simple scoring logic + biological context for LLM reasoning

---

## Development Workflow

**CRITICAL: Always use devtools::load_all() when testing**

We are developing the package and potatoes simultaneously. The installed version is ALWAYS out of date.

```r
# CORRECT - loads latest code from R/ directory
devtools::load_all()
pot <- load_potato("inst/potatoes/my_potato.json")

# WRONG - uses old installed version
library(potato)
pot <- load_potato("inst/potatoes/my_potato.json")
```

When testing any R functions, prefix your Rscript commands with `devtools::load_all()`:

```bash
# Example testing pattern
Rscript -e "devtools::load_all(); pot <- load_potato('inst/potatoes/test.json'); validate_potato(pot)"
```

## Technology Stack

- **R** - Primary interface, analysis, visualization (igraph for DAGs)
- **Python** (via reticulate) - Tool orchestration, file I/O, LLM calls
- **jakomics** - Python utilities for tool execution (blast, kofam, hmmer)
- **Conda** - Bioinformatics tools (kofamscan, hmmer3, blast)
- **JSON** - Potato definitions, tool config
- **Claude API** - LLM agents (builder, converter, analysis)

### Dependency Architecture

```
potato (R package)
    ↓
reticulate bridge
    ↓
inst/python/potato.py (new potato-specific code)
    ↓
jakomics module (tool runners, file utilities)
    ↓
bioinformatics tools (kofamscan, hmmer3, blastp)
```

**Key separation:**
- **jakomics** handles tool execution, file I/O, parsing outputs
- **potato Python code** handles potato JSON, DAG logic, scoring
- **potato R code** provides user interface, visualization, analysis

---

## Potato JSON Structure

Each potato is a self-contained pathway definition with genes and topology:

```json
{
  "id": "pathway_id",
  "name": "Human Readable Name",
  "source": "KEGG M00123 / custom",
  "verified": false,
  "tags": ["metabolism", "energy"],
  "notes": "Brief description",
  
  "input": {
    "compound": "substrate name",
    "kegg_compound": "C00001",
    "targets": ["geneA_1"]
  },
  
  "output": {
    "compound": "product name",
    "kegg_compound": "C00002",
    "sources": ["geneZ_5"]
  },
  
  "nodes": [
    {
      "id": "geneSymbol",
      "step": 1,
      "nodes": ["geneSymbol_1"],
      "type": "enzyme",
      "name": "enzyme name",
      "databases": {
        "kofam": ["K00001"],              // KEGG Orthology IDs
        "blast": ["ref_seq_id"],          // BLAST reference sequence IDs
        "hmm": ["PF00001", "custom_hmm"]  // HMM profile NAMEs (includes PFAM)
      },
      "thresholds": {                     // OPTIONAL: per-gene overrides
        "blast_evalue": 1e-20,
        "hmm_evalue": 1e-15
      },
      "ec": ["1.1.1.1"],
      "required": true,
      "marker": false,
      "notes": "Biological context"
    }
  ],
  
  "edges": [
    {
      "from": "geneA_1",
      "to": "geneB_2",
      "compound": "metabolite name",
      "kegg_compound": "C00031"
    }
  ],
  
  "scoring": {
    "min_fraction": 0.75,
    "marker_mode": "any"
  }
}
```

**Key fields:**
- `verified`: **CRITICAL** - Always set to `false` for new/unverified potatoes. **NEVER set to true**. Only the user manually sets this after validation.
- `input`/`output`: Starting substrate and final product (optional but recommended)
  - For metabolic pathways: actual metabolites (e.g., "D-glucose", "pyruvate")
  - For transporters: location-qualified (e.g., "NH4_external", "NH4_internal")
  - `targets`/`sources`: which DAG nodes connect to input/output
- `databases`: Detection methods using standard types (kofam, blast, hmm, pfam)
- `step`: Sequential step number (or array for bifunctional enzymes)
- `nodes`: DAG node IDs in `id_step` format
- `marker`: Diagnostic gene for this pathway
- `required`: Must be present for pathway completion

**Important notes:**
- **Verified field:** Agents and Claude should ALWAYS set `"verified": false` and NEVER change it to true. Only humans verify potatoes.
- **Active field:** Optional `"active": false` marks deprecated potatoes that shouldn't be loaded by default
- **Standard database types only:** `kofam`, `blast`, `hmm` (no custom names like `kofam118`, `gator_blast`)
- **PFAM profiles:** Go in `hmm` field (PFAM is a type of HMM database), NOT separate `pfam` field
- **HMM profile names:** Use NAME from HMM file header (e.g., `NAME mlrA`), not filename
- **Per-gene thresholds:** Optional `thresholds` field for overriding global defaults (use sparingly)
- **No legacy fields:** Don't use `ko`, `blast_terms`, `hmm_path`, etc. Use `databases` only
- **Notes fields:** Use liberally to document biological context, alternatives, caveats, marker rationale
- **Input/output:** For transporters moving compounds across membranes, use location qualifiers like "_external", "_internal", "_periplasm"

---

## Multi-Pathway Schema (New in v0.9.0)

**Design Philosophy:** Related pathways that share metabolic context should be consolidated into network potatoes with multiple sub-pathways. This provides biological context that isolated potatoes cannot capture.

### When to Use Multi-Pathway Schema

**Use multi-pathway networks when:**
- Pathways are **alternative routes** to the same biological outcome (e.g., Mo vs V nitrogenase)
- Pathways share **metabolic intermediates** and spatial context (e.g., ED pathway variants)
- Finding one variant **explains the absence** of another (e.g., "ED classic absent because ED semi-phos present")
- Pathways **overlay spatially** but serve different purposes (e.g., TCA + glyoxylate shunt)

**Keep as separate potatoes when:**
- Pathways are functionally unrelated
- No shared genes or metabolic context
- Independent detection is more informative

### Multi-Pathway Structure

```json
{
  "id": "entner_doudoroff_network",
  "name": "Entner-Doudoroff Pathway Network",
  "source": "KEGG M00006, M00309, M00308, M00633",
  "active": true,
  "verified": false,
  "tags": ["metabolism", "carbohydrate"],
  "notes": "Four ED variants with different phosphorylation strategies.",
  
  "nodes": [
    {
      "id": "gnaD",
      "name": "gluconate dehydratase",
      "databases": {"kofam": ["K05308"]},
      "ec": ["4.2.1.140"],
      "notes": "Shared by 3 ED variants"
    },
    {
      "id": "edd",
      "name": "phosphogluconate dehydratase",
      "databases": {"kofam": ["K01690"]},
      "ec": ["4.2.1.12"]
    }
    // ... all 20 unique genes (detection methods only)
  ],
  
  "pathways": {
    "classic": {
      "name": "Classic ED (Phosphorylative)",
      "type": "variant",
      "kegg_module": "M00006",
      "notes": "Fully phosphorylated pathway",
      
      "nodes": {
        "glk": {"step": 1, "required": true, "marker": false},
        "zwf": {"step": 2, "required": true, "marker": false},
        "pgl": {"step": 3, "required": false, "marker": false},
        "ybhE": {"step": 3, "required": false, "marker": false},
        "edd": {"step": 4, "required": true, "marker": true},
        "eda": {"step": 5, "required": true, "marker": true}
      },
      
      "edges": [
        {"from": "glk", "to": "zwf", "compound": "D-glucose-6P"},
        {"from": "zwf", "to": "pgl", "compound": "6-phospho-D-glucono-1,5-lactone"},
        {"from": "zwf", "to": "ybhE", "compound": "6-phospho-D-glucono-1,5-lactone"},
        {"from": "pgl", "to": "edd", "compound": "6-phospho-D-gluconate"},
        {"from": "ybhE", "to": "edd", "compound": "6-phospho-D-gluconate"},
        {"from": "edd", "to": "eda", "compound": "KDPG"}
      ],
      
      "input": {
        "compound": "D-glucose",
        "kegg_compound": "C00031",
        "targets": ["glk"]
      },
      "output": {
        "compound": "pyruvate + glyceraldehyde-3P",
        "sources": ["eda"]
      },
      "scoring": {
        "min_fraction": 0.8,
        "marker_mode": "any"
      }
    },
    
    "non_phosphorylative": {
      "name": "Non-Phosphorylative ED",
      "type": "variant",
      "kegg_module": "M00309",
      // ... similar structure
    }
  }
}
```

### Key Principles

**1. Genes Defined Once**
- Global `nodes` array: Detection methods, EC numbers, enzyme names
- Pathway-specific `nodes`: Step number, required/marker status (contextual)
- Same gene can be step 1 in one pathway, step 3 in another
- Same gene can be marker in one pathway, not in another

**2. Pathway Types**
- `"type": "variant"` - Alternative routes to same outcome (Mo vs V nitrogenase, ED variants)
- `"type": "independent"` - Different purpose, shares metabolic space (TCA vs glyoxylate shunt)

**3. Pathway-Specific Attributes**
- `step`, `required`, `marker` - All pathway-specific, defined per pathway
- `edges` - Each pathway has its own topology
- `input`/`output` - Each pathway has its own substrates/products
- `scoring` - Each pathway scored independently

**4. Shared Genes**
- Gene `gnaD` appears in 3 ED pathways at step 1
- Detecting `gnaD` once satisfies all 3 pathways
- Each pathway evaluates completion independently

### Scoring Multi-Pathway Networks

**Current implementation (to be updated):**
Each pathway scored independently, results show per-pathway status:

```
potato: entner_doudoroff_network
  pathway: classic (variant) - present (6/6 genes, 1.0)
  pathway: non_phosphorylative (variant) - absent (2/6 genes, 0.33)
  pathway: semi_phosphorylative (variant) - present (7/8 genes, 0.875)
  pathway: semi_phosphorylative_alt (variant) - absent (1/4 genes, 0.25)

interpretation: ED capability via 2 variants (classic + semi_phos)
                Glucose AND gluconate utilization
```

**Benefits over separate potatoes:**
- Context: "ED classic absent" + "ED semi present" → explains the gap
- Efficiency: Shared gene detected once, counts for multiple pathways
- Biology: Shows metabolic flexibility vs. single strategy

### Visualization Implications

**Network plots** can show:
- All pathways overlaid with color-coding
- Shared nodes highlighted
- Active pathways in bold, inactive grayed out
- Multiple substrate entry points visible

**Planned features:**
- Toggle individual pathways on/off in visualization
- Highlight path-specific edges
- Show which genes satisfy multiple pathways

### Example: TCA + Glyoxylate Network

```json
{
  "id": "tca_network",
  "pathways": {
    "tca_forward": {
      "type": "variant",
      "notes": "Oxidative TCA for respiration",
      "nodes": {
        "citrate_synthase": {"step": 1, "marker": true},
        "aconitase": {"step": 2, "marker": false}
      }
    },
    "tca_reverse": {
      "type": "variant",
      "notes": "Reductive TCA for carbon fixation",
      "nodes": {
        "atp_citrate_lyase": {"step": 1, "marker": true},
        "aconitase": {"step": 2, "marker": false}
      }
    },
    "glyoxylate_shunt": {
      "type": "independent",
      "notes": "C2 assimilation bypassing CO2 loss",
      "nodes": {
        "aceA": {"step": 1, "marker": true},
        "aceB": {"step": 2, "marker": true}
      }
    }
  }
}
```

**Biological interpretation:**
- Forward TCA + glyoxylate → can grow on acetate/fatty acids
- Reverse TCA alone → carbon fixation (chemolithoautotroph)
- Forward + reverse → flexible metabolism

---

## Working Guidelines

### Code Style

**R:**
- Use tidyverse style (snake_case, pipes)
- Document functions with roxygen2
- Keep functions focused and testable

**Python:**
- PEP 8 style
- Type hints for function signatures
- Dataclasses for structured data

**General:**
- Clear variable names (no abbreviations unless standard)
- Comments explain *why*, not *what*
- Keep functions under 50 lines where possible

### File Organization

```
R/                      # R package functions
  ├── io.R             # File reading, potato loading
  ├── tools.R          # Tool execution wrappers  
  ├── scoring.R        # Pathway scoring, DAG traversal
  ├── agents.R         # LLM agent interfaces
  └── visualization.R  # igraph plotting

inst/
  ├── python/          # Python backend
  │   ├── potato.py   # Potato class, DAG handling
  │   ├── tools.py    # Tool runners (kofam, hmmer, blast)
  │   ├── scoring.py  # Scoring engine
  │   └── agents.py   # LLM interactions (Claude API)
  ├── potatoes/        # Potato JSON files
  │   ├── leucine_biosynthesis.json
  │   ├── nitrogen_fixation.json
  │   └── ...
  └── extdata/
      └── tools.json   # Example tool configuration

tests/                 # Test suite
  ├── testthat/
  └── test_data/       # Test genomes, potatoes
```

### Testing Strategy

- **Unit tests** - Each function, edge cases
- **Integration tests** - Full genome → annotation → scoring workflow
- **Validation tests** - Compare against known results (GATOR v1 output)
- **Test data** - Small genomes (<100 genes), simple potatoes (3-5 genes)

### Development Workflow

1. **Implement in Python first** - Core logic (Potato class, tool runners, scoring)
2. **Wrap in R** - Thin wrappers via reticulate
3. **Test incrementally** - Don't wait for full phase completion
4. **Validate early** - Compare against GATOR v1 for known pathways

---

## Recent Simplification Decisions (2026-07-28)

### Database Configuration - Use Standard Types

**Decision:** Database keys in `potato_config.yaml` should be simple type names (`kofam`, `blast`, `hmm`, `patric`), NOT custom names like `kofam113`, `gator_blast`, `nifH_hmm`.

**Note:** `patric` is valid (not `pfam`). PFAM is a type of HMM profile, so it falls under the `hmm` type.

**Rationale:** 
- Portability - potato JSONs reference database types, not user-specific names
- Simplicity - fewer places for discrepancies between config and potatoes
- One database per type is sufficient for most use cases

**Implementation:**
- Config template uses: `kofam:`, `blast:`, `hmm:`, `patric:`
- Potato JSONs reference same: `"databases": {"kofam": [...], "blast": [...]}`
- Validation strictly enforces only these types (errors on any other database name)
- If user needs multiple versions, they manage externally (merge files, use hierarchy)

### Genome Management - Register Don't Copy

**Decision:** Genomes stay where they are on disk, registered with `add_genomes()`. Only accepts protein FASTA files (.faa or .fasta). GenBank conversion removed and deferred to future enhancement.

**Workflow:**
```r
sack <- create_sack()
sack <- add_genomes(sack, "~/data/mags/*.faa")  # Registers .faa/.fasta files, doesn't copy
saveRDS(sack, "sack.rds")  # Standard R save
sack <- readRDS("sack.rds")  # Standard R load
```

**Rationale:** Don't force users to reorganize their existing data structures. Keep foundation simple - GenBank conversion can be added later (see ROADMAP.md).

### S7 Object vs Folder - Both Called "Sack"

**Decision:** Accept that "sack" refers to both the folder structure and the S7 object. Context makes it clear.

- **Folder sack:** Created by `initialize_potato_sack()` - directory with config, potatoes/, etc.
- **S7 PotatoSack:** In-memory object created by `create_sack()` from folder files

`create_sack()` constructs the S7 object from directory files. Use standard `saveRDS()`/`readRDS()` for persistence.

### No Examples in Documentation

**Decision:** Don't add `@examples` to roxygen documentation. Extra noise, not helpful yet.

### Removed Features (Deferred to Future)

**Removed during foundation cleanup to reduce complexity before core workflows exist:**

- **Legacy schema support:** All backward compatibility code removed
- **Database name customization:** Simplified to type-based keys (kofam, blast, hmm, patric)
- **Potato `type` field:** Removed (pathway vs. structural_complex) - not used for any logic
- **GenBank conversion:** Removed from `add_genomes()` - deferred to future (see ROADMAP)
- **Message logging system:** Removed (R/messages.R) - deferred to future (see ROADMAP)
- **Provenance export functions:** Removed - future feature, not needed yet
- **combine_sacks():** Removed - complex function for combining multiple annotation runs
- **Custom save/load functions:** Removed `save_potato_sack()`, `load_saved_sack()`, `inspect_saved_sack()` - use standard `saveRDS()`/`readRDS()` instead
- **RDS auto-detection in load:** Simplified - `create_sack()` only constructs from directory, doesn't load RDS files

---

## Design Rationale (Why These Decisions)

### Why DAG Structure Works Now

**The Breakthrough:** The bipartite DAG wasn't the problem - Excel was!

**Reality Check:**
- DAG structure is sound - matches how biochemists think
- Excel cell limits made it unusable (string syntax was a workaround)
- JSON removes all constraints → DAG becomes viable and elegant
- Visually graphable, handles KEGG module complexity, extensible

### Why Self-Contained Potatoes

**Old thinking:** Central gene registry (like Excel gene sheet)

**New thinking:** Each potato defines its own genes

**Benefits:**
- Portability - share a single JSON file, it works
- Flexibility - can override definitions when needed
- Version control - each potato evolves independently
- No external dependencies

**Trade-off:** Possible inconsistencies across potatoes (addressed via validation warnings)

### Why Gene Specificity Matters

**Biological intuition:** Not all genes are equally informative

**Example:**
- `gapA` (glyceraldehyde-3-P dehydrogenase) - found in 50+ pathways, low information
- `mlrA` (microcystin linearase) - only in microcystin degradation, high information

**Impact:** Transforms scoring from simplistic (5/8 genes) to contextual (high-confidence vs. low-confidence finds based on which genes)

**Implementation:** Pre-compute specificity = 1 / (number of potatoes containing gene)

### Why Compounds Are Edge Metadata

**Old v2 attempt:** Bipartite graph with compound nodes structurally required

**New v1:** Compounds as optional edge metadata

**Benefits:**
- Scoring logic stays enzyme-centric (simple)
- Preserves biological context for LLM reasoning
- Optional - can omit when not relevant

### Why Tool Availability Varies

**Reality:** Not all users have all tools (KEGG is expensive/restricted)

**Solution:** Graceful degradation with transparent reporting

**Implementation:**
- Potatoes specify what they need (`databases: {kofam: [...], blast: [...]}`)
- Runtime uses what's available in user's tools.json
- Reports which tools were used vs. requested
- Warns about genes only detectable via unavailable tools

**Impact:** Potatoes are shareable across users with different tool setups

### Why LLMs Are Assistants, Not Oracles

**What LLMs ARE good for:**
- Pre-computing specificity scores (reasoning baked in)
- Building potato drafts from KEGG modules
- Interpreting results in natural language
- Suggesting alternatives (with caveats)

**What LLMs ARE NOT good for:**
- Gene-to-function annotation (use established tools)
- Runtime scoring decisions (too slow, too variable)
- Replacing expert validation

**Approach:** Hybrid - LLM-enhanced but deterministic at runtime

---

## Complete Workflow Status (v0.9.3)

### Fully Implemented and Tested ✓

**Foundation (v0.5.x)**
- ✅ Project initialization and sack creation
- ✅ Config loading and validation
- ✅ Potato loading and validation
- ✅ GenomeFile S7 class for serialization-safe genome storage
- ✅ Standard R save/load (saveRDS/readRDS)

**Annotation Tools (v0.6.x - v0.7.x)**
- ✅ Kofam annotation with parallel execution (v0.6.0)
  - Per-gene thresholds from KEGG stored in results
- ✅ BLAST annotation with filtered databases (v0.6.1)
  - Global e-value and bitscore thresholds
- ✅ HMM annotation with profile extraction (v0.6.2)
  - Extracts trusted cutoffs (TC) from HMM profiles
  - TC values stored in results (`tc_threshold` column) for per-profile thresholding
  - Profiles without TC fall back to global e-value threshold
- ✅ All tools use consistent architecture:
  - Parallel workers execute shell commands
  - Sequential parsing with jakomics
  - Filtered databases (extract only needed terms)
  - File provenance (raw outputs + command logs)
  - Nested tibble results structure
  - Progress bars (progressr/cli)
- ✅ Conda path detection (v0.7.1)
  - `find_conda()` helper searches PATH, CONDA_EXE, common locations
  - Works out-of-box when conda is shell function (not in R's PATH)

**Scoring (v0.7.0, v0.9.3)**
- ✅ `score_pathways()` - Apply quality thresholds to annotation hits
- ✅ Per-gene thresholding:
  - Kofam: uses KEGG per-gene threshold (column: `threshold`)
  - HMM: uses per-profile TC when available, else e-value (column: `tc_threshold`)
  - BLAST: global e-value and bitscore thresholds
- ✅ Handle OR branches (alternative genes at same step)
- ✅ Calculate pathway completion fractions
- ✅ Determine presence/absence based on min_fraction threshold
- ✅ Store results in `sack@scores` tibble
- ✅ Multi-pathway network scoring (v0.9.3)
  - Each pathway scored independently
  - Results include `pathway` and `pathway_name` columns
  - Handles shared genes across pathways
- ✅ Essential-only scoring metrics (v0.9.3)
  - `steps_detected_essential`, `steps_total_essential`, `fraction_essential`, `present_essential`
  - Tracks required genes separately from all genes
  - Returns NA when no required genes defined
- ✅ Per-pathway thresholds included in results (`min_fraction` column)

**Visualization (v0.7.0 - v0.9.0)**
- ✅ `plot_potato()` - Dual-mode visualization system
  - **Static mode** (`interactive = FALSE`): ggraph network plots
    - Force-directed layouts for multi-pathway networks
    - Step-based layouts for single-pathway potatoes
    - Pathway convex hulls with ggforce
    - Publication-quality PDF/PNG export
  - **Interactive mode** (`interactive = TRUE`): visNetwork plots
    - Drag-and-drop node positioning
    - Full viewport display (100vh)
    - Export coordinates to JSON
    - Zoom, pan, highlight on hover
    - Navigation controls
  - **Curated layouts**: Supports embedded x,y coordinates in potato JSON
    - Dual coordinate systems: `x,y` (enzyme-only) and `x_compounds,y_compounds` (with compounds)
    - Hybrid approach: uses curated coords where available, layout algorithm for new nodes
    - Auto-detection of coordinate type when importing
  - **Pathway filtering**: `pathway` parameter to show single pathway from multi-pathway network
  - **Compound display**: `show_compounds` toggle for bipartite graphs
    - Deduplicated compounds (same compound appears once)
    - Compounds positioned between connected enzymes
    - Different shapes: triangles for compounds, circles for genes
- ✅ `plot_pathway_heatmap()` - ggplot2 tile heatmap across genomes
- ✅ `plot_genome_pathways()` - ggplot2 horizontal bars for one genome
  - Shows per-pathway threshold markers (red vertical bars)
  - Each pathway displays its actual min_fraction threshold
- ✅ `plot_pathway_summary()` - ggplot2 stacked bars (pathways per genome)
- ✅ `export_potato_dot()` - Export to graphviz format
- ✅ `update_potato_coordinates()` - Import visNetwork coordinates to potato JSON
- ✅ All plots use tidyverse/ggplot2 (no base graphics)
- ✅ `potato_theme()` - Consistent theming with transparent backgrounds

**Analysis Functions (v0.7.1, v0.9.3)**
- ✅ `summarize_missing_genes()` - Identify genes systematically missing across genomes
- ✅ `find_near_miss_pathways()` - Find pathways just below detection threshold
- ✅ `plot_near_miss_pathways()` - Visualize near-miss status with color coding
- ✅ `get_gene_results()` - Export gene-level annotation results (v0.9.3)
  - Returns tibble with all hits across tools (kofam, blast, hmm)
  - Includes `passed` column for threshold filtering
  - Kofam uses per-gene threshold, BLAST uses global e-value/bitscore, HMM uses TC or e-value
- ✅ `get_pathway_scores()` - Export pathway-level scores (v0.9.3)
  - Returns tibble with all scoring metrics
  - Includes `potato_hash` for version tracking
  - Shows both all-steps and essential-only metrics
  - Includes `min_fraction` threshold per pathway

**Multi-Pathway Networks (v0.9.0)**
- ✅ Multi-pathway JSON schema implemented
  - Global `nodes` array: gene definitions with detection methods
  - `pathways` object: pathway-specific attributes (step, required, marker, edges)
  - Pathway types: "variant" (alternatives) vs "independent" (different purposes)
- ✅ Gene-based graph structure
  - Same gene appears once, shared across pathways
  - Edges deduplicated automatically
  - No step-based node IDs (uses bare gene IDs)
- ✅ `print_potato()` handles multi-pathway networks
  - Prints each pathway with compact notation
  - Shows alternatives, markers, optional genes
- ✅ `validate_potato()` works with multi-pathway networks
  - Validates global nodes + pathway-specific attributes
  - Checks for cycles, missing references
- ✅ Potato hashing for version tracking
  - `potato_hash` column in annotation results
  - Only functional fields hashed (excludes metadata)
  - Tracks which potato version was used for annotation
- ✅ Example: `entner_doudoroff_network.json`
  - 4 ED pathway variants + gluconate transport + galactose catabolism
  - 23 unique genes, 6 pathways
  - Demonstrates pathway overlap and metabolic flexibility

### Ready for Testing
All features above have been implemented and basic testing done on:
- 22 genomes (marine isolates)
- 7 potatoes (glyoxylate cycle, entner-doudoroff, etc.)
- All three annotation tools
- Scoring with default thresholds
- All visualization types

### Not Yet Implemented
- Gene specificity weighting in scoring
- Marker gene emphasis in scoring
- Multi-pathway scoring (currently score single-pathway potatoes only)
- Threshold sensitivity analysis (high-priority - see ROADMAP Phase 3.4)
- LLM agents (builder, converter, analysis)
- Multi-sack comparisons
- Advanced DAG traversal (currently simple step-based, uses step counting)
- Genbank → FAA conversion (deferred from v0.5)

### Known Limitations
- **Scoring:** Simple step counting (fraction detected) rather than true DAG traversal
  - Doesn't verify connectivity from input to output
  - A pathway with 80% of genes but a broken middle step still scores 0.8
  - Good enough for most presence/absence calls, especially with incomplete MAGs
- **HMM TC values:** Only available when profiles have embedded TC lines
  - PFAM profiles have TC values
  - Custom HMM profiles (e.g., mlr) typically don't have TC
  - Falls back to global e-value threshold when TC absent

## Current Working Functions (v0.9.3)

### User Workflow Functions
- `initialize_potato_sack(path)` - Create project folder structure
- `create_sack(path = NULL)` - Construct PotatoSack S7 object from directory
- `validate_potato(potato)` - Validate potato structure and schema (⚠️ needs update for multi-pathway)
- `add_genomes(sack, path)` - Register .faa/.fasta files with sack
- `run_kofam(sack, potato_names, conda_env, workers, overwrite)` - Run kofam annotation
- Standard R: `saveRDS(sack, "file.rds")` and `readRDS("file.rds")`

### Potato Functions
- `load_potato(path)` - Load single potato JSON (multi-pathway networks supported)
- `load_potatoes(dir, tags = NULL)` - Load all potatoes from directory (filters active=false)
- `load_test_potato()` - Load example test potato
- `print_potato(potato, compact, show_compounds, show_ko, show_ec)` - Text-based pathway view
  - Handles both single-pathway and multi-pathway networks
  - Compact notation: `*` = marker, `^` = optional, `{n}` = complex, `(A|B)` = alternatives
  - Example: `[D-glucose-6-P] -> zwf -> (pgl^ | ybhE^) -> edd* -> eda*`
  - Multi-pathway: prints each pathway separately with its topology
- `get_enzyme_nodes(potato)` - Extract enzyme nodes
- `get_detection_terms(potato, database_name)` - Extract KO/blast/hmm terms
- `get_marker_genes(potato)` - Extract marker gene nodes
- `build_potato_graph(potato)` - Build igraph DAG (gene-based for multi-pathway)
- `build_bipartite_graph(potato)` - Build bipartite graph with compound nodes
- `print_validation(validation_result)` - Pretty print validation
- `compute_potato_hash(potato)` - Generate MD5 hash of functional fields
- `update_potato_coordinates(potato_path, coords_path, with_compounds)` - Import coordinates
  - Auto-detects if coordinates include compounds
  - Saves to x/y or x_compounds/y_compounds fields

### Annotation Functions
- `run_kofam(sack, potato_names, conda_env, workers, overwrite)` - Kofam annotation
- `create_kofam_hal(sack, potato_names)` - Create .hal file for kofam
- `run_blast(sack, potato_names, conda_env, workers, overwrite)` - BLAST annotation
- `create_blast_db(sack, potato_names)` - Extract BLAST reference sequences
- `run_hmm(sack, potato_names, conda_env, workers, overwrite)` - HMM annotation
- `create_hmm_profile(sack, potato_names)` - Extract HMM profiles by NAME

### Scoring Functions
- `score_pathways(sack, kofam_threshold, blast_evalue, blast_bitscore, hmm_evalue)` - Score all pathways

### Visualization Functions
- `plot_potato(potato, sack, genome_name, layout)` - Network plot with detection status
- `plot_pathway_heatmap(sack, cluster_rows, cluster_cols)` - Presence/absence heatmap
- `plot_genome_pathways(sack, genome_name, show_thresholds)` - Completion bars with per-pathway thresholds
- `plot_pathway_summary(sack)` - Pathways detected per genome
- `export_potato_dot(potato, file)` - Export to graphviz DOT format
- `potato_theme()` - Consistent theme for all plots (transparent backgrounds)

### Analysis Functions
- `summarize_missing_genes(sack, potato_name, min_genomes)` - Find systematically missing genes
- `find_near_miss_pathways(sack, buffer)` - Identify pathways just below threshold
- `plot_near_miss_pathways(sack, genome_name, buffer)` - Visualize near-miss status
- `get_gene_results(sack)` - Export gene-level annotation results with threshold pass/fail
- `get_pathway_scores(sack)` - Export pathway scores with essential metrics and potato_hash
- `get_node_status(potato, sack, genome_name)` - Get detection status for plotting (internal)

### Config Functions
- `load_potato_config(config_path = NULL)` - Load and validate config YAML
- `find_potato_sack(path = NULL)` - Search upward for potato_config.yaml

### Internal/Helper Functions
- `GenomeFile` - S7 class for serialization-safe genome metadata
- `jakomics_to_genome_file()` - Convert Python FILE objects to R S7 objects
- `find_conda()` - Locate conda executable (PATH, CONDA_EXE, common locations)

### Test Suite
- 91 passing tests (v0.9.3)
- Test coverage: 27% (up from 11.35%)
- Tests in: test-potato-class.R, test-potato-sack.R, test-save-load.R, test-config.R, test-export.R, test-scoring.R, test-multi-pathway.R

---

## Kofam Annotation Implementation (v0.6.0) ✓

Successfully implemented parallel kofam annotation workflow. Key architectural decisions:

### Serialization-Safe Genome Storage

**Problem:** Python objects from jakomics don't survive `saveRDS()`/`readRDS()` - they become invalid pointers.

**Solution:** Created `GenomeFile` S7 class (pure R) to store genome metadata:
```r
GenomeFile <- S7::new_class(
  properties = list(
    short_name = class_character,
    file_path = class_character,
    name = class_character,
    file_type = class_character,
    md5 = class_character
  )
)
```

**Workflow:**
1. `add_genomes()` uses jakomics to discover/validate files
2. Converts Python FILE objects → GenomeFile S7 objects via `jakomics_to_genome_file()`
3. Stores in `sack@genomes` slot (list of GenomeFile objects)
4. Safe to serialize with `saveRDS()` - no Python pointers

### Parallel Execution Architecture

**Pattern:** Workers execute commands only, sequential parsing after completion

**Why:** Python objects (jakomics modules, potato data) can't serialize to parallel workers

**Implementation:**
```r
# Worker function - NO Python dependencies
run_kofam_cmd <- function(genome_path, genome_name, hal_path, ko_list, conda_env) {
  cmd <- sprintf("conda run -n %s exec_annotation ...", conda_env, ...)
  output <- system(cmd, intern = TRUE)
  list(output = output, command = cmd)
}

# Parallel execution (furrr)
results <- furrr::future_map(seq_along(genome_paths), function(i) {
  run_kofam_cmd(genome_paths[i], genome_names[i], ...)
})

# Sequential parsing (purrr) - AFTER parallel work done
kofam_results <- purrr::map(raw_outputs, function(output_lines) {
  # Parse with jakomics (Python)
  hits <- list()
  for (line in output_lines) {
    clean_line <- sub("^[*\\s]+", "", line)  # Strip leading * or space
    hits[[length(hits) + 1]] <- kegg$KOFAM(clean_line, ...)
  }
  parsed <- kegg$parse_kofam_hits(hits)
  kofam_hits_to_tibble(parsed, potato_data)
})
```

### Conda Environment Support

**Problem:** Parallel workers can't access pre-activated conda environment

**Solution:** Added `conda_env` parameter throughout jakomics:
- Centralized in `jakomics.utilities.system_call()`
- If `conda_env` provided, wraps command: `conda run -n {env} {command}`
- Each worker independently activates environment for subprocess

**jakomics changes:**
```python
def system_call(call, echo=False, run=True, return_type='err', conda_env=None):
    if conda_env:
        call = f"conda run -n {conda_env} {call}"
    # ... rest of function
```

### File Provenance and Logging

**Timestamped directories:** `results/annotations/{timestamp}/`
- Timestamp created once per annotation session in `sack@metadata$annotation_session`
- All tools for same run share same timestamp directory

**Files saved per run:**
- `.hal` file - list of HMM profile paths (MD5 hash filename)
- `{genome}.kofam.txt` - raw kofam output for each genome
- `kofam.log` - TSV with columns: `genome`, `command` (full shell command used)

**Benefits:**
- Reproducibility - can re-run exact command
- Debugging - inspect raw tool outputs
- Provenance - know what was run when

### Results Data Structure

**Nested tibble pattern:**
```r
sack@results <- tibble(
  genome = character(),         # Genome short_name
  kofam = list(),              # List column: one tibble per genome
  blast = list(),              # Future: one tibble per genome
  hmm = list()                 # Future: one tibble per genome
)

# Each kofam[[i]] tibble has:
# potato, node_id, step, gene_id, ko, score, evalue, threshold, passed
```

**Usage:**
```r
# View all results
sack@results %>% unnest(cols = kofam)

# Filter to specific pathway
sack@results %>% 
  unnest(cols = kofam) %>% 
  filter(potato == "glyoxylate_cycle")
```

### Progress Bars

**Parallel (furrr):** Uses progressr package with cli styling
```r
progressr::handlers(progressr::handler_cli(...))
progressr::with_progress({
  p <- progressr::progressor(along = genome_paths)
  results <- furrr::future_map(..., function(i) {
    result <- run_cmd(...)
    p()  # Update progress
    result
  })
})
```

**Sequential (purrr):** Uses cli progress bar
```r
cli::cli_progress_bar("Running kofam", total = length(genome_paths))
results <- purrr::map(..., function(i) {
  cli::cli_progress_update()
  run_cmd(...)
})
cli::cli_progress_done()
```

### Messaging with cli Package

All user-facing messages use cli for consistent styling:
```r
cli::cli_alert_info("Preparing kofam annotation...")
cli::cli_alert_success("Created {.file {basename(hal_path)}}")
cli::cli_alert_warning("Overwriting existing results")
cli::cli_abort(c(
  "kofam results already exist",
  "i" = "Use {.code overwrite = TRUE} to replace"
))
```

### Configuration Defaults

**Pattern:** User params > config > hardcoded defaults

```r
# Get conda_env from config if not provided
if (is.null(conda_env)) {
  conda_env <- sack@config$annotation$conda_env
}

# Get workers from config if not provided
if (is.null(workers)) {
  workers <- sack@config$annotation$workers
  if (is.null(workers)) workers <- 1
}
```

**Config structure:**
```yaml
annotation:
  parallel: true
  workers: 4
  conda_env: potato
```

### Annotation vs Scoring Separation

**Critical design principle:** Annotation = data collection, Scoring = interpretation

**Annotation phase:**
- Collect ALL hits regardless of score/evalue thresholds
- Store raw scores, evalues, thresholds as metadata
- Do NOT include `passed` or similar boolean columns
- Do NOT filter results based on thresholds

**Scoring phase (to be implemented):**
- Apply thresholds in pathway context
- Consider gene specificity weighting
- Handle required vs optional genes
- Generate pathway-level confidence scores

**Why separate?**
- Users can adjust thresholds later without re-running annotation
- Different pathways may need different stringency
- Scoring can consider pathway context (specificity, completeness)

**Implementation notes:**
- **Kofam:** Removed `-T` flag from exec_annotation, removed `passed` column, keep all hits
- **BLAST (future):** Don't filter by evalue/bitscore/pident - return all hits
- **HMM (future):** Return all hits. Note: HMM profiles may have "trusted cutoff" (TC) lines - we'll handle these in scoring phase, not annotation phase

### Lessons for HMM and BLAST Implementation

Apply the same patterns:

**1. Worker function pattern:**
```r
run_hmm_cmd <- function(genome_path, genome_name, hmm_profile, conda_env) {
  cmd <- sprintf("conda run -n %s hmmsearch ...", conda_env, hmm_profile, genome_path)
  output <- system(cmd, intern = TRUE)
  list(output = output, command = cmd)
}
```

**2. Sequential parsing:**
```r
hmm_results <- purrr::map(raw_outputs, function(output_lines) {
  # Parse with jakomics.hmm
  hits <- hmm$parse_hmmsearch_output(output_lines)
  hmm_hits_to_tibble(hits, potato_data)
})
```

**3. File saving:**
```r
# In same timestamp directory as kofam
annotation_dir <- file.path(sack@sack_root, "results", "annotations", 
                           sack@metadata$annotation_session)

# Save raw outputs
output_file <- file.path(annotation_dir, paste0(genome_name, ".hmm.txt"))
writeLines(raw_output, output_file)

# Save command log
log_file <- file.path(annotation_dir, "hmm.log")
log_lines <- paste0(genome_names, "\t", commands)
writeLines(c("genome\tcommand", log_lines), log_file)
```

**4. Nested tibble in results:**
```r
sack@results$hmm <- hmm_results  # Add hmm column
sack@results$blast <- blast_results  # Add blast column
```

**5. HMM-specific considerations:**
- HMM profiles may be concatenated (multiple profiles in one file)
- Extract profile NAME from HMM file (not filename) for detection terms
- May need to create temporary concatenated HMM file like .hal for kofam

**6. BLAST-specific considerations:**
- May need to create BLAST database from reference sequences
- BLAST database can be pre-built or built on-the-fly
- Reference sequences stored in `databases.blast.files` (can be multiple)

## Phase 1 Priorities (MVP)

The goal is a working prototype: **one potato, one genome, no LLM features**.

### Immediate Next Steps: Minimal Viable Test

**Philosophy:** Prove the architecture with the absolute simplest case, then add complexity.

#### Step 0: Minimal Test Case Setup

Create the simplest possible test:
- **One genome** - single FAA file (~50-100 proteins, E. coli subset or test genome)
- **One potato** - 3 genes, linear pathway, no branches
- **One tool** - start with BLAST (simplest) or kofam (most realistic)

```
tests/fixtures/
├── genomes/
│   └── test_genome.faa           # Small test genome
├── potatoes/
│   └── test_pathway.json         # 3 genes: geneA -> geneB -> geneC
└── tools/
    └── test_tools.json           # Just BLAST or kofam configured
```

**Test pathway example:**
```json
{
  "id": "test_pathway",
  "name": "Test Pathway (Linear)",
  "tags": ["test"],
  "nodes": [
    {
      "id": "geneA",
      "step": 1,
      "nodes": ["geneA_1"],
      "type": "enzyme",
      "name": "Gene A",
      "databases": {"kofam": ["K00001"]},
      "required": true
    },
    {
      "id": "geneB",
      "step": 2,
      "nodes": ["geneB_2"],
      "type": "enzyme",
      "name": "Gene B",
      "databases": {"kofam": ["K00002"]},
      "required": true
    },
    {
      "id": "geneC",
      "step": 3,
      "nodes": ["geneC_3"],
      "type": "enzyme",
      "name": "Gene C",
      "databases": {"kofam": ["K00003"]},
      "required": true
    }
  ],
  "edges": [
    {"from": "geneA_1", "to": "geneB_2"},
    {"from": "geneB_2", "to": "geneC_3"}
  ],
  "scoring": {
    "min_fraction": 1.0
  }
}
```

#### Step 1: Python Backend (Minimal)

**File: `inst/python/potato_minimal.py`**

```python
from jakomics import kegg, blast, hmm
from jakomics.genome import GENOME
import json

class Potato:
    """Loads and represents a potato JSON file"""
    def __init__(self, json_path):
        with open(json_path) as f:
            self.data = json.load(f)
        self.id = self.data['id']
        self.nodes = self.data['nodes']
        self.edges = self.data['edges']
    
    def get_ko_list(self):
        """Extract all KO IDs for kofam searching"""
        kos = []
        for node in self.nodes:
            if 'databases' in node and 'kofam' in node['databases']:
                kos.extend(node['databases']['kofam'])
        return kos

class ToolRunner:
    """Runs annotation tools via jakomics"""
    def __init__(self, tools_config_path):
        with open(tools_config_path) as f:
            self.config = json.load(f)
    
    def run_kofam(self, faa_path, ko_list):
        """Run kofam via jakomics"""
        cfg = self.config['kofam']
        hits = kegg.run_kofam(
            faa_path=faa_path,
            hal_path=None,  # Will need to generate
            temp_dir="temp_ko",
            ko_list=cfg['ko_list'],
            score_as_ratio=False,
            echo=False
        )
        return kegg.parse_kofam_hits(hits)

def score_pathway_simple(potato, found_genes):
    """Dead simple scoring: count found/total"""
    total = len(potato.nodes)
    found = len([g for g in found_genes if g in [n['id'] for n in potato.nodes]])
    return {
        'pathway': potato.id,
        'found': found,
        'total': total,
        'present': found == total
    }
```

#### Step 2: R Wrapper (Minimal)

**File: `R/minimal.R`**

```r
#' Minimal test: one genome, one potato, one tool
#' @export
test_minimal <- function(faa_path, potato_path, tools_path) {
  
  # Load Python backend
  potato_py <- reticulate::import_from_path("potato_minimal", "inst/python")
  
  # Load potato
  potato <- potato_py$Potato(potato_path)
  
  # Load tools
  runner <- potato_py$ToolRunner(tools_path)
  
  # Run kofam
  ko_list <- potato$get_ko_list()
  results <- runner$run_kofam(faa_path, ko_list)
  
  # Score pathway (simplistic)
  found_genes <- names(results)  # Gene IDs that had hits
  score <- potato_py$score_pathway_simple(potato, found_genes)
  
  return(score)
}
```

#### Step 3: Manual Test

```r
# In R console
library(potato)

result <- test_minimal(
  faa_path = "tests/fixtures/genomes/test_genome.faa",
  potato_path = "tests/fixtures/potatoes/test_pathway.json",
  tools_path = "tests/fixtures/tools/test_tools.json"
)

print(result)
# Expected: list(pathway="test_pathway", found=3, total=3, present=TRUE)
```

#### Step 4: Incrementally Add Complexity

Once minimal case works:

1. ✅ **One genome, one potato, one tool** (START HERE)
2. Add: Multiple tools (kofam + blast)
3. Add: OR branches in pathway (ilvH | ilvM)
4. Add: Proper DAG traversal scoring
5. Add: Multiple potatoes
6. Add: Multiple genomes (parallelization)
7. Add: Genbank → FAA conversion
8. Add: Output file writing
9. Add: Gene specificity weighting
10. Add: LLM agents

---

### Full Implementation Steps (After Minimal Works)

1. **Define JSON schemas** formally (use JSON Schema spec)
   - potato.schema.json
   - tools.schema.json

2. **Python Potato class** (enhanced)
   ```python
   class Potato:
       def __init__(self, json_path):
           # Load, validate JSON
           # Build DAG representation (networkx?)
           
       def validate(self):
           # Check: no cycles, required nodes exist, etc.
           
       def get_genes(self):
           # Return list of gene definitions
   ```

3. **Python tool runners** (enhanced)
   - Abstract interface: `Tool.run(faa_path, gene_defs) -> results`
   - Support kofam, blast, hmmer, pfam

4. **Scoring engine** (proper DAG)
   - DAG traversal with networkx
   - Handle OR branches
   - Output: gene-level + pathway-level dataframes

### Testing Plan

Create 3 test potatoes manually:
1. **Simple linear** - test_pathway (3 genes, sequential) ← START HERE
2. **With OR branch** - leucine biosynthesis (ilvH | ilvM alternative)
3. **Complex** - nitrogen fixation (multiple gene clusters, optional genes)

Test genomes:
- Minimal test genome (50-100 genes with known KOs) ← START HERE
- Known E. coli genome (complete, well-annotated)
- Known MAG (70% complete, test missing genes)
- Negative control (no pathway genes present)

---

## LLM Agent Guidelines

### When Building the Agent Features (Phase 4)

**Builder Agent:**
- Use Claude API (Anthropic SDK)
- Prompt template: Include potato JSON schema in system prompt
- For KEGG modules: Parse API response carefully (nested parentheses)
- Always generate draft + ask user to review

**Converter Agent:**
- Read Excel carefully (pandas)
- Parse v1 syntax: `->` = sequential, `|` = OR, `+` = AND within step
- Handle edge cases (empty cells, malformed syntax)
- Generate validation report (what converted successfully, what needs review)

**Analysis Agent:**
- Input: scoring results + optional MAG completeness
- Output: Natural language interpretation
- Be conservative: Flag uncertainty, don't overclaim

**General LLM Principles:**
- Pre-compute where possible (specificity scores baked in)
- LLM for interpretation/assistance, not runtime decisions
- Always validate LLM outputs programmatically
- Graceful failure (if API down, core tool still works)

---

## Common Tasks

### Adding a New Potato

1. Create JSON file in `inst/potatoes/`
2. Define nodes (genes with detection methods)
3. Define edges (pathway structure)
4. Set scoring parameters
5. Validate: `potato <- load_potato("my_potato.json")`
6. Visualize: `plot_potato_dag(potato)`

### Adding a New Tool Type

1. Add to `tools.json` schema
2. Implement Python runner in `inst/python/tools.py`
3. Add R wrapper in `R/tools.R`
4. Update potato JSON schema (allow new tool type in gene definitions)
5. Test with known results

### Running Full Pipeline

```r
# Load package
library(potato)

# Configure tools
tools_config <- load_tools_config("tools.json")

# Load potatoes
potatoes <- load_potatoes(dir = "inst/potatoes/", tags = c("amino_acid"))

# Annotate genome
results <- annotate_genome(
  faa_path = "test_genome.faa",
  potatoes = potatoes,
  tools_config = tools_config
)

# View results
results$genes      # Gene-level hits
results$pathways   # Pathway-level confidence scores
```

---

## Important Context

### Biological Domain Knowledge

**MAGs (Metagenome-Assembled Genomes):**
- Typically 50-100% complete
- Real genes may be missing (not sequenced or not assembled)
- Completeness affects interpretation (not all absences are real)

**Annotation Tool Disagreement:**
- Tools disagree more for phylogenetically distant organisms
- Multiple tools = redundancy, confidence
- No single tool is ground truth

**Gene Specificity:**
- Central metabolism genes (gapA, eno, pgi) = ubiquitous, low information
- Pathway-specific genes (nifH, mlrA, pmoA) = rare, high information
- Weighting by specificity improves confidence scoring

**Functional Analogs:**
- Some enzymes have substrate promiscuity
- EC number helps (EC 1.1.1.* = oxidoreductases on CH-OH)
- BUT: substrate specificity matters, can't just swap any enzyme
- LLM suggestions must be flagged as unvalidated

### User Needs

**Primary users:**
- Microbial ecologists analyzing hundreds of MAGs
- Domain experts with specialized pathway knowledge
- Not bioinformaticians (need simple interface)

**Key use cases:**
1. "Does this microbe fix nitrogen?" (nifH pathway)
2. "Can this microbe degrade microcystin?" (mlrA pathway)
3. "What metabolic capabilities are in this community?" (all potatoes)
4. "I discovered a new pathway, can I add it?" (custom potato)

**Output requirements:**
- TSV files (compatible with Excel, R, Python)
- Both gene-level detail and pathway-level summary
- Confidence scores, not just binary present/absent
- Explanations (why high/low confidence?)

---

## Migration from GATOR v1

### Backward Compatibility

**Output format:**
- Keep TSV structure where possible
- Add new columns for confidence, tools_used, specificity
- User scripts expecting v1 output should mostly work

**Database conversion:**
- Converter agent will migrate `gator_db.xlsx` → potato JSONs
- One-time operation per user
- Validation report shows what needs manual review

**Workflow:**
```r
# Old (GATOR v1)
gator --in_dir mags/ --gator_db gator_db.xlsx

# New (POTATO v1)
library(potato)
convert_gator_db("gator_db.xlsx", output_dir = "potatoes/")
annotate_genomes(genome_dir = "mags/", potato_dir = "potatoes/")
```

---

## Pitfalls to Avoid

### 1. Over-Engineering Early
- **Don't:** Build full LLM agent system before core scoring works
- **Do:** Phase 1 MVP first, then add features

### 2. Ambiguous Specificity
- **Don't:** Compute specificity naively (gene name string matching)
- **Do:** Consider gene definition (same name, different KOs = different genes?)

### 3. Tool Execution Safety
- **Don't:** Pass unsanitized paths to system commands
- **Do:** Validate paths, use subprocess safely, handle errors gracefully

### 4. DAG Cycles
- **Don't:** Assume user-provided JSON is valid
- **Do:** Validate DAG has no cycles, required nodes exist, edges connect real nodes

### 5. LLM Over-Reliance
- **Don't:** Make core functionality depend on LLM availability
- **Do:** LLM is enhancement layer, tool works without it

### 6. Hardcoded Paths
- **Don't:** Assume database locations
- **Do:** Everything in tools.json, validate on load

### 7. Silent Failures
- **Don't:** Continue silently when tools fail
- **Do:** Warn loudly, report in output what succeeded vs. failed

---

## Open Questions to Resolve During Implementation

1. **Python DAG library:** networkx vs. custom? (networkx is heavier but battle-tested)

2. **Specificity computation:** Pre-compute on package load, or lazily on first use?

3. **Canonical genes:** Optional feature or skip entirely for v1?

4. **JSON Schema validation:** Strict (fail on unknown fields) or permissive (warn)?

5. **Parallel execution:** Multiprocessing for multiple genomes? (Yes, but keep it simple)

6. **Output format:** Add JSON output option, or TSV sufficient?

7. **Reaction metadata:** Required in potato JSON, or optional? (Optional for v1)

8. **Edge compound info:** Always include, or allow omission? (Allow omission)

---

## User Preferences (Jeff-Specific)

- **R preference:** Use R wherever possible, Python backend only when necessary
- **Naming:** "Potato" not "pathway" in user-facing functions (`load_potato()` not `load_pathway()`)
- **Validation:** Warn by default, not strict fail (users want flexibility)
- **Documentation:** Biological context in docs (explain *why*, not just *how*)
- **LLMs:** Embrace for assistance, but don't make core functionality depend on them

---

## Getting Started (For Future Sessions)

When you return to work on POTATO:

1. **Read ROADMAP.md** - Understand full plan
2. **Check this file** - Refresh on context and guidelines  
3. **Ask:** "What phase are we on? What's the next concrete task?"
4. **Implement incrementally** - Small PRs, test as you go
5. **Validate frequently** - Does output match expectations?

---

## Contact & Collaboration

**Primary Author:** Jeff Kimbrel  
**Original Tool:** GATOR (Python, `/Users/kimbrel1/Github/gator`)  
**New Tool:** POTATO (R package, this repo)

**Conda Environment:** `potato` (see environment.yaml)

**Related Projects:**
- **jakomics** - Python utility library (https://github.com/jeffkimbrel/jakomics)
- **gator** - Original v1 tool (reference for validation)

---

## Version History

- **v0.9.3** - Essential-only scoring metrics, result export functions (get_gene_results, get_pathway_scores), multi-line pathway hover text, test coverage 27% (91 tests)
- **v0.9.2** - Multi-pathway scoring implementation (score each pathway independently)
- **v0.9.1** - Input/output compound visualization and static plot improvements  
- **v0.9.0** - Dual-mode visualization (visNetwork + ggraph), curated coordinate system, multi-pathway networks
- **v0.7.0** - Scoring and visualization complete. All plots use ggplot2/ggraph
- **v0.6.2** - HMM annotation with profile extraction
- **v0.6.1** - BLAST annotation with filtered databases
- **v0.6.0** - Kofam annotation fully implemented with parallel execution, file provenance, GenomeFile serialization
- **v0.5.1** - Simplified save/load (standard R), removed orphan functions
- **v0.5.0** - Stripped to bare bones foundation, 35 tests passing
- **v0.4.0** - Bug fixes and UX improvements (pre-cleanup)
- **v1.0.0** - Target: gene specificity weighting + LLM agents + polished documentation

## Complete Example Workflow (v0.7.0)

```r
library(potato)

# 1. Initialize project
initialize_potato_sack("~/my_project")
sack <- create_sack("~/my_project")

# 2. Add genomes
sack <- add_genomes(sack, "~/genomes/*.faa")

# 3. Run all annotation tools
sack <- run_kofam(sack, conda_env = "potato", workers = 8)
sack <- run_blast(sack, conda_env = "potato", workers = 8)
sack <- run_hmm(sack, conda_env = "potato", workers = 8)

# 4. Score pathways
sack <- score_pathways(sack)

# 5. Save results
saveRDS(sack, "sack.rds")

# 6. Visualize results
plot_pathway_heatmap(sack)
plot_genome_pathways(sack, "genome_name")
plot_potato(sack@potatoes$glyoxylate_cycle, sack, "genome_name")
plot_pathway_summary(sack)

# 7. Export results
gene_results <- get_gene_results(sack)
pathway_scores <- get_pathway_scores(sack)
write.csv(gene_results, "gene_results.csv")
write.csv(pathway_scores, "pathway_scores.csv")
```

## What to Work on Next

**Immediate (v0.7.1):**
- Test with larger dataset (100+ genomes)
- Add gene specificity weighting to scoring
- Export functions for results (write TSV/CSV)

**Near-term (v0.8.0):**
- Marker gene emphasis in scoring
- Required vs optional gene handling
- Pathway-level confidence scores (not just fraction)

**Completed (v0.9.0 - v0.9.3):**
- ✅ Multi-pathway schema implementation
  - ✅ Schema design finalized (pathways field, variant/independent types)
  - ✅ ED network potato created (entner_doudoroff_network.json)
  - ✅ Active flag for deprecating old potatoes
  - ✅ Validation system updated
  - ✅ Scoring system updated (score pathways independently)
  - ✅ Visualization supports multi-pathway networks
  - ✅ Essential-only scoring metrics
  - ✅ Result export functions

**Future (v0.10.0+):**
- LLM agent for potato building
- LLM agent for result interpretation
- Advanced DAG traversal algorithms
- Multi-sack comparisons

---

## Current Potato Inventory

### Active Potatoes (loaded by default)

**Single-pathway potatoes:**
- `bhac.json` - Glyoxylate BHAC pathway
- `glyoxylate_cycle.json` - Glyoxylate cycle
- `microcystin_degradation.json` - Microcystin degradation (mlrABCD)
- `nitrogen_fixation_mo.json` - Nitrogen fixation (Mo-dependent)
- `nitrogen_fixation_v.json` - Nitrogen fixation (V-dependent)
- `test_glycolysis.json` - Test pathway (3 steps)

**Multi-pathway networks:**
- `entner_doudoroff_network.json` - ED pathway network (4 variants: classic, non_phosphorylative, semi_phosphorylative, semi_phosphorylative_alt)

### Inactive Potatoes (active: false, kept for reference)

- `entner_doudoroff_classic.json` - DEPRECATED: Consolidated into network
- `entner_doudoroff_np.json` - DEPRECATED: Consolidated into network
- `entner_doudoroff_semi_phos.json` - DEPRECATED: Consolidated into network
- `entner_doudoroff_semi_phos_alt.json` - DEPRECATED: Consolidated into network

**Note:** Inactive potatoes can still be loaded explicitly with `load_potato()`, but won't be included in `load_potatoes()` by default.

### Future Consolidations

Candidates for multi-pathway networks:
- TCA + glyoxylate shunt + reverse TCA
- Nitrogen fixation (already split into Mo/V variants, could be unified network)

---

Last updated: 2026-07-30
