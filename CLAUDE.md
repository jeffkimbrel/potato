# POTATO - Project Context for Claude

## Project Overview

**POTATO** (Pathway annOTATOr) is an R package for annotating MAGs (metagenome-assembled genomes) against curated metabolic pathways. It's the successor to GATOR (Genome annotATOR), redesigned around self-contained pathway definitions (potatoes) as DAG structures in JSON.

**Current Status:** v0.10.2-dev (2026-08-13) - V2 schema migration complete. Provenance tracking implemented for all annotation and scoring steps. Full reproducibility with command templates and tool versions.

**Key Innovation:** Each "potato" (pathway) is a self-contained JSON file defining:
- Genes with multi-tool detection methods (KEGG, PFAM, BLAST, HMM)
- Pathway structure as a DAG (handles complex branching)
- Confidence scoring via gene specificity weighting
- LLM-assisted building and interpretation

---

## Important Files

- **[ROADMAP.md](ROADMAP.md)** - Complete implementation plan, phases, file formats, future enhancements
- **[WORKFLOW.md](WORKFLOW.md)** - Function call architecture and workflow
- **[PROVENANCE.md](PROVENANCE.md)** - Provenance tracking system documentation
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

**CRITICAL: NEVER commit changes without explicit user approval**

When making changes to files:
1. Make the changes
2. Stage them with `git add`
3. Show the diff to the user
4. **WAIT for explicit approval to commit**
5. Only commit when user says to commit

**NEVER run `git commit` on your own.** The user needs to review diffs before committing.

## Technology Stack

**R** (tidyverse, igraph) + **Python** (via reticulate, jakomics utilities) + **Conda** (bioinformatics tools: kofamscan, hmmer3, blastp) + **JSON** (potato definitions)

**Architecture:** R interface → Python orchestration (jakomics) → bioinformatics tools

---

## Potato JSON Structure (V2 Schema)

**CRITICAL:** POTATO now uses V2 schema exclusively. All potatoes use this network structure.

Each potato is a self-contained network with genes, compounds, and pathways:

```json
{
  "schema_version": "v2",
  "id": "pathway_network_id",
  "name": "Human Readable Network Name",
  "source": "KEGG M00123, M00456 / custom",
  "tags": ["metabolism", "energy"],
  "notes": "Brief description",
  
  "genes": [
    {
      "id": "geneSymbol",
      "name": "enzyme name",
      "type": "enzyme",
      "databases": {
        "kofam": ["K00001"],              // KEGG Orthology IDs
        "blast": ["ref_seq_id"],          // BLAST reference sequence IDs
        "hmm": ["PF00001", "custom_hmm"]  // HMM profile NAMEs (includes PFAM)
      },
      "ec": ["1.1.1.1"],
      "reactions": ["R00001"],
      "notes": "Biological context"
    },
    {
      "id": "complexABC",
      "name": "multi-subunit protein complex",
      "type": "complex",
      "components": ["geneA", "geneB", "geneC"],
      "notes": "All components required for function"
    }
  ],
  
  "compounds": [
    {
      "id": "C00031",
      "name": "D-glucose",
      "kegg_compound": "C00031"
    }
  ],
  
  "pathways": {
    "pathway_id": {
      "name": "Pathway Display Name",
      "type": "variant",  // or "independent"
      "kegg_module": "M00123",
      "verified": false,  // ALWAYS false unless manually validated
      "notes": "Pathway-specific notes",
      
      "genes": ["geneA", "geneB", "chaperone1"],  // OPTIONAL: explicit gene list (includes non-catalytic genes)
      "input": ["C00031"],
      "output": ["C00002"],
      
      "edges": [
        {
          "from": "C00031",
          "to": "geneA",
          "required": true,
          "marker": true,
          "reaction": "R00001"
        },
        {
          "from": "geneA",
          "to": "C00118",
          "required": true,
          "marker": true,
          "reaction": "R00001"
        },
        {
          "from": "C00118",
          "to": "geneB",
          "required": false,
          "marker": false,
          "reaction": "R00002"
        }
      ],
      
      "scoring": {
        "min_fraction": 0.75,
        "max_gaps": 1,
        "marker_mode": "any"
      }
    }
  }
}
```

**Key V2 Schema Features:**
- `schema_version`: **Required** - Must be "v2" (schema is now frozen - no breaking changes)
- `genes`: Global gene definitions with detection methods OR complex definitions
- `type`: Gene types: "enzyme", "transport", "chaperone", or "complex"
- **Protein complexes:** Multi-component units where all components required
  - Complex entries use `"type": "complex"` and `"components": ["geneA", "geneB", "geneC"]`
  - Components must exist as separate genes with detection methods
  - Complex appears as single node in network plot
  - Scoring checks ALL components detected (not yet implemented)
  - Example: `ureA + ureB + ureC` from GATOR becomes ureABC complex
- `compounds`: Global compound definitions
- `pathways`: One or more pathways, each with its own topology
- `edges`: Connect genes/compounds (or complexes), carry `required`/`marker` attributes
- `verified`: **CRITICAL** - Per-pathway field, always `false` for new pathways. **NEVER set to true**. Only humans verify.
  - **NEVER EDIT VERIFIED PATHWAYS**: If `"verified": true`, DO NOT make edits without explicit user approval.

**Important Notes:**
- **No step numbers** - Genes are counted directly, no sequential steps
- **No node IDs** - Edges reference gene/compound IDs directly (or complex IDs for multi-component units)
- **required/marker on edges** - Not on genes, because same gene can have different roles in different pathways (different substrates, different pathway context)
- **Empty edges allowed** - Transport pathways and other special cases can have `edges: []` (valid as of v0.11.1)
- **Optional genes array** - Pathways can explicitly list genes (including non-catalytic genes like chaperones, regulators) via `"genes": ["geneA", "geneB"]`. If omitted, genes are inferred from edges. For complexes, list components (not complex ID).
- **Standard database types only:** `kofam`, `blast`, `hmm` (no custom names)
- **PFAM profiles:** Go in `hmm` field (PFAM is a type of HMM database)
- **Input/output:** Arrays of compound IDs (e.g., `["C00048"]` or `["C00031", "C00024"]`) that define core pathway boundaries for gap-based scoring. IDs must exist in the `compounds` array. For transporters, use location qualifiers in compound IDs ("NH4_external", "NH4_internal"). **Mandatory for gap-based scoring, optional for fraction-based.**
- **Pathway types:** 
  - `variant` - Alternative routes to same outcome (Mo-nif vs V-nif)
  - `independent` - Different functions, shared metabolic space (TCA vs glyoxylate shunt)
- **Scoring methods:**
  - `min_fraction`: Fraction-based scoring (detected genes / total genes ≥ threshold)
  - `max_gaps`: Gap-based scoring (# missing genes in best path from input→output ≤ threshold, requires input/output)
  - Both can be present - pathway passes if EITHER method succeeds
  - Neither present - defaults to `min_fraction: 0.67`

**Schema Format Notes (v0.11.0+):**
- **Current format (v0.11.0+):** `genes` array, simplified `input`/`output` as arrays of compound IDs
- **Legacy formats may exist:** Older V2 potatoes used `nodes` array and verbose input/output objects
- When updating old potatoes, convert to current format for consistency

---

## Gap-Based Scoring (v0.11.0+)

**Motivation:** Fraction-based scoring has limitations:
- Float precision issues (0.67 threshold fails 2/3 genes at 0.6666...)
- Ignores pathway topology (80% genes with broken middle step still scores 0.8)
- Doesn't distinguish complete vs incomplete routes in multi-pathway networks

**Gap-based scoring** addresses this by checking DAG connectivity: "Can you traverse from input→output with ≤N missing genes?"

### How It Works

1. **Define core pathway boundaries** with `input` and `output`:
```json
"input": ["C00048"],      // Array of compound IDs from compounds array
"output": ["C00631"]      // Single or multiple compounds
```

Input/output are **arrays of compound IDs** (not compound objects). The IDs must match entries in the potato's `compounds` array. Multiple inputs/outputs are supported (e.g., `["C00024", "C00036"]` for TCA cycle).

2. **Set max allowed gaps** in scoring block:
```json
"scoring": {
  "max_gaps": 1,        // Allow 1 missing gene in best path
  "min_fraction": 0.67  // Optional: also run fraction-based scoring
}
```

3. **Scoring logic:**
   - Find all paths from `input` → `output` in the DAG
   - For each path, count missing (non-detected) genes
   - If any path has ≤ `max_gaps` missing genes → pathway PRESENT
   - Ignores `required`/`marker` attributes (purely topological)

### Upstream Context

Pathways can have edges upstream of `input` or downstream of `output` for biological context:

```json
"edges": [
  // Upstream context (not scored)
  {"from": "C00031", "to": "glk", "required": true, "marker": false},
  {"from": "glk", "to": "C00668", "required": true, "marker": false},
  
  // CORE PATHWAY: C00668 (input) → C00022 (output)
  {"from": "C00668", "to": "zwf", "required": true, "marker": false},
  // ... rest of core pathway
]
```

Gap scoring only counts genes between `input` and `output`. Upstream/downstream genes are documentation only.

### Dual Scoring

If both `min_fraction` and `max_gaps` are present:
- Calculate both metrics independently
- Pathway passes if **EITHER** method succeeds
- Results include separate columns: `present_fraction`, `present_gaps`

### When to Use

- **Gap-based (`max_gaps`):** Pathways where connectivity matters (linear, branched)
- **Fraction-based (`min_fraction`):** Pathways where gene count matters more than topology
- **Both:** Maximum flexibility - pathway passes by either criterion

### Cyclic Pathways

For true metabolic cycles (TCA, BHAC), gap-based scoring requires defining artificial start/end points via `input`/`output`. Alternatively, use fraction-based scoring only for cycles. (Virtual edges for cycle documentation: future enhancement, see ROADMAP.md)

---

## V2 Schema Design Philosophy

**All potatoes use V2 schema** - Related pathways that share metabolic context are grouped into network potatoes with multiple sub-pathways. This provides biological context that isolated pathways cannot capture.

### When to Group Multiple Pathways

**Group pathways together when:**
- Pathways are **alternative routes** to the same biological outcome (e.g., Mo vs V nitrogenase)
- Pathways share **metabolic intermediates** and spatial context (e.g., ED pathway variants)
- Finding one variant **explains the absence** of another (e.g., "ED classic absent because ED semi-phos present")
- Pathways **overlay spatially** but serve different purposes (e.g., TCA + glyoxylate shunt)

**Keep as separate potato files when:**
- Pathways are functionally unrelated
- No shared genes or metabolic context
- Independent detection is more informative

**Schema structure:** See "Potato JSON Structure (V2 Schema)" section above for complete example.

### Key Principles

**1. Genes Defined Once**
- Global `genes` array: Detection methods, EC numbers, enzyme names
- No step numbers, no required/marker on genes
- Same gene can have different roles in different pathways

**2. Pathway Types**
- `"type": "variant"` - Alternative routes to same outcome (Mo vs V nitrogenase, ED variants)
- `"type": "independent"` - Different purpose, shares metabolic space (TCA vs glyoxylate shunt)

**3. Edges Carry Context**
- `required`, `marker` on edges (not genes)
- Edges connect genes and compounds
- Same gene can be required/marker in one pathway, optional/non-marker in another
- **Why pathway-specific:** Same gene can act on different substrates in different pathways
  - Example: bifunctional enzymes have different roles depending on metabolic context
  - Edge-level attributes allow this flexibility

**4. No Step Numbers**
- V2 scoring counts genes detected, not sequential steps
- Simpler than V1, still effective for presence/absence

**5. Shared Genes**
- Gene `gnaD` appears in 3 ED pathways
- Detecting `gnaD` once satisfies all 3 pathways
- Each pathway evaluates completion independently

**6. Transport Pathways with Empty Edges**
- Transport pathways can have empty `edges: []` arrays
- Valid pattern: no internal reactions, just substrate relocation
- Example: gluconate_transport moves compound across membrane without enzymatic steps
- Validation handles this correctly (fixed in v0.11.1)

**7. Per-Pathway Verified Flag**
- `verified` is per-pathway field, not per-potato
- Each pathway in multi-pathway network has independent verification status
- Default: `"verified": false` for all new pathways
- Only humans set to true after manual validation

**8. Protein Complexes**
- Multi-component complexes where all subunits required for function
- Complex defined in `genes` array with `"type": "complex"` and `"components": ["subA", "subB", "subC"]`
- Each component is a separate gene entry with detection methods
- **Edges reference complex ID** (appears as single node in network plot)
- Pathway's `genes` array lists components (not complex ID) for scoring
- **GATOR translation:** `geneA + geneB + geneC` becomes a complex; `geneA | geneB` remains alternatives (separate edges)
- **Alternative complexes:** Multiple complexes can be alternatives (e.g., `ureABC` OR `ureAB_ureC` complex)
- Validation checks all components exist as genes

**Why complexes in edges (not individual components)?**

With alternative complexes like `ureA + ureB + ureC | ureAB + ureC`, individual gene "required" flags are ambiguous:
- ureA required? Only if taking ureABC route
- ureB required? Only if taking ureABC route
- ureC required? Always (shared between alternatives)
- ureAB required? Only if taking ureAB_ureC route

By putting complexes in edges:
```json
"edges": [
  {"from": "C00086", "to": "ureABC", "required": true, "marker": true},
  {"from": "C00086", "to": "ureAB_ureC", "required": true, "marker": true}
]
```

The logic is unambiguous: "Either ureABC required OR ureAB_ureC required (alternatives)."

**Scoring resolves complexes:**
1. Check if ALL components of ureABC detected → ureABC "present"
2. Check if ALL components of ureAB_ureC detected → ureAB_ureC "present"
3. Pathway passes if any required alternative is complete

This avoids conditional requirements at gene level and matches biological reality (complex acts as functional unit).

### Scoring V2 Potatoes

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

## Complete Workflow Status (v0.10.2)

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
- ✅ `update_potato_coordinates()` - Import visNetwork coordinates to potato JSON
- ✅ All plots use tidyverse/ggplot2 (no base graphics)
- ✅ `potato_theme()` - Consistent theming with transparent backgrounds

**Analysis Functions (v0.7.1, v0.9.3, v0.10.1.9005)**
- ✅ `summarize_missing_genes()` - Identify genes systematically missing across genomes
  - Distinguishes three states: no hits, below threshold, passing threshold
  - Uses internal helper functions for threshold checks
- ✅ `find_near_miss_pathways()` - Find pathways just below detection threshold
  - Returns filtered tibble of near-misses with distance from threshold
- ✅ `inspect_gene_thresholds()` - Detailed view of annotation hits vs thresholds
  - Shows score, threshold, margin for each hit
  - Helps debug threshold issues and identify borderline genes
- ✅ `get_gene_results()` - Export gene-level annotation results (v0.9.3)
  - Returns tibble with all hits across tools (kofam, blast, hmm)
  - Includes `passed` column for threshold filtering
  - Kofam uses per-gene threshold, BLAST uses global e-value/bitscore, HMM uses TC or e-value
- ✅ `get_pathway_scores()` - Export pathway-level scores (v0.9.3)
  - Returns tibble with all scoring metrics
  - Includes `potato_hash` for version tracking
  - Shows both all-steps and essential-only metrics
  - Includes `min_fraction` threshold per pathway

**Provenance Reporting Functions (v0.10.2)**
- ✅ Two-tier function pattern for qmd/RMarkdown rendering
  - `summarize_*()` functions display execution messages from provenance
  - `get_*_details()` functions return detailed results with plots
- ✅ Message storage in provenance
  - Annotation functions collect messages during execution
  - Messages stored as tibbles in `sack@provenance$*$messages`
  - Enables rendering in eval=FALSE chunks
- ✅ `summarize_add_genomes()` - Display genome addition history
  - Returns list with $summary (per-call statistics), $messages (type, message tibble), $status
  - Reconstructs messages from provenance data
- ✅ `get_genome_details()` - Genome QC metrics
  - Returns $summary (file paths, sizes, protein counts, added_in_call), $plot (protein count bar chart)
- ✅ `summarize_kofam()` - Display kofam annotation messages
  - Returns $summary (timestamp, n_genomes, n_potatoes, n_kos, tool_version), $messages, $status
  - Backward compatible for sacks created before message tracking
- ✅ `get_kofam_details()` - Kofam hit statistics
  - Returns $summary (per-genome hits), $per_potato (aggregated stats), $plot (hits per genome)
- ✅ `summarize_blast()` - Display BLAST annotation messages
  - Returns $summary (timestamp, n_genomes, n_subjects), $messages, $status
- ✅ `get_blast_details()` - BLAST hit statistics
  - Returns $summary, $per_potato, $plot
- ✅ `summarize_hmm()` - Display HMM annotation messages
  - Returns $summary (timestamp, n_genomes, n_profiles), $messages, $status
- ✅ `get_hmm_details()` - HMM hit statistics
  - Returns $summary, $per_potato, $plot
- ✅ All summarize functions default to `verbose = FALSE` (return data without printing)

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

**Current status:** Full annotation pipeline working (kofam/BLAST/HMM), V2 schema complete, visualization, provenance tracking. See ROADMAP.md for planned features.

## Package Functions

**Key workflow:** Use `devtools::load_all()` when testing (not `library(potato)`) to load latest code. All functions are documented with roxygen2 - see R package documentation for details.

**Main workflow functions:** `initialize_potato_sack()`, `create_sack()`, `add_genomes()`, `run_kofam()`, `run_blast()`, `run_hmm()`, `score_pathways()`, visualization functions (`plot_potato()`, `view_pathway_detail()`, `plot_pathway_heatmap()`), analysis functions (`get_gene_results()`, `get_pathway_scores()`).

---

## Key Architectural Patterns

### Annotation Tool Implementation

**Key patterns used across all annotation tools (kofam, BLAST, HMM):**

1. **Serialization-safe storage:** `GenomeFile` S7 class stores genome metadata (no Python pointers survive `saveRDS()`)
2. **Parallel execution pattern:** Workers execute shell commands only, sequential parsing after completion (Python objects can't serialize to parallel workers)
3. **Conda environment support:** `conda_env` parameter wraps commands with `conda run -n {env}`
4. **File provenance:** Timestamped directories with raw outputs + command logs for reproducibility
5. **Nested tibble results:** `sack@results` tibble with list columns (kofam, blast, hmm) - one tibble per genome
6. **Annotation vs Scoring separation:** Annotation collects ALL hits (no filtering), scoring applies thresholds in pathway context
7. **Progress bars:** progressr (parallel) + cli (sequential)
8. **User messages:** cli package for consistent styling

**Critical: Annotation = data collection, Scoring = interpretation.** Store raw scores/evalues, apply thresholds later so users can adjust without re-running.

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

## Current Potato Inventory

Example potatoes in `inst/potatoes/`: glyoxylate cycle, BHAC, microcystin degradation, nitrogen fixation (Mo/V variants), Entner-Doudoroff network (multi-pathway example with 6 pathways). Use `Glob` to discover available potatoes.
