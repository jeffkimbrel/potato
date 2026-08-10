# POTATO v1 - Development Roadmap

**Current Version:** v0.9.3 (2026-08-10)

**Status:** Multi-pathway networks fully implemented and production-ready. ED network verified and tested. Enhanced visualization and print functions. Test coverage at 27%.

**What works:** Annotate genomes with kofam/BLAST/HMM → score pathways (all + essential genes) → interactive/static visualization with curated layouts → export results as tibbles → analyze near-misses → print pathways with compound details

**Recent additions (2026-08-10):**
- Compound name normalization (merge "A + B" and "B + A")
- Enhanced `print_potato()` with `show_databases` parameter and angle-bracket compound display
- Transporter pathway support (empty edges arrays)
- Verification system for multi-pathway networks (per-pathway `verified` field)
- Entner-Doudoroff network completed with 6 pathways (4 verified, 2 independent)

**What's next:** Gene specificity weighting, threshold sensitivity analysis, annotation workflow optimization

---

## Vision

POTATO (Pathway annOTATOr) is the successor to GATOR, designed to rapidly annotate collections of genomes (MAGs) against curated metabolic pathways. The core innovation is **self-contained potato files** (pathway definitions as DAG structures in JSON) with **LLM-assisted analysis** for handling MAG incompleteness and functional analogs.

### Key Design Principles

1. **Self-contained potatoes** - each pathway is an independent, portable unit
2. **DAG-based pathway logic** - intuitive, graphable, handles complex branching
3. **Tool-agnostic detection** - works with whatever annotation tools the user has
4. **Confidence-aware scoring** - interprets results in context of gene specificity and MAG completeness
5. **LLM augmentation** - assists in database building and result interpretation

---

## Architecture

### Technology Stack

- **R package** - primary user interface
- **Python backend** (via reticulate) - tool orchestration, file I/O, LLM integration
- **Conda environment** - bioinformatics tools (kofamscan, hmmer3, blast)
- **JSON-based data** - potato definitions, tool configurations

### Core Components

```
potato/
├── R/                      # R package code
│   ├── io.R               # File reading, potato loading
│   ├── tools.R            # Tool execution wrappers
│   ├── scoring.R          # Pathway scoring logic
│   └── agents.R           # LLM agent interfaces
├── inst/
│   ├── python/            # Python backend
│   │   ├── tools.py      # Tool execution
│   │   ├── scoring.py    # DAG traversal, scoring
│   │   └── agents.py     # LLM interactions
│   ├── potatoes/         # Potato JSON files
│   └── extdata/
│       └── tools.json    # User tool configuration
├── tools.json            # User-specific tool paths
└── environment.yaml      # Conda environment spec
```

---

## File Formats

### 1. Potato JSON Structure (v0.8.0 - v0.9.0)

**Two schema versions:**
1. **Single-pathway potatoes** (current, v0.8.0) - One pathway per file
2. **Multi-pathway networks** (NEW, v0.9.0) - Related pathways in one file

#### Single-Pathway Schema (v0.8.0)

```json
{
  "id": "microcystin_degradation",
  "name": "Microcystin Degradation",
  "source": "Literature (Bourne et al. 1996)",
  "active": true,
  "verified": false,
  "tags": ["toxin_degradation"],
  "notes": "mlrABCD gene cluster",
  
  "input": {
    "compound": "microcystin-LR (cyclic)",
    "kegg_compound": "C05371",
    "targets": ["mlrA_1"]
  },
  
  "output": {
    "compound": "ADDA + small peptides",
    "sources": ["mlrC_3"]
  },
  
  "nodes": [
    {
      "id": "mlrA",
      "step": 1,
      "nodes": ["mlrA_1"],
      "type": "enzyme",
      "name": "microcystinase",
      "databases": {
        "kofam": ["K20071"],
        "hmm": ["mlrA_aligned"],
        "blast": ["mlrA_AF411068_partial", "QVQ24103.1_mlrA"]
      },
      "ec": ["3.4.99.-"],
      "required": true,
      "marker": true,
      "notes": "Diagnostic enzyme - cleaves cyclic MC-LR"
    },
    {
      "id": "mlrB",
      "step": 2,
      "nodes": ["mlrB_2"],
      "type": "enzyme",
      "name": "serine hydrolase",
      "databases": {
        "hmm": ["mlrB_aligned"],
        "blast": ["mlrB_AF411069_partial"]
      },
      "required": true,
      "marker": false
    }
  ],
  
  "edges": [
    {"from": "mlrA_1", "to": "mlrB_2", "compound": "linear microcystin-LR"}
  ],
  
  "scoring": {
    "min_fraction": 0.75,
    "marker_mode": "any"
  }
}
```

#### Multi-Pathway Network Schema (v0.9.0 - NEW)

For related pathways sharing metabolic context:

```json
{
  "id": "entner_doudoroff_network",
  "name": "Entner-Doudoroff Pathway Network",
  "source": "KEGG M00006, M00309, M00308, M00633",
  "active": true,
  "verified": false,
  "tags": ["metabolism", "carbohydrate"],
  "notes": "Four ED variants with different phosphorylation strategies",
  
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
        "edd": {"step": 4, "required": true, "marker": true},
        "eda": {"step": 5, "required": true, "marker": true}
      },
      
      "edges": [
        {"from": "glk", "to": "zwf", "compound": "D-glucose-6P"},
        {"from": "zwf", "to": "edd", "compound": "6-phosphogluconate"},
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
      "nodes": {
        "gnaD": {"step": 1, "required": true, "marker": true},
        "kdgA": {"step": 2, "required": true, "marker": true}
      },
      "edges": [
        {"from": "gnaD", "to": "kdgA", "compound": "KDG"}
      ],
      "scoring": {"min_fraction": 0.75}
    }
  }
}
```

**Key Features:**
- **Active flag**: `"active": false` deprecates old potatoes (kept for reference, not loaded by default)
- **Genes defined once**: Global `nodes` array has detection methods
- **Pathway-specific attributes**: Each pathway defines `step`, `required`, `marker` for its own context
- **Pathway types**: `"variant"` (alternative routes) vs `"independent"` (different purposes, shared space)
- **Standard database types**: `kofam`, `blast`, `hmm` only (no custom names)
- **Compounds on edges**: Optional metabolite information
- **Per-pathway scoring**: Each pathway scored independently

**Threshold Philosophy (v0.8.0):**
- **Kofam**: Uses per-gene thresholds from KEGG automatically
- **HMM**: Uses per-profile TC (trusted cutoff) when available, else global e-value
- **BLAST**: Uses global e-value and bitscore thresholds
- **Optional per-gene overrides**: `thresholds` field in nodes (use sparingly)

### 2. Tool Configuration (`tools.json`)

User-specific paths to annotation databases and **default thresholds**. Lives in project root or `~/.potato/tools.json`.

```json
{
  "kofam": {
    "enabled": true,
    "ko_profiles_dir": "/Users/kimbrel1/databases/kofam/profiles/",
    "ko_list": "/Users/kimbrel1/databases/kofam/ko_list",
    "executable": "exec_annotation",
    "default_thresholds": {
      "score": 1.0,
      "evalue": 1e-10
    }
  },
  "hmmer": {
    "enabled": true,
    "executable": "hmmsearch",
    "default_thresholds": {
      "evalue": 1e-15,
      "bitscore": 50,
      "domain_evalue": 1e-10
    }
  },
  "blast": {
    "enabled": true,
    "executable": "blastp",
    "default_thresholds": {
      "evalue": 1e-15,
      "bitscore": 100,
      "pident": 30,
      "qcovs": 50
    }
  },
  "pfam": {
    "enabled": false,
    "hmm_database": "/path/to/Pfam-A.hmm",
    "default_thresholds": {
      "evalue": 1e-10
    }
  },
  "custom_hmm": {
    "enabled": true,
    "hmm_dir": "/Users/kimbrel1/databases/custom_hmms/",
    "default_thresholds": {
      "evalue": 1e-15
    }
  }
}
```

**Threshold Precedence (most specific wins):**
1. Per-gene threshold in potato JSON (`nodes[i].thresholds.blast_evalue`)
2. Default threshold in tools.json (`blast.default_thresholds.evalue`)
3. Built-in fallback (hardcoded in potato code)

**Behavior:**
- If a potato specifies a disabled tool, continue with available tools + warn
- Track which tools were used vs. requested in output
- Flag genes only detectable via unavailable tools

### 3. Gene Registry (Optional)

For enforcing consistency across potatoes. Lives at `inst/canonical_genes.json`.

```json
{
  "ilvB": {
    "name": "acetolactate synthase catalytic subunit",
    "ko": ["K01652", "K01653"],
    "pfam": ["PF00205"],
    "ec": ["2.2.1.6"]
  }
}
```

**Usage:**
- Potatoes can reference `"canonical": "ilvB"` instead of redefining
- Validation warns if potato's definition conflicts with canonical
- Supports flexible overrides when needed

---

## Implementation Phases

### Phase 0: Foundation ✓ (Complete)

**Status:** v0.5.0 - v0.5.1 completed foundational cleanup

- [x] R package skeleton
- [x] Reticulate bridge to Python
- [x] Conda environment with bioinformatics tools
- [x] S7 classes (Potato, PotatoSack, GenomeFile)
- [x] Core workflow functions (initialize, create, validate, add_genomes)
- [x] Config loading and validation
- [x] Potato loading and validation
- [x] Test suite (35 tests passing)
- [x] Standard save/load (saveRDS/readRDS)
- [x] Serialization-safe genome storage (GenomeFile S7 class)
- [x] 7 example potatoes with new schema
- [ ] Define potato JSON schema formally (JSON Schema file) - deferred
- [ ] Define tools.json schema - deferred

### Phase 0.5: Kofam Annotation Implementation ✓ (COMPLETE)

**Status:** v0.6.0 - Kofam annotation fully working

**Completed:**
- [x] Kofam execution with parallel workers via furrr
- [x] Conda environment support (centralized in jakomics)
- [x] File provenance (raw outputs + command logs in timestamped directories)
- [x] GenomeFile S7 class for serialization-safe genome storage
- [x] Nested tibble results structure (one row per genome, list column per tool)
- [x] Progress bars (cli for sequential, progressr for parallel)
- [x] All messaging via cli package
- [x] Config defaults (conda_env and workers from config)
- [x] .hal file creation for kofam profiles
- [x] Sequential parsing after parallel execution
- [x] Overwrite protection with user confirmation

**Architecture proven:**
- Workers execute shell commands only (no Python serialization issues)
- Sequential post-processing handles parsing with jakomics
- Python FILE objects converted to R S7 GenomeFile objects for safe serialization
- All annotation tools will follow same pattern

**Next:** Implement run_hmm() and run_blast() using same architecture

---

### Phase 1: Core Functionality (MVP)

**Goal:** Run annotation tools against genomes, store results.

#### 1.1 Data Structures ✓
- [x] R S7 classes (Potato, PotatoSack, GenomeFile) - fully implemented
- [x] Python integration via reticulate
- [x] Config loading and validation
- [x] Potato loading and validation
- [x] DAG validation (no cycles, required nodes exist, etc.)

#### 1.2 Tool Execution ✓ (COMPLETE)
- [x] **Kofam** - R/run_kofam.R (v0.6.0)
  - [x] Parallel execution with furrr
  - [x] Conda environment support
  - [x] File provenance (raw outputs + logs)
  - [x] Nested tibble results
  - [x] Progress bars
  
- [x] **HMM** - R/run_hmm.R (v0.6.2, enhanced v0.7.1)
  - [x] Parallel execution following kofam pattern
  - [x] Extract HMM profile NAMEs from files
  - [x] Extract trusted cutoffs (TC) from profiles (v0.7.1)
  - [x] TC values stored in results for per-profile thresholding (v0.7.1)
  - [x] Concatenated profile support
  - [x] File provenance and command logging
  
- [x] **BLAST** - R/run_blast.R (v0.6.1)
  - [x] Parallel execution following kofam pattern
  - [x] Create filtered BLAST databases from reference sequences
  - [x] Extract sequences from configured files
  - [x] File provenance and command logging

- [x] **Conda path detection** (v0.7.1)
  - [x] `find_conda()` helper searches PATH, CONDA_EXE, common locations
  - [x] Applied to all three annotation tools
  - [x] Works when conda is shell function (not in R's PATH)

- [ ] R `annotate_all()` function (FUTURE)
  - Wrapper that runs all configured tools (kofam, hmm, blast)
  - Uses same timestamp directory for all tools
  - Returns sack with all tool results populated

#### 1.3 Scoring Engine ✓ (COMPLETE - v0.7.0)
- [x] Pathway scoring with quality thresholds
  - [x] Handle OR branches (alternative genes at same step)
  - [x] Calculate min_fraction thresholds (per-potato)
  - [x] Per-gene thresholding:
    - Kofam: uses KEGG per-gene threshold
    - HMM: uses per-profile TC when available, else e-value
    - BLAST: global e-value and bitscore thresholds
- [x] Results stored in `sack@scores` tibble
- [ ] Output formats (FUTURE):
  - [ ] Gene-level: `genome_gator.tsv` (gene, product, tool, score, locus_tag)
  - [ ] Pathway-level: `genome_potato.tsv` (pathway, present, steps, confidence)

**Note:** Current scoring is simple step counting (fraction detected), not true DAG traversal. Doesn't verify connectivity from input to output. Sufficient for most presence/absence calls.

#### 1.4 Visualization ✓ (COMPLETE - v0.7.0, enhanced v0.7.1)
- [x] `plot_potato()` - Network visualization of pathway DAG
  - [x] Accepts potato path or object
  - [x] Optional genome detection overlay
  - [x] Fixed node sizing and labeling (v0.7.1)
  - [x] Support for input/output compound nodes (v0.7.1)
- [x] `plot_pathway_heatmap()` - Presence/absence across genomes
- [x] `plot_genome_pathways()` - Completion bars for single genome
  - [x] Per-pathway threshold markers (v0.7.1)
- [x] `plot_pathway_summary()` - Stacked bars, pathways per genome
- [x] `print_potato()` - Text-based pathway view (v0.7.1)
  - [x] Compact notation: `*` = marker, `^` = optional, `{n}` = complex, `(A|B)` = alternatives
  - [x] Shows input/output compounds
  - [x] Optional EC numbers and KO IDs
  - [x] Example: `[D-glucose-6-P] -> zwf -> (pgl^ | ybhE^) -> edd* -> eda*`
- [x] `potato_theme()` - Consistent theming with transparent backgrounds (v0.7.1)
- [x] All visualizations use ggplot2/ggraph

#### 1.5 Analysis Functions ✓ (COMPLETE - v0.7.1)

**Goal:** Diagnostic tools for understanding why pathways are/aren't detected and identifying threshold issues.

- [x] `summarize_missing_genes(sack, potato_name, min_genomes)` - Identify systematically missing genes
  - Shows which genes are missing most often across genomes
  - Returns: tibble with `gene_id`, `times_missing`, `fraction_missing`
  - Use case: "Gene X is missing in 90% of genomes → database/threshold issue?"
  
  ```r
  missing <- summarize_missing_genes(sack)
  missing %>% filter(fraction_missing > 0.8)  # genes missing in >80% of genomes
  ```

- [x] `find_near_miss_pathways(sack, buffer)` - Find pathways just below threshold
  - Identifies pathways within `buffer` distance of their threshold
  - Returns: tibble with `distance_from_threshold`, `total_steps_detected`, `total_steps`
  - Use case: "Which pathways would flip to present with slight threshold adjustment?"
  
  ```r
  near_miss <- find_near_miss_pathways(sack, buffer = 0.1)
  # Shows pathways at 0.65-0.75 when threshold is 0.75
  ```

- [x] `plot_near_miss_pathways(sack, genome_name, buffer)` - Visualize near-miss status
  - Color codes: green = present, orange = near miss, gray = absent
  - Shows which pathways are "almost there"
  
  ```r
  plot_near_miss_pathways(sack, genome_name = "Muricauda_sp_ARW7G5W", buffer = 0.1)
  ```

**Why these matter:** Help users diagnose threshold tuning issues and understand why pathways aren't detected. Especially important when using database-level thresholds (BLAST, HMM) rather than per-gene thresholds (kofam).

#### 1.6 Basic Testing & Quality Control
- [x] Create test potatoes (11 potatoes: glyoxylate cycle, entner-doudoroff variants, nitrogen fixation, etc.)
- [x] Test genomes (22 marine isolate genomes)
- [x] Integration testing: genome → annotation → scoring → visualization
- [x] **Potato verification system** (v0.7.1)
  - [x] All potatoes have `"verified": false` field
  - [x] Agents and LLMs instructed to NEVER set verified to true
  - [x] Manual verification workflow established
  - [x] Text-based review with `print_potato()` for quick validation
- [x] **Quality issues discovered** (v0.7.1)
  - [x] Build-potato agent made substrate specificity errors (EC 3.1.1.31 vs 3.1.1.17)
  - [x] Fixed Entner-Doudoroff pathway errors
  - [x] Split nitrogen fixation into Mo-dependent and V-dependent potatoes
  - [x] Updated agent instructions to validate EC numbers against substrate chemistry
- [ ] Unit tests for scoring logic
- [ ] **HMM TC testing**: Add bacteriorhodopsin potato (uses PFAM with TC) to test per-profile trusted cutoff thresholds

### Phase 0.9: Multi-Pathway Schema (IN PROGRESS - v0.9.0)

**Goal:** Support related pathways in single network potatoes for biological context.

#### 0.9.1 Schema Design ✓ (COMPLETE)
- [x] Design multi-pathway JSON structure
  - [x] Global `nodes` array (genes defined once with detection methods)
  - [x] `pathways` field (each pathway has its own step/required/marker context)
  - [x] Pathway types: `"variant"` vs `"independent"`
  - [x] Per-pathway edges, input/output, scoring
- [x] Add `active` flag to deprecate consolidated potatoes
- [x] Document schema in CLAUDE.md
- [x] Create example: `entner_doudoroff_network.json` (4 ED variants)
- [x] Mark old ED potatoes as `active: false`

#### 0.9.2 Code Updates ✓ (COMPLETE)
- [x] Update `load_potato()` to handle both schemas
  - [x] Detect single-pathway vs multi-pathway
  - [x] Parse `pathways` field into internal representation
  - [x] Maintain backward compatibility
- [x] Update `load_potatoes()` to filter by `active` flag
  - [x] Default: only load active potatoes
  - [x] Optional `include_inactive = TRUE` parameter
- [x] Update `validate_potato()` for multi-pathway schema
  - [x] Validate pathway-specific node references
  - [x] Check edges reference valid nodes within pathway
  - [x] Verify step numbers within each pathway
  - [x] Check marker genes exist per pathway
- [ ] Update `score_pathways()` to score independently (TODO)
  - [ ] For single-pathway: score as before
  - [ ] For multi-pathway: score each pathway separately
  - [ ] Results show per-pathway scores
  - [ ] Aggregate network-level interpretation
- [x] Update `print_potato()` for multi-pathway
  - [x] Show network summary
  - [x] Option to print individual pathways
  - [x] Compact view for each pathway
- [x] Add potato hashing for version tracking
  - [x] `compute_potato_hash()` generates MD5 of functional fields
  - [x] `potato_hash` column added to annotation results

#### 0.9.2.1 Multi-Pathway Scoring ✓ (COMPLETE - v0.9.3)
- [x] Update `score_pathways()` to score multi-pathway networks
  - [x] Each pathway scored independently
  - [x] Results include `pathway` and `pathway_name` columns
  - [x] Single-pathway potatoes work as before
  - [x] Handles shared genes across pathways
- [x] Essential-only scoring metrics
  - [x] Add `steps_detected_essential`, `steps_total_essential`, `fraction_essential`, `present_essential`
  - [x] Tracks required genes separately from all genes
  - [x] NA when no required genes defined
- [x] Add `min_fraction` column to scores
  - [x] Shows per-pathway threshold (from potato JSON or default 0.75)
  - [x] Makes threshold explicit in results
- [x] Result export functions
  - [x] `get_gene_results()` - Gene-level hits with `passed` column
  - [x] `get_pathway_scores()` - Pathway scores with potato_hash
  - [x] Both return tibbles (no write_csv parameter)
- [x] Interactive plot improvements
  - [x] Multi-line pathway lists in hover text
  - [x] Better formatting for genes in multiple pathways
- [x] Test coverage improvements
  - [x] 91 passing tests (up from ~35)
  - [x] Test coverage: 27% (up from 11.35%)
  - [x] New test files: test-export.R, test-scoring.R, test-multi-pathway.R

#### 0.9.3 Visualization Updates ✓ (COMPLETE - v0.9.0/v0.9.1)
- [x] Replace plotly with visNetwork for interactive mode
  - [x] Full viewport display (100vh)
  - [x] Drag-and-drop node positioning
  - [x] Export coordinates to JSON
  - [x] JavaScript function: `exportCoordinates()`
  - [x] Zoom, pan, highlight on hover
  - [x] Navigation controls
- [x] Update `plot_potato()` for multi-pathway networks
  - [x] Option to show all pathways overlaid
  - [x] Option to show single pathway (`pathway` parameter)
  - [x] Gene-based graph structure (not step-based)
  - [x] Shared nodes appear once
  - [x] Force-directed layout for multi-pathway (KK, FR, sugiyama)
- [x] Curated layout system
  - [x] Support embedded x,y coordinates in potato JSON
  - [x] `update_potato_coordinates()` imports coordinates
  - [x] Auto-detection of compound vs enzyme-only layouts
  - [x] Dual coordinate systems: `x,y` and `x_compounds,y_compounds`
  - [x] Hybrid approach: curated coords + layout for new nodes
- [x] Compound visualization improvements
  - [x] Deduplicate compounds (same compound appears once)
  - [x] Position compounds between connected enzymes
  - [x] Different shapes: triangles for compounds, circles for genes
  - [x] Bipartite graph support (`show_compounds` parameter)
- [x] Static mode improvements
  - [x] Pathway convex hulls with ggforce (multi-pathway)
  - [x] Force-directed layouts for multi-pathway
  - [x] No warnings when switching between interactive/static
- [ ] Update pathway heatmap for networks (TODO)
  - [ ] Rows = pathways (not potatoes)
  - [ ] Show variant completion status
- [ ] Add network-level summary plots (TODO)

#### 0.9.4 Build-Potato Agent Updates
- [ ] Update agent instructions for building networks
  - [ ] When to use multi-pathway vs single
  - [ ] How to structure variants
  - [ ] Biological context for network decisions
- [ ] Add validation step for networks
- [ ] Update examples

#### 0.9.5 Consolidate Existing Potatoes (IN PROGRESS)
- [x] Expand `entner_doudoroff_network.json`
  - [x] 4 ED pathway variants
  - [x] Gluconate transport (independent)
  - [x] Galactose catabolism (independent)
  - [x] 23 unique genes, 6 pathways
  - [x] Mark old ED potatoes as `active: false`
- [ ] Create `tca_network.json`
  - [ ] TCA forward (variant)
  - [ ] TCA reverse (variant)
  - [ ] Glyoxylate shunt (independent)
- [ ] Review nitrogen fixation potatoes
  - [ ] Consider network: Mo vs V variants
  - [ ] Or keep separate (already minimal overlap)
- [ ] Identify other candidates

### Phase 0.9 Complete When:
- [x] Can load and validate both single and multi-pathway potatoes
- [ ] Scoring works correctly for networks (per-pathway results) - **NEXT PRIORITY**
- [x] Visualization shows network context
  - [x] Interactive visNetwork with coordinate export
  - [x] Static ggraph with pathway hulls
  - [x] Pathway filtering parameter
- [x] At least 1 network potato created and tested (ED network with 6 pathways)
- [ ] Build-potato agent can create networks (deferred)
- [x] Documentation updated (CLAUDE.md, ROADMAP.md)

### Phase 2: Input/Output Validation & DAG Connectivity

**Goal:** Elevate input/output from plotting metadata to core validation requirements. Pathways must define biologically complete units with reachable inputs and connected outputs.

**Key Insight:** A pathway isn't "complete" just because you have 75% of genes. It's complete when you can transform **available inputs → desired outputs** through detected genes. This fundamentally changes how we score pathways.

#### 2.1 Input/Output as First-Class Citizens

**Current status:** Input/output fields are optional, treated as metadata for visualization

**New requirement:** Input/output are **mandatory** for pathway validation

- [ ] **Schema enforcement**: Input/output required for all pathways
  - Validation fails if missing (with clear error message)
  - Build-potato agent MUST collect input/output (not optional)
  - Update all existing potatoes to add missing input/output fields
  
- [ ] **Biological completeness validation**:
  - Input must be realistic (can be imported/synthesized/is environmental)
  - Output must connect to something (next pathway step, export, storage)
  - Isolated pathways with unreachable inputs flagged as incomplete
  
- [ ] **Compound standardization**:
  - Use KEGG compound IDs for matching (C#####)
  - Support location qualifiers for membrane transport (_exterior, _periplasm, _cytoplasm)
  - Handle compound name variants (glucose vs D-glucose vs α-D-glucose)

**Why this matters:**

```
Example: npED pathway in isolation
- Genes: 5/6 detected (83% complete)
- Input: D-gluconate (but where does it come from?)
- Traditional scoring: PRESENT (>75% threshold)
- Reality: NON-FUNCTIONAL (no gluconate source)

Example: ED network (gluconate_transport + npED)
- gluconate_transport: imports D-gluconate_exterior → D-gluconate_cytoplasm
- npED: transforms D-gluconate_cytoplasm → pyruvate + 2-PG
- Network forms COMPLETE BIOLOGICAL UNIT
- Can transform environmental gluconate → central metabolites
```

#### 2.2 Within-Network DAG Connectivity

**Goal:** Multi-pathway networks must form complete biological units with connected inputs/outputs across pathways.

- [ ] **Build meta-graph across pathways in network**:
  - Link pathway outputs to pathway inputs via compound matching
  - `gluconate_transport.output.compound` → `non_phosphorylative.input.compound`
  - Track compound flow through network
  
- [ ] **Check input reachability**:
  - For pathway with input requirement, verify:
    - Can input be imported? (transport pathway detected)
    - Can input be synthesized? (upstream pathway/genes detected)
    - Is input environmental? (always available)
  - Report: REACHABLE / UNREACHABLE / OPTIONAL
  
- [ ] **DAG traversal from input to output**:
  - Given detected genes, simulate flow through pathway DAG
  - Check: Can you reach output from input with detected genes?
  - Report: COMPLETE PATH / BROKEN (at step X) / INPUT UNAVAILABLE
  
- [ ] **Network-level status**:
  - Evaluate network as functional unit (not just individual pathways)
  - "gluconate_transport + npED = complete unit"
  - "npED alone (no gluconate source) = incomplete"

**Scoring output includes:**
```
Pathway: non_phosphorylative (in ED network)
Gene completion: 5/6 genes = 0.83
Input status: REACHABLE (via gluconate_transport)
DAG connectivity: COMPLETE PATH (input → output traversable)
Network status: COMPLETE BIOLOGICAL UNIT
Overall: FUNCTIONAL
```

vs

```
Pathway: non_phosphorylative (isolated)
Gene completion: 5/6 genes = 0.83
Input status: UNREACHABLE (no gluconate source detected)
DAG connectivity: CANNOT START (input unavailable)
Overall: NON-FUNCTIONAL (genes present but pathway inactive)
```

#### 2.3 Multi-Potato Support
- [ ] Load all potatoes from `inst/potatoes/` directory
- [ ] Filter by tags: `run_potatoes(tags = c("amino_acid_metabolism"))`
- [ ] Batch scoring across multiple potatoes

#### 2.4 Consistency Validation
- [ ] Scan all potatoes for gene definitions
- [ ] Build gene usage matrix (gene → list of potatoes)
- [ ] Compute specificity scores for all genes
- [ ] Validate: flag inconsistent definitions of same gene
- [ ] Validation modes:
  - `strict`: fail on conflicts
  - `warn`: flag conflicts but continue
  - `silent`: no checks

#### 2.3 Visualization
- [x] R function to plot potato DAG with `igraph` (v0.7.0: `plot_potato()`)
- [x] Highlight found vs. missing nodes per genome
- [x] Export to graphviz DOT format
- [x] Generate pathway summary reports (heatmap, barplots)
- [ ] **Faceting for `plot_potato()`** - Show multiple genomes/potatoes in one plot
  - `plot_potato(..., genome_names = c("g1", "g2"), facet_by = "genome")` - One potato, multiple genomes
  - `plot_potato(..., potato_list, genome_name = "g1", facet_by = "potato")` - Multiple potatoes, one genome  
  - `plot_potato(..., potato_list, genome_names, facet_by = "both")` - Grid of potatoes × genomes
  - Technical challenge: ggraph faceting requires careful layout coordination across subplots

### Phase 3: Confidence Scoring & Metabolic Network Analysis

**Goal:** Interpret results in context of gene specificity, MAG completeness, and genome-wide metabolic capability.

#### 3.1 Specificity Weighting
- [ ] Pre-compute specificity scores across potato database
- [ ] Implement weighted scoring:
  ```
  confidence = sum(specificity of found genes) / sum(specificity of all genes)
  ```
- [ ] Output includes:
  - Raw score (5/8 genes found)
  - Weighted score (0.85 confidence accounting for specificity)
  - List of high-value genes found vs. missing

#### 3.2 MAG Completeness Integration
- [ ] Accept optional genome completeness estimate (CheckM, BUSCO)
- [ ] Adjust confidence: "3 genes missing, but MAG is 60% complete"
- [ ] Bayesian-style adjustment (or simpler heuristic)
- [ ] Flag: likely present vs. likely absent vs. uncertain

#### 3.3 Enhanced Output
- [ ] Confidence levels in pathway output
- [ ] Explanation text: "Found leuB, leuC, leuD (pathway-specific). Missing ilvE (also in valine biosynthesis). Likely present given 78% MAG completeness."

#### 3.4 Cross-Potato Metabolic Network Analysis

**Goal:** Evaluate pathway completeness in context of genome-wide metabolic capability. A pathway is complete if its input can be satisfied by another detected pathway's output.

**Biological Reality:** Microbial metabolism is non-compartmentalized
- Eukaryotes: Compartmentalized metabolism (glycolysis in cytoplasm, TCA in mitochondria, Calvin in chloroplast)
- Bacteria: One cytoplasmic soup - all pathways co-localized
- **Implication:** Intermediates from one pathway are immediately available for another
- Pyruvate from ED pathway → available for TCA, fermentation, amino acid synthesis simultaneously
- Only membrane barriers matter (exterior ↔ periplasm ↔ cytoplasm)

**The Transporter Gatekeeper Problem:**

Traditional metabolic models assume "everything extracellular is available intracellularly" (black box transport). This breaks down in real genomes:
- Gluconate abundant in environment but genome lacks transporter → effectively absent
- Transporters are gatekeepers - without them, substrates unavailable
- Must explicitly model transport pathways (e.g., `gluconate_transport`)

**Cross-Potato Linking Example:**

```
Potato 1 (ED network):
  npED pathway: DETECTED
  Output: pyruvate + D-glyceraldehyde (C00577)
  Status: COMPLETE

Potato 2 (serine biosynthesis):
  Input: D-glyceraldehyde (C00577)
  Genes: 75% detected
  Traditional scoring: INCOMPLETE (where does glyceraldehyde come from?)
  
Cross-potato reasoning:
  Check genome for C00577 sources → npED produces it!
  Input: SATISFIED (from npED)
  Serine biosynthesis: COMPLETE (uses npED intermediate)
```

**Implementation:**

- [ ] **Build genome-wide metabolite pool map**:
  - Available from environment: glucose, NH4+, O2, etc.
  - Available after transport: compounds with detected transporters
  - Produced by metabolism: track outputs from all detected pathways
  - Consumed by metabolism: track inputs to all detected pathways
  
- [ ] **Compound matching via KEGG IDs**:
  - Standardize on C##### identifiers for cross-potato matching
  - Handle location qualifiers (_exterior, _cytoplasm)
  - "D-glyceraldehyde" = C00577 regardless of pathway naming
  
- [ ] **Genome-level metabolic graph**:
  - Nodes: compounds (available pools)
  - Edges: pathways (transformations)
  - Build per-genome, per-annotation-run
  - Cache for reuse
  
- [ ] **Contextual pathway scoring**:
  ```r
  # Old: score pathways independently
  score_pathways(sack)
  
  # New: score in metabolic context
  score_metabolic_network(sack, genome_name)
  # Builds genome-wide graph
  # Evaluates pathway completeness considering all detected pathways
  # Returns: which pathways are functionally complete given available inputs
  ```

**Output format:**

```
Genome: Marine_Bacterium_X

Functional Metabolic Units:
✓ Glucose → Pyruvate (classic ED complete)
✓ Gluconate import → Pyruvate (gluconate_transport + npED)
✓ D-glyceraldehyde → Serine (npED intermediate feeds serine biosynthesis)
✓ Acetate import → TCA (acetate transporter + TCA forward)

Incomplete Pathways:
⚠ Leucine biosynthesis: 80% genes but missing α-ketoisovalerate source
  - Check: Valine biosynthesis (produces α-ketoisovalerate) - NOT DETECTED
  - Pathway cannot function without intermediate source
  
⚠ Aromatic amino acid biosynthesis: 70% genes but missing chorismate
  - Shikimate pathway NOT DETECTED
  - No alternative chorismate source
  - Likely imports aromatic amino acids (auxotroph)
```

**Statistical Validation:**

Once 50-100 potatoes exist, analyze compound sharing:
- Build compound co-occurrence matrix (which compounds appear in which pathways)
- Calculate average: Z = compounds per pathway
- If Z > 1: compounds are inherently shared (validates cross-potato reasoning)
- Cluster pathways by shared metabolites:
  - Do clusters match functional domains? (glycolysis variants)
  - Or extensive cross-domain linkage? (central metabolism ↔ amino acids ↔ toxins)

**Research Questions:**
- How often do compounds link completely different domains? (EMP ↔ toxin production)
- Is microbial metabolism modular or highly interconnected?
- Can we predict pathway co-occurrence based on shared intermediates?

**Implementation Challenges:**

1. **Compound name standardization** - Always use KEGG C##### for matching
2. **Computational cost** - Build graph per-genome on demand, cache results
3. **Circular dependencies** - Detect cycles, report separately
4. **Multiple sources** - If 3 pathways produce X, which feeds pathway Y? (report all)
5. **Intermediate vs terminal** - Track compound pool availability

**Priority:** Medium-High - Implement after Phase 2 (input/output validation) proves the concept within networks

**This makes POTATO unique:** Most tools score pathways independently. POTATO reasons about metabolic networks: "You have 60% of pathway X, but intermediate Y comes from pathway Z → pathway X is complete."

#### 3.5 Threshold Sensitivity Analysis ⚠️ HIGH-PRIORITY, HIGH-LOAD

**Goal:** Identify "gate-keeper genes" whose thresholds block pathway detection and quantify threshold sensitivity.

**Motivation:** BLAST and HMM use database-level thresholds (not per-gene like kofam), so we're not fine-tuning individual genes' needs. This analysis helps users understand which genes would benefit from relaxed thresholds.

**Implementation:**
- [ ] `analyze_threshold_sensitivity(sack, relaxation_steps = c(0.1, 0.2, 0.5))`
  - For each gene in each pathway:
    - Test relaxing threshold by 10%, 20%, 50%
    - Track which pathways become detected at each step
    - Identify "gate-keeper genes" (single gene blocking pathway)
  - Return: tibble with gene, pathway, current_status, relaxed_status, new_hits_gained
  
- [ ] `plot_threshold_sensitivity(sack, pathway_name)`
  - Visualize pathway completion vs. threshold relaxation
  - Highlight which genes are most sensitive
  - Show "tipping points" where pathway flips to present
  
- [ ] Integration with scoring:
  - Flag pathways "close to threshold" (within buffer)
  - Suggest per-gene threshold adjustments for potato JSON
  - Generate recommended `thresholds:` blocks for potato nodes

**Example output:**
```
Gene: mlrA (microcystin_degradation)
Current: 0/22 genomes detected (e-value 1e-10)
Relaxed 50%: 3/22 genomes detected (e-value 1e-5)
Impact: pathway detected in 3 additional genomes
Recommendation: Add to potato JSON:
  "thresholds": {"hmm_evalue": 1e-5}
```

**Why high-load:**
- Requires re-running hit detection at multiple thresholds (computationally expensive)
- Need to track and compare results across multiple threshold scenarios
- Complex data structure to represent multi-dimensional sensitivity
- UI/visualization needs careful design to be interpretable

### Phase 4: LLM Integration

**Goal:** LLM agents for database building and result interpretation.

#### 4.1 Builder Agent ✓ (IMPLEMENTED - v0.7.1)

Interactive skill for creating new potato JSON files. Invoked with `/build-potato` command.

**Status: Operational with quality control measures**

**Capabilities:**
- [x] **KEGG module import** - Fetch and parse KEGG module definitions
  - [x] Parse module syntax (parentheses for OR, + for AND)
  - [x] Fetch KO details from KEGG API
  - [x] Validate EC numbers against substrate chemistry
  - [x] Generate input/output compound fields
  - [x] Set `"verified": false` (NEVER true)
- [x] **Custom pathway creation** - Conversational pathway building
  - [x] Push back on vague requests ("nitrogen fixation" → "which type?")
  - [x] Suggest PFAM domains and BLAST references
  - [x] Ask about input/output compounds
  - [x] Handle transporters with location qualifiers (_external, _internal)
- [x] **Potato migration** - Update old potatoes to new schema
  - [x] Map custom database names to standard types
  - [x] Add verified field if missing
  - [x] Add input/output fields if missing
- [x] **GATOR Excel conversion** - Parse old GATOR v1 spreadsheets
  - [x] Parse string syntax (->  | + operators)
  - [x] Map to potato JSON DAG structure

**Quality control lessons learned (v0.7.1):**
- ⚠️ Agent made substrate specificity errors (wrong EC numbers)
- ⚠️ Example: Used EC 3.1.1.31 (6-phosphogluconolactonase) instead of EC 3.1.1.17 (gluconolactonase) in non-phosphorylative pathway
- ✓ Solution: Enhanced validation instructions to cross-check EC substrate specificity
- ✓ Added `print_potato()` for quick text-based verification
- ✓ Mandatory `"verified": false` field - only humans verify potatoes
- ✓ Agent now validates KEGG module definitions more carefully

**Known limitations:**
- Agent can make biochemical errors - human verification required
- Does not guarantee substrate specificity correctness
- User must manually set `verified: true` after thorough review
- [ ] LLM asks clarifying questions
- [ ] Generates potato JSON

#### 4.2 Converter Agent

Migrate existing gator_db.xlsx to potato JSON files.

```r
convert_gator_db(excel_path = "gator_db.xlsx", 
                 output_dir = "inst/potatoes/")
```
- [ ] Reads Excel sheets (db, gene, pathway)
- [ ] Parses v1 pathway syntax (`->`, `|`, `+`)
- [ ] Converts to DAG structure
- [ ] One JSON per pathway
- [ ] Validation report

#### 4.3 Analysis Agent - LLM-Enhanced Pathway Interpretation

**Goal:** Use LLM to interpret pathway scoring results with biological reasoning that goes beyond algorithmic scoring.

**Key Innovation:** LLM provides context-aware interpretation of gene presence patterns, accounting for:
1. **Gene specificity** - Distinguish diagnostic genes (mlrA) from promiscuous housekeeping genes (gapA)
2. **Cross-pathway context** - Recognize when detected genes actually belong to a different pathway (PPP vs ED overlap)
3. **Biological feasibility** - Assess whether partial pathways are likely functional

**Architecture:**
- **Algorithmic scoring first** (fast, deterministic): `sack <- score_pathways(sack)`
- **LLM interpretation second** (optional, contextual): `interpret_pathway_llm(sack, potato_name, genome_name)`
- **R interface via elmer** - Simple LLM access (Anthropic, OpenAI, local models)
- **Config in potato_config.yaml**:
  ```yaml
  llm:
    provider: "anthropic"  # or "openai", "ollama"
    api_key_env: "ANTHROPIC_API_KEY"
    model: "claude-sonnet-4-5"
  ```

**Implementation:**

##### 4.3.1 Gene Specificity Interpretation
- [ ] Pre-compute gene specificity scores across all potatoes (1 / num_potatoes containing gene)
- [ ] LLM prompt includes:
  - Detected genes with their specificity scores
  - Pathway completion percentage
  - Which genes are required vs optional
  - Which genes are markers
- [ ] LLM output: "60% complete with mlrABC (highly specific) = strong evidence" vs "60% complete with gapA/pgk (promiscuous) = weak evidence"

**Example:**
```r
interpret_pathway_llm(sack, 
  potato_name = "microcystin_degradation", 
  genome_name = "Haloferax_sp"
)

# LLM response:
# "Microcystin degradation pathway shows STRONG evidence (60% complete):
#  - mlrA, mlrB, mlrC detected (all pathway-specific, avg specificity 0.98)
#  - Missing mlrD (required but sometimes absent in partial operons)
#  - High confidence despite incompleteness given marker gene presence"
```

##### 4.3.2 Cross-Pathway Context Awareness
- [ ] LLM aware of major central metabolism pathways (PPP, EMP, TCA, ED)
- [ ] Optional: Fetch KEGG context via API (`include_kegg_context = TRUE`)
  - For each detected KO, query which KEGG modules it appears in
  - Format: "K00036: appears in M00004 (PPP), M00008 (ED)"
- [ ] Pre-loaded pathway context (curated definitions of common pathways)
- [ ] LLM checks: "Are detected genes ALSO in PPP/EMP/TCA? Could this be false positive from overlap?"

**Example:**
```r
interpret_pathway_llm(sack,
  potato_name = "entner_doudoroff_classic",
  genome_name = "E_coli_K12",
  include_kegg_context = TRUE
)

# LLM response:
# "ED pathway shows WEAK evidence (60% complete):
#  - Detected: zwf, pgl, gnd (all shared with oxidative PPP)
#  - Missing: edd, eda (ED-specific markers, CRITICAL for ED)
#  - KEGG context confirms zwf/pgl in both M00004 (PPP) and M00008 (ED)
#  - Interpretation: This is likely oxidative PPP, NOT ED pathway"
```

**Prompt Structure:**
```r
CENTRAL_METABOLISM_CONTEXT <- "
Common metabolic pathways and their diagnostic genes:

Pentose Phosphate Pathway (PPP):
- Oxidative: zwf (K00036), pgl (K01057), gnd (K00033) - NOT ED-specific
- Non-oxidative: rpe, rpi, tkt, tal

Entner-Doudoroff (ED):
- Classic: glk, zwf, pgl, edd (K01690)*, eda (K01625)* - edd/eda are diagnostic
- Non-phosphorylative: gnaD, kdgA/kdpgA, gadh/cutABC/aor, gck

Glycolysis (EMP):
- Lower: gapA, pgk, pgm, eno, pyk - highly promiscuous
- Upper: pgi, pfk, fba

TCA cycle:
- Forward: citrate synthase, aconitase, isocitrate DH, etc.
- Reverse: ATP citrate lyase, etc.
"

interpret_pathway_llm <- function(sack, potato_name, genome_name, 
                                   include_kegg_context = FALSE) {
  # Get scoring results
  scores <- sack@scores %>% filter(potato == potato_name, genome == genome_name)
  gene_hits <- get_gene_results(sack) %>% 
    filter(potato == potato_name, genome == genome_name, passed == TRUE)
  
  # Build prompt
  prompt <- glue::glue("
    Pathway: {potato_name}
    Completion: {scores$fraction} ({scores$steps_detected}/{scores$steps_total})
    Essential completion: {scores$fraction_essential}
    
    Detected genes: {paste(gene_hits$node_id, collapse=', ')}
    Gene specificity scores: {paste(gene_hits$specificity, collapse=', ')}
    Marker genes detected: {paste(gene_hits$node_id[gene_hits$marker], collapse=', ')}
    
    Context: {CENTRAL_METABOLISM_CONTEXT}
  ")
  
  # Optionally add KEGG API context
  if (include_kegg_context) {
    kegg_info <- fetch_kegg_context_for_genes(gene_hits$ko)
    prompt <- paste(prompt, "\n\nKEGG module lookup:\n", kegg_info)
  }
  
  # Critical questions for LLM
  prompt <- paste(prompt, "
    CRITICAL: Before concluding {potato_name} is present, verify:
    1. Are detected genes ALSO found in PPP/EMP/TCA? (check context above)
    2. Are pathway-SPECIFIC genes present? (not just shared housekeeping)
    3. Could this be a false positive from metabolic overlap?
    4. Do the detected genes have high or low specificity scores?
    
    Provide biological interpretation: Is this STRONG or WEAK evidence for pathway presence?
  ")
  
  # Query LLM via elmer
  chat <- elmer::chat_anthropic(
    model = sack@config$llm$model,
    api_key = Sys.getenv(sack@config$llm$api_key_env)
  )
  chat$chat(prompt)
}
```

##### 4.3.3 Functional Analog Suggestions
- [ ] LLM uses reaction/EC metadata to suggest alternatives
- [ ] "Missing ketol-acid reductoisomerase (EC 1.1.1.86). Found hypothetical protein with EC 1.1.1.* domain. Possible functional analog?"
- [ ] Flag as unvalidated, requires expert review
- [ ] Consider phylogenetic distance: "Archaeal genome missing bacterial homolog - check for archaeal-specific isoforms"

##### 4.3.4 MAG Completeness Integration
- [ ] Accept optional genome completeness estimate (CheckM, BUSCO)
- [ ] Adjust confidence: "3 genes missing, but MAG is 60% complete - absence may be due to incompleteness"
- [ ] Bayesian-style adjustment or simpler heuristic

**Priority:** Medium-High - Provides major value-add over pure algorithmic scoring. Implement after Phase 3.1 (Specificity Weighting) is complete so specificity scores are available

#### 4.4 Enhancement Agent - Potato Validation & Expansion

**Goal:** LLM agent suggests additional detection methods, missing genes, and structural improvements for existing potatoes.

**Use Case:** User creates a potato and wants to know: "What am I missing? Are there better detection methods? Should I add regulatory genes?"

**What the agent is good at:**

✅ **Cross-database mapping**
- Suggest PFAM domains for genes currently only detected by KO
- Map to PATRIC/RAST annotations
- LLMs have broad training on database cross-references

✅ **Suggesting alternative/missing genes**
- "You have nifH but missing nifD/nifK (nitrogenase complex subunits)"
- "Consider vnfH/vnfG for V-dependent variant"
- Pattern recognition across pathway definitions

✅ **Regulatory context**
- "This pathway typically regulated by ntrC"
- "FNR regulates ED pathway genes in E. coli"
- Literature knowledge about operons/regulons

✅ **Structural gap identification**
- "Step 2 produces compound X but no step consumes it"
- "Missing transporter for substrate uptake"
- DAG completeness reasoning

**What would be risky (requires mitigation):**

❌ **Database ID hallucination**
- LLM might invent PFAM/PATRIC IDs that don't exist
- Training data is old but databases update constantly
- **Mitigation:** API verification required before adding to potato

❌ **EC number / substrate specificity errors**
- Already seen build-potato make EC chemistry mistakes
- **Mitigation:** User must verify EC substrate chemistry, not trust blindly

❌ **Rare/novel pathway bias**
- LLM training biased toward well-studied pathways
- Custom/niche pathways may get bad suggestions
- **Mitigation:** Lower confidence when source is "custom" not "KEGG"

**Implementation Design:**

```r
suggest_potato_enhancements(potato, check_databases = TRUE)
```

**Returns structured suggestions:**
```
Potato: entner_doudoroff_network (classic pathway)

DETECTION METHOD SUGGESTIONS:
✓ G6PD (K00036) - VERIFIED
  - PFAM PF00479 (G6PD C-terminal) ✓ verified via API
  - PATRIC PLF_1240_00001700 ✓ verified
  - Confidence: HIGH
  
⚠ edd (K01690) - CHECK NEEDED  
  - PFAM PF03446 (edd domain) - API lookup failed, verify manually
  - Confidence: MEDIUM

MISSING PATHWAY COMPONENTS:
! No glucokinase regulator detected
  - Consider: crp (catabolite repression, K02529)
  - Reasoning: ED pathway subject to glucose catabolite repression
  - Confidence: MEDIUM (regulatory, not structural)

STRUCTURAL CHECKS:
✓ All steps have downstream connections
✓ No orphan compounds
! Input compound (D-glucose) - no transporter defined
  - Consider: glucose PTS system or porins for uptake
  - Confidence: LOW (organism-dependent)
```

**Key Features:**
1. **API verification built-in** - `check_databases = TRUE` queries KEGG/PFAM/PATRIC APIs to verify suggested IDs
2. **Confidence scores** - HIGH (API verified) / MEDIUM (plausible) / LOW (speculative)
3. **Reasoning included** - explains *why* each suggestion
4. **Categorized output** - detection methods vs pathway components vs structure
5. **User validates** - returns checklist, doesn't auto-edit JSON (user decides what to add)

**Implementation Steps:**

- [ ] LLM analyzes potato structure
  - Input: `print_potato()` compact summary
  - Prompt: "Suggest improvements: additional detection methods, missing genes, regulatory components, structural issues"
  - Output: Structured JSON with suggestions + reasoning + confidence
  
- [ ] API verification layer (`check_databases = TRUE`)
  - `verify_pfam_id()` - Query PFAM API
  - `verify_ko_id()` - Query KEGG API
  - `verify_patric_id()` - Query PATRIC/RAST if available
  - Mark each suggestion: ✓ verified / ⚠ check needed / ✗ not found
  
- [ ] Structure suggestions into categories
  - Detection methods (add to existing genes)
  - Missing pathway components (new genes)
  - Regulatory genes (optional additions)
  - Structural issues (gaps, orphans, missing transporters)
  
- [ ] Interactive workflow
  - User reviews suggestions
  - Manually adds approved items to potato JSON
  - Re-run `validate_potato()` after changes

**Example Prompt Structure:**
```r
prompt <- glue::glue("
  Analyze this metabolic pathway and suggest improvements:
  
  Pathway: {potato@name}
  Source: {potato@source}
  {print_potato(potato, compact = TRUE)}
  
  Suggest:
  1. Additional detection methods (PFAM, PATRIC, HMM) for existing genes
  2. Missing alternative genes (isoforms, paralogs, variants)
  3. Regulatory genes commonly associated with this pathway
  4. Structural issues (missing transporters, cofactor dependencies, orphan compounds)
  
  For each suggestion provide:
  - Specific database ID (PFAM/PATRIC/KO)
  - Confidence level (high/medium/low)
  - Reasoning (why this is relevant)
  
  Be conservative - only suggest well-established associations.
  Do NOT suggest genes unless you're confident they belong in this pathway.
")
```

**Use as Research Assistant:**
- Agent finds leads, user verifies before adding
- Think: grad student literature search, not oracle
- Especially useful for:
  - Adding PFAM to KEGG-only pathways
  - Checking for missing regulatory genes
  - Identifying structural gaps

**Priority:** Medium - Implement after Phase 4.3 (Analysis Agent) when LLM infrastructure exists

#### 4.5 LLM Infrastructure
- [ ] Claude API integration (via Anthropic SDK) or elmer wrapper
- [ ] Prompt templates for each agent mode (builder, analyzer, enhancer)
- [ ] Token management, caching
- [ ] Error handling, retry logic
- [ ] API verification helpers (KEGG, PFAM, PATRIC lookups)

### Phase 5: Advanced Features (Future)

#### 5.0 Ghost Edges for Cyclic Pathways

**Problem:** Some metabolic pathways are truly cyclic (e.g., BHAC cycle, TCA cycle), but validation rejects cycles because pathways must be DAGs for scoring logic.

**Current limitation:**
- BHAC has glycine regenerated by bhcA and fed back to bhcC (biological reality)
- Adding edge `bhcA → bhcC` creates cycle
- Validation fails: "Pathway 'bhac_cycle': contains cycles (must be DAG)"

**Proposed solution: Ghost edges**
- Special edge type `"ghost": true` that visualization renders but validation/scoring ignores
- Allows biological accuracy in visualizations without breaking DAG logic

**Example:**
```json
{
  "edges": [
    {"from": "bhcC", "to": "bhcB", "compound": "β-Hydroxyaspartate"},
    {"from": "bhcB", "to": "bhcD", "compound": "Iminosuccinate"},
    {"from": "bhcD", "to": "bhcA", "compound": "L-Aspartate"},
    {
      "from": "bhcA", 
      "to": "bhcC", 
      "compound": "Glycine",
      "ghost": true,
      "notes": "Glycine regeneration - closes the cycle"
    }
  ]
}
```

**Implementation:**
- Graph builders check for `"ghost": true` and skip when building igraph for validation/scoring
- Visualization includes ghost edges with distinct styling (dashed lines?)
- Print functions show ghost edges with annotation

**Benefits:**
- Validates as DAG (scoring works)
- Visualizes as cycle (biological accuracy)
- Documents cyclic nature explicitly

---

#### 5.1 Structured Input/Output Compound Storage

**Problem:** Input/output compounds are currently stored as single strings (e.g., `"compound": "pyruvate + glyceraldehyde-3P"`), which makes parsing difficult and error-prone.

**Current limitation:**
- Single string field can contain multiple compounds separated by " + " or " and "
- KEGG compound IDs may not align with compound order in string
- Splitting logic is fragile and prone to edge cases
- Visualizations currently show full string as single node (simplified but not ideal)

**Proposed solution: Structured arrays**

Store input/output compounds as arrays with explicit compound/KEGG ID pairs:

```json
{
  "input": {
    "compounds": [
      {"name": "D-glucose", "kegg_compound": "C00031"}
    ],
    "targets": ["glk"]
  },
  "output": {
    "compounds": [
      {"name": "pyruvate", "kegg_compound": "C00022"},
      {"name": "glyceraldehyde-3P", "kegg_compound": "C00118"}
    ],
    "sources": ["eda"]
  }
}
```

**Benefits:**
- Unambiguous compound boundaries
- Each compound has its own KEGG ID
- Easy to create separate nodes for each compound in visualizations
- Simpler deduplication logic (compare KEGG IDs directly)
- Better for programmatic pathway construction

**Implementation notes:**
- Maintain backward compatibility with legacy string format
- Auto-upgrade on load: detect string format, parse to structured format
- Validation warns about unparseable strings
- Build-potato agent generates structured format by default

**Status:** Planned for future enhancement

---

#### 5.1 Self-Contained Database Embedding

**Goal:** Embed HMM/BLAST sequences directly in potato JSON for true portability.

**Current approach:**
- Potato references "K00001" (external identifier)
- User must have kofam profiles installed locally
- tools.json points to database location

**Enhanced approach:**
- Potato embeds the actual HMM text or BLAST sequences
- No external database installation needed
- Truly portable, shareable potatoes

**Example:**
```json
{
  "nodes": [
    {
      "id": "ilvB",
      "name": "acetolactate synthase",
      "ko": ["K01652"],
      "embedded_data": {
        "hmm": "HMMER3/f [3.3.2 | Nov 2020]\nNAME  K01652\n...",
        "blast_seqs": ">K01652_ref1\nMKLITVGAAG...\n>K01652_ref2\nMRLVSVGAAG..."
      }
    }
  ]
}
```

**Use Cases:**
- ✅ Custom HMMs (user-built, no licensing issues)
- ✅ BLAST reference sequences (user-curated from public databases)
- ❌ KEGG HMMs (no explicit license = all rights reserved, do NOT embed/redistribute)
- ✅ PFAM HMMs (open license, CC0/public domain)
- ✅ TIGRFAMs (open license)
- ✅ User-generated data (sequences, alignments, custom HMMs)

**Trade-offs:**

*Pros:*
- Ultimate portability (email a potato, it works)
- Version control (exact HMM used is preserved)
- No setup required (no tools.json needed)
- Perfect for custom/niche pathways

*Cons:*
- File size (HMMs are text but can be large, 5-50KB each)
- Licensing concerns (KEGG is restrictive about redistribution)
- Redundancy (same HMM in multiple potatoes)
- Update propagation (improved HMM = update all potatoes)

**Implementation Strategy:**

1. **Optional, not required** - potatoes can reference external DBs OR embed data
2. **Hybrid mode** - embed custom HMMs, reference standard databases
3. **Compression** - gzip embedded HMM text (reduce size ~70%)
4. **Licensing check** - warn if embedding KEGG data
5. **Fallback** - if embedded data exists, use it; else fall back to tools.json

**JSON Structure:**
```json
{
  "nodes": [
    {
      "id": "mlrA",
      "name": "microcystinase",
      "detection_methods": [
        {
          "type": "hmm",
          "source": "embedded",
          "data": "HMMER3/f...",
          "compressed": true
        },
        {
          "type": "blast",
          "source": "embedded",
          "sequences": [
            {"id": "mlrA_ref1", "seq": "MKLITV..."},
            {"id": "mlrA_ref2", "seq": "MRLVSV..."}
          ]
        },
        {
          "type": "ko",
          "source": "external",
          "ids": ["K01234"]
        }
      ]
    }
  ]
}
```

**Licensing Status:**
- ✅ **PFAM** - CC0 public domain, freely redistributable
- ✅ **TIGRFAMs** - Open license (now merged with PFAM)
- ✅ **UniProt sequences** - CC BY 4.0, can use with attribution
- ❌ **KEGG HMMs** - No explicit license on profiles (KofamScan software is MIT, but HMM data is not). Publicly available ≠ redistributable. **Do not embed KEGG HMMs in distributed potatoes.**
- ✅ **Custom HMMs** - User-generated, no restrictions

**Practical Guidance:**
- **For personal use:** Embed whatever you want (KEGG, PFAM, custom)
- **For public distribution:** Only embed openly licensed data (PFAM, custom, UniProt with attribution)
- **For KEGG-dependent potatoes:** Distribute potato JSON with KO references only, let users provide their own KEGG installation via tools.json

**Builder Agent Integration:**
When building a potato, agent asks:
- "Do you want to embed the HMM/sequence data for portability?"
- "Embedding will increase file size by ~XXX KB but make the potato fully portable"
- **Licensing check:** "Detected KEGG HMMs. These cannot be embedded in publicly distributed potatoes (no redistribution license). Embed for personal use only, or omit for public distribution."
- Safe to embed: PFAM, TIGRFAMs, custom HMMs, UniProt sequences (with attribution)

**Tools.json becomes optional:**
- If all potatoes use embedded data, tools.json not needed
- If mixing embedded + external, tools.json provides fallback paths

---

#### 5.2 Pathway Evolution Tracking
- [ ] Track potato changes over time (git-based versioning)
- [ ] Compare old vs new pathway definitions
- [ ] Re-score genomes when potatoes are updated
- [ ] Report: "pathway X now detected in 5 additional genomes due to updated thresholds"

#### 5.3 Cross-Potato Queries
- [ ] "Which potatoes produce acetyl-CoA?"
- [ ] "Show me all pathways downstream of pyruvate"
- [ ] Network visualization across potatoes

#### 5.4 ESM-2 Embeddings (Protein LM)
- [ ] Sequence-similarity rescue for missed annotations
- [ ] BacPT/Bacformer as pre-computed embedding database
- [ ] ANN search for functional annotation

#### 5.5 LLM-Assisted Potato Refinement
- [ ] LLM suggests pathway improvements based on detected patterns
- [ ] "I found bifunctional enzyme X that shortcuts step 2-3, should we add it?"
- [ ] Suggests new edges/nodes to existing potatoes
- [ ] User reviews, accepts/rejects
- [ ] Updates potato JSON with user approval

#### 5.6 Alternative Visualization Libraries (Future Consideration)
- [ ] **g6R** - Explore as alternative to visNetwork for pathway-based clustering
  - Automatic pathway-based node coloring via combos
  - More modern graph rendering (AntV G6)
  - May provide better visual grouping for multi-pathway networks
  - Note: visNetwork works well currently, this is a "nice to have" for future enhancement
  - Reference: https://www.cynkra.com/blog/2025-07-10-g6R/

---

## Key Decisions

### 1. Gene Definition Strategy

**Decision:** Self-contained potatoes with optional canonical registry

**Rationale:**
- Portability (each potato is independent)
- Flexibility (can override when needed)
- Validation warns of inconsistencies without blocking

**Implementation:**
- Default: define genes inline in potato JSON
- Optional: reference canonical genes
- Validation mode configurable (strict/warn/silent)

### 2. DAG vs. String Syntax

**Decision:** DAG structure in JSON

**Rationale:**
- Handles KEGG module complexity
- Visually graphable (human-reviewable)
- Not constrained by Excel cell limits
- LLM can reason about structure
- Easy to add alternative paths

**Trade-off:** More verbose than string, but that's acceptable in JSON

### 3. Compound Nodes

**Decision:** Compounds as edge metadata, not structural nodes

**Rationale:**
- Scoring logic stays enzyme-centric
- Preserves biological context
- Enables LLM reasoning about chemistry
- Optional (can omit if not relevant)

**Future:** Could promote to full nodes for metabolic network analysis

### 4. Tool Configuration

**Decision:** Separate tools.json, continue with unavailable tools

**Rationale:**
- User-specific paths (not hardcoded)
- Shareable potatoes work for users with different tool setups
- Graceful degradation (use what's available)
- Transparent reporting (tools used vs. requested)

### 5. Specificity Weighting

**Decision:** Pre-computed using LLM, baked into potato database

**Rationale:**
- Fast at runtime (no per-genome LLM calls)
- Contextual reasoning during pre-computation
- Recompute when database changes

**Implementation:**
- Scan all potatoes on load/update
- Build gene usage matrix
- Compute specificity (1 / num_potatoes)
- Optional LLM enhancement (considers biochemical context)

---

## Success Criteria

### Phase 1 MVP Complete When:
- ✅ Can run a potato against a genome
- ✅ Produces gene + pathway output files
- ✅ Results match expected behavior (validated against test cases)

### Phase 2 Complete When:
- ✅ Can work with 20+ potatoes simultaneously
- ✅ Tag-based filtering works
- ✅ Validation catches inconsistencies
- ✅ Can visualize potato DAGs

### Phase 3 Complete When:
- ✅ Confidence scores incorporate specificity weighting
- ✅ MAG completeness affects interpretation
- ✅ Output explains why pathway is likely present/absent

### Phase 4 Complete When:
- ✅ Can build potato from KEGG module M00570
- ✅ Can convert all gator_db.xlsx pathways to JSON
- ✅ Analysis agent produces useful interpretations

### Phase 5 Complete When:
- ✅ LLM suggests functional analogs for missing enzymes
- ✅ Can query cross-potato relationships
- ✅ DAG evolution suggestions are useful

---

## Open Questions

1. **Specificity computation**: Simple count-based or LLM-enhanced with biochemical context?
2. **Canonical gene registry**: Required, optional, or skip entirely?
3. **Output formats**: TSV only, or also JSON/HDF5 for large datasets?
4. **LLM provider**: Claude only, or support multiple (GPT-4, local models)?
5. **Validation strictness**: Default to strict or warn mode?
6. **Potato versioning**: How to handle updates to potato definitions over time?
7. **Reaction metadata**: Required, optional, or generate on-demand from EC numbers?
8. **Edge compound info**: Always required, or optional?
9. **Embedded database data**: Should this be Phase 5 or earlier? KEGG licensing resolved: do not redistribute KEGG HMMs.
10. **File size limits**: Max potato JSON size with embedded HMMs? Should we compress by default?
11. **Licensing enforcement**: Should builder agent prevent KEGG embedding, or just warn?
12. **Database consensus/agreement**: How to track and report when multiple databases agree or disagree on gene function assignments?

---

## Future Enhancements (Post v1.0)

### Visualization Enhancements

- [ ] **Gradient heatmap**: `plot_pathway_heatmap()` show completion fraction (0.0-1.0) as color gradient instead of binary present/absent
  - Current: green = present, gray = absent
  - Proposed: color gradient from red (0%) → yellow (50%) → green (100%)
  - Shows "almost there" pathways more clearly
  
- [ ] **Faceting support**: Add faceting to `plot_potato()` for comparing detection across multiple genomes side-by-side

### Threshold Sensitivity Analysis

**Goal**: Identify which gene thresholds are blocking pathway detection and test optimal threshold settings.

**Problem**: Some pathways may be absent only because one or two genes have overly strict thresholds. Hard to know which genes to relax without systematic testing.

**Proposed Feature**: Per-gene threshold relaxation analysis

```r
# Test what happens if each gene's threshold is relaxed
sensitivity <- threshold_sensitivity_analysis(
  sack, 
  relaxation_steps = c(0.1, 0.2, 0.3),  # Try 10%, 20%, 30% relaxation
  pathways = c("glyoxylate_cycle"),     # Optional: specific pathways
  genomes = c("genome1", "genome2")     # Optional: specific genomes
)

# Returns tibble showing:
# - pathway, genome, gene, original_threshold, relaxation_level
# - pathway_detected (TRUE/FALSE at each relaxation)
# - gate_keeper (TRUE if this gene blocks detection)
```

**Use Cases:**
- "aceA has blast e-value 1e-20, but relaxing to 1e-18 would call glyoxylate cycle present in 12 genomes"
- "This pathway is robustly absent - even relaxing all thresholds 50% doesn't call it present"
- "nifH threshold is too strict - it's a gate-keeper gene preventing nitrogen fixation detection"

**Visualization:**
```r
plot_threshold_sensitivity(sensitivity)
# Heatmap: genes × relaxation levels, colored by pathway status change
```

**Priority**: Medium - implement after basic scoring analysis functions are tested

---

### GenBank Conversion

**Current Status**: v1.0 only accepts FAA (protein FASTA) files as input. GenBank conversion logic was removed during foundation cleanup.

**Needed**: Function to convert GenBank files (.gb, .gbk, .gbff) to FAA format for annotation.

**Implementation Options**:
1. Use jakomics `GENOME.from_genbank()` to parse GenBank and extract proteins
2. Call Biopython SeqIO directly
3. Use external tool (e.g., `genbank_to_fasta.py`)

**Priority**: Medium - Many users have GenBank files from NCBI/IMG, but FAA is increasingly common for MAGs.

---

### Workflow Message Logging

**Current Status**: Message collection system removed during foundation cleanup to reduce complexity before core workflows exist.

**Needed**: System to collect and retrieve warnings, errors, and info messages during annotation workflows. Useful for reviewing what happened during batch runs on hundreds of MAGs.

**Implementation**:
1. `add_sack_message()` - Internal function to append messages to `sack@messages` list
2. `sack_msg()` - Helper that prints (if verbose) AND stores messages
3. `sack_messages()` - User-facing retrieval/filtering (by level, stage, genome)
4. `sack_message_summary()` - Formatted summary showing errors/warnings/info counts

**Use Cases**:
- "Which genomes had kofam annotation failures?"
- "Show me all warnings from the scoring stage"
- "Did any genomes trigger low-completeness warnings?"

**Design Considerations**:
- Should messages be stored in sack object or written to log file?
- How to handle message volume in large batch runs (1000+ genomes)?
- Integration with R's message()/warning() system vs. separate logging

**Priority**: Low - Implement after core annotation workflows are stable and generating messages worth collecting.

---

### Database Agreement & Consensus Tracking

**Problem**: When a potato gene can be detected by multiple databases (e.g., KOfam + BLAST), it's valuable to know:
1. **Agreement** - Multiple databases confirm the same gene has the expected function (high confidence)
2. **Conflict** - Different genes in the genome are assigned the same function by different databases (ambiguity)
3. **Complementarity** - One database finds the gene, others don't (need to investigate why)

**Use Cases**:
- "KOfam says gene A is nifH (K02588), BLAST confirms gene A hits nifH reference → HIGH CONFIDENCE"
- "KOfam says gene A is aceA, but BLAST says gene B is aceA → CONFLICT/AMBIGUITY"
- "Only BLAST found this gene, KOfam missed it → Check HMM quality or score threshold"

**Current Status**: 
- v1.0 stores results by database but doesn't analyze cross-database agreement
- No current potatoes use multiple databases for the same gene (roadmapped for future)

**Proposed Implementation (v1.5+)**:

1. **Potato schema extension**: Allow genes to specify multiple databases
   ```json
   {
     "id": "nifH",
     "databases": {
       "kofam113": ["K02588"],
       "gator_blast": ["nifH_ref1", "nifH_ref2"]
     },
     "consensus_required": "any"  // or "all", "majority"
   }
   ```

2. **Scoring enhancements**: Track which databases contributed to detection
   ```r
   sack@scores:
     detected_genes: ["nifH"]
     detected_by_database: list(nifH = c("kofam113", "gator_blast"))  # consensus
     consensus_level: "high"  # high/medium/low based on agreement
   ```

3. **Conflict detection**: Flag when different genes match different databases
   ```r
   get_conflicts(sack) %>%
     filter(potato == "glyoxylate_cycle", gene == "aceA")
   # Shows: gene_X from kofam, gene_Y from blast → investigate
   ```

4. **Summary functions**: 
   ```r
   summarize_database_agreement(sack)
   # Returns: % of genes with multi-DB confirmation, conflicts found, etc.
   ```

**Design Considerations**:
- Should consensus be **required** (all DBs must agree) or **additive** (any DB finding counts)?
- How to weight databases? (BLAST may be more promiscuous than KOfam)
- Should conflicts block pathway detection or just flag for review?
- How to handle DB-specific thresholds? (KOfam has per-KO thresholds, BLAST uses global e-value)

**Priority**: Medium (post-v1.0) - Wait for real use cases with multi-DB potatoes before implementing

---

## Migration Path

### For Existing GATOR Users

1. **Install potato R package** + conda environment
2. **Configure tools.json** with local database paths
3. **Convert gator_db.xlsx**: 
   ```r
   convert_gator_db("gator_db.xlsx", output_dir = "my_potatoes/")
   ```
4. **Run on genomes**:
   ```r
   initialize_potato_sack("my_project")
   sack <- create_sack("my_project")
   sack <- add_genomes(sack, "my_mags/*.faa")
   # Annotation workflow (to be implemented)
   ```
5. **Inspect results** (same output format as v1)

### Backward Compatibility

- Output file formats match gator v1 (TSV structure)
- Column names preserved where possible
- Additional columns for new features (confidence, tools_used, etc.)
- Scripts expecting v1 output should mostly work

---

## Next Steps

### Completed (v0.5.0 - v0.5.1)
1. ✅ Complete this roadmap
2. ✅ Simplify foundation - remove premature features
3. ✅ Establish core workflow (initialize, create, validate, add_genomes)
4. ✅ Create test suite (35 tests passing)
5. ✅ Standardize save/load (use R's saveRDS/readRDS)
6. ✅ Update all 7 example potatoes to new schema
7. ✅ Strict validation (databases field, standard types only)

### Next Priorities
1. [ ] Implement annotation workflow (annotate_simple.R skeleton exists)
2. [ ] Connect jakomics tool runners (kofam, blast, hmm)
3. [ ] Design output format (tibble with nested genes)
4. [ ] Scoring logic (simple fraction first, then DAG-aware)
5. [ ] Integration test with real genome
6. [ ] Iterate based on real usage

---

## Contributors

- Jeff Kimbrel (primary author)

Last updated: 2026-07-29
