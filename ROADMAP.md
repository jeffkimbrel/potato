# POTATO v1 - Development Roadmap

**Current Version:** v0.8.0 (2026-07-30)

**Status:** Complete annotation pipeline with scoring, visualization, text-based pathway views, and LLM-assisted potato building. Ready for real-world testing and potato database expansion.

**What works:** Annotate genomes with kofam/BLAST/HMM → score pathways → visualize results → analyze near-misses

**What's next:** Gene specificity weighting, threshold sensitivity analysis, more verified potatoes

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

### 1. Potato JSON Structure

Each potato is a self-contained pathway definition stored as a DAG.

```json
{
  "id": "leucine_biosynthesis",
  "name": "Leucine Biosynthesis",
  "source": "KEGG M00570",
  "verified": false,
  "tags": ["amino_acid_metabolism", "branched_chain"],
  "notes": "Canonical leucine biosynthesis pathway",
  
  "input": {
    "compound": "pyruvate",
    "kegg_compound": "C00022",
    "targets": ["ilvB_1"]
  },
  
  "output": {
    "compound": "L-leucine",
    "kegg_compound": "C00123",
    "sources": ["leuA_5"]
  },
  
  "nodes": [
    {
      "id": "ilvB",
      "type": "enzyme",
      "name": "acetolactate synthase catalytic subunit",
      "ko": ["K01652", "K01653"],
      "pfam": ["PF00205"],
      "ec": ["2.2.1.6"],
      "hmm": ["ilvB.hmm"],
      "blast_db": "ilvB_refs.faa",
      "required": true,
      "specificity": 0.75,
      "appears_in": ["leucine_biosynthesis", "valine_biosynthesis", "isoleucine_biosynthesis"],
      "thresholds": {
        "kofam_score": 1.0,
        "kofam_evalue": 1e-10,
        "blast_evalue": 1e-20,
        "blast_bitscore": 100,
        "blast_pident": 30,
        "hmm_evalue": 1e-15,
        "hmm_bitscore": 50,
        "hmm_domain_evalue": 1e-10
      },
      "reaction": {
        "type": "C-C bond formation",
        "substrate": "pyruvate",
        "product": "2-acetolactate"
      }
    },
    {
      "id": "ilvH",
      "type": "enzyme",
      "name": "acetolactate synthase small subunit",
      "ko": ["K01654"],
      "required": false,
      "specificity": 0.80,
      "thresholds": {
        "kofam_score": 1.0,
        "blast_evalue": 1e-15
      }
    },
    {
      "id": "ilvM",
      "type": "enzyme", 
      "name": "acetolactate synthase small subunit (alternative)",
      "ko": ["K11258"],
      "required": false,
      "specificity": 0.80
    },
    {
      "id": "ilvC",
      "type": "enzyme",
      "name": "ketol-acid reductoisomerase",
      "ko": ["K00053"],
      "pfam": ["PF07991"],
      "ec": ["1.1.1.86"],
      "required": true,
      "specificity": 0.85
    }
  ],
  
  "edges": [
    {
      "from": "ilvB",
      "to": "ilvH",
      "compound": "complex formation",
      "kegg_compound": null
    },
    {
      "from": "ilvB", 
      "to": "ilvM",
      "compound": "complex formation (alternative)",
      "kegg_compound": null
    },
    {
      "from": "ilvH",
      "to": "ilvC",
      "compound": "2-acetolactate",
      "kegg_compound": "C06010"
    },
    {
      "from": "ilvM",
      "to": "ilvC", 
      "compound": "2-acetolactate",
      "kegg_compound": "C06010"
    }
  ],
  
  "scoring": {
    "min_fraction": 0.75,
    "required_nodes": ["ilvB", "ilvC"],
    "use_specificity_weighting": true
  }
}
```

**Key Features:**
- **Nodes**: Enzyme definitions with detection methods (KO, PFAM, EC, HMM, BLAST)
- **Edges**: Connections between enzymes, decorated with metabolite info
- **Scoring**: Pathway-specific thresholds and requirements
- **Specificity**: Pre-computed across potato database (how unique is this gene?)
- **Thresholds**: Per-gene score cutoffs (e-value, bitscore, percent identity, etc.)
- **Reaction metadata**: Enables LLM reasoning about functional analogs (optional)

**Threshold Philosophy:**
- Each gene can specify custom thresholds for each detection method
- If not specified, use sensible defaults from tools.json or built-in
- Allows stricter thresholds for high-confidence genes (e.g., `blast_evalue: 1e-50`)
- Allows looser thresholds for divergent genes (e.g., `blast_pident: 25`)
- Tool-specific: kofam score, blast e-value/bitscore/pident, HMM e-value/bitscore

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
  - Returns: tibble with `distance_from_threshold`, `steps_detected`, `steps_total`
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

### Phase 2: Database Management

**Goal:** Work with collections of potatoes, validate consistency.

#### 2.1 Multi-Potato Support
- [ ] Load all potatoes from `inst/potatoes/` directory
- [ ] Filter by tags: `run_potatoes(tags = c("amino_acid_metabolism"))`
- [ ] Batch scoring across multiple potatoes

#### 2.2 Consistency Validation
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

### Phase 3: Confidence Scoring

**Goal:** Interpret results in context of gene specificity and MAG completeness.

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

#### 3.4 Threshold Sensitivity Analysis ⚠️ HIGH-PRIORITY, HIGH-LOAD

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

#### 4.3 Analysis Agent

Interprets results per-genome.

```r
interpret_results(genome_results, genome_completeness = 0.78)
```
- [ ] Reads potato scoring results
- [ ] Considers specificity, MAG completeness
- [ ] Generates natural language interpretation
- [ ] Suggests functional analogs for missing genes (Phase 5)

#### 4.4 LLM Infrastructure
- [ ] Claude API integration (via Anthropic SDK)
- [ ] Prompt templates for each agent mode
- [ ] Token management, caching
- [ ] Error handling, retry logic

### Phase 5: Advanced Features (Future)

#### 5.0 Self-Contained Database Embedding

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

#### 5.1 Functional Analog Suggestions
- [ ] LLM uses reaction metadata to suggest alternatives
- [ ] "Missing ketol-acid reductoisomerase (EC 1.1.1.86). Found hypothetical protein with EC 1.1.1.* domain. Possible functional analog?"
- [ ] Flag as unvalidated, requires expert review

#### 5.2 Cross-Potato Queries
- [ ] "Which potatoes produce acetyl-CoA?"
- [ ] "Show me all pathways downstream of pyruvate"
- [ ] Network visualization across potatoes

#### 5.3 ESM-2 Embeddings (Protein LM)
- [ ] Sequence-similarity rescue for missed annotations
- [ ] BacPT/Bacformer as pre-computed embedding database
- [ ] ANN search for functional annotation

#### 5.4 DAG Evolution Suggestions
- [ ] LLM: "I know of a bifunctional enzyme that shortcuts A → B → C to A → C"
- [ ] Suggest new edges/nodes to existing potatoes
- [ ] User reviews, accepts/rejects

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
