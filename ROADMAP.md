# POTATO v1 - Development Roadmap

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
  "tags": ["amino_acid_metabolism", "branched_chain"],
  "notes": "Canonical leucine biosynthesis pathway",
  
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

### Phase 0: Foundation ✓ (Partially Complete)

- [x] R package skeleton
- [x] Reticulate bridge to Python
- [x] Conda environment with bioinformatics tools
- [ ] Define potato JSON schema formally (JSON Schema file)
- [ ] Define tools.json schema

### Phase 0.5: Minimal Viable Test (NEW)

**Goal:** Prove the architecture with the absolute simplest possible case.

**Test Case:**
- ✅ One genome (FAA file, ~50 proteins)
- ✅ One potato (3 genes, linear pathway: A → B → C)
- ✅ One tool (kofam via jakomics)
- ✅ No LLM, no parallelization, no complexity

**Why Start Here:**
Validate that the core pieces connect:
1. Load potato JSON → Python object
2. Extract KO terms from potato
3. Run kofam via jakomics.kegg.run_kofam()
4. Map results back to potato genes
5. Score pathway (3/3 = present)
6. Return result to R

**Files to Create:**
```
tests/fixtures/
├── genomes/test_genome.faa          # Small test genome
├── potatoes/test_linear.json        # 3 genes: gapA → pgk → eno
└── tools/test_tools.json            # Minimal kofam config

inst/python/potato_minimal.py        # Bare minimum code
R/minimal.R                          # R wrapper for test
```

**Success Criteria:**
```r
result <- test_minimal(
  "tests/fixtures/genomes/test_genome.faa",
  "tests/fixtures/potatoes/test_linear.json",
  "tests/fixtures/tools/test_tools.json"
)
# Returns: list(pathway="test_linear", found=3, total=3, present=TRUE)
```

**Incremental Complexity Ladder:**
Once minimal works, add one feature at a time:
1. ✅ One genome, one potato, one tool ← START HERE
2. Output files (TSV writing)
3. Multiple tools (kofam + blast)
4. OR branches in pathway
5. Proper DAG traversal (networkx)
6. Required vs. optional genes
7. Multiple potatoes
8. Multiple genomes (parallelization)
9. Genbank input (conversion)
10. Gene specificity weighting
11. LLM agents

---

### Phase 1: Core Functionality (MVP)

**Goal:** Run a single potato against a single genome without LLM features.

#### 1.1 Data Structures
- [ ] Python `Potato` class - loads/validates JSON, represents DAG
- [ ] Python `ToolConfig` class - loads/validates tools.json
- [ ] R wrapper functions for loading potatoes
- [ ] DAG validation (no cycles, required nodes exist, etc.)

#### 1.2 Tool Execution
- [ ] Python tool runners:
  - [ ] KofamScan wrapper
  - [ ] HMMER3 wrapper  
  - [ ] BLASTP wrapper
  - [ ] PFAM wrapper (via HMMER)
- [ ] R `annotate_genome()` function
  - Takes: genome FAA file, potato, tools.json
  - Returns: annotation results (which genes detected by which tools)

#### 1.3 Scoring Engine
- [ ] DAG traversal for pathway scoring
  - [ ] Handle OR branches (ilvH | ilvM)
  - [ ] Handle required vs. optional nodes
  - [ ] Calculate min_fraction thresholds
- [ ] Output formats:
  - [ ] Gene-level: `genome_gator.tsv` (gene, product, tool, score, locus_tag)
  - [ ] Pathway-level: `genome_potato.tsv` (pathway, present, steps, confidence)

#### 1.4 Basic Testing
- [ ] Create 3-5 test potatoes (leucine biosynthesis, nitrogen fixation, etc.)
- [ ] Test genomes (known MAGs with varying completeness)
- [ ] Unit tests for scoring logic
- [ ] Integration test: genome → annotation → scoring → output

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
- [ ] R function to plot potato DAG with `igraph`
- [ ] Highlight found vs. missing nodes per genome
- [ ] Export to graphviz DOT format
- [ ] Generate pathway summary reports (HTML/PDF)

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

### Phase 4: LLM Integration

**Goal:** LLM agents for database building and result interpretation.

#### 4.1 Builder Agent

Helps create new potato JSON files.

**Mode 1: From KEGG Module**
```r
build_potato(kegg_module = "M00570")
```
- [ ] Fetch from KEGG REST API
- [ ] Parse KEGG module definition syntax
- [ ] Extract KO IDs, pathway structure, compound info
- [ ] Generate potato JSON draft
- [ ] User reviews/edits, saves

**Mode 2: From Gene List**
```r
build_potato(genes = c("rbdA", "rbdB", "rbdC", "rbdD"), 
             pathway_name = "R-body synthesis")
```
- [ ] LLM searches for KO/PFAM/EC for each gene
- [ ] Suggests pathway logic
- [ ] Generates potato JSON draft

**Mode 3: Interactive**
```r
build_potato(prompt = "I want a shikimate pathway potato")
```
- [ ] Conversational refinement
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
   results <- annotate_genomes(
     genome_dir = "my_mags/",
     potato_dir = "my_potatoes/",
     tools_config = "tools.json"
   )
   ```
5. **Inspect results** (same output format as v1)

### Backward Compatibility

- Output file formats match gator v1 (TSV structure)
- Column names preserved where possible
- Additional columns for new features (confidence, tools_used, etc.)
- Scripts expecting v1 output should mostly work

---

## Next Steps

1. ✅ Complete this roadmap (done!)
2. [ ] Formalize JSON schemas (potato, tools.json)
3. [ ] Implement Phase 1.1 (Python data structures)
4. [ ] Create 3 test potatoes manually
5. [ ] Implement Phase 1.2 (tool execution)
6. [ ] Build scoring engine (Phase 1.3)
7. [ ] Integration test with real genome
8. [ ] Iterate...

---

## Contributors

- Jeff Kimbrel (primary author)
- Claude (design discussions, implementation assistance)

Last updated: 2026-07-22
