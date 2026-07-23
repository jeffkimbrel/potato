# POTATO - Project Context for Claude

## Project Overview

**POTATO** (Pathway annOTATOr) is an R package for annotating MAGs (metagenome-assembled genomes) against curated metabolic pathways. It's the successor to GATOR (Genome annotATOR), redesigned around self-contained pathway definitions (potatoes) as DAG structures in JSON.

**Current Status:** Foundation exists (R package skeleton, reticulate setup, conda env). Starting Phase 1 implementation.

**Key Innovation:** Each "potato" (pathway) is a self-contained JSON file defining:
- Genes with multi-tool detection methods (KEGG, PFAM, BLAST, HMM)
- Pathway structure as a DAG (handles complex branching)
- Confidence scoring via gene specificity weighting
- LLM-assisted building and interpretation

---

## Important Files

- **[ROADMAP.md](ROADMAP.md)** - Complete implementation plan, phases, file formats
- **[BRAINSTORM.md](BRAINSTORM.md)** - Design evolution, key insights, dead ends
- **[environment.yaml](environment.yaml)** - Conda environment with bioinformatics tools
- **R/zzz.R** - Reticulate setup, loads Python backend on package load
- **inst/python/** - Python backend code (will be rewritten for v1)
- **inst/extdata/potato_db.xlsx** - Old GATOR v1 database (will be converted)

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
    {"id": "geneA", "type": "enzyme", "name": "Gene A", "ko": ["K00001"], "required": true},
    {"id": "geneB", "type": "enzyme", "name": "Gene B", "ko": ["K00002"], "required": true},
    {"id": "geneC", "type": "enzyme", "name": "Gene C", "ko": ["K00003"], "required": true}
  ],
  "edges": [
    {"from": "geneA", "to": "geneB"},
    {"from": "geneB", "to": "geneC"}
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
            if 'ko' in node:
                kos.extend(node['ko'])
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

- **v0.0.5** - Current, foundation exists (reticulate, conda env)
- **v1.0.0** - Target for Phase 1 MVP completion

---

Last updated: 2026-07-22
