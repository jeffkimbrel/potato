# POTATO Implementation Notes

## Key Revelations & Decisions

### 2026-07-23: Per-Gene Score Thresholds

**Problem:** GATOR v1 uses hardcoded thresholds (e.g., BLAST e-value 1e-15)

**Limitation:** Different genes have different requirements:
- Highly conserved genes → can use strict thresholds (e-value 1e-50)
- Divergent genes → need relaxed thresholds (e-value 1e-5, pident 25%)
- Context-dependent (MAG quality, phylogenetic distance)

**Solution:** Per-gene thresholds in potato JSON

```json
{
  "nodes": [
    {
      "id": "nifH",
      "ko": ["K02588"],
      "thresholds": {
        "kofam_score": 1.2,
        "blast_evalue": 1e-30,
        "blast_bitscore": 200,
        "blast_pident": 40,
        "hmm_evalue": 1e-20
      }
    }
  ]
}
```

**Precedence (most specific wins):**
1. Per-gene threshold in potato → `nodes[i].thresholds.blast_evalue`
2. Default in tools.json → `blast.default_thresholds.evalue`
3. Built-in fallback → hardcoded (e.g., 1e-15)

**Tool-Specific Thresholds:**
- **kofam**: score (ratio), evalue
- **blast**: evalue, bitscore, pident (% identity), qcovs (query coverage)
- **hmmer**: evalue, bitscore, domain_evalue

**Benefits:**
- ✅ Gene-specific stringency (strict for conserved, relaxed for divergent)
- ✅ Reproducibility (thresholds documented in potato)
- ✅ Flexibility (override per-gene when needed)
- ✅ User control (defaults in tools.json, overrides in potato)

**Implementation:** Tool runners check gene-specific thresholds first, fall back to defaults

---

### 2026-07-23: Self-Contained Database Embedding

**Idea:** Embed HMM/BLAST sequences directly in potato JSON files for true portability.

**Current limitation:**
- Potato references "K01652" (external identifier)
- User must have kofam profiles installed at path in tools.json
- Sharing potato requires recipient to set up databases

**Enhanced approach:**
- Potato embeds actual HMM text or BLAST sequences
- No external database needed
- Email a potato JSON, it just works

**Technical details:**
```json
{
  "nodes": [
    {
      "id": "mlrA",
      "detection_methods": [
        {
          "type": "hmm",
          "source": "embedded",
          "data": "HMMER3/f [3.3.2 | Nov 2020]\nNAME mlrA\nLENG 542\n...",
          "compressed": true
        }
      ]
    }
  ]
}
```

**Benefits:**
- Perfect for custom/niche pathways
- Version control (exact HMM preserved)
- No setup friction for potato recipients
- Works offline

**Challenges:**
- File size (5-50KB per HMM, mitigate with gzip)
- KEGG licensing (redistribution restrictions - need investigation)
- Redundancy (same HMM in multiple potatoes)
- Updates (improved HMM requires updating all potatoes)

**Implementation approach:**
- Optional, not required (hybrid mode: embed custom, reference standard)
- Fallback logic (use embedded if present, else external via tools.json)
- Compression by default (gzip HMM text)
- Licensing warnings in builder agent

**Status:** Proposed for Phase 5, pending licensing investigation

---

### 2026-07-22: Architecture Simplification

**Discovery:** We don't need to reimplement tool execution logic - jakomics already handles it!

**Old approach (rejected):**
```
potato → reticulate → inst/python/gator.py → jakomics
```

**New approach (cleaner):**
```
potato (R) → reticulate → potato.py (DAG/scoring only)
                       → jakomics (tool execution)
```

**What this means:**
- Remove `inst/python/gator.py`, `metadata.py`, `pathway.py` (old GATOR v1 code)
- Write new potato-specific code: `potato.py`, `scoring.py`
- Call jakomics functions directly:
  - `jakomics.kegg.run_kofam()` + `kegg.parse_kofam_hits()`
  - `jakomics.blast.run_blast()`
  - `jakomics.hmm.run_hmmsearch()` + `hmm.parse_hmm_hits()`
  - `jakomics.genome.GENOME` for genbank → FAA conversion
  - `jakomics.utilities.get_files()` for file discovery (already wrapped in R!)

**Benefits:**
- ✅ Don't maintain duplicate tool-running code
- ✅ jakomics handles messy subprocess/parsing details
- ✅ potato code focuses on DAG logic, scoring, LLM agents
- ✅ Updates to jakomics automatically benefit potato

### Logic to Extract from GATOR v1

**File handling (already available via jakomics):**
```python
# Genbank → FAA conversion
from jakomics.genome import GENOME
gbk = GENOME(genbank_path)
gbk.genbank_to_fasta(write_faa=output_faa)

# File discovery (already wrapped in R/io.R!)
from jakomics.utilities import get_files
genomes = get_files(file_list, in_dir, ["faa", "gb", "gbk", "gbff"])
```

**Tool execution patterns (call jakomics directly):**
```python
from jakomics import kegg, blast, hmm

# KofamScan
hits = kegg.run_kofam(faa_path=..., hal_path=..., ko_list=..., echo=False)
results = kegg.parse_kofam_hits(hits)

# BLAST
results = blast.run_blast(type="prot", q=faa_path, db=db_path, e=1e-15, make=False)

# HMMER
hmm.run_hmmsearch(faa_path, log_file, output_file, hmm_db, cut_tc=False, echo=False)
results = hmm.parse_hmm_hits(output_file)
```

**What potato needs to implement (new code):**
- Potato class (load/validate JSON, represent as DAG)
- Scoring engine (traverse DAG, handle OR branches, compute confidence)
- Output formatting (gene-level + pathway-level TSVs)
- LLM agents (builder, converter, analysis)

---

## Minimal Viable Test Strategy

**Philosophy:** Prove the architecture with the absolute simplest case before adding complexity.

### Step 0: The Minimal Test Case

**Setup:**
- ✅ One genome (FAA file, 50-100 proteins)
- ✅ One potato (3 genes, linear, no branches)
- ✅ One tool (kofam or blast)
- ✅ No LLM features
- ✅ No parallelization
- ✅ No file format conversion

**Test structure:**
```
tests/fixtures/
├── genomes/
│   └── test_genome.faa           # Small E. coli subset
├── potatoes/
│   └── test_linear.json          # geneA -> geneB -> geneC
├── tools/
│   └── test_tools.json           # Minimal kofam config
└── expected/
    ├── test_genome_genes.tsv     # Expected gene hits
    └── test_genome_pathway.tsv   # Expected pathway score
```

**Success criteria:**
- ✅ Loads potato JSON
- ✅ Runs kofam via jakomics
- ✅ Maps hits to potato genes
- ✅ Scores pathway (3/3 = present)
- ✅ Returns result

### Incremental Complexity Ladder

Once minimal case works, add one feature at a time:

1. ✅ **One genome, one potato, one tool** ← START HERE
2. ✅ Add: Proper output files (TSV writing)
3. ✅ Add: Multiple tools (kofam + blast)
4. ✅ Add: OR branches in pathway (geneA -> (geneB | geneC) -> geneD)
5. ✅ Add: Proper DAG traversal (networkx-based)
6. ✅ Add: Required vs. optional genes
7. ✅ Add: Multiple potatoes in one run
8. ✅ Add: Multiple genomes (parallelization)
9. ✅ Add: Genbank input files (conversion via jakomics)
10. ✅ Add: Gene specificity weighting
11. ✅ Add: Confidence scoring
12. ✅ Add: LLM agents (converter, builder, analysis)

**Why this order?**
- Core logic first (potato loading, tool execution, scoring)
- Complexity adds gradually (OR branches, then full DAG)
- Output/UX next (files, multiple genomes)
- Advanced features last (specificity, LLM)

---

## File Removal Plan

### Files to Remove (Old GATOR v1 Code)

```
inst/python/gator.py      # Old CLI tool - replaced by new architecture
inst/python/metadata.py   # Excel loading - replaced by JSON loading
inst/python/pathway.py    # String parsing - replaced by DAG scoring
```

### Files to Create (New POTATO v1 Code)

**Phase 0 (Minimal):**
```
inst/python/potato_minimal.py    # Bare minimum: load potato, run tool, score
tests/fixtures/                  # Test data
tests/testthat/test_minimal.R    # R test for minimal case
```

**Phase 1 (Full MVP):**
```
inst/python/potato.py            # Potato class, JSON loading, validation
inst/python/tools.py             # Tool runner wrappers (calls jakomics)
inst/python/scoring.py           # DAG traversal, pathway scoring
inst/schemas/potato.schema.json  # JSON Schema for potato format
inst/schemas/tools.schema.json   # JSON Schema for tools config
```

**Phase 2+:**
```
inst/python/agents.py            # LLM agent code (Claude API)
inst/potatoes/*.json             # Converted potato files
```

---

## Python Code Patterns

### Calling jakomics for Tool Execution

```python
from jakomics import kegg, blast, hmm
from jakomics.genome import GENOME
from jakomics.utilities import get_files

class ToolRunner:
    """Wrapper for jakomics tool functions"""
    
    def __init__(self, tools_config):
        self.config = tools_config
    
    def get_threshold(self, tool_name, threshold_name, gene_thresholds=None):
        """Get threshold with precedence: gene > tools.json > built-in"""
        # Check gene-specific first
        if gene_thresholds and threshold_name in gene_thresholds:
            return gene_thresholds[threshold_name]
        
        # Check tools.json defaults
        if 'default_thresholds' in self.config[tool_name]:
            defaults = self.config[tool_name]['default_thresholds']
            if threshold_name in defaults:
                return defaults[threshold_name]
        
        # Built-in fallbacks
        fallbacks = {
            'kofam_score': 1.0,
            'kofam_evalue': 1e-10,
            'blast_evalue': 1e-15,
            'blast_bitscore': 50,
            'blast_pident': 30,
            'blast_qcovs': 50,
            'hmm_evalue': 1e-15,
            'hmm_bitscore': 50,
            'hmm_domain_evalue': 1e-10
        }
        return fallbacks.get(f"{tool_name}_{threshold_name}", None)
    
    def run_kofam(self, faa_path, ko_list, gene_thresholds=None):
        """Run KofamScan via jakomics"""
        cfg = self.config['kofam']
        if not cfg['enabled']:
            return {}
        
        # Create temporary HAL file (list of HMM paths)
        hal_path = self._create_hal_file(ko_list, cfg['ko_profiles_dir'])
        
        hits = kegg.run_kofam(
            faa_path=faa_path,
            hal_path=hal_path,
            temp_dir="temp_kofam",
            ko_list=cfg['ko_list'],
            score_as_ratio=False,
            echo=False
        )
        
        return kegg.parse_kofam_hits(hits)
    
    def run_blast(self, faa_path, db_path, gene_thresholds=None):
        """Run BLASTP via jakomics"""
        if not self.config['blast']['enabled']:
            return {}
        
        # Get thresholds
        evalue = self.get_threshold('blast', 'evalue', gene_thresholds)
        
        results = blast.run_blast(
            type="prot",
            q=faa_path,
            db=db_path,
            e=evalue,
            make=False,
            return_query_results=False
        )
        
        # Filter by additional thresholds (bitscore, pident, qcovs)
        bitscore_cutoff = self.get_threshold('blast', 'bitscore', gene_thresholds)
        pident_cutoff = self.get_threshold('blast', 'pident', gene_thresholds)
        
        filtered = {}
        for gene_id, hits in results.items():
            filtered[gene_id] = [
                hit for hit in hits
                if hit.bitscore >= bitscore_cutoff and 
                   hit.pident >= pident_cutoff
            ]
        
        return filtered
    
    def run_hmmer(self, faa_path, hmm_path):
        """Run HMMER3 via jakomics"""
        if not self.config['hmmer']['enabled']:
            return {}
        
        import tempfile
        import os
        
        log_file = tempfile.mktemp(suffix=".log")
        output_file = tempfile.mktemp(suffix=".out")
        
        hmm.run_hmmsearch(
            faa_path, log_file, output_file, hmm_path,
            cut_tc=False, echo=False, run=True
        )
        
        results = hmm.parse_hmm_hits(output_file)
        
        # Cleanup
        os.remove(log_file)
        os.remove(output_file)
        
        return results
    
    def convert_genbank(self, gbk_path, output_faa):
        """Convert genbank to FAA via jakomics"""
        gbk = GENOME(gbk_path)
        gbk.genbank_to_fasta(write_faa=output_faa, return_gene_dict=False)
        return output_faa
```

### Potato Class Structure

```python
import json
import networkx as nx

class Potato:
    """Represents a pathway definition from JSON"""
    
    def __init__(self, json_path):
        with open(json_path) as f:
            self.data = json.load(f)
        
        self.id = self.data['id']
        self.name = self.data['name']
        self.tags = self.data.get('tags', [])
        self.nodes = self.data['nodes']
        self.edges = self.data['edges']
        self.scoring = self.data.get('scoring', {})
        
        # Build DAG
        self.graph = self._build_graph()
        
        # Validate
        self.validate()
    
    def _build_graph(self):
        """Convert nodes/edges to networkx DiGraph"""
        G = nx.DiGraph()
        
        # Add nodes
        for node in self.nodes:
            G.add_node(node['id'], **node)
        
        # Add edges
        for edge in self.edges:
            G.add_edge(edge['from'], edge['to'], **edge)
        
        return G
    
    def validate(self):
        """Check for common errors"""
        # Check for cycles
        if not nx.is_directed_acyclic_graph(self.graph):
            raise ValueError(f"Potato {self.id} contains cycles!")
        
        # Check required nodes exist
        required = [n for n in self.nodes if n.get('required', False)]
        for node in required:
            if node['id'] not in self.graph:
                raise ValueError(f"Required node {node['id']} not in graph")
        
        # Check edge endpoints exist
        for edge in self.edges:
            if edge['from'] not in self.graph:
                raise ValueError(f"Edge from {edge['from']} - node doesn't exist")
            if edge['to'] not in self.graph:
                raise ValueError(f"Edge to {edge['to']} - node doesn't exist")
    
    def get_gene_defs(self):
        """Return list of gene definitions for annotation"""
        return [n for n in self.nodes if n['type'] == 'enzyme']
    
    def get_detection_terms(self, tool_type):
        """Get all KO/PFAM/etc terms for a specific tool"""
        terms = []
        for node in self.nodes:
            if node['type'] == 'enzyme' and tool_type in node:
                terms.extend(node[tool_type])
        return list(set(terms))
```

---

## R Code Patterns

### Loading Potatoes

```r
#' Load a potato JSON file
#' @param path Path to potato JSON file
#' @return Potato object (Python via reticulate)
#' @export
load_potato <- function(path) {
  if (!file.exists(path)) {
    stop("Potato file not found: ", path, call. = FALSE)
  }
  
  potato_py <- reticulate::import_from_path("potato", "inst/python")
  potato_py$Potato(path)
}

#' Load all potatoes from a directory
#' @param dir Directory containing potato JSON files
#' @param tags Optional: filter by tags
#' @return List of potato objects
#' @export
load_potatoes <- function(dir, tags = NULL) {
  json_files <- list.files(dir, pattern = "\\.json$", full.names = TRUE)
  
  potatoes <- lapply(json_files, load_potato)
  names(potatoes) <- sapply(potatoes, function(p) p$id)
  
  # Filter by tags if specified
  if (!is.null(tags)) {
    potatoes <- Filter(function(p) {
      any(tags %in% p$tags)
    }, potatoes)
  }
  
  potatoes
}
```

### Running Annotation

```r
#' Annotate a genome against potatoes
#' @param faa_path Path to FAA file (or genbank)
#' @param potatoes List of potato objects
#' @param tools_config Path to tools.json
#' @return List with $genes and $pathways dataframes
#' @export
annotate_genome <- function(faa_path, potatoes, tools_config) {
  
  # Load Python backend
  potato_py <- reticulate::import_from_path("potato", "inst/python")
  tools_py <- reticulate::import_from_path("tools", "inst/python")
  
  # Convert genbank if needed
  if (grepl("\\.(gb|gbk|gbff)$", faa_path)) {
    runner <- tools_py$ToolRunner(tools_config)
    faa_path <- runner$convert_genbank(faa_path, 
                                       tempfile(fileext = ".faa"))
  }
  
  # Run tools
  runner <- tools_py$ToolRunner(tools_config)
  annotation_results <- list()
  
  for (potato in potatoes) {
    # Get what to search for
    ko_list <- potato$get_detection_terms('ko')
    
    # Run kofam
    if (length(ko_list) > 0) {
      annotation_results$kofam <- runner$run_kofam(faa_path, ko_list)
    }
  }
  
  # Score pathways
  pathway_results <- lapply(potatoes, function(p) {
    scoring_py <- reticulate::import_from_path("scoring", "inst/python")
    scoring_py$score_pathway(p, annotation_results)
  })
  
  list(
    genes = format_gene_results(annotation_results),
    pathways = format_pathway_results(pathway_results)
  )
}
```

---

## Next Immediate Actions

### 1. Create Minimal Test Fixtures

```bash
mkdir -p tests/fixtures/{genomes,potatoes,tools}
```

**Create `tests/fixtures/potatoes/test_linear.json`:**
```json
{
  "id": "test_linear",
  "name": "Test Linear Pathway",
  "tags": ["test"],
  "nodes": [
    {
      "id": "gapA",
      "type": "enzyme",
      "name": "glyceraldehyde-3-phosphate dehydrogenase",
      "ko": ["K00134"],
      "required": true
    },
    {
      "id": "pgk",
      "type": "enzyme",
      "name": "phosphoglycerate kinase",
      "ko": ["K00927"],
      "required": true
    },
    {
      "id": "eno",
      "type": "enzyme",
      "name": "enolase",
      "ko": ["K01689"],
      "required": true
    }
  ],
  "edges": [
    {"from": "gapA", "to": "pgk"},
    {"from": "pgk", "to": "eno"}
  ],
  "scoring": {
    "min_fraction": 1.0
  }
}
```

**Create `tests/fixtures/tools/test_tools.json`:**
```json
{
  "kofam": {
    "enabled": true,
    "ko_profiles_dir": "/path/to/kofam/profiles",
    "ko_list": "/path/to/kofam/ko_list"
  }
}
```

### 2. Implement Minimal Python Code

Create `inst/python/potato_minimal.py` with:
- Potato class (load JSON, basic validation)
- ToolRunner class (call jakomics.kegg.run_kofam)
- score_pathway_simple function (count found genes)

### 3. Implement Minimal R Wrapper

Create `R/minimal.R` with:
- `test_minimal()` function
- Loads potato, runs tool, scores pathway
- Returns simple result

### 4. Test Manually

```r
devtools::load_all()
result <- test_minimal(
  "tests/fixtures/genomes/test_genome.faa",
  "tests/fixtures/potatoes/test_linear.json", 
  "tests/fixtures/tools/test_tools.json"
)
```

### 5. Once Working, Add One Feature at a Time

Follow the incremental complexity ladder from Step 0 above.

---

Last updated: 2026-07-22
