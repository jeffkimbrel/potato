# POTATO - Development Roadmap

**Current Version:** v0.10.2-dev (2026-08-14)

**Status:** V2 schema complete. Full annotation pipeline (kofam/BLAST/HMM) → scoring (all/essential genes) → visualization (interactive/static) → provenance tracking

**Up Next:**

1. Potato locking - prevent modification when results exist
2. Incremental annotation - merge results from multiple runs
3. Threshold sensitivity analysis

---

## Vision

POTATO (Pathway annOTATOr) annotates genome collections (MAGs) against curated metabolic pathways. Core innovation: **self-contained potato files** (pathway definitions as DAG structures in JSON) with **LLM-assisted analysis**.

### Design Principles

1. **Self-contained potatoes** - each pathway is independent, portable
2. **DAG-based pathway logic** - intuitive, graphable, handles complex branching
3. **Tool-agnostic detection** - works with whatever annotation tools available
4. **Confidence-aware scoring** - interprets results in context of gene specificity and MAG completeness
5. **LLM augmentation** - assists in database building and result interpretation

---

## V2 Schema Reference

```json
{
  "schema_version": "v2",
  "id": "pathway_id",
  "name": "Human Readable Name",
  "source": "KEGG M00123 / custom",
  "tags": ["metabolism"],
  "notes": "Brief description",
  
  "genes": [
    {
      "id": "geneSymbol",
      "name": "enzyme name",
      "databases": {
        "kofam": ["K00001"],
        "blast": ["ref_seq_id"],
        "hmm": ["PF00001"]
      },
      "ec": ["1.1.1.1"],
      "reactions": ["R00001"]
    }
  ],
  
  "compounds": [
    {"id": "C00031", "name": "D-glucose"}
  ],
  
  "pathways": {
    "pathway_id": {
      "name": "Pathway Display Name",
      "type": "variant",  // or "independent"
      "verified": false,
      "edges": [
        {
          "from": "C00031",
          "to": "geneA",
          "required": true,
          "marker": true,
          "reaction": "R00001"
        }
      ],
      "scoring": {
        "min_fraction": 0.75,
        "marker_mode": "any"
      }
    }
  }
}
```

**Key features:**
- `genes`: Global definitions (detection methods only)
- `pathways`: One or more pathways with topology
- `edges`: Connect genes/compounds, carry `required`/`marker`
- `verified`: Always `false` for new pathways (only humans verify)
- Pathway types: `variant` (alternatives) or `independent` (different functions)

---

## Data Integrity & Workflow

### Gap-Based Scoring System **[HIGH PRIORITY]**

**Problem:** Fraction-based scoring has fundamental limitations:
- Float precision issues (0.67 threshold fails 2/3 genes at 0.6666...)
- Doesn't distinguish complete vs. incomplete routes in multi-pathway networks
- Ignores DAG connectivity (80% with middle gap fails, 60% with full path succeeds)

**Solution:** Add gap-based scoring as alternative to fraction-based

**Phase 1: Schema Update**
- [ ] Add `scoring.method` field: "fraction" (default), "gaps", or "both"
- [ ] Add `scoring.max_gaps` field: integer, number of allowed missing genes in best path
- [ ] Update validation to support both methods
- [ ] Both methods optional; fall back to run-level defaults if not specified
- [ ] Document in schema reference

**Phase 2: Scoring Implementation**
- [ ] Implement DAG traversal to find best path (fewest gaps)
- [ ] Calculate gap count for each pathway
- [ ] When `method: "both"`, calculate and return both scores
- [ ] Update results tibble to include: `gap_count`, `gap_present` (boolean)
- [ ] Maintain backward compatibility with existing fraction-only potatoes

**Example potato scoring block:**
```json
"scoring": {
  "method": "gaps",
  "max_gaps": 1,
  "min_fraction": 0.67  // ignored when method = "gaps"
}
```

**Benefits:**
- Biologically meaningful ("allow 1 missing gene")
- Handles alternative routes correctly
- No float precision issues
- True DAG traversal (checks connectivity)

---

### Potato Locking System **[HIGH PRIORITY]**

**Problem:** User edits potato JSON after annotation → hash mismatch → results meaningless

**Solution:**

- [ ] Detect hash mismatches on sack load
- [ ] Add `potatoes_locked` flag to PotatoSack
- [ ] Check lock status in modification functions
- [ ] Clear error: revert potato or clear results

**Why critical:** Prevents data corruption, results always match definitions

---

### Incremental Annotation **[HIGH PRIORITY]**

**Problem:** Running annotation twice errors unless `overwrite = TRUE` (replaces ALL results). Can't add potatoes incrementally.

**Solution:**

- [ ] Detect conflicts (genome + potato_hash)
- [ ] Merge logic: skip conflicts or replace based on `overwrite`
- [ ] Update provenance for merged runs
- [ ] Apply to kofam, BLAST, HMM

**Use case:**

```r
sack <- run_kofam(sack, potato_names = "pathway_A")
sack <- run_kofam(sack, potato_names = "pathway_B")  # Should merge
```

---

### Threshold Sensitivity Analysis **[MEDIUM PRIORITY]**

**Goal:** Understand how threshold choices affect results

**Implementation:**

- [ ] `analyze_thresholds(sack)` sweeps thresholds
- [ ] Show which pathways flip present/absent
- [ ] Identify threshold-sensitive vs robust pathways
- [ ] Visualization: completion vs threshold
- [ ] Export results as tibble

---

## Scoring & Analysis

### Specificity Weighting **[MEDIUM PRIORITY]**

**Goal:** Weight genes by specificity - finding `nifH` (pathway-specific) is more informative than `gapA` (ubiquitous)

**Implementation:**

- [ ] Pre-compute specificity: `1 / (num potatoes containing gene)`
- [ ] Weighted scoring: `sum(specificity of found) / sum(specificity of all)`
- [ ] Output: raw score + weighted score + high-value genes

**Context-Aware Threshold Adjustment:**

- [ ] Query KEGG API for pathway membership per KO
- [ ] Identify pathway-specific marker genes (only in one pathway)
- [ ] Relax thresholds for borderline hits when marker genes present
- [ ] Example: `gcl` scores 450 (threshold 678), BUT `glxK` (specific marker) present → likely real pathway

**Use case:** Genome has glxR (ubiquitous) + glxK (specific marker only in glycerate pathway) + borderline gcl hit → context suggests gcl is real, not a paralog

---

### Input/Output Validation **[MEDIUM PRIORITY]**

**Goal:** Pathways aren't "complete" just by gene count - must transform **available inputs → desired outputs**

**Current:** Input/output optional metadata  
**New:** Input/output mandatory, validated for reachability

**Implementation:**

- [ ] Schema: require input/output for all pathways
- [ ] Validate biological completeness
- [ ] Compound matching via KEGG IDs
- [ ] Check input reachability (transport? synthesis? environmental?)
- [ ] DAG traversal: can you reach output from input with detected genes?

**Example:**

```
npED alone: 5/6 genes (83%) → but no gluconate source → NON-FUNCTIONAL
ED network: transport + npED → gluconate available → FUNCTIONAL
```

---

### Cross-Potato Metabolic Networks **[MEDIUM PRIORITY]**

**Goal:** Evaluate pathways in context of whole-genome metabolism

**Biological reality:** Bacterial metabolism is one cytoplasmic soup - intermediates from one pathway immediately available for another

**Implementation:**

- [ ] Build genome-wide metabolite pool (available, produced, consumed)
- [ ] Compound matching via KEGG IDs
- [ ] `score_metabolic_network(sack, genome)` evaluates in context
- [ ] "ED produces glyceraldehyde → serine biosynthesis uses it → both COMPLETE"

---

### MAG Completeness Integration **[LOW PRIORITY]**

**Goal:** Adjust confidence based on genome completeness

**Implementation:**

- [ ] Accept CheckM/BUSCO estimates
- [ ] Adjust: "3 genes missing, MAG 60% complete"
- [ ] Flag: likely present / likely absent / uncertain

---

### Normalized Margin Calculation **[LOW PRIORITY]**

**Goal:** Make threshold margins comparable across different metric types in `inspect_gene_thresholds()`

**Problem:** E-value margins (1e-11) vs bitscore margins (10) vs score margins (50) are on completely different scales

**Solution:**

- [ ] Add percent-of-threshold calculation: `(score/threshold) * 100` or `(threshold/evalue) * 100`
- [ ] Values >100% = passing, <100% = failing
- [ ] Allows meaningful sorting across different metric types
- [ ] Alternative: add margin_type column and sort within types only

---

## Visualization

### V2 Coordinate System **[MEDIUM PRIORITY]**

**Problem:** V1 had coordinate export → curated layouts. Not carried to V2. Workflow broken.

**Solution:**

- [ ] Add export button to `plot_v2_interactive()`
- [ ] Update `update_potato_coordinates()` for V2
- [ ] Update `align_coordinates()` for V2
- [ ] Normalize pixel coords to plot coords

**Workflow:**

```r
plot_v2(pot, interactive = TRUE, layout = "fr")
# Arrange nodes → Export → Import
update_potato_coordinates("pathway.json", "coords.json")
```

---

### Near-Miss Pathway Visualization **[LOW PRIORITY]**

**Goal:** Visual way to identify pathways close to detection threshold

**Challenge:** Doesn't scale well - faceting by genome fails with 100+ genomes, faceting by pathway fails with dozens of pathways

**Current workaround:** `find_near_miss_pathways()` returns filtered tibble, `plot_genome_pathways()` shows per-genome view with thresholds

**Future idea:** Interactive filtering or small-multiples view for targeted subsets

---

## LLM Integration

### Build-Potato Agent **[MEDIUM PRIORITY]**

**Goal:** LLM-assisted potato creation from KEGG modules or literature

**Implementation:**

- [ ] Parse KEGG module structure → potato draft
- [ ] Extract gene IDs, EC numbers, compounds
- [ ] Handle complex KEGG syntax (branches, alternatives)
- [ ] Suggest BLAST/HMM when KEGG insufficient
- [ ] User reviews/edits → validates

---

### Result Interpretation Agent **[LOW PRIORITY]**

**Goal:** Natural language summaries of results

**Implementation:**

- [ ] Input: scores + optional MAG completeness
- [ ] Output: "Genome has complete ED (classic variant). Missing edd but MAG 60% complete, likely present."
- [ ] Flag uncertainty, suggest follow-up

---

## Advanced Features

### Embed Detection Data in Potatoes **[LOW PRIORITY]**

**Goal:** Ultimate portability - potato contains all HMM/BLAST data

**Implementation:**

- [ ] Optional fields: `hmm_profile_data`, `blast_sequences`
- [ ] Makes tools.json optional
- [ ] Build custom databases from scratch

---

### Multi-Sack Comparisons **[LOW PRIORITY]**

**Goal:** Compare annotation results across projects/studies

**Implementation:**

- [ ] `compare_sacks(sack_list)`
- [ ] Shared vs unique pathways
- [ ] Batch visualization
- [ ] Meta-analysis helpers

---

### Inter-Potato Networks **[LOW PRIORITY]**

**Goal:** Link potatoes via compound flow (serine uses ED outputs)

**Implementation:**

- [ ] Build inter-potato compound graph
- [ ] Evaluate higher-order metabolic capabilities
- [ ] "Genome has glycolysis + TCA + amino acids = autotroph"

---

## Priority Levels

- **HIGH** = Data integrity, critical bugs, blocking workflows
- **MEDIUM** = Valuable features, improves accuracy/usability  
- **LOW** = Nice to have, future enhancements

---

**Version history:** See git log and CLAUDE.md  
**Implementation details:** See CLAUDE.md for agent context
