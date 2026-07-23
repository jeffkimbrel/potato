# POTATO v1 - Design Brainstorm

## The Journey from GATOR v1 to POTATO v1

This document captures the design evolution from the original GATOR tool to the new POTATO architecture. It preserves the key insights, dead ends, and breakthroughs from our brainstorming session.

---

## Starting Point: GATOR v1

### What Worked Well

**Core Concept:**
- Excel-driven database with 3 sheets:
  - **db**: Tool definitions (kofam, BLAST, HMM, PFAM, PATRIC)
  - **gene**: Gene registry mapping gene names → detection methods
  - **pathway**: Pathway definitions using string syntax

**Key Innovation:**
- Annotation redundancy: a gene is "found" if detected by ANY listed tool
- Handles database disagreement (especially for non-E. coli genomes)
- Domain experts can define novel pathways without rigid KEGG constraints

**String Syntax:**
```
ilvB -> ilvH | ilvM -> ilvC -> ilvD -> ilvE
```
- `->` : sequential steps
- `|` : alternatives (OR)
- `+` : required co-occurrence (AND)

### What Didn't Work

**The Excel Constraint:**
- Pathway definitions crammed into single Excel cell
- Becomes unreadable with nested logic
- Ambiguous operator precedence: `ilvB + (ilvH | ilvM) -> ilvC`
- Can't handle convergent paths (branches that merge downstream)
- Hard to visualize complex pathways

**Missing Features:**
- No confidence scoring (binary pass/fail)
- Can't handle MAG incompleteness probabilistically
- Can't suggest functional analogs
- No gene specificity weighting (central metabolism genes vs. pathway-specific)

**Example Problem:**
Reverse TCA cycle has 2 alternative sub-pathways (rTCA1 vs. rTCA2) using different enzymes for one step, then converging. String syntax can't express this clearly.

---

## Design Exploration: V2 Attempt (2026-07-22)

### Initial V2 Proposal

From `/Users/kimbrel1/JAK_obsidian/Syntheses/20260722 - gator_v2_design_summary.md`:

**Key Ideas:**
1. **LLM assistance** for two problems:
   - MAG incompleteness tolerance
   - Functional analog detection
   
2. **Bipartite DAG structure:**
   - Alternating compound nodes (metabolites) ↔ enzyme nodes (genes)
   - Structurally mirrors biochemistry
   - Handles OR branches and convergence naturally
   
3. **JSON per pathway** instead of Excel sheet 3
   - Enzyme nodes reference gene sheet by ID
   - Visual authoring tool for experts

**Example Structure:**
```json
{
  "nodes": [
    {"id": "c1", "type": "compound", "name": "cyclic-microcystin"},
    {"id": "e1", "type": "enzyme", "gene_sheet_id": 353},
    {"id": "c2", "type": "compound", "name": "linearized-microcystin"}
  ],
  "edges": [
    {"from": "c1", "to": "e1"},
    {"from": "e1", "to": "c2"}
  ]
}
```

### Why V2 Stalled

**Realization:** The bipartite DAG wasn't rejected - the **Excel constraint** was!

With JSON files, we get:
- ✅ Full structural representation (not flattened to string)
- ✅ Visually graphable
- ✅ No cell size limits
- ✅ Can handle KEGG module complexity

**Key insight:** "The DAG is intuitive - you have different paths from A to B, and you can quantify which steps are missing. A string definition gets hard to parse, but a DAG can always be graphed so a human can see how the routes work."

---

## Current Brainstorm Session (2026-07-22)

### Fresh Start Questions

**Q: What do you want in potato v1?**

1. **First-class tool support** - CLI access to kofam, BLAST, HMMER on their terms
2. **Agentic skills** - LLM help building databases and analyzing results
3. **JSON folder structure** - many files, not one Excel
4. **Gene + pathway outputs** - both contextual views
5. **LLM-assisted completion** - handle MAG incompleteness

### Key Design Decisions

#### 1. Self-Contained Potatoes

**Original thought:** Gene definitions separate (like Excel gene sheet)

**New insight:** "I like the idea of self-contained JSON where genes are defined alongside the pathway"

**Benefits:**
- Each potato is portable/shareable
- No external dependencies
- Can have different definitions for same gene name in different contexts
- Easy to version control

**Challenge:** What if gene definitions conflict?

**Resolution:** 
- Validation can warn about inconsistencies
- Optional canonical gene registry for common genes
- Flexibility > strict enforcement

#### 2. Potato Categories/Tags

**Need:** Users want to run subsets of potatoes
- All potatoes (comprehensive scan)
- Specific pathway (`leucine_biosynthesis.json`)
- Category (`tags: ["taxonomy"]` for 16S, rps3, etc.)
- Category (`tags: ["photoheterotrophy"]`)

**LLM bonus:** Can find potatoes matching intent even without explicit tags

#### 3. Tool Configuration

**Problem:** Users have different file locations
- Jeff's kofam HMMs: `/Users/kimbrel1/databases/kofam/`
- Another user's: `/data/kofam/`

**Solution:** Separate `tools.json` config file
```json
{
  "kofam": {
    "enabled": true,
    "ko_profiles_dir": "/path/to/kofam/profiles/"
  }
}
```

**Graceful degradation:** If potato specifies KEGG but user has it disabled:
- ✅ Continue with available tools
- ⚠️ Warn about missing tools
- 📊 Report `tools_requested` vs. `tools_used` in output
- Flag genes only searchable via unavailable tools

#### 4. DAG Structure (The Breakthrough)

**Why DAG works now:**

1. **JSON removes constraints** - no Excel cell size limits
2. **Visually graphable** - always reviewable by humans
3. **Handles KEGG modules** - complex nested logic
4. **Extensible** - LLM can suggest new paths
5. **Intuitive** - matches how biochemists think

**Example:** Bifunctional enzyme
- Traditional: `A -> B -> C`
- With bifunctional shortcut: Add edge `A -> C` directly
- Both paths coexist in DAG

**LLM enhancement:** "I know of a bifunctional enzyme that can go A → C directly, want me to add it to the potato?"

#### 5. Compounds in DAG

**Initial question:** Do we need compound nodes (metabolites)?

**Biological reality:**
- Enzymes have substrates and products
- Different potatoes can connect via shared intermediates
- Metabolite context helps LLM reasoning

**Example use case:** Missing enzyme strips OH group from compound X
- LLM knows substrate/product
- Can search for EC 1.1.1.* (oxidoreductases on CH-OH)
- Suggests possible functional analogs
- BUT: must be cautious (substrate specificity matters)

**Resolution:** Compounds as **edge metadata**, not structural nodes

```json
{
  "edges": [
    {
      "from": "ilvB",
      "to": "ilvC",
      "compound": "2-acetolactate",
      "kegg_compound": "C06010"
    }
  ]
}
```

**Benefits:**
- Scoring logic stays simple (enzyme-centric)
- Preserves biological context
- Enables cross-potato queries
- Optional (can omit if user doesn't care)

#### 6. Gene Specificity Weighting

**Problem:** Not all genes are equally informative

**Example:**
- **Low specificity:** `gapA` (glyceraldehyde-3-P dehydrogenase)
  - Found in 50+ pathways
  - Presence tells you almost nothing about THIS pathway
  
- **High specificity:** `mlrA` (microcystin linearase)  
  - Only in microcystin degradation pathway
  - Presence is strong evidence pathway exists

**Biological intuition:**
"If I have 5/8 genes found, but those 5 are only involved in THIS pathway, that's strong evidence the pathway exists even with 3 missing. But if those 5 are central metabolism genes in 10 other pathways, it tells me nothing."

**Solution:** Pre-computed specificity scores

```
specificity = 1 / (number of potatoes containing this gene)
```

**Weighted scoring:**
```
confidence = sum(specificity of found genes) / sum(specificity of all genes)
```

**Example:**
- Found 5/8 genes, all low-specificity → 35% confidence
- Found 5/8 genes, including 3 high-specificity unique ones → 85% confidence

**LLM interpretation:**
"You found leuB, leuC, leuD (unique to leucine biosynthesis) but missing ilvE (also in valine/isoleucine). High confidence pathway exists - ilvE likely missing due to 78% MAG completeness."

**Implementation:** Pre-computed using LLM reasoning, baked into potato database (hybrid approach)

#### 7. Agentic Features

**Two distinct agent roles:**

**A. Potato Builder Agent** (database curation)

*Mode 1: From KEGG module*
```r
build_potato(kegg_module = "M00012")
```
- Fetch from KEGG REST API
- Parse KEGG module syntax (complex nested parentheses)
- Extract KOs, pathway structure, compounds
- Generate potato JSON draft
- User reviews/edits

**KEGG module example:**
```
(K00030,K00031) (K00164+K00658+K00382,K01902+K01903)
```

*Mode 2: From gene list*
```r
build_potato(genes = c("rbdA", "rbdB", "rbdC", "rbdD"))
```
- LLM searches for KOs, PFAMs, ECs
- Suggests pathway logic
- Generates draft

*Mode 3: Interactive*
```r
build_potato(prompt = "I want a shikimate pathway potato")
```
- Conversational refinement

**B. Genome Analysis Agent** (result interpretation)

```r
interpret_results(genome_results, completeness = 0.78)
```
- Considers gene specificity
- Considers MAG completeness
- Generates explanations
- Suggests functional analogs (Phase 5)

**C. Converter Agent** (migration)

```r
convert_gator_db("gator_db.xlsx")
```
- Reads old Excel format
- Parses v1 pathway syntax
- Generates potato JSON files
- Validation report

#### 8. Functional Analog Detection (Advanced)

**The challenge:**
- MAGs are 50-100% complete (real genes may be missing)
- Some enzymes have promiscuity (can work on related substrates)
- EC classification helps (EC 1.1.1.- = oxidoreductases on CH-OH)
- BUT substrate specificity matters (can't just swap any dehydrogenase)

**Proposed approach:**

Potato genes include reaction metadata:
```json
{
  "id": "ilvC",
  "reaction": {
    "type": "oxidoreductase",
    "substrate": "(S)-2-acetolactate", 
    "product": "(R)-2,3-dihydroxy-3-methylbutanoate",
    "chemistry": "reduction + isomerization"
  }
}
```

**LLM reasoning:**
- See what chemistry is needed
- Know substrate specificity
- Search for EC 1.1.1.* alternatives
- Flag as "possible alternative - unvalidated"

**Output confidence levels:**
- ✅ Canonical gene found (high confidence)
- ⚠️ Possible functional analog (unvalidated, expert review needed)
- ❌ Missing, MAG 85% complete (may be in missing 15%)
- ❌ Missing, MAG 98% complete (likely truly absent)

**Caution:** This is Phase 5 (advanced). Core tool must work without it.

---

## Architecture Evolution

### Initial Concept (Pre-Brainstorm)
```
R package
├── inst/extdata/potato_db.xlsx  (v1 Excel format)
└── inst/python/                  (copied from gator v1)
```

### Final Architecture (Post-Brainstorm)
```
R package (primary interface)
├── Python backend (via reticulate)
│   ├── Tool orchestration
│   ├── DAG scoring
│   └── LLM integration
├── JSON potato files (self-contained pathways)
├── tools.json (user tool configuration)
└── Optional canonical_genes.json (consistency)
```

---

## Key Insights & Breakthroughs

### Insight #1: The DAG Was Good, Excel Was Bad

**Mistake:** Thinking the bipartite DAG was too complex

**Reality:** The DAG structure was sound - the Excel constraint made it unusable

**Breakthrough:** JSON files remove all constraints → DAG becomes viable and elegant

### Insight #2: Self-Contained is Better Than Centralized

**Old thinking:** Central gene registry (like Excel gene sheet)

**New thinking:** Each potato defines its own genes

**Benefit:** Portability, flexibility, version control, independent evolution

### Insight #3: Not All Genes Are Equal

**Biological intuition:** Pathway-specific genes are more informative than ubiquitous ones

**Breakthrough:** Specificity weighting transforms scoring from simplistic (5/8 genes) to contextual (high-confidence vs. low-confidence finds)

**Impact:** Makes probabilistic interpretation possible

### Insight #4: Compounds Inform But Don't Constrain

**Old v2:** Bipartite graph with compound nodes structurally required

**New v1:** Compounds as optional edge metadata

**Benefit:** Simple scoring logic + biological context when needed

### Insight #5: Tool Availability is Variable

**Reality:** Not all users have all tools (KEGG is expensive/restricted)

**Solution:** Graceful degradation with transparent reporting

**Impact:** Potatoes are shareable across users with different setups

### Insight #6: LLMs Are Assistants, Not Oracles

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

## Dead Ends & Lessons

### Dead End #1: String Syntax Extensions

**Tried:** Make string syntax more expressive (nested parentheses, etc.)

**Problem:** Always hits ambiguity wall

**Lesson:** Switch representations (DAG) rather than patching syntax

### Dead End #2: Central Gene Registry Required

**Tried:** Enforce canonical definitions

**Problem:** Too rigid, conflicts with portability

**Lesson:** Make it optional with validation warnings

### Dead End #3: LLM at Runtime for Every Genome

**Tried:** LLM analyzes each genome's results live

**Problem:** Too slow, expensive, variable

**Lesson:** Pre-compute where possible, LLM for interpretation layer only

---

## Open Questions (Still Unresolved)

1. **Specificity computation:** How much LLM reasoning vs. simple counting?

2. **Canonical gene registry:** Required, optional, or skip entirely?

3. **Reaction metadata:** Always required, optional, or generate on-demand?

4. **Validation strictness:** Default to strict enforcement or warnings?

5. **Potato versioning:** How to track updates to potato definitions over time?

6. **Output formats:** TSV sufficient, or also JSON/HDF5/Parquet for big data?

7. **LLM provider:** Claude only, or multi-provider (GPT-4, local models)?

8. **Phase 5 timing:** When is functional analog detection actually ready for production?

---

## What Makes POTATO Different

### vs. GATOR v1

| Aspect | GATOR v1 | POTATO v1 |
|--------|----------|-----------|
| **Pathway format** | String in Excel cell | DAG in JSON |
| **Gene definitions** | Central registry | Self-contained |
| **Scoring** | Binary pass/fail | Confidence-weighted |
| **Tool handling** | All-or-nothing | Graceful degradation |
| **Visualization** | N/A | Graphable DAGs |
| **LLM features** | None | Builder + analysis agents |
| **Specificity** | Not considered | Pre-computed weighting |

### vs. KEGG Modules

| Aspect | KEGG Modules | POTATO |
|--------|--------------|--------|
| **Pathway coverage** | 500+ canonical | User-defined + canonical |
| **Custom pathways** | Not supported | First-class |
| **Syntax** | KEGG-specific | JSON DAG (human + machine) |
| **Visualization** | KEGG pathway maps | igraph, graphviz |
| **Tool flexibility** | KEGG only | Any annotation tool |
| **Confidence** | Binary | Weighted + MAG-aware |

### vs. Other Tools

**vs. AntiSMASH / BiG-SCAPE:**
- Focus: General metabolism vs. specialized (BGCs)
- Scope: Any pathway vs. secondary metabolites only

**vs. MetaCyc / Pathway Tools:**
- Curation: User-driven vs. centralized database
- Flexibility: Custom pathways easy vs. rigid ontology

**vs. DRAM:**
- Approach: Pathway-centric vs. gene-centric annotation summary
- Scoring: Confidence-weighted DAG vs. presence/absence matrix

---

## The Name

**GATOR:** Genome annotATOR

**POTATO:** Pathway annOTATOr

**Evolution:** Gene-centric → Pathway-centric

**Why it matters:** The focus shifted from annotating individual genes to interpreting metabolic capabilities through pathway context. The name change reflects this conceptual shift.

---

## What's Next

1. **Formalize schemas** (JSON Schema for potato, tools.json)
2. **Build Phase 1 MVP** (single potato, single genome, no LLM)
3. **Create test potatoes** (manually author 3-5 for validation)
4. **Implement converter agent** (migrate gator_db.xlsx)
5. **Add specificity weighting** (Phase 3)
6. **Integrate LLM agents** (Phase 4)
7. **Explore functional analogs** (Phase 5)

---

## Success Metrics

**Technical success:**
- Can run 20+ potatoes on 100+ MAGs in reasonable time
- Confidence scores correlate with expert validation
- Converter agent successfully migrates gator_db.xlsx

**User success:**
- Domain experts can build custom potatoes easily
- Results are more interpretable than GATOR v1
- Reduces false positives/negatives vs. binary scoring

**Community success:**
- Shareable potato collections emerge (diazotroph specialists, etc.)
- Published database of validated potatoes
- Adopted beyond original lab

---

## Reflection

This brainstorm revealed that:

1. **The v2 design was 80% right** - we just needed to remove the Excel constraint
2. **Self-contained potatoes** solve more problems than we initially realized
3. **Specificity weighting** is the killer feature for MAG analysis
4. **LLMs are best as assistants** - pre-compute where possible, interpret at the end
5. **The DAG is the right structure** - it just needed breathing room (JSON)

The hardest part was recognizing that the Excel limitation was driving architectural decisions, not biological requirements. Once we removed that constraint, the design opened up.

**Key lesson:** Sometimes the right architecture is blocked by the wrong tooling. Switch tools, unlock design.

---

## Contributors

- Jeff Kimbrel (domain expertise, biological intuition, design direction)
- Claude (design exploration, architectural synthesis, brainstorm facilitation)

Session date: 2026-07-22

Last updated: 2026-07-22
