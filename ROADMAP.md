# POTATO - Development Roadmap

**Current Version:** v0.10.2-dev (2026-08-14)

**Status:** V2 schema complete. Full annotation pipeline
(kofam/BLAST/HMM) → scoring (all/essential genes) → visualization
(interactive/static) → provenance tracking

**Up Next:**

1.  Gap-based scoring implementation in
    [`score_pathways()`](https://jeffkimbrel.github.io/potato/reference/score_pathways.md) -
    DAG traversal logic
2.  Potato splitting and merging functions - split multi-pathway
    networks into separate files, merge separate potatoes into networks
3.  Potato locking - prevent modification when results exist
4.  Incremental annotation - merge results from multiple runs

------------------------------------------------------------------------

## Vision

POTATO (Pathway annOTATOr) annotates genome collections (MAGs) against
curated metabolic pathways. Core innovation: **self-contained potato
files** (pathway definitions as DAG structures in JSON) with
**LLM-assisted analysis**.

### Design Principles

1.  **Self-contained potatoes** - each pathway is independent, portable
2.  **DAG-based pathway logic** - intuitive, graphable, handles complex
    branching
3.  **Tool-agnostic detection** - works with whatever annotation tools
    available
4.  **Confidence-aware scoring** - interprets results in context of gene
    specificity and MAG completeness
5.  **LLM augmentation** - assists in database building and result
    interpretation

------------------------------------------------------------------------

## V2 Schema Reference

**Schema Status:** V2 is frozen (stable). Package versions (v0.11.x,
v0.12.x) track implementation milestones, not schema changes. Future
breaking changes would require V3.

**Complete schema documentation:** See
[inst/POTATO_SCHEMA_V2.md](https://jeffkimbrel.github.io/potato/inst/POTATO_SCHEMA_V2.md)

**Key features:** - `genes`: Global definitions with detection methods
OR complex definitions - `type`: “enzyme”, “transport”, “chaperone”, or
“complex” - Complexes: use `components` array (not `databases`), all
components must exist as genes - `pathways`: One or more pathways with
topology - `edges`: Connect genes/compounds (or complexes), carry
`required`/`marker` - `verified`: Always `false` for new pathways (only
humans verify) - Pathway types: `variant` (alternatives) or
`independent` (different functions) - `input`/`output`: Mandatory for
gap-based scoring, optional for fraction-based - Marker genes:
Pathway-level `markers` field (v0.14.4+) OR edge-level `marker`
attributes - Scoring methods: Fraction-based (`min_fraction`) and/or
gap-based (`max_gaps`)

------------------------------------------------------------------------

## Data Integrity & Workflow

### Potato Splitting and Merging **\[HIGH PRIORITY\]**

**Problem:** Users need to reorganize potatoes as they learn about
pathway relationships: - Multi-pathway networks may need to be split
into separate files for independent analysis - Separate pathway potatoes
may need to be merged into networks when relationships discovered -
Current workflow: manual JSON editing (error-prone, tedious)

**Solution:** Add utility functions for potato reorganization

**Functions needed:**

**1. `write_potato()` - Save potato object to JSON file**

``` r

write_potato(
  potato,                                    # Potato S7 object (from sack@potatoes or load_potato())
  path = "inst/potatoes/my_potato.json",    # Output path
  indent = 2,                                # JSON indentation (default: 2 spaces)
  validate = TRUE                            # Validate before writing (default: TRUE)
)
```

Writes potato S7 object to JSON file with: - Proper field ordering (id,
name, source, schema_version, tags, notes, genes, compounds, pathways) -
Pretty-printed JSON with configurable indentation - Optional validation
before writing (catches structural errors) - Preserves all metadata,
coordinates, etc.

**Use cases:** - User loaded potato from sack, made edits in R, wants to
save: `write_potato(sack@potatoes$ed_network, "updated.json")` - User
built potato programmatically in R, needs to export to file - Batch
operations: load all potatoes, modify, write back out - Symmetric with
[`load_potato()`](https://jeffkimbrel.github.io/potato/reference/load_potato.md) -
what you load, you should be able to write

**2. `split_potato()` - Extract pathways from multi-pathway network**

``` r

split_potato(
  potato_path = "inst/potatoes/ed_network.json",
  output_dir = "inst/potatoes/",
  pathways = c("classic", "non_phosphorylative"),  # Which pathways to extract
  name_pattern = "{potato_id}_{pathway_id}"        # Output filename pattern
)
```

Creates separate potato files, each with: - Only genes used in that
pathway - Single pathway in `pathways` field - Original metadata
(source, notes) preserved - New `notes` field: “Extracted from
{original_potato} on {date}”

**2. `merge_potatoes()` - Combine separate potatoes into network**

``` r

merge_potatoes(
  potato_paths = c("inst/potatoes/ed_classic.json", 
                   "inst/potatoes/ed_np.json"),
  output_path = "inst/potatoes/ed_network.json",
  network_id = "ed_network",
  network_name = "Entner-Doudoroff Network",
  pathway_type = "variant"  # or "independent"
)
```

Creates multi-pathway network with: - Deduplicated genes (same gene ID →
merged into one definition) - All pathways in `pathways` field -
Validates no conflicting gene definitions (same ID, different
databases) - New `notes` field: “Merged from {potato1}, {potato2} on
{date}” - Uses `write_potato()` internally to save output

**Implementation details:** - All three functions should use common
internal helpers for JSON serialization - `write_potato()` is the
foundation - `split_potato()` and `merge_potatoes()` call it - Preserve
`verified` status per pathway - Handle coordinate systems (x/y fields) -
keep if consistent, drop if conflicting - Validate schema compatibility
before merging - Check for duplicate pathway IDs - Warn about genes with
different detection methods (potential conflicts)

**Use cases:** - User built ED network, wants to analyze classic variant
separately → split - User has Mo-nif and V-nif separate, realizes
they’re variants → merge - Testing: split network, run analysis, compare
to network results - Distribution: ship separate potatoes for simple
cases, networks for advanced users

------------------------------------------------------------------------

### Potato Extension & Overlay System **\[MEDIUM PRIORITY\]**

**Problem:** Users need to customize bundled potatoes without modifying
original files: - Adjust scoring thresholds for their specific use
case - Add custom pathway variants to existing networks - Maintain
separation between bundled (package) and custom (user) potatoes - Update
package without losing customizations

**Solution:** Extension/overlay system with `extends` field

**Use cases:**

**1. Load from multiple directories:**

``` r

sack <- create_sack()
sack <- add_potatoes(sack, "inst/potatoes")      # Bundled with package
sack <- add_potatoes(sack, "~/custom_potatoes")  # User custom potatoes
```

**2. Extend existing potato** (add pathways to bundled network):

``` json
{
  "schema_version": "v2",
  "extends": "nitrogen_fixation",
  "pathways": {
    "my_custom_nif": {
      "name": "My Custom Nif Variant",
      "type": "variant",
      "notes": "Custom variant for my specific organism group",
      "edges": [...],
      "scoring": {...}
    }
  }
}
```

Result: Loaded potato has all original pathways PLUS `my_custom_nif`

**3. Override settings** (adjust thresholds without copying whole file):

``` json
{
  "schema_version": "v2",
  "extends": "nitrogen_fixation",
  "pathways": {
    "mo_nif": {
      "scoring": {
        "min_fraction": 0.9  // Override bundled 0.75
      }
    }
  }
}
```

Result: `mo_nif` pathway uses 0.9 threshold, all other pathways
unchanged

**4. Add custom genes to existing network:**

``` json
{
  "schema_version": "v2",
  "extends": "urea_metabolism",
  "genes": [
    {
      "id": "my_custom_urease",
      "name": "novel urease variant",
      "databases": {"blast": ["my_ref_seq"]},
      "type": "enzyme"
    }
  ],
  "pathways": {
    "urease": {
      "genes": ["ureA", "ureB", "ureC", "my_custom_urease"]
    }
  }
}
```

**Implementation details:**

Add `extends` field to schema (optional, references potato `id`)

Merge logic in
[`load_potato()`](https://jeffkimbrel.github.io/potato/reference/load_potato.md):

- Load base potato first
- Deep merge: genes/compounds/pathways are additive
- Pathway-level overrides: if pathway ID exists, merge/override fields
- Gene-level overrides: if gene ID exists, merge databases/notes

Config option for potato search paths (package, user, project-local)

Load order / priority: project-local \> user \> package

Validation: check extended potato exists and is compatible

[`print()`](https://rdrr.io/r/base/print.html) method shows extension
chain: “nitrogen_fixation (extended by custom_nif.json)”

**Design decisions:** - Extensions don’t modify original files
(non-destructive) - Multiple extension levels allowed: A extends B
extends C - Circular extends detection (error if detected) - Extensions
are **additive** - can’t delete pathways/genes from base - Field-level
overrides for scoring, notes, edges - Version compatibility: extension
must match base `schema_version`

**Benefits:** - Package updates don’t overwrite user customizations -
Share custom variants without duplicating entire potato - Test parameter
sensitivity (create extensions with different thresholds) -
Domain-specific tuning (marine microbes, gut microbes, etc.)

------------------------------------------------------------------------

## Data Integrity & Workflow

### Gap-Based Scoring System **\[HIGH PRIORITY\]**

**Problem:** Fraction-based scoring has fundamental limitations: - Float
precision issues (0.67 threshold fails 2/3 genes at 0.6666…) - Doesn’t
distinguish complete vs. incomplete routes in multi-pathway networks -
Ignores DAG connectivity (80% with middle gap fails, 60% with full path
succeeds)

**Solution:** Add gap-based scoring as alternative to fraction-based

**Phase 1: Schema Design ✅ COMPLETE (v0.11.0-dev)** - ✅ Schema
designed: `max_gaps` field in scoring block - ✅ No `method` field
needed - presence of `max_gaps` + `input`/`output` triggers gap
scoring - ✅ Both `min_fraction` and `max_gaps` can coexist - pathway
passes if EITHER succeeds - ✅ Simplified input/output format: arrays of
compound IDs (`["C00048"]`) instead of verbose objects - ✅ Validation
updated: converts JSON arrays, validates compound IDs exist, checks
connectivity - ✅ Visualization updated:
[`view_pathway_detail()`](https://jeffkimbrel.github.io/potato/reference/view_pathway_detail.md)
displays input/output and max_gaps - ✅ Implemented in example potatoes:
`glycerate_pathway_v2.json`, `glyoxylate_shunt_v2.json` - ✅ Style guide
documented in CLAUDE.md “Gap-Based Scoring” section - ✅ Schema
reference updated with examples

**Phase 2: Scoring Implementation** **\[HIGH PRIORITY - NEXT TASK\]** -
\[ \] Update
[`score_pathways()`](https://jeffkimbrel.github.io/potato/reference/score_pathways.md)
to detect gap-based scoring parameters - \[ \] Implement DAG traversal
from `input` → `output` compounds - \[ \] Find best path (fewest missing
genes) considering OR branches - \[ \] Calculate gap count for each
pathway - \[ \] Update results tibble to include: `gap_count`,
`present_gaps` (boolean) - \[ \] When both methods present, calculate
both and report separately - \[ \] Maintain backward compatibility with
fraction-only potatoes

**Example potato scoring block (implemented):**

``` json
"input": ["C00048"],
"output": ["C00631"],
"scoring": {
  "min_fraction": 1.0,  // fraction-based threshold
  "max_gaps": 0,        // gap-based threshold (both can be present)
  "marker_mode": "any"
}
```

**Benefits:** - Biologically meaningful (“allow 1 missing gene”) -
Handles alternative routes correctly - No float precision issues - True
DAG traversal (checks connectivity) - Dual scoring provides flexibility

**Implementation notes:** - Ignore `required`/`marker` attributes for
gap counting (purely topological) - Only traverse edges between `input`
and `output` (ignore upstream/downstream context) - For cyclic pathways,
gap-based scoring requires artificial start/end points via
`input`/`output`

**Phase 3: Scoring Coverage Utility** **\[LOW PRIORITY\]** - \[ \] Add
function to summarize scoring type coverage across potatoes - \[ \]
Report: which pathways use fraction-only, gaps-only, dual, or neither -
\[ \] Identify potatoes missing `input`/`output` (candidates for
gap-based scoring) - \[ \] Could be part of
[`validate_potato()`](https://jeffkimbrel.github.io/potato/reference/validate_potato.md)
output or standalone `summarize_scoring_types()` - \[ \] Useful for
auditing potato collection and tracking migration to gap-based scoring

------------------------------------------------------------------------

### Potato Filtering at Sack Creation **\[MEDIUM PRIORITY\]**

**Problem:** Users may want to create multiple sacks from the same
project with different potato subsets. Currently
[`create_sack()`](https://jeffkimbrel.github.io/potato/reference/create_sack.md)
loads ALL potatoes from the `potatoes/` directory with no filtering
options.

**Use cases:** - Focused analysis: nitrogen metabolism only vs carbon
metabolism only - Performance: load subset for quick testing, all
potatoes for full run - Multiple sacks: `sack_nitrogen`, `sack_carbon`
from one project directory

**Solution:** Add filtering parameters to
[`create_sack()`](https://jeffkimbrel.github.io/potato/reference/create_sack.md)

**Proposed signature:**

``` r

create_sack(
  path = NULL,
  tags = NULL,              # Filter potatoes by tags (e.g., c("nitrogen", "metabolism"))
  pattern = NULL,           # Regex pattern for filenames (e.g., "nitrogen.*\\.json")
  include_inactive = FALSE  # Load potatoes marked active: false
)
```

**Implementation:** - Pass filtering parameters through to
[`load_potatoes()`](https://jeffkimbrel.github.io/potato/reference/load_potatoes.md)
call (line 53 in R/create_sack.R) -
[`load_potatoes()`](https://jeffkimbrel.github.io/potato/reference/load_potatoes.md)
already supports `tags` and `include_inactive` parameters - Add
`pattern` parameter to
[`load_potatoes()`](https://jeffkimbrel.github.io/potato/reference/load_potatoes.md)
for filename filtering

**Example workflows:**

``` r

# Load all nitrogen-related pathways
sack_n <- create_sack(tags = "nitrogen")

# Load specific pathways by filename pattern
sack_denit <- create_sack(pattern = "denitrification|nitrate")

# Load everything including inactive/deprecated potatoes
sack_full <- create_sack(include_inactive = TRUE)

# Combine filters
sack_custom <- create_sack(tags = c("metabolism", "energy"), pattern = "^amino")
```

**Benefits:** - Multiple specialized sacks from one project - Faster
iteration during development (load subset) - Organized analysis
workflows (group related pathways)

------------------------------------------------------------------------

### Potato Locking System **\[HIGH PRIORITY\]**

**Problem:** User edits potato JSON after annotation → hash mismatch →
results meaningless

**Solution:**

Detect hash mismatches on sack load

Add `potatoes_locked` flag to PotatoSack

Check lock status in modification functions

Clear error: revert potato or clear results

**Why critical:** Prevents data corruption, results always match
definitions

------------------------------------------------------------------------

### Potato Version Upgrades **\[LOW PRIORITY\]**

**Problem:** Bundled potatoes in package may get bug fixes or
improvements, but projects use static snapshots from initialization. No
way to “upgrade” a potato mid-analysis without losing existing results.

**Current behavior (keep for now):** - Projects snapshot potatoes at
initialization → reproducible, but no updates - Editing potato JSON
breaks hash → results invalidated - **Friction is good:** forces users
to start over if changing pathway definitions

**Future feature - Versioned potato upgrades:**

Allow multiple versions of the same potato to coexist in a sack,
distinguished by hash suffix:

    glycolate_cycle_a3f4b2c  (v1 - original definition)
    glycolate_cycle_d8e1f9a  (v2 - bug fix, revised thresholds)

**Implementation:**

`upgrade_potato(sack, potato_id)` - loads latest package version
alongside existing

New potato gets hash-suffixed ID: `{original_id}_{hash_short}`

Both potatoes retained with their annotation results

User can compare: did bug fix change pathway calls?

Provenance tracks upgrade event: “glycolate_cycle upgraded from a3f4b2c
→ d8e1f9a”

**Use case:**

``` r

sack <- run_kofam(sack)  # Uses glycolate_cycle (hash a3f4b2c)
# Package releases bug fix for glycolate_cycle
sack <- upgrade_potato(sack, "glycolate_cycle")  # Adds glycolate_cycle_d8e1f9a
sack <- run_kofam(sack)  # Runs on new version only
# Now have results for both versions - can compare
```

**Benefits:** - Compare old vs new definitions side-by-side - Preserve
historical results (reproducibility) - Gradual migration (test new
version before removing old) - Clear audit trail via hash tracking

**Related:** Potato locking system (above) prevents accidental edits.
This feature is for **intentional** upgrades.

------------------------------------------------------------------------

### Update Claude Skills in Projects **\[LOW PRIORITY\]**

**Problem:** Claude agent skill files (`.claude/commands/*.md`) may
improve in package, but projects have static copies from initialization.

**Solution:** Add `update_skills()` function

``` r

update_skills()  # Updates .claude/ directory from package to current project

# Or specify path
update_skills("path/to/project")
```

**Implementation:**

Copy `.claude/` directory from package `inst/.claude/` to project

Warn if local modifications detected (file hash mismatch)

Option: `force = TRUE` to overwrite local changes

Report which files updated, which skipped

Safe operation: agents have no associated data (unlike potatoes)

**Benefits:** - Get bug fixes and improvements to agents - Low risk:
agents don’t have associated results - Users can customize agents
locally and decline updates

**Note:** Unlike potatoes (which need version tracking for
reproducibility), agent updates are safe because they don’t affect
annotation results.

------------------------------------------------------------------------

### Incremental Annotation **\[HIGH PRIORITY\]**

**Problem:** Running annotation twice errors unless `overwrite = TRUE`
(replaces ALL results). Can’t add potatoes incrementally.

**Solution:**

Detect conflicts (genome + potato_hash)

Merge logic: skip conflicts or replace based on `overwrite`

Update provenance for merged runs

Apply to kofam, BLAST, HMM

**Use case:**

``` r

sack <- run_kofam(sack, potato_names = "pathway_A")
sack <- run_kofam(sack, potato_names = "pathway_B")  # Should merge
```

------------------------------------------------------------------------

### Pathway-Level Marker Genes **\[MEDIUM PRIORITY\]**

**Problem:** Current schema ties `marker` attribute to edges, but
biologically, marker genes are properties of the pathway independent of
topology. When genes are in complexes, you can’t mark individual
components.

**Current schema:**

``` json
"edges": [
  {"from": "substrate", "to": "geneA", "marker": true}
]
```

**Issue:** - If geneA is part of complex geneABC, edges reference the
complex, not individual genes - Can’t mark geneA as diagnostic when it’s
buried in a complex - Edge-level markers don’t make sense for
non-catalytic genes (chaperones, regulators)

**Proposed schema enhancement:**

Add pathway-level `markers` field:

``` json
"pathways": {
  "pathway_id": {
    "markers": ["nifH", "vnfH"],  // Gene IDs that are diagnostic
    "edges": [...],
    "scoring": {
      "marker_mode": "any"  // Already exists, just not implemented
    }
  }
}
```

**Implementation:**

Add `markers` field to V2 schema (optional, array of gene IDs)

Update
[`score_pathways()`](https://jeffkimbrel.github.io/potato/reference/score_pathways.md)
to use pathway-level markers for `marker_mode` scoring

`marker_mode: "any"` = pathway present if ANY marker detected

`marker_mode: "all"` = pathway present only if ALL markers detected

Validate: marker gene IDs must exist in pathway’s `genes` array

Visualization:
[`view_pathway_detail()`](https://jeffkimbrel.github.io/potato/reference/view_pathway_detail.md)
shows checkmark in marker column for listed genes

Maintain backward compatibility: if no pathway-level `markers` field,
fall back to edge-level `marker` attributes

**Benefits:** - Clearer biological meaning (these genes ARE markers,
period) - Works with complexes (mark components directly) - Enables
marking non-catalytic genes (electron carriers, chaperones) - Decouples
diagnostic value from topology

**Example:** Nitrogen fixation - nifH and vnfH are diagnostic regardless
of whether they’re in complexes or sequential steps.

------------------------------------------------------------------------

### Threshold Sensitivity Analysis **\[MEDIUM PRIORITY\]**

**Goal:** Understand how threshold choices affect results

**Implementation:**

`analyze_thresholds(sack)` sweeps thresholds

Show which pathways flip present/absent

Identify threshold-sensitive vs robust pathways

Visualization: completion vs threshold

Export results as tibble

------------------------------------------------------------------------

## Scoring & Analysis

### Protein Complex Scoring **\[HIGH PRIORITY - IN PROGRESS\]**

**Goal:** Handle multi-component protein complexes where all components
are required for function

**Current status:** - ✅ Schema defined and frozen as part of V2 - ✅
Visualization complete (complex appears as single node) - ✅ Validation
complete (checks components exist as genes) - \[ \] Scoring logic not
yet implemented

**Schema structure:**

``` json
{
  "genes": [
    // Component genes with detection methods
    {"id": "urtA", "type": "enzyme", "databases": {"kofam": ["K11959"]}},
    {"id": "urtB", "type": "enzyme", "databases": {"kofam": ["K11960"]}},
    {"id": "urtC", "type": "enzyme", "databases": {"kofam": ["K11961"]}},
    {"id": "urtD", "type": "enzyme", "databases": {"kofam": ["K11962"]}},
    {"id": "urtE", "type": "enzyme", "databases": {"kofam": ["K11963"]}},
    
    // Complex definition
    {
      "id": "urtABCDE",
      "name": "urea ABC transporter complex",
      "type": "complex",
      "components": ["urtA", "urtB", "urtC", "urtD", "urtE"],
      "notes": "All components required for function"
    }
  ],
  "pathways": {
    "pathway_id": {
      "genes": ["urtA", "urtB", "urtC", "urtD", "urtE"],  // List components for scoring
      "edges": [
        {"from": "substrate", "to": "urtABCDE", "required": true, "marker": true}  // Complex in edges
      ]
    }
  }
}
```

**Key design decisions:** - Components are regular genes with detection
methods (kofam/blast/hmm) - Complex is additional entry with
`type: "complex"` and `components` array - Complex uses `components`
(not `databases`) - detected via component annotation - Multiple
complexes = alternatives (e.g., `ureABC` OR `ureAB_ureC` like GATOR `|`
operator) - Pathway’s `genes` array lists components (not complex IDs)
for scoring - Edges reference complex IDs for clean network
visualization

**Implementation needed:**

Scoring logic: detect `type: "complex"` entries in genes array

Check if ALL `components` have annotation hits passing thresholds

Mark complex as detected only if all components present

Handle complexes with alternative subunits (e.g., ureABC vs ureAB+ureC)

**Use cases:** - ABC transporters (urtABCDE - 5 components all
required) - Multi-subunit enzymes (ureABC, glcDEF - all subunits
required) - Alternative complex stoichiometries (ureABC vs ureAB+ureC)

**Validation:** - ✅ Checks all components exist as genes with their own
IDs - ✅ Warns if complex has `databases` field (should only have
`components`) - ✅ Warns if pathway uses complex in edges but doesn’t
list components in `genes` array

### Specificity Weighting **\[MEDIUM PRIORITY\]**

**Goal:** Weight genes by specificity - finding `nifH`
(pathway-specific) is more informative than `gapA` (ubiquitous)

**Implementation:**

Pre-compute specificity: `1 / (num potatoes containing gene)`

Weighted scoring: `sum(specificity of found) / sum(specificity of all)`

Output: raw score + weighted score + high-value genes

**Context-Aware Threshold Adjustment:**

Query KEGG API for pathway membership per KO

Identify pathway-specific marker genes (only in one pathway)

Relax thresholds for borderline hits when marker genes present

Example: `gcl` scores 450 (threshold 678), BUT `glxK` (specific marker)
present → likely real pathway

**Use case:** Genome has glxR (ubiquitous) + glxK (specific marker only
in glycerate pathway) + borderline gcl hit → context suggests gcl is
real, not a paralog

------------------------------------------------------------------------

### Required Edge Attributes and OR Branch Scoring **\[MEDIUM PRIORITY\]**

**Goal:** Clarify and implement `required` edge attribute logic for
alternative pathway routes

**Current situation:** - `required` attribute exists on edges but
scoring doesn’t fully utilize it - Unclear how to mark OR branches
(Mo-nif OR V-nif alternatives) - Example: nitrogen fixation has two
routes, all edges marked `required: false`

**Question to resolve:** For OR branches (either Mo-nif OR V-nif
needed): - Should all alternatives be `required: true` (meaning “at
least one needed”)? - Or should all be `required: false` (meaning
“optional” because alternatives exist)?

**Implementation needed:**

Define clear semantics for `required` in OR branches

Update schema documentation with examples

Implement scoring logic that respects `required` attributes

Handle cases where:

- All edges in route A are required IF taking route A
- All edges in route B are required IF taking route B
- At least ONE complete route must be present

**Related to:** Gap-based scoring, protein complex scoring (complexes as
alternatives)

------------------------------------------------------------------------

### Input/Output Validation **\[MEDIUM PRIORITY\]**

**Goal:** Pathways aren’t “complete” just by gene count - must transform
**available inputs → desired outputs**

**Current:** Input/output required for gap-based scoring, optional for
fraction-based  
**Future:** Validate reachability and biological completeness

**Implementation:**

Validate biological completeness for gap-scored pathways

Compound matching via KEGG IDs

Check input reachability (transport? synthesis? environmental?)

Warning: pathway with input/output defined but no path between them

Cross-pathway validation: can other pathways provide missing inputs?

**Example:**

    npED alone: 5/6 genes (83%) → but no gluconate source → NON-FUNCTIONAL
    ED network: transport + npED → gluconate available → FUNCTIONAL

**Note:** Input/output will remain **optional** for fraction-based
scoring (not all pathways need them).

------------------------------------------------------------------------

### Cross-Potato Metabolic Networks **\[MEDIUM PRIORITY\]**

**Goal:** Evaluate pathways in context of whole-genome metabolism

**Biological reality:** Bacterial metabolism is one cytoplasmic soup -
intermediates from one pathway immediately available for another

**Implementation:**

Build genome-wide metabolite pool (available, produced, consumed)

Compound matching via KEGG IDs

`score_metabolic_network(sack, genome)` evaluates in context

“ED produces glyceraldehyde → serine biosynthesis uses it → both
COMPLETE”

------------------------------------------------------------------------

### MAG Completeness Integration **\[LOW PRIORITY\]**

**Goal:** Adjust confidence based on genome completeness

**Implementation:**

Accept CheckM/BUSCO estimates

Adjust: “3 genes missing, MAG 60% complete”

Flag: likely present / likely absent / uncertain

------------------------------------------------------------------------

### Normalized Margin Calculation **\[LOW PRIORITY\]**

**Goal:** Make threshold margins comparable across different metric
types in
[`inspect_gene_thresholds()`](https://jeffkimbrel.github.io/potato/reference/inspect_gene_thresholds.md)

**Problem:** E-value margins (1e-11) vs bitscore margins (10) vs score
margins (50) are on completely different scales

**Solution:**

Add percent-of-threshold calculation: `(score/threshold) * 100` or
`(threshold/evalue) * 100`

Values \>100% = passing, \<100% = failing

Allows meaningful sorting across different metric types

Alternative: add margin_type column and sort within types only

------------------------------------------------------------------------

## Visualization

### Visualize Non-Catalytic Genes **\[MEDIUM PRIORITY\]**

**Problem:** Pathways include non-catalytic genes (chaperones,
regulators, electron carriers, maturation factors) that are listed in
the `genes` array but don’t appear in network plots because they aren’t
in edges.

**Current behavior:** - Only genes with edges (catalytic enzymes,
complexes) are visualized - Non-catalytic genes like nifF (electron
carrier), regulatory genes, chaperones are invisible in plots - They ARE
tracked in scoring and appear in results tables

**Why it matters:** - Accessory genes provide biological context (e.g.,
nifF shows electron transfer pathway) - Users want to see “full picture”
of pathway genes, not just catalytic core - Currently appears incomplete
(“where’s nifF? it’s in the genes list!”)

**Proposed solutions:**

**Option 1: Ghost nodes** - Non-catalytic genes appear as
dimmed/outlined nodes - Not connected by edges (float near relevant
enzymes) - Tooltip shows role: “electron carrier”, “chaperone”, etc. -
Visual distinction from catalytic genes

**Option 2: Sidebar panel** - Network plot shows catalytic genes only -
Separate “Accessory Genes” panel lists non-catalytic genes with
detection status - Less cluttered but separates information

**Option 3: Optional layer toggle** - Default: catalytic genes only -
Toggle “Show accessory genes” adds ghost nodes - User controls
complexity

**Implementation:**

Identify non-catalytic genes: in `genes` array but not referenced in
pathway `edges`

Add visual styling for non-catalytic nodes (dimmed, dashed border)

Position near related complexes (e.g., nifF near nifHDK)

Update legend to distinguish catalytic vs. accessory genes

Consider interactive toggle for complexity control

**Example:** Nitrogen fixation with nifF (electron carrier) would show
nifF as ghost node near nifHDK complex.

------------------------------------------------------------------------

### V2 Coordinate System **\[MEDIUM PRIORITY\]**

**Problem:** V1 had coordinate export → curated layouts. Not carried to
V2. Workflow broken.

**Solution:**

Add export button to
[`plot_potato_interactive()`](https://jeffkimbrel.github.io/potato/reference/plot_potato_interactive.md)

Update
[`update_potato_coordinates()`](https://jeffkimbrel.github.io/potato/reference/update_potato_coordinates.md)
for V2

Update
[`align_coordinates()`](https://jeffkimbrel.github.io/potato/reference/align_coordinates.md)
for V2

Normalize pixel coords to plot coords

**Workflow:**

``` r

plot_potato(pot, interactive = TRUE, layout = "fr")
# Arrange nodes → Export → Import
update_potato_coordinates("pathway.json", "coords.json")
```

------------------------------------------------------------------------

### Near-Miss Pathway Visualization **\[LOW PRIORITY\]**

**Goal:** Visual way to identify pathways close to detection threshold

**Challenge:** Doesn’t scale well - faceting by genome fails with 100+
genomes, faceting by pathway fails with dozens of pathways

**Current workaround:**
[`find_near_miss_pathways()`](https://jeffkimbrel.github.io/potato/reference/find_near_miss_pathways.md)
returns filtered tibble, `plot_genome_pathways()` shows per-genome view
with thresholds

**Future idea:** Interactive filtering or small-multiples view for
targeted subsets

------------------------------------------------------------------------

## LLM Integration

### Build-Potato Agent **\[MEDIUM PRIORITY\]**

**Goal:** LLM-assisted potato creation from KEGG modules or literature

**Implementation:**

Parse KEGG module structure → potato draft

Extract gene IDs, EC numbers, compounds

Handle complex KEGG syntax (branches, alternatives)

Suggest BLAST/HMM when KEGG insufficient

User reviews/edits → validates

------------------------------------------------------------------------

### Result Interpretation Agent **\[LOW PRIORITY\]**

**Goal:** Natural language summaries of results

**Implementation:**

Input: scores + optional MAG completeness

Output: “Genome has complete ED (classic variant). Missing edd but MAG
60% complete, likely present.”

Flag uncertainty, suggest follow-up

------------------------------------------------------------------------

## Advanced Features

### Embed Detection Data in Potatoes **\[LOW PRIORITY\]**

**Goal:** Ultimate portability - potato contains all HMM/BLAST data

**Implementation:**

Optional fields: `hmm_profile_data`, `blast_sequences`

Makes tools.json optional

Build custom databases from scratch

------------------------------------------------------------------------

### Multi-Sack Comparisons **\[LOW PRIORITY\]**

**Goal:** Compare annotation results across projects/studies

**Implementation:**

`compare_sacks(sack_list)`

Shared vs unique pathways

Batch visualization

Meta-analysis helpers

------------------------------------------------------------------------

### Inter-Potato Networks **\[LOW PRIORITY\]**

**Goal:** Link potatoes via compound flow (serine uses ED outputs)

**Implementation:**

Build inter-potato compound graph

Evaluate higher-order metabolic capabilities

“Genome has glycolysis + TCA + amino acids = autotroph”

------------------------------------------------------------------------

## Priority Levels

- **HIGH** = Data integrity, critical bugs, blocking workflows
- **MEDIUM** = Valuable features, improves accuracy/usability  
- **LOW** = Nice to have, future enhancements

------------------------------------------------------------------------

**Version history:** See git log and CLAUDE.md  
**Implementation details:** See CLAUDE.md for agent context
