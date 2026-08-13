# POTATO Provenance Tracking

**Version:** v0.10.2  
**Date:** 2026-08-13

## Overview

POTATO now tracks complete provenance for all annotation and analysis steps. This allows users to:
- Reproduce analyses exactly
- Identify gaps in annotation coverage
- Track which tool versions were used
- Understand what was run when

---

## Provenance Structure

Provenance data is stored in `sack@provenance` with the following structure:

### Genomes (`sack@provenance$genomes`)

Tracks each time genomes were added. List of entries (appended on each `add_genomes()` call):

```r
list(
  timestamp = "2026-08-13 10:23:45",
  n_added = 22,
  genome_names = c("genome1", "genome2", ...),
  genome_paths = c("/full/path/to/genome1.faa", "/full/path/to/genome2.faa", ...)
)
```

### Annotation Tools (`sack@provenance$kofam`, `$blast`, `$hmm`)

Tracks the latest annotation run for each tool (overwrites on re-run):

```r
list(
  timestamp = "2026-08-13 10:25:12",
  tool_version = "kofamscan 1.3.0",
  conda_env = "potato",
  workers = 8,
  potatoes_requested = c("bhac_cycle", "nitrogen_fixation"),  # What user asked for
  potatoes_with_genes = c("bhac_cycle", "nitrogen_fixation"), # Which had genes for this tool
  n_genomes = 22,
  n_kos = 45,  # or n_subjects, n_profiles
  
  commands = list(
    # Reproducible command templates with {placeholders}
    exec_annotation_template = "conda run -n potato exec_annotation ... {genome_path} ...",
    genome_paths = c("/path/to/genome1.faa", ...),
    hal_path = "/path/to/hash.hal",
    hal_content = c("/path/to/K00001.hmm", ...)  # Actual content for reproduction
  )
)
```

**Key fields:**
- `potatoes_requested`: Potatoes the user specified (or all if NULL)
- `potatoes_with_genes`: Potatoes that actually contributed genes for this database
  - If in `requested` but NOT in `with_genes` → potato has no genes for this tool
  - If NOT in `requested` → potato wasn't searched

### Scoring (`sack@provenance$scoring`)

Tracks the latest scoring run (overwrites on re-run):

```r
list(
  timestamp = "2026-08-13 10:30:00",
  thresholds = list(
    kofam_threshold = NULL,  # NULL = per-gene default
    blast_evalue = "1e-10",
    blast_bitscore = 50,
    hmm_evalue = "1e-10"
  ),
  n_pathways = 15,
  n_genomes = 22
)
```

---

## Using Provenance

### View Provenance Summary

```r
print_provenance(sack)
```

Shows formatted summary of all provenance data including:
- When each step was run
- Tool versions used
- Which potatoes were checked
- Command templates for reproduction

### Check Annotation Coverage

```r
plot_annotation_coverage(sack)
```

Creates a heatmap showing:
- **Rows:** All pathways (potato:pathway for multi-pathway networks)
- **Columns:** kofam, blast, hmm databases
- **Colors:**
  - 🟢 **Green**: Pathway was checked (potato in `potatoes_requested` AND has genes)
  - ⚫ **Gray**: Pathway has no genes for this database (N/A)
  - 🔴 **Red**: Pathway HAS genes but was NOT checked (annotation gap!)

This immediately shows:
- What's fully annotated
- Where you have gaps
- What can't be checked (no genes)

### Reproduce Annotations

The command templates can be extracted and run outside POTATO:

```r
# Extract kofam command template
kofam_cmd <- sack@provenance$kofam$commands$exec_annotation_template
genome_paths <- sack@provenance$kofam$commands$genome_paths

# Run for each genome (example)
for (genome_path in genome_paths) {
  cmd <- gsub("\\{genome_path\\}", genome_path, kofam_cmd)
  # Also substitute other placeholders: {hal_path}, {ko_list}, {tmp_dir}
  system(cmd)
}
```

### Check Tool Versions

```r
sack@provenance$kofam$tool_version
sack@provenance$blast$tool_version
sack@provenance$hmm$tool_version
```

### Identify Annotation Gaps

```r
# Which potatoes were requested but have no kofam genes?
requested <- sack@provenance$kofam$potatoes_requested
with_genes <- sack@provenance$kofam$potatoes_with_genes
no_kofam <- setdiff(requested, with_genes)

# Which potatoes were never checked?
all_potatoes <- names(sack@potatoes)
not_requested <- setdiff(all_potatoes, requested)
```

---

## Example Workflow with Provenance

```r
# 1. Initialize and add genomes
sack <- create_sack("~/project")
sack <- add_genomes(sack, "~/genomes/*.faa")
# Provenance: genomes[[1]] records 22 genomes added

# 2. Run kofam on subset
sack <- run_kofam(sack, potato_names = c("nitrogen_fixation", "bhac_cycle"))
# Provenance: kofam records which potatoes requested, which had KOs, tool version

# 3. Check coverage - identify gaps
plot_annotation_coverage(sack)
# Shows: nitrogen_fixation needs blast, bhac_cycle fully covered with kofam

# 4. Fill gaps
sack <- run_blast(sack, potato_names = "nitrogen_fixation")
# Provenance: blast records this run

# 5. Run remaining tools on all potatoes
sack <- run_hmm(sack)
# Provenance: hmm records all potatoes checked

# 6. Score and check what thresholds were used
sack <- score_pathways(sack)
print_provenance(sack)
# Shows exact thresholds used for scoring

# 7. Review complete provenance
print_provenance(sack)
# Outputs:
#   - Genomes: 22 added from ~/genomes/*.faa
#   - Kofam: 2 potatoes, version 1.3.0, 45 KOs
#   - BLAST: 1 potato, version 2.12.0, 12 subjects
#   - HMM: 7 potatoes, version 3.3.2, 89 profiles
#   - Scoring: thresholds used
```

---

## Benefits of Provenance Tracking

### Reproducibility
- Command templates with placeholders can be re-run outside POTATO
- Tool versions recorded for exact reproduction
- Database contents saved (hal_content, fasta_content, profile_content)

### Quality Control
- `plot_annotation_coverage()` immediately shows gaps
- Know which potatoes can't be checked (no genes for database)
- Distinguish "not found" from "not searched"

### Auditing
- Know exactly what was run when
- Track changes over time (genomes appended, tools re-run)
- Share sacks with collaborators - they can see full history

### Debugging
- If results look wrong, check provenance
- Verify correct thresholds were used
- Confirm expected potatoes were searched

---

## Future Enhancements

Potential additions to provenance:
- **Database versions:** Track kofam profiles version, BLAST database version
- **File hashes:** MD5 checksums of input files for change detection
- **Git integration:** Track git commit when analysis was run
- **User tracking:** Record who ran the analysis
- **Comparison:** Functions to diff provenance between sacks

See ROADMAP.md for details.

---

## Technical Notes

### Why Overwrite vs. Append?

**Annotation tools (kofam, blast, hmm):** Provenance is **overwritten** on each run because:
- Only the latest run matters for current results
- Results in `sack@results` are replaced, not appended
- Keeps provenance in sync with results

**Genomes:** Provenance is **appended** because:
- `add_genomes()` can be called multiple times
- Genomes accumulate in `sack@genomes`
- History of additions is useful

### Storage Location

All provenance is in `sack@provenance`, which is part of the PotatoSack S7 object. When you `saveRDS(sack)`, provenance is saved with it.

### Command Template Placeholders

Templates use `{placeholder}` syntax:
- `{genome_path}` - path to genome file
- `{hal_path}` - path to .hal file
- `{blast_db}` - path to BLAST database
- `{hmm_profile}` - path to HMM profile
- `{ko_list}` - path to ko_list file
- `{tmp_dir}` - temporary directory
- `{tblout}` - output table file

Values for placeholders are in `commands$genome_paths`, `commands$hal_path`, etc.

---

Last updated: 2026-08-13
