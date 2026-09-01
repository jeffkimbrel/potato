# Score pathway presence/absence across all genomes

Applies quality thresholds to annotation hits and calculates
pathway-level completion scores. For multi-pathway networks, scores each
pathway independently.

## Usage

``` r
score_pathways(
  sack,
  kofam_threshold = NULL,
  blast_evalue = 1e-10,
  blast_bitscore = 50,
  hmm_evalue = 1e-10
)
```

## Arguments

- sack:

  PotatoSack object with annotation results

- kofam_threshold:

  Score threshold for kofam hits (NULL = use per-gene threshold)

- blast_evalue:

  E-value threshold for BLAST hits (default: 1e-10)

- blast_bitscore:

  Bitscore threshold for BLAST hits (default: 50)

- hmm_evalue:

  E-value threshold for HMM hits without TC (default: 1e-10)

## Value

Modified PotatoSack with scores in @scores. For multi-pathway networks,
scores tibble includes 'pathway' and 'pathway_name' columns with one row
per pathway per genome. Scoring includes both all-steps metrics
(total_steps_detected, total_steps, fraction, present) and required-only
metrics (essential_total_steps_detected, essential_steps,
essential_fraction, essential_pathway_present).

## Details

Handles OR branches (alternative genes), required vs optional genes, and
multi-pathway networks where genes are shared across pathways.

Threshold priority:

- Kofam: Uses per-gene threshold from KEGG (can override with
  kofam_threshold)

- HMM: Uses per-profile TC (trusted cutoff) if available, otherwise
  hmm_evalue

- BLAST: Uses global blast_evalue and blast_bitscore
