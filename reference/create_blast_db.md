# Create filtered BLAST database from potato reference sequences

Extracts only the reference sequences needed by potatoes from the
configured BLAST files, creates a filtered FASTA, and builds a BLAST
database.

## Usage

``` r
create_blast_db(sack, potato_names = NULL)
```

## Arguments

- sack:

  PotatoSack object

- potato_names:

  Vector of potato names (NULL = all)

## Value

List with blast_db path and modified sack
