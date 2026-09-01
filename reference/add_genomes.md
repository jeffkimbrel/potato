# Add genomes to a potato sack

Register genome files with a potato sack without copying them. Can be
called multiple times to add genomes from different locations.

## Usage

``` r
add_genomes(sack, path, validate = TRUE, recursive = FALSE)
```

## Arguments

- sack:

  PotatoSack object

- path:

  Character. Path to protein FASTA (.faa or .fasta) files. Can be:

  - A directory: adds all .faa/.fasta files in directory

  - A wildcard pattern: "~/data/mags/\*.faa"

  - A single file: "~/data/genome1.faa"

  - A vector of files

- validate:

  Logical. If TRUE, validates .faa files (default: TRUE)

- recursive:

  Logical. If path is a directory, search recursively (default: FALSE)

## Value

Modified PotatoSack object with genomes added
