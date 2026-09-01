# Genome File S7 Class

An S7 class representing a genome file with metadata extracted from
jakomics FILE objects. This class is serialization-safe (unlike Python
objects from jakomics).

## Usage

``` r
GenomeFile(
  short_name = character(0),
  file_path = character(0),
  name = character(0),
  file_type = character(0),
  md5 = ""
)
```

## Fields

- `short_name`:

  Short genome name (e.g., "GCA_001234")

- `file_path`:

  Absolute path to genome file

- `name`:

  Full file name with extension

- `file_type`:

  File extension (e.g., "faa", "fasta")

- `md5`:

  MD5 hash of file contents (optional, for change tracking)
