# Repeat-aware mapping pipeline

## Usage

```bash

snakemake --use-conda --use-singularity --cores 4 --notemp -kpn
```

## Dependencies

### CSEM

The trickiest is CSEM v2.4, which must be compiled and installed separately. By default, the program expects that the CSEM directory is located as follows:

```
$HOME/bin/csem-2.4
```

This can be adjusted in `config/config.yaml`.