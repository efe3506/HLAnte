# Validating Input Files

Use `hlante validate` to check whether an input file or directory is
parseable — without running the full annotation pipeline.

## Basic usage

```bash
hlante validate -i sample.genotype.json -t arcashla
```

On success, the command exits with code 0 and prints a summary:

```
✓ Input is valid
  Tool: arcashla
  Files: 1
  Samples: 1
  Locus calls: 10
  Unique loci: HLA-A, HLA-B, HLA-C, HLA-DQB1, HLA-DRB1
```

On failure, the command exits with code 1 and prints the parse error:

```
ERROR: Invalid allele format: 'B*invalid'
```

## Validate a directory

When `--input` is a directory, all tool-appropriate files inside it are
validated:

```bash
hlante validate -i arcashla_outputs/ -t arcashla
```

## Supported tools

Pass the same `-t` / `--tool` argument as for `annotate`:

```bash
hlante validate -i sample.genotype.json   -t arcashla
hlante validate -i sample_genotype.tsv    -t t1k
hlante validate -i sample_final.result.txt -t hlahd
hlante validate -i sample_result.tsv      -t optitype
```

## When to use validate

- Before a large batch run, to detect malformed files early.
- After receiving typing results from an external lab, to confirm
  the format is compatible with HLAnte.
- In CI/CD pipelines to gate on input quality before annotation.
