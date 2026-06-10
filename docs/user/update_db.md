# Updating Databases

HLAnte queries four databases. Local copies must be present for offline or
cached operation. Use `hlante db-update` to download or refresh them.

## Download all databases

```bash
hlante db-update
```

This downloads IPD-IMGT/HLA, PharmGKB, GWAS Catalog, and AFND in sequence.
Requires internet access. Total download size is approximately 200 MB.

## Update a single database

```bash
hlante db-update --db imgt       # IPD-IMGT/HLA allele table
hlante db-update --db pharmgkb   # PharmGKB clinical annotations
hlante db-update --db gwas       # GWAS Catalog bulk TSV
hlante db-update --db afnd       # AFND allele frequency table
```

## Force re-download

By default HLAnte skips databases that are already up-to-date. Use `--force`
to download regardless of age:

```bash
hlante db-update --force
hlante db-update --db gwas --force
```

## Custom storage locations

By default databases are stored in `~/.hlante/`. To override:

```bash
hlante db-update \
    --imgt-dir     /data/hlante/imgt/ \
    --pharmgkb-dir /data/hlante/pharmgkb/ \
    --gwas-cache-dir /data/hlante/gwas/ \
    --afnd-dir     /data/hlante/afnd/
```

Pass the same path flags to `hlante annotate` when running offline:

```bash
hlante annotate -i sample.json -t arcashla \
    --offline \
    --imgt-db-path /data/hlante/imgt/ \
    --pharmgkb-dir /data/hlante/pharmgkb/
```

## Check installed versions

```bash
hlante version
```

Example output:

```
HLAnte v0.1.0
  IPD-IMGT/HLA : 3.64.0  (/home/user/.hlante/imgt_hla)
  PharmGKB     : installed  (/home/user/.hlante/pharmgkb)
  GWAS cache   : 12 file(s) (/home/user/.hlante/cache/gwas)
  AFND         : installed  (/home/user/.hlante/afnd)
```

## Database notes

### IPD-IMGT/HLA

The IPD-IMGT/HLA database releases quarterly. HLAnte will warn if the local
copy is more than 6 months old. Update with:

```bash
hlante db-update --db imgt
```

### AFND

AFND (allelefrequencies.net) does not publish a machine-readable bulk
download. HLAnte ships a curated built-in frequency table (7 loci × 8
population groups) derived from Gonzalez-Galarza et al. 2020. This table
is embedded in the package and does not require downloading.

If you have a custom AFND-format TSV, point to it with:

```bash
hlante db-update --db afnd --afnd-url https://your-server/afnd_custom.tsv
# or for a local file:
hlante db-update --db afnd --afnd-url file:///data/afnd_local.tsv
```
