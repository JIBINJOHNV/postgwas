# PostGWAS

PostGWAS is a command-line toolkit that provides a collection of modules commonly used in post-GWAS analyses.
It supports tasks such as LD-block annotation, fine-mapping, summary-statistic filtering, imputation, MAGMA gene analysis,
PoPS scoring, and more. Each module is available as a standalone command, and a unified **pipeline** interface allows
chaining multiple steps together.

## Available Modules

```
┏━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Command        ┃ Description                               ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ annot_ldblock  │ LD-block annotation                       │
│ finemap        │ Fine-mapping (SuSiE / FINEMAP)            │
│ flames         │ FLAMES integrative scoring                │
│ formatter      │ Convert VCF to tool-specific formats      │
│ harmonisation  │ Harmonisation of input summary statistics │
│ heritability   │ LDSC-based heritability estimation        │
│ imputation     │ Summary-statistic imputation              │
│ ld_clump       │ LD clumping                               │
│ magma          │ MAGMA gene/pathway analysis               │
│ magmacovar     │ MAGMA gene-property model                 │
│ manhattan      │ Manhattan/QQ plot generation              │
│ pipeline       │ Multi-step summary-statistic pipeline     │
│ pops           │ PoPS gene-prioritisation                  │
│ qc             │ QC summary reporting                      │
│ sumstat_filter │ Summary-statistic filtering               │
└────────────────┴───────────────────────────────────────────┘
```

## Basic Usage

Display list of all availble module:

```sh
postgwas --help
docker run --platform=linux/amd64 -it jibinjv/postgwas:1.0 postgwas --help
```

Display help for any module:

```sh
postgwas <module> --help

docker run --platform=linux/amd64 -it jibinjv/postgwas:1.0 postgwas finemap --help
```

Example:

```sh
postgwas finemap --help
```

## Pipeline Usage

```sh
postgwas pipeline --modules finemap
docker run --platform=linux/amd64 -it jibinjv/postgwas:1.0 postgwas pipeline --modules flames --help

```

Optional flags:

- --apply-filter
- --apply-imputation
- --apply-manhattan
- --heritability

Example:

```sh
postgwas pipeline --modules finemap --apply-filter --apply-imputation
docker run --platform=linux/amd64 -it jibinjv/postgwas:1.0 postgwas pipeline --modules flames --apply-filter --apply-imputation --help
```

## Input Requirements

- Harmonised GWAS summary-statistics **VCF**
- Module-specific reference files (LD blocks, LD panels, MAGMA resources, etc.)

Intermediate files are generated automatically by modules.

## Installation (development mode)

```sh
git clone https://github.com/JIBINJOHNV/postgwas.git
cd postgwas
pip install -e .
```

Or install directly from GitHub:

```sh
pip install -U git+https://github.com/JIBINJOHNV/postgwas.git
```

## Project Structure

```
postgwas/
    annot_ldblock/
    finemap/
    formatter/
    harmonisation/
    pipeline/
    pops/
    magma/
    ld_clump/
    qc_summary/
```

## License

Specify your desired license (MIT, Apache-2.0, etc.)

## Notes

- For research use only.
- Ensure genome builds match reference data.
- Some tools (MAGMA, LDSC, PLINK) require separate installation.