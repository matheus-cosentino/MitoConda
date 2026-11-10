# MitoConda - Mitochondrial Amplicon Nanopore Pipeline

![MitoConda Logo](resources/logo/MitoConda_Logo.jpeg)

[![Language](https://img.shields.io/badge/language-Python%20%26%20Shell-blue)](https://www.python.org/)
[![Workflow](https://img.shields.io/badge/workflow-Snakemake-orange)](https://snakemake.readthedocs.io/)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)

---

## Overview

**MitoConda** is a modular, reproducible **Snakemake** pipeline designed for the analysis of mitochondrial amplicon (e.g., **12S**, **MTCYB**) **Oxford Nanopore long-read** sequencing data. It integrates quality control, consensus clustering, taxonomic assignment via Lowest Common Ancestor (LCA), reference database building, multiple sequence alignment, and phylogenetic reconstruction—all controlled through a single wrapper script (`MitoConda.sh`).

---

## Table of Contents
- [Features & Workflow Steps](#features--workflow-steps)
- [Installation](#installation)
- [Preparing Your Data & References](#preparing-your-data--references)
- [Usage Guide](#usage-guide)
  - [CLI Options Reference](#cli-options-reference)
  - [Common Workflow Recipes](#common-workflow-recipes)
  - [Important: Matching `--taxid`, `phylo_names`, and `--lca-rank`](#important-matching---taxid-phylo_names-and---lca-rank)
- [Configuration (`config/config.yaml`)](#configuration-configconfigyaml)
- [Output Structure](#output-structure)
- [Contributing](#contributing)

---

## Features & Workflow Steps

1. **Quality Control (`--qc`)**
   - Filters raw Nanopore reads by quality and length using `Chopper`.
   - Generates read statistics before and after filtering using `NanoStat`.
   - Aggregates QC and taxonomic summaries into an interactive `MultiQC` report.
2. **Clustering & Chimera Removal (`--diversity` / `--phylo`)**
   - Dereplicates reads, clusters consensus sequences by identity (`VSEARCH`), and removes *de novo* chimeras (`uchime_denovo`).
3. **Homology Search & LCA Taxonomic Assignment (`--diversity` / `--phylo`)**
   - Queries non-chimeric consensus sequences against the NCBI mitochondrial database (`blastn`) with automatic download/setup of NCBI `mito` and `taxdb` databases.
   - Computes Lowest Common Ancestor (LCA) lineages using `TaxonKit` (`lca` → `lineage` → `reformat`).
   - Splits and trims sample amplicon sequences (removing Nanopore adapters/barcodes and orienting strands) into taxon-specific FASTA files (`Fasta_by_LCA_<pident>_<rank>/`).
4. **Reference Database & Phylogeny (`--build-db-only` / `--phylo`)**
   - Uses a reference gene bait (`resources/genes/<gene>.fasta`) to extract homologous regions for the target TaxID from the NCBI mitochondrial BLAST database (`blastdbcmd`).
   - Builds a reference multiple sequence alignment (`MAFFT --auto`).
   - Adds sample LCA-classified fragments into the reference alignment (`MAFFT --addfragments --adjustdirection --keeplength`), trims alignment columns (`trimAl`), and reconstructs a Maximum-Likelihood phylogenetic tree (`FastTree`).

---

## Installation

Clone the repository:

```bash
git clone https://github.com/matheus-cosentino/MitoConda.git
cd MitoConda
```

> **Automatic Conda Setup**: You do **not** need to manually create Conda environments before running `MitoConda.sh`. When you run `bash MitoConda.sh ...`, the script automatically checks for (and creates if missing) the base `mitoconda` environment from `workflow/envs/mitoconda.yaml`, and Snakemake automatically builds all tool-specific environments (`chopper`, `vsearch`, `blast`, `taxonkit`, `mafft`, `trimal`, `fasttree`, `multiqc`) inside `.snakemake/conda/`.

If you prefer to create the base environment manually beforehand:

```bash
conda env create -n mitoconda -f workflow/envs/mitoconda.yaml
conda activate mitoconda
```

---

## Preparing Your Data & References

Before running the pipeline, check these two requirements:

### 1. Input FASTQ Files
Place your raw Nanopore reads in your input folder (default: `data/` or any directory passed to `--input`).
- Files **must** end with `.fastq.gz` (e.g., `data/CaviCor/barcode01.fastq.gz`, `data/CaviCor/barcode02.fastq.gz`).

### 2. Reference Gene Bait (`resources/genes/<gene>.fasta`)
To build the reference alignment for phylogeny (`--phylo` or `--build-db-only`), place a FASTA file containing a reference sequence of your target gene in `resources/genes/`:
- The filename **must** match the `--gene` argument:
  - `--gene 12s` → requires `resources/genes/12s.fasta`
  - `--gene MTCYB` → requires `resources/genes/MTCYB.fasta`

---

## Usage Guide

All executions use `MitoConda.sh`, which generates a dynamic `run_overrides.yaml` file on top of `config/config.yaml` and runs Snakemake.

```bash
bash MitoConda.sh --output <OUTPUT_DIR> [OPTIONS]
```

### CLI Options Reference

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--output <DIR>` | **Yes** | — | Output directory where all results and reports will be stored. |
| `--input <DIR>` | No | `data` | Directory containing input `.fastq.gz` files. |
| `--qc` | No | `false` | Run Quality Control (`Chopper`, `NanoStat`, `MultiQC`). |
| `--diversity` | No | `false` | Run BLAST & TaxonKit LCA abundance tables and MultiQC summaries. |
| `--phylo` | No | `false` | Run full phylogenetic workflow (QC/VSEARCH/BLAST/LCA + Reference DB + MAFFT + trimAl + FastTree). |
| `--build-db-only` | No | `false` | Only extract and build the reference gene FASTA (`resources/blast_db/<gene>/<taxid>.fasta`). |
| `--taxid <INT>` | No | `config.yaml` | Target NCBI TaxID for phylogeny (e.g., `40674` for Mammalia, `9397` for Chiroptera). |
| `--gene <STR>` | No | `config.yaml` | Target mitochondrial gene name matching `resources/genes/<gene>.fasta` (e.g., `12s`, `MTCYB`). |
| `--lca-rank <STR>` | No | `config.yaml` | Taxonomic rank used to split LCA FASTAs (`phylum`, `class`, `order`, `family`, `genus`, `species`). |
| `--jobs <INT>` | No | `15` | Maximum number of parallel Snakemake jobs/cores. |
| `--temp-dir <DIR>` | No | `/tmp/${USER}_cavicor` | Temporary directory for intermediate and shadow files. |
| `-h`, `--help` | No | — | Display the help message. |

---

### Common Workflow Recipes

#### 1. Phylogenetic Analysis for a Specific Taxonomic Group (`--phylo`)
Runs read filtering, VSEARCH clustering, BLAST against mitochondrial DB, LCA extraction, reference database alignment, fragment insertion (`MAFFT`), trimming (`trimAl`), and tree building (`FastTree`):

```bash
bash MitoConda.sh \
    --phylo \
    --output cavicor_output \
    --input data/CaviCor \
    --gene 12s \
    --taxid 40674 \
    --lca-rank class
```

#### 2. Complete Pipeline (QC + Diversity + Phylogeny)
Runs all modules and generates the complete `MultiQC` report alongside phylogenetic trees:

```bash
bash MitoConda.sh \
    --qc \
    --diversity \
    --phylo \
    --output cavicor_output \
    --input data/CaviCor \
    --gene 12s \
    --taxid 40674 \
    --lca-rank class
```

#### 3. Only Build the Reference Gene Database (`--build-db-only`)
Extracts and formats homologous sequences for a gene and TaxID from the NCBI mitochondrial BLAST database into `resources/blast_db/<gene>/<taxid>.fasta` without processing sample reads:

```bash
bash MitoConda.sh \
    --build-db-only \
    --output cavicor_output \
    --gene 12s \
    --taxid 40674
```

#### 4. Only Quality Control & Taxonomic Diversity (`--qc --diversity`)
Generates read QC metrics, BLAST/LCA abundance tables, and the MultiQC report without building phylogenetic trees:

```bash
bash MitoConda.sh \
    --qc \
    --diversity \
    --output cavicor_output \
    --input data/CaviCor \
    --lca-rank order
```

---

### Important: Matching `--taxid`, `phylo_names`, and `--lca-rank`

When running `--phylo`, the pipeline selects sample sequences assigned to your target taxon from `<output>/<sample>/Fasta_by_LCA_<pident>_<lca_rank>/<TaxonName>.fasta`.

Make sure that:
1. The TaxID is mapped to its scientific name under `phylo_names` in `config/config.yaml`:
   ```yaml
   phylo_names:
     40674: "Mammalia"      # Rank: class
     9397: "Chiroptera"     # Rank: order
     9398: "Pteropodidae"   # Rank: family
   ```
2. Your `--lca-rank` matches the taxonomic rank of that target taxon:
   - For `40674` (`Mammalia`), use `--lca-rank class`
   - For `9397` (`Chiroptera`), use `--lca-rank order`
   - For `9398` (`Pteropodidae`), use `--lca-rank family`

If you add a new TaxID, add its ID and scientific name to `phylo_names` in `config/config.yaml`.

---

## Configuration (`config/config.yaml`)

Default pipeline parameters are defined in `config/config.yaml` and can be customized directly:

| Section / Parameter | Default | Description |
| :--- | :--- | :--- |
| `chopper.min_quality` | `[10]` | Minimum mean read Phred quality score. |
| `chopper.min_length` | `[300]` | Minimum read length (bp). |
| `chopper.max_length` | `[800]` | Maximum read length (bp). |
| `vsearch.identity` | `[0.9]` | Clustering identity threshold (`0.9` = 90%). |
| `vsearch.minsize` | `[10]` | Minimum cluster abundance size retained after clustering. |
| `blast.perc_identity` | `[70]` | Minimum BLAST percentage identity (`pident`). |
| `blast.evalue` | `[0.001]` | Maximum BLAST E-value threshold. |
| `blast.max_target_seqs` | `[5]` | Maximum number of aligned reference targets kept per consensus read. |
| `phylo_target` | `["40674"]` | Default TaxID(s) for phylogeny (overridden by `--taxid`). |
| `phylo_names` | `{40674: "Mammalia", ...}` | Mapping of TaxID to scientific name used to pick the LCA FASTA file. |
| `phylo_genes` | `["12s"]` | Default gene(s) for phylogeny (overridden by `--gene`). |
| `taxonkit.lca_rank` | `["order"]` | Default taxonomic rank for LCA splitting (overridden by `--lca-rank`). |

---

## Output Structure

```text
<output_dir>/
├── Phylo/
│   └── <taxid>_<gene>_ref_aln.fasta                       # Reference database alignment (MAFFT)
├── <sample>/
│   ├── QC/                                                # Chopper filtered FASTQs & NanoStat reports
│   ├── Vsearch/                                           # Dereplicated, clustered, and non-chimeric FASTAs
│   ├── Blast/                                             # Raw BLASTn output tables
│   ├── LCA/                                               # TaxonKit LCA + formatted lineage tables
│   ├── Abundance/                                         # Read abundance summaries per taxonomic rank
│   ├── Fasta_by_LCA_<pident>_<rank>/                      # Amplicon sequences trimmed & split by taxon (e.g., Mammalia.fasta)
│   └── Phylo/
│       ├── <sample>_<pident>_<taxid>_renamed.fasta        # Sample sequences for target taxon with clean headers
│       ├── <sample>_<taxid>_<pident>_<gene>_Aligned.fasta # Sample fragments aligned against reference DB
│       ├── <sample>_<taxid>_<pident>_<gene>_Aligned_cleaned.fasta # trimAl cleaned alignment
│       └── <sample>_<taxid>_<pident>_<gene>_aligned_tree.nwk      # FastTree phylogenetic tree (Newick format)
├── multiqc_all/
│   └── <pident>_<rank>_multiqc_report.html                # Interactive MultiQC report (when --qc or --diversity is enabled)
└── mitoconda_report.html                                  # Snakemake execution report & workflow provenance
```

Generated `.nwk` tree files can be inspected directly in [iTOL](https://itol.embl.de/), [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), or `ggtree` in R.

---

## Contributing

Contributions are welcome! Please fork the repository, create a feature branch, and submit a pull request. For major changes, open an issue first to discuss your ideas.
