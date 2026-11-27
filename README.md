# MitoConda — Mitochondrial Amplicon Nanopore Pipeline

![MitoConda Logo](resources/logo/MitoConda_Logo.jpeg)

[![Language](https://img.shields.io/badge/language-Python%20%26%20Bash-blue)](https://www.python.org/)
[![Workflow](https://img.shields.io/badge/workflow-Snakemake-orange)](https://snakemake.readthedocs.io/)
[![Report](https://img.shields.io/badge/report-MultiQC%20%2B%20Plotly-brightgreen)](https://multiqc.info/)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)

---

## Overview

**MitoConda** is a modular, reproducible **Snakemake** pipeline designed for the analysis of mitochondrial amplicon (e.g., **12S rRNA**, **MTCYB**) **Oxford Nanopore long-read** sequencing data to characterize broad vertebrate and animal biodiversity (mammals, bats, rodents, primates, carnivores, birds, reptiles, amphibians, and bony/cartilaginous fishes).

It integrates read quality control, consensus clustering and chimera removal, sequence retention/filtering audits, taxonomic assignment via Lowest Common Ancestor (LCA), analytical Hurlbert–Heck rarefaction curves, deep order- and family-stratified reference database construction (>210 vertebrate lineages), multiple sequence alignment, Maximum-Likelihood phylogenetic reconstruction, sequence-by-sequence method comparison (**BLAST Top Hit vs. TaxonKit LCA vs. Phylogeny**), and whole-database pairwise identity/genetic distance profiling—all controlled through a single wrapper script (`MitoConda.sh`) and summarized in a 5-part interactive **MultiQC** report.

---

## Table of Contents
- [Features & Workflow Modules](#features--workflow-modules)
- [MultiQC Report Structure (5 Main Sections)](#multiqc-report-structure-5-main-sections)
- [Project Architecture](#project-architecture)
- [Installation](#installation)
- [Preparing Your Data & References](#preparing-your-data--references)
- [Usage Guide](#usage-guide)
  - [CLI Options Reference](#cli-options-reference)
  - [Common Workflow Recipes](#common-workflow-recipes)
- [Configuration (`config/config.yaml`)](#configuration-configconfigyaml)
- [Output Directory Structure (`all_samples/` & `individual_samples/`)](#output-directory-structure-all_samples--individual_samples)
- [References](#references)
- [Contributing](#contributing)

---

## Features & Workflow Modules

### 1. Unified Quality Control, Clustering, Sequence Audit & Reference Database Build *(Always Active)*
- **Read Filtering (`Chopper`)**: Filters raw Nanopore `.fastq.gz` reads by Phred quality score (`-q`) and amplicon length range (`-l` to `--maxlength`).
- **Long-Read Metrics (`NanoStat`)**: Computes total sequence yield (Mb), mean/median read length, N50, length standard deviation, mean/median Phred quality scores, and `% reads >Q10` and `>Q15`.
- **Consensus Clustering & Chimera Removal (`VSEARCH`)**: Clusters filtered reads by identity (`--cluster_fast`), filters low-abundance clusters/singletons (`--fastx_filter --minsize`), and removes *de novo* PCR chimeras (`--uchime_denovo`).
- **Unified QC & Sequence Filtering Audit Tables**:
  - Consolidates all **Chopper**, **NanoStat**, and **VSEARCH** metrics into a single MultiQC table (`1.1 General Sequencing, Quality Control & Clustering Summary`).
  - Tracks and audits every consensus sequence across downstream checkpoints (`1.2 Consensus Sequence Filtering & Retention Audit` and `1.3 Detailed List of Removed / Filtered Sequences`), distinguishing sequences removed due to **No BLAST Hit** (`<pident%` or non-mitochondrial), **Non-Target Clade** (e.g., plant mitochondria outside `Vertebrata`), or **trimAl / Alignment Coverage** (`< 60%` non-gap coverage).
- **Hierarchically Stratified Reference Database (`db_phylo.smk` + `build_ref_db.sh`)**:
  - Automatically downloads the NCBI mitochondrial (`mito`) and taxonomy (`taxdump`) databases and extracts a taxonomically balanced reference gene database (`resources/blast_db/<gene>/<taxid>.fasta`).
  - Performs a high-capacity BLASTn search using a 19-sequence multi-class 12S bait panel and climbs the local NCBI SQLite taxonomy tree (`taxonomy4blast.sqlite3` / `nodes.dmp`) bottom-up (`Species -> Genus -> Subfamily -> Family -> Suborder -> Order -> Class`) to allocate mitogenomes across **>210 dedicated vertebrate lineage buckets** (including 33 Rodentia suborders/families/subfamilies, all 20 Chiroptera families, 15 Carnivora families, Primates, Marsupialia, Afrotheria, Xenarthra, Eulipotyphla, 24 Aves orders, Squamata/Testudines/Crocodylia families, 33 Amphibia families, 25 Actinopteri orders, Chondrichthyes, Sarcopterygii, and Cyclostomata).
  - Automatically prioritizes any reference mitogenome accession (`sseqid`) that was a direct BLAST hit in the analyzed samples.

### 2. Homology Search, LCA Taxonomic Diversity & Rarefaction (`--diversity`)
- **Mitochondrial Homology Search (`BLASTn`)**: Queries non-chimeric consensus clusters against the full NCBI mitochondrial database (`-task blastn -word_size 11`).
- **Lowest Common Ancestor Assignment (`TaxonKit`)**: Aggregates all top hits per cluster and resolves full NCBI lineages (`taxonkit lca` → `taxonkit lineage` → `taxonkit reformat`), including automatic class resolution for `Sauropsida` (Testudines and Crocodylia).
- **Analytical Hurlbert–Heck Rarefaction Curves**: Computes exact hypergeometric rarefaction curves ($E[S_m]$) across sequencing depths for every barcode across `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`, and `VSEARCH Clusters (OTUs)`, with interactive Linear and $\log_{10}$ X-axis toggles.
- **Abundance Tables & Strand-Normalized FASTAs**: Produces per-sample abundance tables (`Abundance/`) and splits sample sequences (trimmed to the aligned amplicon coordinates and reverse-complemented when on the minus strand) into taxon-specific FASTA files (`Fasta_by_LCA_<pident>_<rank>/`).

### 3. Phylogenetic Reconstruction & Group FASTA Splitting (`--phylo`)
- **Clade Extraction & Header Summarization**: Automatically extracts all consensus clusters matching the target clade (`phylo_target`, e.g., `7742` for `Vertebrata`) directly from the full LCA lineage—independent of `--lca-rank`—and formats clean headers sorted by read abundance (`barcode01_C01_size2518`).
- **Multiple Sequence Alignment & Barcode-Guided Trimming (`MAFFT` + `trimAl`)**:
  - Aligns the reference database (`mafft --auto`) and inserts sample fragments (`mafft --addfragments --adjustdirection --reorder --op 2.5 --keeplength`) both per sample and integrated across all barcodes (`All_Barcodes`).
  - Crops flanking reference-only columns not covered by sample amplicons, trims gappy columns (`trimal -gt 0.7`), enforces `>= 60%` non-gap alignment coverage to preserve alignment quality (recording any shorter dropped fragments in the audit tables), and formats reference headers as `Species_Accession`.
- **Maximum-Likelihood Phylogeny (`FastTree`)**: Reconstructs ML trees under the `GTR+GAMMA` nucleotide model (`fasttree -nt -gtr -gamma`) for each individual barcode and for all barcodes combined.
- **Interactive Midpoint-Rooted Tree Visualization (`Plotly.js` + `SVG`)**:
  - Performs exact $O(N)$ midpoint rooting on the full integrated tree (excluding rogue artifact branches `> 0.50` subs/site from pulling the root), and prunes distant reference branches to display reference taxa within **5 internal nodes of topological distance** from sample clusters while preserving exact root-to-tip branch lengths.
  - Embeds an interactive **Plotly.js** phylogeny in MultiQC (and standalone HTML) with box-zoom, scroll-zoom, patristic distance hover tooltips, and clickable sample legend filtering.
- **Patristic Distance Placement & Putative Taxonomic Level Inference**:
  - Computes exact patristic distances on the integrated ML tree, attaches full taxonomic lineages (`Class > Order > Family > Genus`), and classifies each cluster/taxon into distance-calibrated confidence levels (`Same Species`, `Same Genus`, `Same Family`, `Same Order`, or `Divergent / Putative NUMT`).
- **Pooled FASTAs Split by Taxonomic Group (`fasta_by_group_<taxid>_<pident>_<gene>_<rank>/`)**:
  - Automatically separates all consensus sequences across all samples into group-specific `.fasta` files according to the requested taxonomic rank (`--lca-rank`, e.g., `Family`, `Order`, `Class`), alongside subfolders separated by closest phylogenetic species (`by_species/`) and by BLAST+LCA assignment (`by_lca_<rank>/`).

### 4. Method Comparison (`BLAST Top Hit vs. LCA vs. Phylogeny`) & Whole-Database Profiling
- **Per-FASTA Sequence-by-Sequence Comparison (`Section 4.1`)**: Compares every 12S consensus FASTA sequence side by side across **BLAST Top Hit** (species, accession, `% identity`), **TaxonKit LCA** (assigned taxon and full lineage), and **FastTree Phylogeny** (closest species, patristic distance, putative level, and automated method concordance classification).
- **Richness, Shannon Diversity ($H'$), Pielou Evenness ($J'$) & Jaccard Overlap (`Sections 4.2–4.4`)**: Compares per-sample and pooled dataset taxonomic richness and read abundances at both the user-defined `--lca-rank` and species level.
- **Whole-Database Pairwise Identity & Genetic Distance Profiling (`Sections 5.1–5.2`)**:
  - Computes vectorized pairwise **Sequence Identity (%)** and **Genetic Distance ($p$-distance, substitutions/site)** across 100% of the aligned reference database (`>= 200 bp` overlap).
  - Provides interactive **Plotly.js** histograms (also openable directly as standalone `.html` files in any browser) across `Hierarchy (All Pairs)`, `Hierarchy (Best Hit / BLAST Top Hit)`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, and `Species` (partitioning the entire database globally and across all vertebrate classes).

---

## MultiQC Report Structure (5 Main Sections)

The aggregated report (`<output_dir>/all_samples/<pident>_multiqc_report.html`) is organized into **5 main analytical sections** configured in [config/multiqc_config.yaml](file:///d:/CaviCor12s/config/multiqc_config.yaml):

1. **General Sequencing & Clustering Quality Control**
   - `1.1 General Sequencing, Quality Control & Clustering Summary (Chopper + NanoStat + VSEARCH)`
   - `1.2 Consensus Sequence Filtering & Retention Audit (VSEARCH -> BLAST -> Target Clade -> trimAl)`
   - `1.3 Detailed List of Removed / Filtered Sequences (Non-Target Clade, trimAl & Top Unmatched Clusters)`
2. **BLAST + LCA Taxonomic Diversity & Rarefaction**
   - `2.1 BLAST Top Hits Metrics`
   - `2.2 Taxonomic Abundance by <Rank> (BLAST + LCA)` *(Stacked Barplot)*
   - `2.3 Taxonomic Rarefaction Curves (BLAST + LCA)` *(Interactive Hurlbert–Heck Plotly.js Curves across Phylum, Class, Order, Family, Genus, Species & VSEARCH Clusters)*
3. **Phylogenetic Analysis**
   - `3.1 Integrated Phylogenetic Tree (All Barcodes vs. Closest References <= 5 Nodes)` *(Interactive Plotly.js)*
   - `3.2 Phylogenetic Abundance by Closest Reference Species` *(Stacked Barplot)*
   - `3.3 Phylogenetic Distance Summary (Per Sample)`
   - `3.4 Phylogenetic Placement: Barcodes vs. All Closest Reference Taxa` *(with Full Lineage & Putative Taxonomic Level)*
   - `3.5 Phylogenetic Distance per Cluster (Closest Reference)` *(with Full Lineage & Putative Taxonomic Level)*
4. **Summary Comparison: BLAST + LCA vs. Phylogeny**
   - `4.1 Per-FASTA Sequence Comparison (BLAST Top Hit vs. LCA vs. Phylogeny)`
   - `4.2 Taxonomic Richness & Diversity per Sample (BLAST+LCA vs. Phylogeny)`
   - `4.3 Overall Taxonomic Abundance by <Rank> (BLAST+LCA vs. Phylogeny)`
   - `4.4 Overall Species-Level Abundance & Richness (BLAST+LCA vs. Phylogeny)`
5. **Reference Database Genetic Distance & Identity Summary**
   - `5.1 Reference DB Genetic Distance & Identity Distribution` *(Interactive Dual-Panel Whole-Database Plotly.js Histogram)*
   - `5.2 Reference DB Genetic Distance & Identity Summary (Phylum, Class, Order, Family, Genus, Species)`
   - `Pipeline References`

---

## Project Architecture

The pipeline separates declarative Snakemake rules (`workflow/rules/`) from modular Python and Bash execution scripts (`workflow/scripts/`):

```text
MitoConda/
├── MitoConda.sh                           # Main CLI wrapper script (environment setup & Snakemake runner)
├── config/
│   ├── config.yaml                        # Default pipeline parameters & taxonomic mappings
│   └── multiqc_config.yaml                # MultiQC 5-topic section ordering (report_section_order & custom_content)
├── resources/
│   ├── genes/
│   │   ├── 12s.fasta                      # 19-sequence multi-class vertebrate 12S rRNA bait panel (Mammalia, Aves, Reptilia, Amphibia, Actinopteri, Chondrichthyes)
│   │   └── MTCYB.fasta                    # Reference MTCYB bait sequence
│   └── references.txt                     # Tool citations rendered in MultiQC
└── workflow/
    ├── Snakefile                          # Entry point loading config, environments, and rules
    ├── envs/                              # Conda environment YAML specifications per tool
    ├── rules/
    │   ├── chopper.smk                    # Read quality filtering (Chopper) & QC stats (NanoStat)
    │   ├── vsearch.smk                    # Consensus clustering, minsize filtering & chimera removal (VSEARCH)
    │   ├── db_phylo.smk                   # Stratified reference gene extraction from NCBI mito DB
    │   ├── blast.smk                      # Mitochondrial BLASTn, TaxonKit LCA lineages & per-sample taxon FASTA splitting
    │   ├── phylo.smk                      # Per-sample & integrated MAFFT alignments, trimAl cleaning, FastTree ML trees & group FASTA splitting
    │   └── multiqc.smk                    # MultiQC custom content generation & final report aggregation
    └── scripts/
        ├── functions.py                   # Snakemake target & MultiQC input resolution helpers
        ├── build_ref_db.sh                # Bottom-up SQLite taxonomy stratification (>210 vertebrate lineages) & blastdbcmd extraction
        ├── lca_utils.py                   # Per-sample LCA abundance summarization & amplicon trimming/splitting by taxon
        ├── phylo_utils.py                 # Clade filtering from LCA lineage, header summarization, barcode pooling & group FASTA splitting
        ├── clean_alignment.sh             # Barcode-coverage cropping, trimAl (-gt), >=60% coverage filtering & header formatting
        ├── resolve_ref_taxonomy.sh        # TaxonKit lineage resolution for aligned reference database sequences
        ├── mqc_qc_diversity.py            # MultiQC unified QC table (Chopper+NanoStat+VSEARCH), BLAST, LCA, Rarefaction & References
        ├── mqc_phylo_distance.py          # Patristic distances, sequence filtering audit tables (1.2-1.3) & BLAST+LCA vs. Phylogeny tables (4.1-4.4)
        ├── mqc_phylo_tree.py              # Fast O(N) midpoint rooting, <=5-node topological pruning & standalone/MultiQC Plotly/SVG tree plot
        └── mqc_ref_identity.py            # Vectorized whole-DB pairwise identity/distance matrix & standalone/MultiQC Plotly histograms
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/matheus-cosentino/MitoConda.git
cd MitoConda
```

> **Automatic Conda Setup**: You do **not** need to manually create Conda environments before running `MitoConda.sh`. When you execute `bash MitoConda.sh ...`, the wrapper automatically checks for (and creates if missing) the base `mitoconda` environment from `workflow/envs/mitoconda.yaml`, and Snakemake automatically builds all tool-specific environments (`chopper`, `nanostat`, `vsearch`, `blast`, `taxonkit`, `mafft`, `trimal`, `fasttree`, `multiqc`) inside `.snakemake/conda/`.

---

## Preparing Your Data & References

1. **Input FASTQ Files**:
   Place your demultiplexed Oxford Nanopore reads in your input directory (default: `data/` or any folder passed via `--input`). Files must end with `.fastq.gz` (e.g., `data/CaviCor/barcode01_junto.fastq.gz`).
2. **Reference Gene Bait (`resources/genes/<gene>.fasta`)**:
   - `resources/genes/12s.fasta` includes a 19-sequence multi-class vertebrate 12S rRNA bait panel covering mammals (primates, bats, carnivores, lagomorphs, eulipotyphlans, and 5 rodent suborders/families: Muridae, Cricetidae, Caviomorpha/Octodontidae, Sciuridae, Anomaluridae), birds, reptiles (squamates and turtles), amphibians (anurans and salamanders), ray-finned fishes (`Actinopteri`: Salmoniformes and Siluriformes/Otophysi), and cartilaginous fishes (`Chondrichthyes`).
   - To analyze another mitochondrial marker, place `<gene>.fasta` in `resources/genes/` and pass `--gene <gene>`.

---

## Usage Guide

```bash
bash MitoConda.sh --output <OUTPUT_DIR> [OPTIONS]
```

### CLI Options Reference

| Argument | Required | Default | Description |
| :--- | :---: | :---: | :--- |
| `--output <DIR>` | **Yes** | — | Output directory where `all_samples/` and `individual_samples/` will be saved. |
| `--input <DIR>` | No | `data` | Directory containing raw `.fastq.gz` files. |
| `--diversity` | No | `false` | Enable Taxonomic Diversity Analysis (BLAST + TaxonKit LCA abundance tables, barplots & rarefaction curves). |
| `--phylo` | No | `false` | Enable Phylogenetic Analysis (per-sample & integrated trees, sequence filtering audit, method comparison & reference DB profiling). |
| `--taxid <INT>` | No | `7742` | Target NCBI TaxID for phylogeny (e.g., `7742` for Vertebrata, `40674` for Mammalia, `9989` for Rodentia, `9397` for Chiroptera, `8782` for Aves). |
| `--gene <STR>` | No | `12s` | Target mitochondrial gene matching `resources/genes/<gene>.fasta` (e.g., `12s`, `MTCYB`). |
| `--lca-rank <STR>` | No | `order` | Taxonomic rank for LCA abundance, group FASTAs, and reference DB intra-group histograms (`phylum`, `class`, `order`, `family`, `genus`, `species`). |
| `--jobs <INT>` | No | `15` | Maximum number of parallel Snakemake jobs/cores. |
| `--temp-dir <DIR>` | No | `/tmp/${USER}_cavicor` | Temporary directory for intermediate and shadow files. |
| `-h`, `--help` | No | — | Display help message. |
| `-v`, `--version` | No | — | Display pipeline version. |

> **Note**: Quality Control (`Chopper`, `NanoStat`, `VSEARCH`), Reference Database construction, and `MultiQC` report generation always run automatically.

---

### Common Workflow Recipes

#### 1. Complete Vertebrate Biodiversity & Phylogeny Pipeline (`--diversity --phylo`)
Runs unified QC, VSEARCH clustering, BLASTn + TaxonKit LCA diversity, analytical rarefaction curves, stratified Vertebrata (`7742`) 12S reference alignment, per-sample and integrated `FastTree` phylogenies, group-separated FASTAs, per-FASTA method comparison, Reference DB identity/distance histograms, and the full 5-part `MultiQC` report:

```bash
bash MitoConda.sh \
    --diversity \
    --phylo \
    --output output/cavicor_output \
    --input data/CaviCor \
    --gene 12s \
    --taxid 7742 \
    --lca-rank family
```

#### 2. Clade-Focused Phylogenetic Analysis (`--phylo`)
Reconstructs phylogenies and genetic distance tables focused on a specific vertebrate class or order (e.g., `40674` for `Mammalia`, `9989` for `Rodentia`, or `9397` for `Chiroptera`):

```bash
bash MitoConda.sh \
    --phylo \
    --output output/cavicor_mammalia \
    --input data/CaviCor \
    --gene 12s \
    --taxid 40674 \
    --lca-rank family
```

#### 3. Rapid Taxonomic Diversity Survey Only (`--diversity`)
Runs read QC, VSEARCH clustering, BLASTn homology search, TaxonKit LCA abundance plots, and analytical rarefaction curves without phylogenetic alignment or tree building:

```bash
bash MitoConda.sh \
    --diversity \
    --output output/cavicor_diversity \
    --input data/CaviCor \
    --lca-rank genus
```

---

## Configuration (`config/config.yaml`)

Default parameters in `config/config.yaml` can be customized directly:

| Section / Parameter | Default | Description |
| :--- | :--- | :--- |
| `chopper.min_quality` | `[10]` | Minimum mean read Phred quality score. |
| `chopper.min_length` | `[300]` | Minimum read length (bp). |
| `chopper.max_length` | `[800]` | Maximum read length (bp). |
| `vsearch.identity` | `[0.98]` | VSEARCH clustering identity threshold (`0.98` = 98%). |
| `vsearch.minsize` | `[10]` | Minimum cluster read abundance retained after clustering. |
| `blast.perc_identity` | `[90]` | Minimum BLAST percentage identity (`pident`). |
| `blast.evalue` | `[0.001]` | Maximum BLAST E-value threshold. |
| `blast.max_target_seqs` | `[5]` | Maximum aligned reference targets kept per consensus cluster. |
| `phylo_target` | `["7742"]` | Default NCBI TaxID(s) for phylogeny (`7742` = Vertebrata; overridden by `--taxid`). |
| `phylo_names` | `{7742: "Vertebrata", ...}` | Mapping of TaxID to scientific clade name for LCA lineage filtering. |
| `phylo_genes` | `["12s"]` | Default mitochondrial gene(s) for phylogeny (overridden by `--gene`). |
| `phylo_db.max_per_species` | `1` | Maximum representative mitogenomes kept per species in the reference DB. |
| `phylo_db.max_per_genus` | `5` | Maximum species kept per genus for rich genus-level representation. |
| `phylo_db.max_per_group` | `600` | Default fallback group cap when targeting a single unstratified TaxID. |
| `trimal.gt` | `0.7` | `trimAl` gap threshold (applied after removing reference-only flanking columns). |
| `taxonkit.lca_rank` | `["order"]` | Default taxonomic rank for LCA summaries, group FASTAs, and reference DB histograms (overridden by `--lca-rank`). |

---

## Output Directory Structure (`all_samples/` & `individual_samples/`)

All results inside `<output_dir>/` are organized into **`all_samples/`** (project-wide MultiQC report, pipeline logs/benchmarks, MultiQC plots/tables, categorized integrated phylogeny outputs, and reference DB parameter distributions) and **`individual_samples/`** (per-barcode outputs):

```text
<output_dir>/
├── all_samples/
│   ├── <pident>_multiqc_report.html                                       # Main interactive MultiQC report (5 ordered analytical sections)
│   │
│   ├── pipeline_logs/
│   │   ├── mitoconda_report.html                                          # Snakemake workflow execution report, rule DAG & provenance
│   │   ├── logs/                                                          # Tool execution logs (Chopper, VSEARCH, BLAST, MAFFT, trimAl, FastTree, MultiQC)
│   │   └── benchmarks/                                                    # Per-rule runtime & memory consumption benchmarks (.tsv)
│   │
│   ├── multiqc_data/
│   │   ├── <pident>_multiqc_report_data/                                  # Parsed MultiQC data tables & exported plots
│   │   ├── LCA_Rarefaction_<pident>_<rank>.tsv                            # Tabular Hurlbert-Heck analytical rarefaction curves across ranks
│   │   ├── LCA_Rarefaction_<pident>_<rank>.svg                            # Vector SVG of analytical rarefaction curves
│   │   ├── LCA_Rarefaction_<pident>_<rank>_mqc.html                       # Standalone & MultiQC interactive Plotly.js rarefaction curves
│   │   ├── All_Barcodes_<taxid>_<pident>_<gene>_tree_mqc.html             # Standalone & MultiQC interactive Plotly.js midpoint-rooted tree (<= 5 nodes)
│   │   ├── Sequence_Filtering_Summary_<taxid>_<pident>_<gene>_<rank>_mqc.tsv # Section 1.2 sequence retention & filtering audit per sample
│   │   ├── Removed_Sequences_Detail_<taxid>_<pident>_<gene>_<rank>_mqc.tsv   # Section 1.3 detailed table of sequences removed by Non-Target Clade / trimAl
│   │   ├── FASTA_Method_Comparison_<taxid>_<pident>_<gene>_<rank>_mqc.tsv    # Section 4.1 per-FASTA comparison (BLAST Top Hit vs. LCA vs. Phylogeny)
│   │   ├── Sample_Richness_Comparison_<taxid>_<pident>_<gene>_<rank>_mqc.tsv # Section 4.2 per-sample richness, Shannon H', Pielou J' & Jaccard overlap
│   │   ├── Overall_Rank_Abundance_<taxid>_<pident>_<gene>_<rank>_mqc.tsv     # Section 4.3 overall dataset abundance at the configured rank
│   │   └── Overall_Species_Abundance_<taxid>_<pident>_<gene>_<rank>_mqc.tsv  # Section 4.4 overall dataset abundance at the species level
│   │
│   ├── phylogeny/
│   │   ├── reference_alignment/
│   │   │   ├── <taxid>_<gene>_ref_aln.fasta                               # MAFFT alignment of the hierarchically stratified reference database
│   │   │   └── <taxid>_<gene>_ref_taxonomy.tsv                            # Resolved NCBI taxonomy (Phylum to Species) for all aligned references
│   │   ├── alignments/
│   │   │   ├── All_Barcodes_<taxid>_<pident>_renamed.fasta                # Combined target-clade sample clusters across all barcodes
│   │   │   ├── All_Barcodes_<taxid>_<pident>_<gene>_Aligned.fasta         # Integrated alignment of all barcodes + reference DB (MAFFT --addfragments)
│   │   │   └── All_Barcodes_<taxid>_<pident>_<gene>_Aligned_cleaned.fasta # Barcode-cropped, trimAl (-gt 0.7) & >=60% coverage cleaned alignment
│   │   ├── trees/
│   │   │   ├── All_Barcodes_<taxid>_<pident>_<gene>_aligned_tree.nwk      # Full integrated FastTree ML phylogenetic tree (Newick)
│   │   │   └── All_Barcodes_<taxid>_<pident>_<gene>_tree.svg              # Vector SVG of midpoint-rooted tree (<= 5 internal nodes from samples)
│   │   ├── summaries/
│   │   │   ├── Overall_Abundance_Richness_<taxid>_<pident>_<gene>_<rank>.tsv     # Complete dataset inventory comparing BLAST+LCA vs. Phylogeny
│   │   │   └── Filtered_and_Removed_Sequences_<taxid>_<pident>_<gene>_<rank>.tsv # Complete audit of all consensus sequences (retained & removed)
│   │   └── fasta_by_group_<taxid>_<pident>_<gene>_<rank>/
│   │       ├── <Taxon_at_Requested_Rank>.fasta                            # Consensus sequences across all samples split by identified group at <rank>
│   │       ├── by_species/<Closest_Species>.fasta                         # Consensus sequences split by closest phylogenetic reference species
│   │       └── by_lca_<rank>/<LCA_Taxon>.fasta                            # Consensus sequences split by BLAST+LCA assignment at <rank>
│   │
│   └── parameters_values/
│       ├── <taxid>_<gene>_<rank>_ref_identity_distribution.tsv            # Whole-DB pairwise identity (%) & genetic distance summary across all ranks & taxa
│       ├── <taxid>_<gene>_<rank>_ref_identity_histogram.svg               # Standalone vector SVG of whole-DB identity distributions
│       ├── <taxid>_<gene>_<rank>_ref_identity_hist_mqc.html               # Standalone & MultiQC interactive Plotly.js whole-DB identity/distance histogram (Section 5.1)
│       └── <taxid>_<gene>_<rank>_ref_identity_summary_mqc.tsv             # MultiQC summary table of whole-DB genetic distance & identity (Section 5.2)
│
└── individual_samples/
    └── <sample>/
        ├── Chopper/                                                       # Chopper quality/length filtered FASTQ (.fastq.gz)
        ├── QC/
        │   ├── NanoStat/                                                  # Raw NanoStat long-read quality report (.txt)
        │   └── <sample>_qc_summary_mqc.tsv                                # Unified Chopper + NanoStat + VSEARCH QC summary (Section 1.1)
        ├── Vsearch/                                                       # Dereplicated centroids, non-chimeric consensus FASTAs & clustering log
        ├── Blast/                                                         # Mitochondrial BLASTn hits table & MultiQC top-hits summary
        ├── LCA/                                                           # TaxonKit LCA + reformatted lineage table & MultiQC abundance data
        ├── Abundance/                                                     # Read abundance summary table at the configured rank
        ├── Fasta_by_LCA_<pident>_<rank>/                                  # Trimmed & strand-normalized per-sample FASTAs split by LCA taxon
        └── Phylo/
            ├── <sample>_<pident>_<taxid>_renamed.fasta                    # Sample clusters matching target clade with summarized headers
            ├── <sample>_<taxid>_<pident>_<gene>_Aligned.fasta             # Sample fragments aligned against reference DB
            ├── <sample>_<taxid>_<pident>_<gene>_Aligned_cleaned.fasta     # Barcode-cropped & trimAl (-gt 0.7) cleaned alignment
            ├── <sample>_<taxid>_<pident>_<gene>_aligned_tree.nwk          # Per-sample FastTree ML phylogenetic tree (Newick)
            └── <sample>_<taxid>_<pident>_<gene>_genetic_distance.tsv      # Patristic genetic distances, lineages & putative taxonomic levels
```

Generated `.nwk` tree files contain the complete unpruned reference + sample phylogenies and can be opened directly in [iTOL](https://itol.embl.de/), [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), or `ggtree` in R. All `*_mqc.html` files inside `all_samples/parameters_values/` and `all_samples/multiqc_data/` include a built-in Plotly.js loader and can be opened directly in any web browser as standalone interactive visualizations.

---

## References

- **Chopper & NanoStat (NanoPack)**: De Coster W, Rademakers R. NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics*. 2023;39(5):btad311. doi:10.1093/bioinformatics/btad311
- **VSEARCH**: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: a versatile open source tool for metagenomics. *PeerJ*. 2016;4:e2584. doi:10.7717/peerj.2584
- **BLAST+**: Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. doi:10.1186/1471-2105-10-421
- **TaxonKit**: Shen W, Ren H. TaxonKit: a practical and efficient NCBI Taxonomy toolkit. *Journal of Genetics and Genomics*. 2021;48(9):844-850. doi:10.1016/j.jgg.2021.03.006
- **MAFFT**: Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*. 2013;30(4):772-780. doi:10.1093/molbev/mst010
- **trimAl**: Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 2009;25(15):1972-1973. doi:10.1093/bioinformatics/btp348
- **FastTree 2**: Price MN, Dehal PS, Arkin AP. FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE*. 2010;5(3):e9490. doi:10.1371/journal.pone.0009490
- **Rarefaction (Hurlbert–Heck)**: Hurlbert SH. The nonconcept of species diversity: a critique and alternative parameters. *Ecology*. 1971;52(4):577-586; Heck KL, van Belle G, Simberloff D. Explicit calculation of the rarefaction diversity measurement and the determination of sufficient sample size. *Ecology*. 1975;56(6):1459-1461.
- **MultiQC**: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016;32(19):3047-3048. doi:10.1093/bioinformatics/btw354

---

## Contributing

Contributions are welcome! Please fork the repository, create a feature branch, and submit a pull request. For major changes, open an issue first to discuss your ideas.
