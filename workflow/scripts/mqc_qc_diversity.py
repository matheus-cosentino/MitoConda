############################################################################
#                     MultiQC QC & BLAST/LCA Helpers                       #
#                          MSc. Matheus Cosentino                          #
############################################################################

import os
import re
import gzip
from collections import defaultdict


def run_blast_summary_mqc(blast_out, mqc_file, sample):
    """Generate MultiQC table summary of BLAST top hits for a sample."""
    import pandas as pd
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "staxids", "sscinames", "scomnames", "slineage",
    ]
    try:
        df = pd.read_csv(blast_out, sep="\t", names=cols, header=None)
        total_hits = len(df)
        avg_ident = df["pident"].mean() if total_hits > 0 else 0.0
        max_bit = df["bitscore"].mean() if total_hits > 0 else 0.0
        clusters_id = len(df["qseqid"].unique())
    except Exception:
        clusters_id, avg_ident, max_bit = 0, 0.0, 0.0

    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(mqc_file, "w") as f:
        f.write("# id: blast_metrics\n")
        f.write("# section_name: '2.1 BLAST Top Hits Metrics'\n")
        f.write(
            "# description: '<b>Methodology:</b> Non-chimeric VSEARCH consensus sequences from each barcode are queried against the complete NCBI Mitochondrial RefSeq database (<code>mito</code>) using <code>blastn -task blastn -word_size 11</code>. "
            "<b>Columns &amp; Interpretation:</b> <b>Total Hits</b> is the number of consensus clusters that matched at least one mitochondrial reference above the configured minimum percentage identity and E-value thresholds; "
            "<b>Avg Identity (%)</b> is the mean sequence similarity across all retained top hits (values &ge;97&ndash;99% typically indicate exact species matches in the database, whereas 85&ndash;95% indicates congeneric/confamilial relatives or nuclear mitochondrial pseudogenes [NUMTs]); "
            "<b>Max Bitscore</b> reflects the mean alignment bitscore (higher values indicate longer, higher-confidence alignments across the target amplicon).'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("Sample\tTotal Hits\tAvg Identity\tMax Bitscore\n")
        f.write(f"{sample}\t{clusters_id}\t{avg_ident:.2f}\t{max_bit:.2f}\n")


def run_lca_abundance_mqc(lca_out, mqc_file, sample, rank):
    """Generate MultiQC stacked bargraph data for LCA taxonomic abundance."""
    rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    rank_label = rank.title()
    rank_index = rank_map.get(rank, 4)

    taxon_counts = defaultdict(int)
    processed_ids = set()

    if os.path.exists(lca_out):
        with open(lca_out, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue

                query_id = parts[0]
                if query_id in processed_ids:
                    continue
                processed_ids.add(query_id)

                size_match = re.search(r"size=(\d+)", query_id)
                count = int(size_match.group(1)) if size_match else 1

                full_lineage = parts[-2] if len(parts) >= 2 else ""
                lineage = parts[-1]
                taxa = [t.strip() for t in lineage.split(";")]
                if len(taxa) > 1 and not taxa[1] and "Sauropsida" in full_lineage:
                    taxa[1] = "Sauropsida"
                non_empty_taxa = [t for t in taxa if t]
                if len(taxa) > rank_index and taxa[rank_index]:
                    taxon = taxa[rank_index]
                elif non_empty_taxa:
                    taxon = non_empty_taxa[-1]
                else:
                    taxon = "Unclassified"

                taxon_counts[taxon] += count

    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(mqc_file, "w") as f:
        f.write(f"# id: {rank}_abundance_plot\n")
        f.write(f"# section_name: '2.2 Taxonomic Abundance by {rank_label} (BLAST + LCA)'\n")
        f.write(
            f"# description: '<b>Methodology:</b> For each consensus cluster, all top BLASTn hits are passed to <code>TaxonKit lca</code> to compute the Lowest Common Ancestor (LCA) TaxID, followed by <code>TaxonKit lineage</code> and <code>reformat</code>. "
            f"The cluster read abundance (<code>size=N</code> from VSEARCH) is then summed per assigned <b>{rank_label}</b>. "
            f"<b>How to Interpret:</b> This stacked barplot shows the absolute or relative read composition across barcodes at the <b>{rank_label}</b> level. "
            f"When a cluster matches multiple closely related species equally well, the LCA algorithm conservatively assigns it to their lowest shared ancestor; if that ancestor is above {rank_label}, the higher-level taxon name is displayed.'\n"
        )
        f.write("# plot_type: 'bargraph'\n")
        f.write("# pconfig:\n")
        f.write(f"#    title: 'Reads per {rank_label} (BLAST + TaxonKit LCA)'\n")
        f.write("#    ylab: 'Number of Reads'\n")

        taxa = sorted(taxon_counts.keys())
        f.write("Sample\t" + "\t".join(taxa) + "\n")
        counts_str = "\t".join(str(taxon_counts[t]) for t in taxa)
        f.write(f"{sample}\t{counts_str}\n")


def run_qc_summary_mqc(raw_fastq, filtered_fastq, nanostat_txt, vsearch_log, consenso_fasta, mqc_file, sample):
    """
    Generate a unified MultiQC table combining Chopper read filtering,
    NanoStat long-read quality/length metrics, and VSEARCH clustering & chimera removal.
    """
    def count_fastq_reads(fastq_path):
        try:
            with gzip.open(fastq_path, "rb") as f:
                return sum(1 for _ in f) // 4
        except Exception:
            return 0

    raw_reads = count_fastq_reads(raw_fastq)

    # Parse NanoStat output
    ns_reads = 0
    total_bases = 0.0
    mean_len = 0.0
    median_len = 0.0
    read_n50 = 0
    stdev_len = 0.0
    mean_q = 0.0
    median_q = 0.0
    pct_q10 = 0.0
    pct_q15 = 0.0

    if os.path.exists(nanostat_txt):
        with open(nanostat_txt, "r") as f:
            for line in f:
                s = line.strip()
                if s.startswith("Mean read length:"):
                    mean_len = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith("Mean read quality:"):
                    mean_q = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith("Median read length:"):
                    median_len = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith("Median read quality:"):
                    median_q = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith("Number of reads:"):
                    ns_reads = int(float(s.split(":", 1)[1].strip().replace(",", "")))
                elif s.startswith("Read length N50:"):
                    read_n50 = int(float(s.split(":", 1)[1].strip().replace(",", "")))
                elif s.startswith("STDEV read length:"):
                    stdev_len = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith("Total bases:"):
                    total_bases = float(s.split(":", 1)[1].strip().replace(",", ""))
                elif s.startswith(">Q10:"):
                    m = re.search(r"\(([\d.]+)%\)", s)
                    if m:
                        pct_q10 = float(m.group(1))
                elif s.startswith(">Q15:"):
                    m = re.search(r"\(([\d.]+)%\)", s)
                    if m:
                        pct_q15 = float(m.group(1))

    # Parse VSEARCH log
    vs_reads = 0
    initial_clusters = 0
    singletons = 0
    pct_singletons = 0.0
    clusters_minsize = 0
    chimeras_removed = 0
    final_clusters = 0
    clustered_reads = 0

    if os.path.exists(vsearch_log):
        with open(vsearch_log, "r") as f:
            content = f.read()

        reads_match = re.search(r"in (\d+) seqs", content)
        clusters_match = re.search(r"Clusters: (\d+)", content)
        singletons_match = re.search(r"Singletons: (\d+), ([\d.]+)% of seqs", content)
        kept_matches = re.findall(r"(\d+) sequences kept", content)
        chim_matches = re.findall(
            r"Found (\d+) \([\d.]+%\) chimeras, (\d+) \([\d.]+%\) non-chimeras",
            content,
        )
        abund_matches = re.findall(
            r"(\d+) \([\d.]+%\) chimeras, (\d+) \([\d.]+%\) non-chimeras",
            content,
        )

        if reads_match:
            vs_reads = int(reads_match.group(1))
        if clusters_match:
            initial_clusters = int(clusters_match.group(1))
        if singletons_match:
            singletons = int(singletons_match.group(1))
            pct_singletons = float(singletons_match.group(2))
        if kept_matches:
            clusters_minsize = int(kept_matches[-1])
        if chim_matches:
            chimeras_removed = int(chim_matches[-1][0])
            final_clusters = int(chim_matches[-1][1])
        if len(abund_matches) >= len(chim_matches) * 2 and abund_matches:
            clustered_reads = int(abund_matches[-1][1])

    # Cross-check final consensus FASTA directly for exact cluster & read counts
    if consenso_fasta and os.path.exists(consenso_fasta):
        fasta_clusters = 0
        fasta_reads = 0
        with open(consenso_fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    fasta_clusters += 1
                    m = re.search(r"size=(\d+)", line)
                    fasta_reads += int(m.group(1)) if m else 1
        final_clusters = fasta_clusters
        clustered_reads = fasta_reads

    filtered_reads = ns_reads or vs_reads or count_fastq_reads(filtered_fastq)
    pct_retained = (filtered_reads / raw_reads * 100.0) if raw_reads > 0 else 0.0
    total_bases_mb = total_bases / 1e6
    pct_clustered_reads = (clustered_reads / filtered_reads * 100.0) if filtered_reads > 0 else 0.0

    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(mqc_file, "w") as f:
        f.write("# id: sequencing_qc_summary\n")
        f.write("# section_name: '1.1 General Sequencing, Quality Control & Clustering Summary (Chopper + NanoStat + VSEARCH)'\n")
        f.write(
            "# description: '<b>Methodology:</b> Consolidates the entire read preprocessing and clustering quality control workflow per sample into a single table: "
            "<b>(1) Chopper Quality &amp; Length Filtering</b> filters raw Oxford Nanopore reads by minimum Phred score (<code>-q</code>) and amplicon length window (<code>-l</code> to <code>--maxlength</code>); "
            "<b>(2) NanoStat Long-Read Profiling</b> characterizes the filtered read length distribution (mean, median, N50, SD, total yield in Mb) and basecall accuracy (mean/median Phred Q-scores, % reads &gt;Q10 and &gt;Q15); and "
            "<b>(3) VSEARCH Consensus Clustering &amp; Chimera Removal</b> clusters filtered reads (<code>--cluster_fast</code>), removes low-abundance noise/singletons (<code>--fastx_filter --minsize</code>), and eliminates PCR chimeras (<code>--uchime_denovo</code>) to produce the final consensus sequences. "
            "<b>Columns &amp; Interpretation:</b> Compare <b>Raw Reads</b> vs. <b>Filtered Reads (% Retained)</b> to assess amplicon specificity; verify that <b>Median Length (bp)</b> matches the target 12S rRNA amplicon (~350&ndash;450 bp) and <b>Median Quality (Q)</b> exceeds Q10&ndash;Q12; and check <b>Final Clusters</b> and <b>% Reads in Final Clusters</b> to evaluate how effectively filtered reads were consolidated into high-confidence, non-chimeric consensus sequences.'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("# pconfig:\n")
        f.write("#    id: 'sequencing_qc_summary_table'\n")
        f.write("#    title: 'Sequencing, Read Quality & VSEARCH Clustering Metrics'\n")
        f.write("#    col1_header: 'Sample'\n")
        f.write("# headers:\n")
        f.write("#    Raw_Reads:\n")
        f.write("#        title: 'Raw Reads'\n")
        f.write("#        description: 'Chopper: Total initial reads in the raw demultiplexed FASTQ file'\n")
        f.write("#        scale: 'Blues'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Filtered_Reads:\n")
        f.write("#        title: 'Filtered Reads'\n")
        f.write("#        description: 'Chopper: High-quality reads retained after Phred quality and length filtering'\n")
        f.write("#        scale: 'Blues'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Pct_Retained:\n")
        f.write("#        title: '% Retained'\n")
        f.write("#        description: 'Chopper: Percentage of raw reads surviving quality and length filtering'\n")
        f.write("#        suffix: '%'\n")
        f.write("#        scale: 'YlGn'\n")
        f.write("#        format: '{:,.2f}'\n")
        f.write("#    Total_Bases_Mb:\n")
        f.write("#        title: 'Yield (Mb)'\n")
        f.write("#        description: 'NanoStat: Total filtered sequence yield in megabases (Mb)'\n")
        f.write("#        suffix: ' Mb'\n")
        f.write("#        scale: 'Purples'\n")
        f.write("#        format: '{:,.2f}'\n")
        f.write("#    Mean_Length:\n")
        f.write("#        title: 'Mean Len (bp)'\n")
        f.write("#        description: 'NanoStat: Mean length of filtered amplicon reads (bp)'\n")
        f.write("#        scale: 'Teal'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Median_Length:\n")
        f.write("#        title: 'Median Len (bp)'\n")
        f.write("#        description: 'NanoStat: Median length of filtered amplicon reads (bp)'\n")
        f.write("#        scale: 'Teal'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Read_N50:\n")
        f.write("#        title: 'N50 (bp)'\n")
        f.write("#        description: 'NanoStat: Read length N50 of filtered reads (bp)'\n")
        f.write("#        scale: 'Teal'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Stdev_Length:\n")
        f.write("#        title: 'SD Len (bp)'\n")
        f.write("#        description: 'NanoStat: Standard deviation of filtered read length (bp)'\n")
        f.write("#        scale: 'Greys'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Mean_Q:\n")
        f.write("#        title: 'Mean Q'\n")
        f.write("#        description: 'NanoStat: Mean basecall Phred quality score of filtered reads'\n")
        f.write("#        scale: 'RdYlGn'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Median_Q:\n")
        f.write("#        title: 'Median Q'\n")
        f.write("#        description: 'NanoStat: Median basecall Phred quality score of filtered reads'\n")
        f.write("#        scale: 'RdYlGn'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Pct_Q10:\n")
        f.write("#        title: '% >Q10'\n")
        f.write("#        description: 'NanoStat: Percentage of filtered reads with mean Phred quality > Q10 (>=90% accuracy)'\n")
        f.write("#        suffix: '%'\n")
        f.write("#        scale: 'YlGn'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Pct_Q15:\n")
        f.write("#        title: '% >Q15'\n")
        f.write("#        description: 'NanoStat: Percentage of filtered reads with mean Phred quality > Q15 (>=96.8% accuracy)'\n")
        f.write("#        suffix: '%'\n")
        f.write("#        scale: 'YlGn'\n")
        f.write("#        format: '{:,.1f}'\n")
        f.write("#    Initial_Clusters:\n")
        f.write("#        title: 'Initial Clusters'\n")
        f.write("#        description: 'VSEARCH: Total initial clusters formed by --cluster_fast before abundance filtering'\n")
        f.write("#        scale: 'Oranges'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Singletons:\n")
        f.write("#        title: 'Singletons'\n")
        f.write("#        description: 'VSEARCH: Number of clusters containing only 1 read (size=1)'\n")
        f.write("#        scale: 'OrRd'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Pct_Singletons:\n")
        f.write("#        title: '% Singletons'\n")
        f.write("#        description: 'VSEARCH: Percentage of filtered reads falling into singleton clusters'\n")
        f.write("#        suffix: '%'\n")
        f.write("#        scale: 'OrRd'\n")
        f.write("#        format: '{:,.2f}'\n")
        f.write("#    Clusters_MinSize:\n")
        f.write("#        title: 'Clusters (>=MinSize)'\n")
        f.write("#        description: 'VSEARCH: Consensus clusters retained after minimum read abundance filtering (--minsize)'\n")
        f.write("#        scale: 'Greens'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Chimeras_Removed:\n")
        f.write("#        title: 'Chimeras'\n")
        f.write("#        description: 'VSEARCH: Consensus clusters identified and removed as de novo PCR chimeras (--uchime_denovo)'\n")
        f.write("#        scale: 'Reds'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Final_Clusters:\n")
        f.write("#        title: 'Final Clusters'\n")
        f.write("#        description: 'VSEARCH: Final non-chimeric consensus sequences retained for BLAST/LCA and Phylogeny'\n")
        f.write("#        scale: 'YlGnBu'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Clustered_Reads:\n")
        f.write("#        title: 'Clustered Reads'\n")
        f.write("#        description: 'VSEARCH: Total filtered reads represented in final non-chimeric consensus clusters (sum of size=N)'\n")
        f.write("#        scale: 'Blues'\n")
        f.write("#        format: '{:,.0f}'\n")
        f.write("#    Pct_Clustered_Reads:\n")
        f.write("#        title: '% Clustered Reads'\n")
        f.write("#        description: 'VSEARCH: Percentage of filtered reads accounted for by the final non-chimeric consensus clusters'\n")
        f.write("#        suffix: '%'\n")
        f.write("#        scale: 'YlGn'\n")
        f.write("#        format: '{:,.2f}'\n")
        f.write(
            "Sample\tRaw_Reads\tFiltered_Reads\tPct_Retained\tTotal_Bases_Mb\t"
            "Mean_Length\tMedian_Length\tRead_N50\tStdev_Length\tMean_Q\tMedian_Q\tPct_Q10\tPct_Q15\t"
            "Initial_Clusters\tSingletons\tPct_Singletons\tClusters_MinSize\tChimeras_Removed\tFinal_Clusters\tClustered_Reads\tPct_Clustered_Reads\n"
        )
        f.write(
            f"{sample}\t{raw_reads}\t{filtered_reads}\t{pct_retained:.2f}\t{total_bases_mb:.2f}\t"
            f"{mean_len:.1f}\t{median_len:.1f}\t{read_n50}\t{stdev_len:.1f}\t{mean_q:.1f}\t{median_q:.1f}\t{pct_q10:.1f}\t{pct_q15:.1f}\t"
            f"{initial_clusters}\t{singletons}\t{pct_singletons:.2f}\t{clusters_minsize}\t{chimeras_removed}\t{final_clusters}\t{clustered_reads}\t{pct_clustered_reads:.2f}\n"
        )


def run_vsearch_summary_mqc(log_path, mqc_file, sample):
    """Parse VSEARCH log and write MultiQC summary table."""
    with open(log_path, "r") as f:
        content = f.read()

    reads_match = re.search(r"in (\d+) seqs", content)
    clusters_match = re.search(r"Clusters: (\d+)", content)
    singletons_match = re.search(r"Singletons: (\d+), ([\d.]+)% of seqs", content)

    reads = reads_match.group(1) if reads_match else "0"
    clusters = clusters_match.group(1) if clusters_match else "0"
    singletons = singletons_match.group(1) if singletons_match else "0"
    perc_sing = singletons_match.group(2) if singletons_match else "0"

    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(mqc_file, "w") as f:
        f.write("# id: vsearch_stats\n")
        f.write("# section_name: '1.3 VSEARCH Clustering Summary'\n")
        f.write(
            "# description: '<b>Methodology:</b> Filtered Nanopore reads are clustered into consensus sequences using <code>vsearch --cluster_fast</code> at the configured sequence similarity threshold, filtered by minimum cluster abundance (<code>--fastx_filter --minsize</code>) to suppress sequencing noise, and screened for PCR chimeras via <code>vsearch --uchime_denovo</code>. "
            "<b>Columns &amp; Interpretation:</b> <b>Total Reads</b> is the number of quality-filtered input reads entering VSEARCH; "
            "<b>Clusters</b> is the total number of initial clusters formed before abundance filtering; "
            "<b>Singletons</b> and <b>% Singletons</b> indicate clusters composed of a single read (typically random Nanopore sequencing error outliers, which are automatically discarded by the minimum cluster size filter).'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("Sample\tTotal Reads\tClusters\tSingletons\t% Singletons\n")
        f.write(f"{sample}\t{reads}\t{clusters}\t{singletons}\t{perc_sing}\n")


def run_chopper_summary_mqc(raw_fastq, filtered_fastq, mqc_file, sample):
    """Count reads before and after Chopper quality filtering for MultiQC."""
    def count_reads(fastq_path):
        try:
            with gzip.open(fastq_path, "rb") as f:
                return sum(1 for _ in f) // 4
        except FileNotFoundError:
            return 0

    raw_reads = count_reads(raw_fastq)
    filtered_reads = count_reads(filtered_fastq)
    survival_rate = (filtered_reads / raw_reads * 100.0) if raw_reads > 0 else 0.0

    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(mqc_file, "w") as f:
        f.write("# id: chopper_stats\n")
        f.write("# section_name: '1.1 Chopper Filtering Summary'\n")
        f.write(
            "# description: '<b>Methodology:</b> Raw demultiplexed Oxford Nanopore reads are filtered using <code>Chopper</code> by minimum mean Phred quality score (<code>-q</code>) and expected amplicon length window (<code>-l</code> min length to <code>--maxlength</code> max length) to remove primer dimers, short fragments, and concatenated chimeric reads. "
            "<b>Columns &amp; Interpretation:</b> <b>Raw Reads</b> is the total initial read count in the input FASTQ file; "
            "<b>Filtered Reads</b> is the number of high-quality full-length amplicon reads retained for downstream clustering; "
            "<b>% Retained</b> is the survival rate after filtering (low retention usually indicates off-target short fragments or degraded template DNA).'\n"
        )
        f.write("# plot_type: 'table'\n")
        f.write("Sample\tRaw Reads\tFiltered Reads\t% Retained\n")
        f.write(f"{sample}\t{raw_reads}\t{filtered_reads}\t{survival_rate:.2f}\n")


def run_references_section_mqc(refs_txt, mqc_file):
    """Build MultiQC HTML references section from resources/references.txt."""
    os.makedirs(os.path.dirname(mqc_file), exist_ok=True)
    with open(refs_txt, "r") as f:
        ref_text = f.read()

    html_content = """<!--
id: 'references_section'
section_name: 'Pipeline References'
description: 'Bibliographic citations for the bioinformatics tools and algorithms integrated into the MitoConda workflow.'
plot_type: 'html'
-->
<div class="well">
    <h3>Pipeline References &amp; Citations</h3>
    <p>Please cite the following tools and methodologies when publishing results generated with <b>MitoConda</b>:</p>
    <ul style="line-height: 1.8;">
"""
    for line in ref_text.strip().split("\n"):
        if line.strip():
            html_content += f"        <li><small>{line}</small></li>\n"

    html_content += """    </ul>
</div>
"""
    with open(mqc_file, "w") as f:
        f.write(html_content)


def run_lca_rarefaction_mqc(lca_files, tsv_out, svg_out, mqc_html, samples, pident, defined_rank):
    """
    Compute exact analytical Hurlbert-Heck hypergeometric rarefaction curves across
    sequencing depths for all samples from BLAST+LCA results (across all taxonomic ranks
    and VSEARCH clusters), and output TSV, standalone SVG, and interactive Plotly MultiQC HTML.
    """
    import math
    import json

    os.makedirs(os.path.dirname(tsv_out), exist_ok=True)
    rank_keys = ["phylum", "class", "order", "family", "genus", "species", "clusters"]
    rank_idx_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    defined_rank = str(defined_rank).lower() if defined_rank else "order"
    if defined_rank not in rank_keys:
        defined_rank = "order"

    samples_list = list(samples)
    # sample_counts[rank_key][sample] -> list of integer read abundances per unique taxon/cluster
    sample_counts = {rk: {} for rk in rank_keys}

    for s_name, lca_path in zip(samples_list, lca_files):
        short_s = re.sub(r"_junto$", "", s_name)
        taxa_maps = {rk: defaultdict(int) for rk in rank_keys}
        if os.path.exists(lca_path):
            seen_q = set()
            with open(lca_path, "r") as lf:
                for line in lf:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    qid = parts[0]
                    if qid in seen_q:
                        continue
                    seen_q.add(qid)

                    sz_m = re.search(r"size=(\d+)", qid)
                    reads = int(sz_m.group(1)) if sz_m else 1

                    full_lin = parts[-2] if len(parts) >= 2 else ""
                    ref_lin = parts[-1]
                    taxa = [t.strip() for t in ref_lin.split(";")]
                    while len(taxa) < 6:
                        taxa.append("")
                    if not taxa[1] and "Sauropsida" in full_lin:
                        taxa[1] = "Sauropsida"
                    non_empty = [t for t in taxa if t]

                    for rk, r_idx in rank_idx_map.items():
                        if taxa[r_idx]:
                            t_label = taxa[r_idx]
                        elif non_empty:
                            t_label = f"Unclassified ({non_empty[-1]})"
                        else:
                            t_label = "Unclassified"
                        taxa_maps[rk][t_label] += reads

                    taxa_maps["clusters"][qid] += reads

        for rk in rank_keys:
            counts_vec = [c for c in taxa_maps[rk].values() if c > 0]
            sample_counts[rk][short_s] = counts_vec

    def expected_richness(counts_vec, m):
        """Exact Hurlbert (1971) / Heck et al. (1975) hypergeometric expected richness at depth m."""
        n_tot = sum(counts_vec)
        if m <= 0 or n_tot <= 0:
            return 0.0
        if m >= n_tot:
            return float(len(counts_vec))
        lg_n_tot = math.lgamma(n_tot + 1)
        lg_n_minus_m = math.lgamma(n_tot - m + 1)
        exp_s = 0.0
        for ni in counts_vec:
            if n_tot - ni < m:
                exp_s += 1.0
            else:
                log_prob_absent = (
                    math.lgamma(n_tot - ni + 1)
                    + lg_n_minus_m
                    - math.lgamma(n_tot - ni - m + 1)
                    - lg_n_tot
                )
                exp_s += 1.0 - math.exp(max(-700.0, min(0.0, log_prob_absent)))
        return exp_s

    palette = [
        "#e11d48", "#2563eb", "#16a34a", "#d97706", "#9333ea",
        "#0891b2", "#ea580c", "#4f46e5", "#0d9488", "#c026d3",
    ]
    sorted_short_samples = sorted(sample_counts[defined_rank].keys())
    sample_color_map = {s: palette[i % len(palette)] for i, s in enumerate(sorted_short_samples)}

    curves_by_rank = {}
    with open(tsv_out, "w") as tf:
        tf.write("Rank\tSample\tTotal_Reads\tTotal_Richness\tSubsample_Depth\tExpected_Richness\n")
        for rk in rank_keys:
            rk_series = []
            for s in sorted_short_samples:
                c_vec = sample_counts[rk].get(s, [])
                n_tot = sum(c_vec)
                s_tot = len(c_vec)
                if n_tot <= 0:
                    continue
                # Generate ~65 smooth depth points combining log-spaced early points and linear points up to N
                raw_depths = set([1, n_tot])
                for p in range(1, 35):
                    d_log = int(round(10 ** (math.log10(max(1, n_tot)) * (p / 35.0))))
                    raw_depths.add(max(1, min(n_tot, d_log)))
                for p in range(1, 35):
                    d_lin = int(round(n_tot * (p / 35.0)))
                    raw_depths.add(max(1, min(n_tot, d_lin)))
                depths = sorted(raw_depths)
                rich_vals = [round(expected_richness(c_vec, d), 3) for d in depths]

                if rk == defined_rank:
                    for d_v, r_v in zip(depths, rich_vals):
                        tf.write(f"{rk}\t{s}\t{n_tot}\t{s_tot}\t{d_v}\t{r_v:.3f}\n")

                rk_series.append({
                    "sample": s,
                    "color": sample_color_map[s],
                    "n_reads": n_tot,
                    "richness": s_tot,
                    "x": [0] + depths,
                    "y": [0.0] + rich_vals,
                })
            curves_by_rank[rk] = rk_series

    # Build standalone SVG for defined_rank
    svg_w, svg_h = 1000, 520
    left_x, right_x = 75, svg_w - 250
    top_y, bottom_y = 65, svg_h - 60
    pw, ph = right_x - left_x, bottom_y - top_y
    def_series = curves_by_rank.get(defined_rank, [])
    max_reads = max((s["n_reads"] for s in def_series), default=100)
    max_rich = max((s["richness"] for s in def_series), default=5)
    max_rich_axis = max(5.0, max_rich * 1.12)

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {svg_w} {svg_h}" width="100%" height="100%" style="background:#ffffff; font-family:Arial,sans-serif;">',
        f'<text x="25" y="32" font-size="18" font-weight="bold" fill="#0f172a">Taxonomic Rarefaction Curves — BLAST + TaxonKit LCA ({defined_rank.title()})</text>',
        f'<text x="25" y="50" font-size="12" fill="#475569">Exact Hurlbert-Heck hypergeometric expected taxonomic richness E[S_m] as a function of sequencing depth (reads)</text>',
        f'<rect x="{left_x}" y="{top_y}" width="{pw}" height="{ph}" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1"/>',
    ]

    for step_i in range(6):
        xv = int(round((max_reads / 5.0) * step_i))
        px = left_x + (xv / float(max(1, max_reads))) * pw
        svg_lines.append(f'<line x1="{px:.1f}" y1="{top_y}" x2="{px:.1f}" y2="{bottom_y}" stroke="#e2e8f0" stroke-width="1"/>')
        svg_lines.append(f'<text x="{px:.1f}" y="{bottom_y + 18}" font-size="10" fill="#475569" text-anchor="middle">{xv:,}</text>')

    for step_i in range(6):
        yv = (max_rich_axis / 5.0) * step_i
        py = bottom_y - (yv / max_rich_axis) * ph
        svg_lines.append(f'<line x1="{left_x}" y1="{py:.1f}" x2="{right_x}" y2="{py:.1f}" stroke="#e2e8f0" stroke-width="1"/>')
        svg_lines.append(f'<text x="{left_x - 8}" y="{py + 3:.1f}" font-size="10" fill="#475569" text-anchor="end">{yv:.1f}</text>')

    for idx, s_obj in enumerate(def_series):
        col = s_obj["color"]
        pts = []
        for xv, yv in zip(s_obj["x"], s_obj["y"]):
            px = left_x + (xv / float(max(1, max_reads))) * pw
            py = bottom_y - (yv / max_rich_axis) * ph
            pts.append(f"{px:.1f},{py:.1f}")
        svg_lines.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="2.2"/>')
        end_x = left_x + (s_obj["n_reads"] / float(max(1, max_reads))) * pw
        end_y = bottom_y - (s_obj["richness"] / max_rich_axis) * ph
        svg_lines.append(f'<circle cx="{end_x:.1f}" cy="{end_y:.1f}" r="4" fill="{col}" stroke="#ffffff" stroke-width="1"/>')
        ly = top_y + 18 + idx * 22
        svg_lines.append(f'<line x1="{right_x + 15}" y1="{ly - 4}" x2="{right_x + 33}" y2="{ly - 4}" stroke="{col}" stroke-width="2.5"/>')
        svg_lines.append(f'<circle cx="{right_x + 24}" cy="{ly - 4}" r="3.5" fill="{col}"/>')
        svg_lines.append(
            f'<text x="{right_x + 40}" y="{ly}" font-size="11" fill="#1e293b">{s_obj["sample"]} (S={s_obj["richness"]}, N={s_obj["n_reads"]:,})</text>'
        )

    svg_lines.append(f'<text x="{left_x + pw / 2:.1f}" y="{svg_h - 15}" font-size="12" font-weight="bold" fill="#1e293b" text-anchor="middle">Sequencing Depth (Number of Classified Reads)</text>')
    svg_lines.append("</svg>")
    with open(svg_out, "w") as sf:
        sf.write("\n".join(svg_lines))

    div_id = f"mqc_lca_rarefaction_{pident}_{defined_rank}"
    payload_json = json.dumps({
        "defined_rank": defined_rank,
        "curves": curves_by_rank,
    })

    html_out = f"""<!--
id: 'lca_rarefaction_curve'
section_name: '2.3 Taxonomic Rarefaction Curves (BLAST + LCA)'
description: '<b>Methodology:</b> Analytical Hurlbert (1971) &amp; Heck et al. (1975) hypergeometric rarefaction curves showing expected taxonomic richness <i>E[S<sub>m</sub>]</i> as a function of sequencing depth (<i>m</i> subsampled reads) for each sample/barcode. <b>How to Interpret:</b> A curve that reaches an horizontal plateau (asymptote) indicates that sequencing depth was sufficient to saturate the detectable diversity in that sample, whereas a steeply rising curve indicates that additional sequencing depth would likely recover more taxa. <b>Interactive Controls:</b> Switch the taxonomic resolution across <b>Class</b>, <b>Order</b>, <b>Family</b>, <b>Genus</b>, <b>Species</b>, or <b>VSEARCH Clusters (OTUs)</b>, and toggle the X-axis between <b>Linear Scale</b> and <b>Log<sub>10</sub> Scale</b> to inspect both low-depth and high-depth samples simultaneously.'
plot_type: 'html'
-->
<style>
.mqc-custom-content-plot {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }}
.mqc-custom-content-plot .btn {{ cursor: pointer; border-radius: 4px; padding: 3px 8px; font-size: 12px; line-height: 1.4; border: 1px solid #cbd5e1; background: #f8fafc; color: #334155; transition: all 0.15s; }}
.mqc-custom-content-plot .btn:hover {{ background: #e2e8f0; }}
.mqc-custom-content-plot .btn-primary {{ background: #2563eb !important; color: #ffffff !important; border-color: #1d4ed8 !important; }}
</style>
<div class="mqc-custom-content-plot" style="border:1px solid #e2e8f0; border-radius:8px; padding:14px; background:#ffffff;">
    <div style="display:flex; flex-wrap:wrap; gap:8px; align-items:center; justify-content:space-between; margin-bottom:12px; padding-bottom:10px; border-bottom:1px solid #f1f5f9;">
        <div style="display:flex; flex-wrap:wrap; gap:6px; align-items:center;">
            <span style="font-size:12px; font-weight:700; color:#334155; margin-right:4px;">Taxonomic Level:</span>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="phylum">Phylum</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="class">Class</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="order">Order</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="family">Family</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="genus">Genus</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="species">Species</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-rank-btn_{div_id}" data-rank="clusters">VSEARCH Clusters (OTUs)</button>
        </div>
        <div style="display:flex; gap:6px; align-items:center;">
            <span style="font-size:12px; font-weight:700; color:#334155;">X-Axis Scale:</span>
            <button type="button" class="btn btn-xs btn-primary lca-rare-scale-btn_{div_id}" data-scale="linear">Linear (Reads)</button>
            <button type="button" class="btn btn-xs btn-default lca-rare-scale-btn_{div_id}" data-scale="log">Log10 (Reads)</button>
        </div>
    </div>
    <div id="{div_id}" style="width:100%; height:480px;"></div>
</div>
<script type="text/javascript">
(function initMitoCondaLcaRarefaction() {{
    if (typeof Plotly === "undefined") {{
        if (!window._mitoconda_loading_plotly) {{
            window._mitoconda_loading_plotly = true;
            var s = document.createElement("script");
            s.src = "https://cdn.plot.ly/plotly-2.35.2.min.js";
            document.head.appendChild(s);
        }}
        setTimeout(initMitoCondaLcaRarefaction, 100);
        return;
    }}
    var payload = {payload_json};
    var targetEl = document.getElementById("{div_id}");
    if (!targetEl) return;

    var currentRank = payload.defined_rank || "order";
    var currentScale = "linear";

    function renderRarefaction() {{
        var seriesList = payload.curves[currentRank] || [];
        var isLog = (currentScale === "log");
        var traces = [];

        for (var i = 0; i < seriesList.length; i++) {{
            var s = seriesList[i];
            var xVals = isLog ? s.x.slice(1) : s.x;
            var yVals = isLog ? s.y.slice(1) : s.y;
            traces.push({{
                type: "scatter",
                mode: "lines+markers",
                name: s.sample + " (S=" + s.richness + ", N=" + s.n_reads.toLocaleString() + ")",
                x: xVals,
                y: yVals,
                line: {{ color: s.color, width: 2.5, shape: "spline" }},
                marker: {{ color: s.color, size: 5 }},
                hovertemplate:
                    "<b>" + s.sample + "</b><br>" +
                    "Subsampled Reads (m): <b>%{{x:,}}</b><br>" +
                    "Expected Richness E[S<sub>m</sub>]: <b>%{{y:.2f}}</b><br>" +
                    "Total Sample Richness: <b>" + s.richness + "</b> (Total Reads: " + s.n_reads.toLocaleString() + ")" +
                    "<extra></extra>"
            }});
        }}

        var rankTitle = currentRank === "clusters" ? "VSEARCH Consensus Clusters (OTUs)" : (currentRank.charAt(0).toUpperCase() + currentRank.slice(1) + " Level");
        var layout = {{
            title: {{
                text: "Analytical Rarefaction Curves — " + rankTitle + " (" + (isLog ? "Log10 Depth Scale" : "Linear Depth Scale") + ")",
                font: {{ size: 14, color: "#0f172a" }},
                x: 0.01
            }},
            dragmode: "zoom",
            hovermode: "closest",
            margin: {{ l: 65, r: 20, t: 45, b: 55 }},
            plot_bgcolor: "#ffffff",
            paper_bgcolor: "#ffffff",
            xaxis: {{
                title: "Sequencing Depth (Number of Classified Reads)",
                type: isLog ? "log" : "linear",
                gridcolor: "#f1f5f9",
                zeroline: true
            }},
            yaxis: {{
                title: "Expected Taxonomic Richness E[Sm] (" + rankTitle + ")",
                rangemode: "tozero",
                gridcolor: "#f1f5f9",
                zeroline: true
            }},
            legend: {{
                orientation: "v",
                x: 1.01,
                y: 1,
                bgcolor: "rgba(255,255,255,0.92)",
                bordercolor: "#e2e8f0",
                borderwidth: 1,
                font: {{ size: 11 }}
            }}
        }};

        Plotly.react(targetEl, traces, layout, {{ responsive: true, displaylogo: false, scrollZoom: true }});

        var rBtns = document.querySelectorAll(".lca-rare-rank-btn_{div_id}");
        for (var j = 0; j < rBtns.length; j++) {{
            var rk = rBtns[j].getAttribute("data-rank");
            rBtns[j].className = "btn btn-xs " + (rk === currentRank ? "btn-primary" : "btn-default") + " lca-rare-rank-btn_{div_id}";
        }}
        var sBtns = document.querySelectorAll(".lca-rare-scale-btn_{div_id}");
        for (var k = 0; k < sBtns.length; k++) {{
            var sc = sBtns[k].getAttribute("data-scale");
            sBtns[k].className = "btn btn-xs " + (sc === currentScale ? "btn-primary" : "btn-default") + " lca-rare-scale-btn_{div_id}";
        }}
    }}

    var rBtns = document.querySelectorAll(".lca-rare-rank-btn_{div_id}");
    for (var j = 0; j < rBtns.length; j++) {{
        rBtns[j].addEventListener("click", function() {{
            currentRank = this.getAttribute("data-rank");
            renderRarefaction();
        }});
    }}
    var sBtns = document.querySelectorAll(".lca-rare-scale-btn_{div_id}");
    for (var k = 0; k < sBtns.length; k++) {{
        sBtns[k].addEventListener("click", function() {{
            currentScale = this.getAttribute("data-scale");
            renderRarefaction();
        }});
    }}

    renderRarefaction();
}})();
</script>
"""
    with open(mqc_html, "w") as hf:
        hf.write(html_out)

