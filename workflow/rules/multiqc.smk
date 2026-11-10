######################################################################
#                            MultiQC Rules                           #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

rule blast_summary_mqc:
  input:
    blast_out = "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
  output:
    mqc_file = "{out_dir}/{sample}/Blast/{sample}_{pident}_blast_summary_mqc.tsv"
  run:
    import pandas as pd
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "sscinames", "scomnames", "slineage"]
    try:
      df = pd.read_csv(input.blast_out, sep='\t', names=cols, header=None)
      total_hits = len(df)
      avg_ident = df['pident'].mean() if total_hits > 0 else 0
      max_bit = df['bitscore'].mean() if total_hits > 0 else 0
      clusters_id = len(df['qseqid'].unique())
    except:
      total_hits, avg_ident, max_bit = 0, 0, 0
    with open(output.mqc_file, 'w') as f:
      f.write("# id: blast_metrics\n")
      f.write("# section_name: 'BLAST Top Hits Metrics'\n")
      f.write("# plot_type: 'table'\n")
      f.write("Sample\tTotal Hits\tAvg Identity\tMax Bitscore\n")
      f.write(f"{wildcards.sample}\t{clusters_id}\t{avg_ident:.2f}\t{max_bit}\n")

rule mqc_lca_abundance:
    input:
        lca_out = "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
    output:
        mqc_file = "{out_dir}/{sample}/LCA/{sample}_{pident}_{rank}_mqc.tsv"
    run:
        import os
        from collections import defaultdict
        import re

        rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
        rank_label = wildcards.rank.title()
        rank_index = rank_map.get(wildcards.rank, 4)

        taxon_counts = defaultdict(int)
        processed_ids = set()

        with open(input.lca_out, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                query_id = parts[0]
                if query_id in processed_ids:
                    continue
                processed_ids.add(query_id)

                size_match = re.search(r"size=(\d+)", query_id)
                count = int(size_match.group(1)) if size_match else 1

                lineage = parts[-1]
                taxa = [t.strip() for t in lineage.split(';')]
                non_empty_taxa = [t for t in taxa if t]
                if len(taxa) > rank_index and taxa[rank_index]:
                    taxon = taxa[rank_index]
                elif non_empty_taxa:
                    taxon = non_empty_taxa[-1]
                else:
                    taxon = "Unclassified"

                taxon_counts[taxon] += count

        os.makedirs(os.path.dirname(output.mqc_file), exist_ok=True)

        with open(output.mqc_file, 'w') as f:
            f.write(f"# id: {wildcards.rank}_abundance_plot\n")
            f.write(f"# section_name: 'Taxonomic Abundance ({rank_label})'\n")
            f.write("# plot_type: 'bargraph'\n")
            f.write("# pconfig:\n")
            f.write(f"#    title: 'Reads per {rank_label}'\n")
            f.write("#    ylab: 'Number of Reads'\n")

            taxa = sorted(taxon_counts.keys())
            f.write("Sample\t" + "\t".join(taxa) + "\n")
            counts_str = "\t".join(str(taxon_counts[t]) for t in taxa)
            f.write(f"{wildcards.sample}\t{counts_str}\n")

rule vsearch_summary_mqc:
    input:
        log = os.path.join(OUT_DIR, "{sample}/Vsearch/{sample}_Vsearch.log")
    output:
        mqc_file = os.path.join(OUT_DIR, "{sample}/Vsearch/{sample}_vsearch_mqc.tsv")
    run:
        import re
        import os

        # Reading the log content
        with open(input.log, 'r') as f:
            content = f.read()

        # Regex to capture data from the log
        # E.g.: "50587495 nt in 89317 seqs"
        reads_match = re.search(r"in (\d+) seqs", content)
        # E.g.: "Clusters: 12945"
        clusters_match = re.search(r"Clusters: (\d+)", content)
        # E.g.: "Singletons: 11506, 12.9% of seqs"
        singletons_match = re.search(r"Singletons: (\d+), ([\d.]+)% of seqs", content)

        reads = reads_match.group(1) if reads_match else "0"
        clusters = clusters_match.group(1) if clusters_match else "0"
        singletons = singletons_match.group(1) if singletons_match else "0"
        perc_sing = singletons_match.group(2) if singletons_match else "0"

        # Writing the file formatted for MultiQC
        with open(output.mqc_file, 'w') as f:
            f.write("# id: vsearch_stats\n")
            f.write("# section_name: 'VSEARCH Clustering Summary'\n")
            f.write("# plot_type: 'table'\n")
            f.write("Sample\tTotal Reads\tClusters\tSingletons\t% Singletons\n")
            f.write(f"{wildcards.sample}\t{reads}\t{clusters}\t{singletons}\t{perc_sing}\n")



rule chopper_summary_mqc:
    message:
        """
        > MultiQC >> Generate Chopper Summary <<
        > Input >> {input.raw} & {input.filtered}
        > Output >> {output.mqc_file}
        """
    input:
        raw = os.path.join(DATA, "{sample}.fastq.gz"),
        filtered = "{out_dir}/{sample}/Chopper/{sample}_filtered_fastq.gz"
    output:
        mqc_file = "{out_dir}/{sample}/Chopper/{sample}_chopper_mqc.tsv"
    threads: 1
    run:
        import gzip
        import os

        def count_reads(fastq_path):
            try:
                with gzip.open(fastq_path, 'rb') as f:
                    return sum(1 for _ in f) // 4
            except FileNotFoundError:
                return 0

        raw_reads = count_reads(input.raw)
        filtered_reads = count_reads(input.filtered)
        survival_rate = (filtered_reads / raw_reads * 100) if raw_reads > 0 else 0

        os.makedirs(os.path.dirname(output.mqc_file), exist_ok=True)

        with open(output.mqc_file, 'w') as f:
            f.write("# id: chopper_stats\n")
            f.write("# section_name: 'Chopper Filtering Summary'\n")
            f.write("# plot_type: 'table'\n")
            f.write("Sample\tRaw Reads\tFiltered Reads\t% Retained\n")
            f.write(f"{wildcards.sample}\t{raw_reads}\t{filtered_reads}\t{survival_rate:.2f}\n")

rule references_section:
    message:
        """
        > Generate References Section for MultiQC Report
        > Input >> resources/references.txt
        > Output >> {output}
        """
    input:
        refs = "resources/references.txt"
    output:
        mqc_file = os.path.join(OUT_DIR, "multiqc_all", "references_mqc.html")
    run:
        os.makedirs(os.path.dirname(output.mqc_file), exist_ok=True)
        
        with open(input.refs, 'r') as f:
            ref_text = f.read()
        
        html_content = f"""<!-- id: references_section -->
<!-- section_name: 'References' -->
<!-- plot_type: 'html' -->
<div class="well">
    <h3>Pipeline References</h3>
    <p>This analysis pipeline incorporates the following tools and methodologies:</p>
    <ul style="line-height: 1.8;">
"""
        
        for line in ref_text.strip().split('\n'):
            if line.strip():
                html_content += f"        <li><small>{line}</small></li>\n"
        
        html_content += """    </ul>
</div>
"""
        
        with open(output.mqc_file, 'w') as f:
            f.write(html_content)

rule multiqc_aggregate:
  message:
    """
    > Generate MultiQC HTML report for all samples
    > Input: {input.files}
    """ 
  conda:
    MULTIQC
  input:
    files = get_multiqc_inputs,
    refs = os.path.join(OUT_DIR, "multiqc_all", "references_mqc.html"),
    config = "config/multiqc_config.yaml"
  output:
    report = os.path.join(OUT_DIR, "multiqc_all", "{pident}_multiqc_report.html"),
    # MultiQC appends '_data' to the filename, not the wildcard.
    data_dir = directory(os.path.join(OUT_DIR, "multiqc_all", "{pident}_multiqc_report_data"))
  params:
    extra = "--title 'MitoConda Aggregate Report'"
  log:
    os.path.join(OUT_DIR, "logs", "{pident}_multiqc_aggregate.log")
  shell:
    """
    multiqc \
      --quiet \
      --export \
      --force  \
      --outdir $(dirname {output.report}) \
      --filename $(basename {output.report}) \
      --config {input.config} \
      {params.extra} \
      {input.files} {input.refs} > {log} 2>&1 
    """

