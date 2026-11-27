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

import os
from scripts.mqc_qc_diversity import (
    run_blast_summary_mqc,
    run_lca_abundance_mqc,
    run_lca_rarefaction_mqc,
    run_qc_summary_mqc,
    run_vsearch_summary_mqc,
    run_chopper_summary_mqc,
    run_references_section_mqc,
)
from scripts.mqc_phylo_distance import (
    run_phylo_distance_mqc,
    run_overall_abundance_richness_mqc,
)
from scripts.mqc_phylo_tree import run_phylo_tree_plot_mqc
from scripts.mqc_ref_identity import run_ref_identity_histogram_mqc


rule blast_summary_mqc:
    message:
        """
        > MultiQC >> Generate BLAST Top Hits Summary <<
        > Input >> {input.blast_out} <<
        > Output >> {output.mqc_file} <<
        """
    input:
        blast_out = "{out_dir}/individual_samples/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
    output:
        mqc_file = "{out_dir}/individual_samples/{sample}/Blast/{sample}_{pident}_blast_summary_mqc.tsv"
    run:
        run_blast_summary_mqc(input.blast_out, output.mqc_file, wildcards.sample)


rule mqc_lca_abundance:
    message:
        """
        > MultiQC >> Generate LCA Taxonomic Abundance ({wildcards.rank}) <<
        > Input >> {input.lca_out} <<
        > Output >> {output.mqc_file} <<
        """
    input:
        lca_out = "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
    output:
        mqc_file = "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_{rank}_mqc.tsv"
    run:
        run_lca_abundance_mqc(input.lca_out, output.mqc_file, wildcards.sample, wildcards.rank)


rule lca_rarefaction_mqc:
    message:
        """
        > MultiQC >> Generate Analytical LCA Taxonomic Rarefaction Curves ({wildcards.rank}, pident: {wildcards.pident}) <<
        > Input >> {input.lca_files} <<
        > Output >> {output.tsv_out}, {output.svg_out}, {output.mqc_html} <<
        """
    input:
        lca_files = expand("{out_dir}/individual_samples/{sample}/LCA/{sample}_{{pident}}_LCA_Lineage.txt", out_dir=OUT_DIR, sample=SAMPLE),
        script_dep = "workflow/scripts/mqc_qc_diversity.py"
    output:
        tsv_out = "{out_dir}/all_samples/multiqc_data/LCA_Rarefaction_{pident}_{rank}.tsv",
        svg_out = "{out_dir}/all_samples/multiqc_data/LCA_Rarefaction_{pident}_{rank}.svg",
        mqc_html = "{out_dir}/all_samples/multiqc_data/LCA_Rarefaction_{pident}_{rank}_mqc.html"
    run:
        run_lca_rarefaction_mqc(
            input.lca_files,
            output.tsv_out,
            output.svg_out,
            output.mqc_html,
            SAMPLE,
            wildcards.pident,
            wildcards.rank,
        )


rule qc_summary_mqc:
    message:
        """
        > MultiQC >> Generate Unified Sequencing, Quality Control & Clustering Summary (Chopper + NanoStat + VSEARCH) <<
        > Input >> {input.raw}, {input.filtered}, {input.nanostat}, {input.vsearch_log} & {input.consenso} <<
        > Output >> {output.mqc_file} <<
        """
    input:
        raw = os.path.join(DATA, "{sample}.fastq.gz"),
        filtered = "{out_dir}/individual_samples/{sample}/Chopper/{sample}_filtered_fastq.gz",
        nanostat = "{out_dir}/individual_samples/{sample}/QC/NanoStat/{sample}_NanoStat.txt",
        vsearch_log = "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_Vsearch.log",
        consenso = "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_consenso.fasta",
        script_dep = "workflow/scripts/mqc_qc_diversity.py"
    output:
        mqc_file = "{out_dir}/individual_samples/{sample}/QC/{sample}_qc_summary_mqc.tsv"
    threads: 1
    run:
        run_qc_summary_mqc(
            input.raw,
            input.filtered,
            input.nanostat,
            input.vsearch_log,
            input.consenso,
            output.mqc_file,
            wildcards.sample,
        )


rule vsearch_summary_mqc:
    message:
        """
        > MultiQC >> Generate VSEARCH Clustering Summary <<
        > Input >> {input.log} <<
        > Output >> {output.mqc_file} <<
        """
    input:
        log = os.path.join(OUT_DIR, "individual_samples/{sample}/Vsearch/{sample}_Vsearch.log")
    output:
        mqc_file = os.path.join(OUT_DIR, "individual_samples/{sample}/Vsearch/{sample}_vsearch_mqc.tsv")
    run:
        run_vsearch_summary_mqc(input.log, output.mqc_file, wildcards.sample)


rule chopper_summary_mqc:
    message:
        """
        > MultiQC >> Generate Chopper Filtering Summary <<
        > Input >> {input.raw} & {input.filtered} <<
        > Output >> {output.mqc_file} <<
        """
    input:
        raw = os.path.join(DATA, "{sample}.fastq.gz"),
        filtered = "{out_dir}/individual_samples/{sample}/Chopper/{sample}_filtered_fastq.gz"
    output:
        mqc_file = "{out_dir}/individual_samples/{sample}/Chopper/{sample}_chopper_mqc.tsv"
    threads: 1
    run:
        run_chopper_summary_mqc(input.raw, input.filtered, output.mqc_file, wildcards.sample)


rule phylo_distance_mqc:
    message:
        """
        > MultiQC >> Calculate Patristic Genetic Distance to Nearest Reference (Integrated FastTree) <<
        > Input >> {input.tree} & {input.ref_tax} <<
        > Output >> {output.dist_tsv}, {output.summary_mqc}, {output.species_mqc}, {output.bargraph_mqc}, {output.clusters_mqc} <<
        """
    input:
        tree = "{out_dir}/all_samples/phylogeny/trees/All_Barcodes_{phylo_target}_{pident}_{gene}_aligned_tree.nwk",
        ref_tax = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_taxonomy.tsv",
        script_dep = "workflow/scripts/mqc_phylo_distance.py"
    output:
        dist_tsv = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_genetic_distance.tsv",
        summary_mqc = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_phylo_summary_mqc.tsv",
        species_mqc = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_phylo_species_mqc.tsv",
        bargraph_mqc = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_phylo_bargraph_mqc.tsv",
        clusters_mqc = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_phylo_clusters_mqc.tsv"
    run:
        run_phylo_distance_mqc(
            input.tree,
            input.ref_tax,
            output.dist_tsv,
            output.summary_mqc,
            output.species_mqc,
            output.bargraph_mqc,
            output.clusters_mqc,
            wildcards.sample,
        )


rule overall_abundance_richness_mqc:
    message:
        """
        > MultiQC >> Generate Per-FASTA & Overall Taxonomic Comparison for Rank '{wildcards.rank}' (BLAST+LCA vs. Phylogeny) <<
        > Input >> {input.blast_files}, {input.lca_files} & {input.phylo_files} <<
        > Output >> {output.full_tsv}, {output.fasta_comp_mqc}, {output.richness_mqc}, {output.rank_abund_mqc}, {output.species_abund_mqc} <<
        """
    input:
        blast_files = expand("{out_dir}/individual_samples/{sample}/Blast/{sample}_{{pident}}_Blastn_12s.txt", out_dir=OUT_DIR, sample=SAMPLE),
        lca_files = expand("{out_dir}/individual_samples/{sample}/LCA/{sample}_{{pident}}_LCA_Lineage.txt", out_dir=OUT_DIR, sample=SAMPLE),
        renamed_fastas = expand("{out_dir}/individual_samples/{sample}/Phylo/{sample}_{{pident}}_{{phylo_target}}_renamed.fasta", out_dir=OUT_DIR, sample=SAMPLE) if config.get("modules", {}).get("phylogeny") else [],
        phylo_files = expand("{out_dir}/individual_samples/{sample}/Phylo/{sample}_{{phylo_target}}_{{pident}}_{{gene}}_genetic_distance.tsv", out_dir=OUT_DIR, sample=SAMPLE) if config.get("modules", {}).get("phylogeny") else [],
        ref_tax = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_taxonomy.tsv" if config.get("modules", {}).get("phylogeny") else [],
        script_dep = "workflow/scripts/mqc_phylo_distance.py"
    output:
        full_tsv = "{out_dir}/all_samples/phylogeny/summaries/Overall_Abundance_Richness_{phylo_target}_{pident}_{gene}_{rank}.tsv",
        filtered_audit_tsv = "{out_dir}/all_samples/phylogeny/summaries/Filtered_and_Removed_Sequences_{phylo_target}_{pident}_{gene}_{rank}.tsv",
        filter_summary_mqc = "{out_dir}/all_samples/multiqc_data/Sequence_Filtering_Summary_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv",
        removed_seqs_mqc = "{out_dir}/all_samples/multiqc_data/Removed_Sequences_Detail_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv",
        fasta_comp_mqc = "{out_dir}/all_samples/multiqc_data/FASTA_Method_Comparison_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv",
        richness_mqc = "{out_dir}/all_samples/multiqc_data/Sample_Richness_Comparison_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv",
        rank_abund_mqc = "{out_dir}/all_samples/multiqc_data/Overall_Rank_Abundance_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv",
        species_abund_mqc = "{out_dir}/all_samples/multiqc_data/Overall_Species_Abundance_{phylo_target}_{pident}_{gene}_{rank}_mqc.tsv"
    run:
        run_overall_abundance_richness_mqc(
            input.blast_files,
            input.lca_files,
            input.renamed_fastas,
            input.phylo_files,
            input.ref_tax,
            output.full_tsv,
            output.fasta_comp_mqc,
            output.richness_mqc,
            output.rank_abund_mqc,
            output.species_abund_mqc,
            SAMPLE,
            config,
            wildcards.rank,
            output.filtered_audit_tsv,
            output.filter_summary_mqc,
            output.removed_seqs_mqc,
        )


rule phylo_tree_plot_mqc:
    message:
        """
        > MultiQC >> Render Integrated Phylogenetic Tree SVG & Plotly HTML <<
        > Input >> {input.tree} <<
        > Output >> {output.svg_file} & {output.mqc_html} <<
        """
    input:
        tree = "{out_dir}/all_samples/phylogeny/trees/All_Barcodes_{phylo_target}_{pident}_{gene}_aligned_tree.nwk",
        script_dep = "workflow/scripts/mqc_phylo_tree.py"
    output:
        svg_file = "{out_dir}/all_samples/phylogeny/trees/All_Barcodes_{phylo_target}_{pident}_{gene}_tree.svg",
        mqc_html = "{out_dir}/all_samples/multiqc_data/All_Barcodes_{phylo_target}_{pident}_{gene}_tree_mqc.html"
    run:
        run_phylo_tree_plot_mqc(
            input.tree,
            output.svg_file,
            output.mqc_html,
            SAMPLE,
            wildcards.phylo_target,
            wildcards.pident,
            wildcards.gene,
        )


rule ref_identity_histogram_mqc:
    message:
        """
        > MultiQC >> Compute Pairwise Sequence Identity & Genetic Distance Distribution ({wildcards.phylo_target}, {wildcards.gene}, rank: {wildcards.rank}) <<
        > Input >> {input.ref_aln} & {input.tax_tsv} <<
        > Output >> {output.dist_tsv}, {output.summary_mqc}, {output.hist_svg}, {output.mqc_html} <<
        """
    input:
        ref_aln = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_aln.fasta",
        tax_tsv = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_taxonomy.tsv",
        script_dep = "workflow/scripts/mqc_ref_identity.py"
    output:
        dist_tsv = "{out_dir}/all_samples/parameters_values/{phylo_target}_{gene}_{rank}_ref_identity_distribution.tsv",
        summary_mqc = "{out_dir}/all_samples/parameters_values/{phylo_target}_{gene}_{rank}_ref_identity_summary_mqc.tsv",
        hist_svg = "{out_dir}/all_samples/parameters_values/{phylo_target}_{gene}_{rank}_ref_identity_histogram.svg",
        mqc_html = "{out_dir}/all_samples/parameters_values/{phylo_target}_{gene}_{rank}_ref_identity_hist_mqc.html"
    run:
        run_ref_identity_histogram_mqc(
            input.ref_aln,
            input.tax_tsv,
            output.dist_tsv,
            output.summary_mqc,
            output.hist_svg,
            output.mqc_html,
            wildcards.phylo_target,
            wildcards.gene,
            config,
            wildcards.rank,
        )


rule multiqc_aggregate:
    message:
        """
        > MultiQC >> Aggregate Final HTML Report <<
        > Input >> {input.files} <<
        > Output >> {output.report} & {output.data_dir} <<
        """
    conda:
        MULTIQC
    input:
        files = get_multiqc_inputs,
        refs = "resources/references.txt",
        config = "config/multiqc_config.yaml"
    output:
        report = os.path.join(OUT_DIR, "all_samples", "{pident}_multiqc_report.html"),
        data_dir = directory(os.path.join(OUT_DIR, "all_samples", "multiqc_data", "{pident}_multiqc_report_data"))
    params:
        extra = "--title 'MitoConda Aggregate Report'"
    log:
        os.path.join(OUT_DIR, "all_samples", "pipeline_logs", "logs", "{pident}_multiqc_aggregate.log")
    benchmark:
        os.path.join(OUT_DIR, "all_samples", "pipeline_logs", "benchmarks", "{pident}_multiqc_aggregate.tsv")
    shell:
        """
        tmp_ref="$(dirname {output.report})/_tmp_references_mqc.html"
        PYTHONPATH=workflow python3 -c 'from scripts.mqc_qc_diversity import run_references_section_mqc; run_references_section_mqc("{input.refs}", "'"$tmp_ref"'")'
        multiqc \
          --quiet \
          --force  \
          --outdir $(dirname {output.report}) \
          --filename $(basename {output.report}) \
          --config {input.config} \
          {params.extra} \
          {input.files} "$tmp_ref" > {log} 2>&1
        rm -f "$tmp_ref"
        mkdir -p $(dirname {output.data_dir})
        rm -rf {output.data_dir}
        if [ -d "$(dirname {output.report})/{wildcards.pident}_multiqc_report_data" ]; then
            mv "$(dirname {output.report})/{wildcards.pident}_multiqc_report_data" {output.data_dir}
        fi
        """
