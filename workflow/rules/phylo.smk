######################################################################
#                            Rules Phylogeny                         #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

import os
from scripts.phylo_utils import (
    rename_sample_fasta_headers,
    combine_renamed_fastas,
    split_all_barcodes_fastas_by_group,
)

wildcard_constraints:
    phylo_target = r"\d+",
    pident = r"\d+"


rule rename_fasta_headers:
    message:
        """
        > Phylo Prep >> Extract Target Clade ({wildcards.phylo_target}) & Summarize Sample Headers <<
        > Input >> {input.lca} & {input.query} <<
        > Output >> {output.renamed_fasta} <<
        """
    input:
        lca = "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt",
        query = "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_consenso.fasta"
    output:
        renamed_fasta = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta"
    run:
        rename_sample_fasta_headers(
            input.lca,
            input.query,
            output.renamed_fasta,
            wildcards.sample,
            wildcards.phylo_target,
            config,
        )


rule aln_ref_taxid:
    message:
        """
        > MAFFT >> Align Reference Database ({wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.ref_fasta} <<
        > Output >> {output.ref_aln} <<
        """
    input:
        ref_fasta = "resources/blast_db/{gene}/{phylo_target}.fasta"
    output:
        ref_aln = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_aln.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/mafft_ref_{phylo_target}_{gene}.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/mafft_ref_{phylo_target}_{gene}.tsv"
    shell:
        """
        mafft --thread {threads} --auto {input.ref_fasta} > {output.ref_aln} 2> {log}
        """


rule ref_aln_taxonomy:
    message:
        """
        > TaxonKit >> Resolve Taxonomic Lineage for Reference Alignment ({wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.ref_aln} & {input.db} <<
        > Output >> {output.tax_tsv} <<
        """
    input:
        ref_aln = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_aln.fasta",
        db = "resources/taxonomy/taxdump"
    output:
        tax_tsv = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_taxonomy.tsv"
    conda:
        TAXONKIT
    threads:
        1
    shell:
        """
        bash workflow/scripts/resolve_ref_taxonomy.sh {input.ref_aln} {input.db} {output.tax_tsv}
        """


rule addfragments_align:
    message:
        """
        > MAFFT >> Add Sample Fragments to Reference Alignment ({wildcards.sample}, {wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.fragments} & {input.ref_aln} <<
        > Output >> {output.aligned} <<
        """
    input:
        fragments = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta",
        ref_aln = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_aln.fasta"
    output:
        aligned = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/mafft_addfragments_{sample}_{phylo_target}_{pident}_{gene}.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/mafft_addfragments_{sample}_{phylo_target}_{pident}_{gene}.tsv"
    shell:
        """
        if [ -s {input.fragments} ]; then
            mafft --addfragments {input.fragments} --adjustdirection --reorder --op 2.5 --keeplength --thread {threads} {input.ref_aln} > {output.aligned} 2> {log}
        else
            cp {input.ref_aln} {output.aligned}
            echo "Input fragments file is empty. Copied reference alignment." > {log}
        fi
        """


rule clean_alignment:
    message:
        """
        > trimAl >> Crop Reference-Only Columns & Trim Alignment (-gt {params.gt}) for {wildcards.sample} <<
        > Input >> {input.aligned} & {input.fragments} <<
        > Output >> {output.cleaned_alignment} <<
        """
    input:
        aligned = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta",
        fragments = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta",
        script_dep = "workflow/scripts/clean_alignment.sh"
    output:
        cleaned_alignment = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    params:
        gt = config.get("trimal", {}).get("gt", 0.7)
    conda:
        TRIMAL
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/trimal_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        bash workflow/scripts/clean_alignment.sh {input.fragments} {input.aligned} {output.cleaned_alignment} {params.gt} {log}
        """


rule build_tree_with_fasttree:
    message:
        """
        > FastTree >> Reconstruct Phylogenetic Tree ({wildcards.sample}, {wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.aligned} <<
        > Output >> {output.tree} <<
        """
    input:
        aligned = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    output:
        tree = "{out_dir}/individual_samples/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_aligned_tree.nwk"
    conda:
        FASTTREE
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/fasttree_{sample}_{phylo_target}_{pident}_{gene}.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/fasttree_{sample}_{phylo_target}_{pident}_{gene}.tsv"
    shell:
        """
        fasttree -nt -gtr -gamma < {input.aligned} > {output.tree} 2> {log}
        """


######################################################################
#        Integrated Phylogeny Across All Samples / Barcodes          #
######################################################################

rule combine_all_barcodes_fasta:
    message:
        """
        > Phylo Prep >> Combine Renamed FASTA Sequences Across All Barcodes ({wildcards.phylo_target}, pident: {wildcards.pident}) <<
        > Input >> {input.fastas} <<
        > Output >> {output.combined_fasta} <<
        """
    input:
        fastas = expand("{out_dir}/individual_samples/{sample}/Phylo/{sample}_{{pident}}_{{phylo_target}}_renamed.fasta", out_dir=OUT_DIR, sample=SAMPLE)
    output:
        combined_fasta = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_renamed.fasta"
    run:
        combine_renamed_fastas(input.fastas, output.combined_fasta)


rule addfragments_align_all_barcodes:
    message:
        """
        > MAFFT >> Align Combined Fragments Across All Barcodes to Reference Alignment ({wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.fragments} & {input.ref_aln} <<
        > Output >> {output.aligned} <<
        """
    input:
        fragments = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_renamed.fasta",
        ref_aln = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_aln.fasta"
    output:
        aligned = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    conda:
        MAFFT
    threads:
        8
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/mafft_addfragments_All_Barcodes_{phylo_target}_{pident}_{gene}.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/mafft_addfragments_All_Barcodes_{phylo_target}_{pident}_{gene}.tsv"
    shell:
        """
        if [ -s {input.fragments} ]; then
            mafft --addfragments {input.fragments} --adjustdirection --reorder --op 2.5 --keeplength --thread {threads} {input.ref_aln} > {output.aligned} 2> {log}
        else
            cp {input.ref_aln} {output.aligned}
            echo "Combined fragments file is empty. Copied reference alignment." > {log}
        fi
        """


rule clean_alignment_all_barcodes:
    message:
        """
        > trimAl >> Crop Reference-Only Columns & Trim Integrated Alignment (-gt {params.gt}) Across All Barcodes <<
        > Input >> {input.aligned} & {input.fragments} <<
        > Output >> {output.cleaned_alignment} <<
        """
    input:
        aligned = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_{gene}_Aligned.fasta",
        fragments = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_renamed.fasta",
        script_dep = "workflow/scripts/clean_alignment.sh"
    output:
        cleaned_alignment = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    params:
        gt = config.get("trimal", {}).get("gt", 0.7)
    conda:
        TRIMAL
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/trimal_All_Barcodes_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        bash workflow/scripts/clean_alignment.sh {input.fragments} {input.aligned} {output.cleaned_alignment} {params.gt} {log}
        """


rule build_integrated_tree_fasttree:
    message:
        """
        > FastTree >> Reconstruct Integrated Phylogenetic Tree Across All Barcodes ({wildcards.phylo_target}, {wildcards.gene}) <<
        > Input >> {input.aligned} <<
        > Output >> {output.tree} <<
        """
    input:
        aligned = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    output:
        tree = "{out_dir}/all_samples/phylogeny/trees/All_Barcodes_{phylo_target}_{pident}_{gene}_aligned_tree.nwk"
    conda:
        FASTTREE
    threads:
        4
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/fasttree_All_Barcodes_{phylo_target}_{pident}_{gene}.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/fasttree_All_Barcodes_{phylo_target}_{pident}_{gene}.tsv"
    shell:
        """
        fasttree -nt -gtr -gamma < {input.aligned} > {output.tree} 2> {log}
        """


rule split_all_samples_fastas_by_group:
    message:
        """
        > Phylo Group FASTAs >> Split Consensus FASTAs Across All Samples by Identified Group ({wildcards.rank}) <<
        > Input >> {input.combined_fasta}, {input.phylo_files}, {input.lca_files} & {input.ref_tax} <<
        > Output >> {output.group_dir} <<
        """
    input:
        combined_fasta = "{out_dir}/all_samples/phylogeny/alignments/All_Barcodes_{phylo_target}_{pident}_renamed.fasta",
        phylo_files = expand("{out_dir}/individual_samples/{sample}/Phylo/{sample}_{{phylo_target}}_{{pident}}_{{gene}}_genetic_distance.tsv", out_dir=OUT_DIR, sample=SAMPLE),
        lca_files = expand("{out_dir}/individual_samples/{sample}/LCA/{sample}_{{pident}}_LCA_Lineage.txt", out_dir=OUT_DIR, sample=SAMPLE),
        ref_tax = "{out_dir}/all_samples/phylogeny/reference_alignment/{phylo_target}_{gene}_ref_taxonomy.tsv",
        script_dep = "workflow/scripts/phylo_utils.py"
    output:
        group_dir = directory("{out_dir}/all_samples/phylogeny/fasta_by_group_{phylo_target}_{pident}_{gene}_{rank}")
    run:
        split_all_barcodes_fastas_by_group(
            input.combined_fasta,
            input.phylo_files,
            input.lca_files,
            input.ref_tax,
            output.group_dir,
            wildcards.rank,
        )
