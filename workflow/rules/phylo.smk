######################################################################
#                            Rules Phylogeny                        #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

# Wildcard constraints for phylogenetic analysis
wildcard_constraints:
    phylo_target = r"\d+"  # phylo_target should be numeric (TaxID)
import os





rule rename_fasta_headers:
    message:
        """
        > Rename FASTA headers with sample name
        > Input >> Fasta by LCA directory
        > Output >> {output.renamed_fasta}
        """
    input:
        fasta_dir = "{out_dir}/{sample}/Fasta_by_LCA_{pident}_" + config["taxonkit"]["lca_rank"][0] + "/"
    output:
        renamed_fasta = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta"
    run:
        import os
        import re
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.renamed_fasta), exist_ok=True)
        
        target_name = config['phylo_names'][int(wildcards.phylo_target)]
        fasta_path = os.path.join(input.fasta_dir, f"{target_name}.fasta")
        
        if not os.path.exists(fasta_path):
            # Create an empty FASTA if the taxon was not found by LCA
            with open(output.renamed_fasta, 'w') as out_fh:
                pass
        else:
            with open(output.renamed_fasta, 'w') as out_fh:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    base_id = record.id.split()[0]
                    
                    # Replace special characters that break Newick parsing (like semicolons from VSEARCH)
                    clean_id = re.sub(r'[;=:,()[\]]', '_', base_id)
                    
                    original_desc = record.description if record.description else base_id
                    record.id = f"{wildcards.sample}_{clean_id}"
                    record.description = f"{wildcards.sample}|{original_desc}"
                    SeqIO.write(record, out_fh, "fasta")


rule aln_ref_taxid:
    message:
        """
        > Align reference database for phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> resources/blast_db/{wildcards.gene}/{wildcards.phylo_target}.fasta
        > Output >> {output.ref_aln}
        """
    input:
        ref_fasta = "resources/blast_db/{gene}/{phylo_target}.fasta"
    output:
        ref_aln = "{out_dir}/Phylo/{phylo_target}_{gene}_ref_aln.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/logs/mafft_ref_{phylo_target}_{gene}.log"
    shell:
        """
        mafft --thread {threads} --auto {input.ref_fasta} > {output.ref_aln} 2> {log}
        """





rule addfragments_align:
    message:
        """
        > Add sample fragments to reference alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.fragments} & {input.ref_aln}
        > Output >> {output.aligned}
        """
    input:
        fragments = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta",
        ref_aln = "{out_dir}/Phylo/{phylo_target}_{gene}_ref_aln.fasta"
    output:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/logs/mafft_addfragments_{sample}_{phylo_target}_{pident}_{gene}.log"
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
        > Clean alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.aligned}
        > Output >> {output.cleaned_alignment}
        """
    input:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    output:
        cleaned_alignment = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    conda:
        TRIMAL
    threads:
        4
    log:
        "{out_dir}/logs/trimal_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        trimal -in {input.aligned} -out {output.cleaned_alignment} -automated1 > {log} 2>&1
        """


rule build_tree_with_fasttree:
    message:
        """
        > Build phylogenetic tree from alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.aligned}
        > Output >> {output.tree}
        """
    input:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    output:
        tree = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_aligned_tree.nwk"
    conda:
        FASTTREE
    threads:
        4
    log:
        "{out_dir}/logs/fasttree_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        fasttree -nt -gtr -gamma < {input.aligned} > {output.tree} 2> {log}
        """

# Note: Generated tree files (.nwk) can be visualized using external tools like:
# - iTOL (Interactive Tree of Life): https://itol.embl.de/
# - FigTree: http://tree.bio.ed.ac.uk/software/figtree/
# - R packages: ape, ggtree, phytools
# - Python packages: dendropy, toytree
