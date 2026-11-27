######################################################################
#                            Rule DB Phylo                           #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################


rule get_fasta_from_db:
    message:
        """
        > BLAST >> Extract Fasta Sequences from Mitochondrial BLAST Database <<
        > Output >> {output.fasta} <<
        """
    input:
        db = "resources/blast_db/mito"
    output:
        fasta = "resources/blast_db/FastaTaxid/{phylo_target}.fasta"
    threads:
        1
    conda:
        BLAST
    log:
        "logs/Databases_log/get_{phylo_target}_fasta.log"
    shell:
        """
        export BLASTDB="resources/blast_db/mito:${{BLASTDB:-}}"
        blastdbcmd -db resources/blast_db/mito/mito -taxids {wildcards.phylo_target} -out {output.fasta}
        """


rule blast_gene_hits:
    message:
        """
        > BLAST >> Find balanced, non-redundant gene hits for {wildcards.gene} in mitochondrial database (taxid: {wildcards.taxid}) <<
        > Output >> {output.blast_results} <<
        """
    input:
        gene_query = "resources/genes/{gene}.fasta",
        db = "resources/blast_db/mito",
        script_dep = "workflow/scripts/build_ref_db.sh"
    output:
        blast_results = "resources/blast_db/{gene}/{gene}_{taxid}_blast_results.tsv"
    params:
        max_per_species = config.get("phylo_db", {}).get("max_per_species", 1),
        max_per_genus = config.get("phylo_db", {}).get("max_per_genus", 5),
        max_per_group = config.get("phylo_db", {}).get("max_per_group", 600)
    threads:
        4
    conda:
        BLAST
    log:
        "logs/Databases_log/blast_{gene}_{taxid}.log"
    shell:
        """
        bash workflow/scripts/build_ref_db.sh blast_hits \
            {input.gene_query} {wildcards.taxid} {output.blast_results} \
            {params.max_per_species} {params.max_per_genus} {params.max_per_group} \
            {threads} {log}
        """


rule extract_gene_raw_fasta:
    message:
        """
        > BLASTDBCMD >> Extract raw gene FASTA for {wildcards.gene} (taxid: {wildcards.taxid}) <<
        > Output >> {output.raw_fasta} <<
        """
    input:
        blast_results = "resources/blast_db/{gene}/{gene}_{taxid}_blast_results.tsv",
        db = "resources/blast_db/mito"
    output:
        raw_fasta = "resources/blast_db/{gene}/{taxid}.raw"
    threads:
        1
    conda:
        BLAST
    log:
        "logs/Databases_log/extract_{gene}_{taxid}.log"
    shell:
        """
        bash workflow/scripts/build_ref_db.sh extract_raw {input.blast_results} {output.raw_fasta} {log}
        """


rule normalize_gene_fasta:
    message:
        """
        > FASTA Normalize >> Format Reference Gene FASTA Headers ({wildcards.gene}, taxid: {wildcards.taxid}) <<
        > Input >> {input.raw_fasta} <<
        > Output >> {output.gene_fasta} <<
        """
    input:
        raw_fasta = "resources/blast_db/{gene}/{taxid}.raw"
    output:
        gene_fasta = "resources/blast_db/{gene}/{taxid}.fasta"
    threads:
        1
    shell:
        """
        bash workflow/scripts/build_ref_db.sh normalize_headers {input.raw_fasta} {output.gene_fasta}
        """
