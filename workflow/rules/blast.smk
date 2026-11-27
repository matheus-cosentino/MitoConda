######################################################################
#                             Rule Blast                             #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

from scripts.lca_utils import summarize_lca_abundance, split_fastas_by_lca


rule get_mito_db:
    message:
        """
        > BLAST >> Download Mitochondrial BLAST Database <<
        > Output >> {output} <<
        """
    output:
        directory("resources/blast_db/mito")
    conda:
        BLAST
    log:
        "logs/Databases_log/get_mito_db.log"
    shell:
        """
        mkdir -p {output}
        echo "Downloading DB..." > {log}
        wget "https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz" -O {output}/mito.tar.gz >> {log} 2>&1
        echo "Downloading taxdb..." >> {log}
        wget "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz" -O {output}/taxdb.tar.gz >> {log} 2>&1
        echo "Extracting Files..." >> {log}
        tar -xzvf {output}/mito.tar.gz -C {output} >> {log} 2>&1
        tar -xzvf {output}/taxdb.tar.gz -C {output} >> {log} 2>&1
        echo "Done!" >> {log}
        """


rule blast_mito:
    message:
        """
        > BLASTn >> Mitochondrial BLASTn <<
        > Input >> {input.query} <<
        > Output >> {output} <<
        > Identity >> {wildcards.pident} <<
        """
    input:
        query = "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_consenso.fasta",
        db_dir = "resources/blast_db/mito"
    output:
        "{out_dir}/individual_samples/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
    threads:
        4
    params:
        max_target_seqs = config["blast"]["max_target_seqs"][0],
        evalue = config["blast"]["evalue"][0],
        qcov_hsp_perc = config["blast"].get("qcov", 50),
        task = "blastn",
        word_size = 11,
        taxid_filter = ("-taxids " + ",".join(str(t) for t in (config["blast"]["taxid"] if isinstance(config["blast"]["taxid"], list) else [config["blast"]["taxid"]]))) if config["blast"].get("taxid") else ""
    conda:
        BLAST
    log:
        "{out_dir}/all_samples/pipeline_logs/logs/blast_{sample}_{pident}_Blastn_12s.log"
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/blast_{sample}_{pident}_Blastn_12s.tsv"
    shell:
        """
        export BLASTDB=$(pwd)/{input.db_dir}

        blastn \
            -task {params.task} \
            -word_size {params.word_size} \
            -query {input.query} \
            -db mito \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -perc_identity {wildcards.pident} \
            -qcov_hsp_perc {params.qcov_hsp_perc} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames slineage"  \
            -out {output} 2>&1 >> {log}
        """


rule get_taxdump:
    message:
        """
        > TaxDump DB >> Download Database <<
        > Output >> {output} <<
        """
    conda:
        TAXONKIT
    threads:
        1
    output:
        directory("resources/taxonomy/taxdump")
    log:
        "logs/Databases_log/get_taxdump.log"
    shell:
        """
        mkdir -p {output}
        wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {output}/taxdump.tar.gz >> {log} 2>&1
        tar -zxvf {output}/taxdump.tar.gz -C {output} >> {log} 2>&1
        rm {output}/taxdump.tar.gz
        """


rule get_lineages:
    message:
        """
        > LCA TaxDump >> LCA Algorithm <<
        > Input >> {input.blast} & {input.db} <<
        > Output >> {output} <<
        """
    input:
        blast = "{out_dir}/individual_samples/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt",
        db = rules.get_taxdump.output
    output:
        "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
    conda:
        TAXONKIT
    threads:
        1
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/lca_{sample}_{pident}.tsv"
    shell:
        """
        if [ -s {input.blast} ]; then
            awk -F'\t' 'BEGIN {{OFS="\t"}} {{
                q = $1
                t = $13
                gsub(/,/, ";", t)
                if (!(q in seen)) {{
                    seen[q] = 1
                    order[++n] = q
                    line[q] = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12
                    taxids[q] = t
                }} else if (t != "") {{
                    taxids[q] = taxids[q] ";" t
                }}
            }} END {{
                for (i = 1; i <= n; i++) {{
                    q = order[i]
                    print line[q], taxids[q]
                }}
            }}' {input.blast} | \
            taxonkit lca -i 13 -s ";" --data-dir {input.db} | \
            taxonkit lineage -i 14 --data-dir {input.db} | \
            taxonkit reformat -i 15 --data-dir {input.db} -f "{{p}};{{c}};{{o}};{{f}};{{g}};{{s}}" > {output}
        else
            > {output}
        fi
        """


rule summarize_blast_lca:
    message:
        """
        > LCA Summarize >> Generate LCA Summary <<
        > Input >> {input.blast_out} <<
        > Output >> {output.abundance_table} <<
        """
    input:
        blast_out = "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
    output:
        abundance_table = "{out_dir}/individual_samples/{sample}/Abundance/{sample}_{pident}_Abundance_{rank}.tsv"
    run:
        summarize_lca_abundance(input.blast_out, output.abundance_table, wildcards.rank)


rule extract_fastas_by_lca:
    message:
        """
        > LCA FASTA >> Split & Trim Consensus FASTAs by LCA Taxon <<
        > Input >> {input.lca} & {input.query} <<
        > Output >> {output.fasta_dir} <<
        """
    input:
        lca = "{out_dir}/individual_samples/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt",
        query = "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_consenso.fasta"
    output:
        fasta_dir = directory("{out_dir}/individual_samples/{sample}/Fasta_by_LCA_{pident}_{rank}/")
    run:
        split_fastas_by_lca(input.lca, input.query, output.fasta_dir, wildcards.rank)
