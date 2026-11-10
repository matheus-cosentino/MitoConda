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
        query="{out_dir}/{sample}/Vsearch/{sample}_consenso.fasta", 
        db_dir="resources/blast_db/mito"
    output:
        "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
    threads: 
        4
    params:
        max_target_seqs=config["blast"]["max_target_seqs"][0],
        evalue=config["blast"]["evalue"][0],
        qcov_hsp_perc=config["blast"].get("qcov", 50), # Permite alinhar mesmo com adaptador/barcode ONT nas pontas
        task="blastn", # blastn padrão é mais sensível que megablast para 12S (70% pident)
        word_size=11, # 11 é padrão do blastn para capturar sequências mais divergentes
        taxid=",".join(str(t) for t in (config["phylo_target"] if isinstance(config["phylo_target"], list) else [config["phylo_target"]]))
    conda:
        BLAST
    log:
        "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.log"
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
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames slineage" -taxids {params.taxid} \
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
    blast = "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt",
    db = rules.get_taxdump.output
  output:
    "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
  conda:
    TAXONKIT
  threads: 
    1
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
    blast_out= "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
  output:
    abundance_table="{out_dir}/{sample}/Abundance/{sample}_{pident}_Abundance_{rank}.tsv"
  run:
    import os
    from collections import defaultdict
    import re

    os.makedirs(os.path.dirname(output.abundance_table), exist_ok=True)

    rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
    rank_index = rank_map.get(wildcards.rank, 4)

    taxon_abundance = defaultdict(int)
    processed_centroids = set()

    with open(input.blast_out, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            
            centroid_id = parts[0]
            if centroid_id in processed_centroids:
                continue
            
            processed_centroids.add(centroid_id)

            # Extract 'size' (abundance of reads in cluster)
            size_match = re.search(r"size=(\d+)", centroid_id)
            count = int(size_match.group(1)) if size_match else 1
            
            # Get the lineage (LCA already calculated by TaxonKit)
            lineage = parts[-1] 
            taxa = [t.strip() for t in lineage.split(';')]
            non_empty_taxa = [t for t in taxa if t]
            
            if len(taxa) > rank_index and taxa[rank_index]:
                taxon = taxa[rank_index]
            elif non_empty_taxa:
                taxon = non_empty_taxa[-1]
            else:
                taxon = "Unclassified"

            taxon_abundance[taxon] += count

    # Save final table
    with open(output.abundance_table, 'w') as out:
        out.write("Taxon\tAbundance\n")
        for taxon, total in sorted(taxon_abundance.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{taxon}\t{total}\n")

rule extract_fastas_by_lca:
    message:
        """
        > Extracting Fasta by LCA taxon
        > Input: {input.lca}
        > Output Directory: {output.fasta_dir}
        """
    input:
        lca = "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt",
        query = "{out_dir}/{sample}/Vsearch/{sample}_consenso.fasta"
    output:
        fasta_dir = directory("{out_dir}/{sample}/Fasta_by_LCA_{pident}_{rank}/")
    run:
        import os
        import glob
        from Bio import SeqIO

        rank_map = {"phylum": 0, "class": 1, "order": 2, "family": 3, "genus": 4, "species": 5}
        rank_index = rank_map.get(wildcards.rank, 4)

        # 1. Map each sequence ID to the taxon assigned by LCA and its BLAST alignment coordinates/strand
        seq_info = {}
        with open(input.lca, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue

                query_id = parts[0]
                if query_id in seq_info:
                    continue

                try:
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    sstart = int(parts[8])
                    send = int(parts[9])
                except ValueError:
                    qstart, qend, sstart, send = None, None, None, None

                lineage = parts[-1]
                taxa = [t.strip() for t in lineage.split(';')]
                non_empty_taxa = [t for t in taxa if t]
                if len(taxa) > rank_index and taxa[rank_index]:
                    taxon = taxa[rank_index]
                elif non_empty_taxa:
                    taxon = non_empty_taxa[-1]
                else:
                    taxon = "Unclassified"

                taxon_clean = taxon.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
                seq_info[query_id] = (taxon_clean, qstart, qend, sstart, send)

        # 2. Ensure output directory exists and remove stale FASTA files
        os.makedirs(output.fasta_dir, exist_ok=True)
        for old_fasta in glob.glob(os.path.join(output.fasta_dir, "*.fasta")):
            os.remove(old_fasta)

        # 3. Read input FASTA, trim to aligned 12S region, orient to forward strand, and distribute sequences
        handles = {}
        try:
            for record in SeqIO.parse(input.query, "fasta"):
                if record.id in seq_info:
                    t_name, qstart, qend, sstart, send = seq_info[record.id]
                    if qstart is not None and qend is not None and 1 <= qstart < qend <= len(record.seq):
                        record.seq = record.seq[qstart - 1 : qend]
                    if sstart is not None and send is not None and sstart > send:
                        record.seq = record.seq.reverse_complement()

                    file_path = os.path.join(output.fasta_dir, f"{t_name}.fasta")

                    if t_name not in handles:
                        handles[t_name] = open(file_path, "w")

                    SeqIO.write(record, handles[t_name], "fasta")
        finally:
            for h in handles.values():
                h.close()
