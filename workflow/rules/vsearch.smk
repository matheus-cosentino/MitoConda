######################################################################
#                             Rules Vsearch.                         #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

rule vsearch:
  message:
    """
    > VSEARCH >> Cluster Amplicon Reads <<
    > Input >> {input.chopper} <<
    > Output >> {output} <<
    """
  input:
    chopper= "{out_dir}/individual_samples/{sample}/Chopper/{sample}_filtered_fastq.gz"
  output:
    centroid= "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_centroids.fasta",
    consenso= "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_consenso.fasta"
  params:
    identity = config["vsearch"]["identity"][0],
    minsize = config["vsearch"]["minsize"][0]
  conda:
    VSEARCH
  threads: 
    2
  log:
    "{out_dir}/individual_samples/{sample}/Vsearch/{sample}_Vsearch.log"
  benchmark:
    "{out_dir}/all_samples/pipeline_logs/benchmarks/vsearch_{sample}.tsv"
  shell:
    """ 
    > {log}
    vsearch --cluster_fast {input.chopper} --id {params.identity} --threads {threads} --centroids {output.centroid}.tmp --consout {output.consenso}.tmp --sizeout >> {log} 2>&1
    vsearch --fastx_filter {output.centroid}.tmp --minsize {params.minsize} --fastaout {output.centroid}.filtered >> {log} 2>&1
    vsearch --fastx_filter {output.consenso}.tmp --minsize {params.minsize} --fastaout {output.consenso}.filtered >> {log} 2>&1
    vsearch --uchime_denovo {output.centroid}.filtered --nonchimeras {output.centroid} >> {log} 2>&1
    vsearch --uchime_denovo {output.consenso}.filtered --nonchimeras {output.consenso} >> {log} 2>&1
    rm {output.centroid}.tmp {output.consenso}.tmp {output.centroid}.filtered {output.consenso}.filtered
    mkdir -p {wildcards.out_dir}/all_samples/pipeline_logs/logs
    cp {log} {wildcards.out_dir}/all_samples/pipeline_logs/logs/vsearch_{wildcards.sample}.log
    """