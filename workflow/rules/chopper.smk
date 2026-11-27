######################################################################
#                             Rule Chopper                           #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

rule clean_reads:
  message:
    """
    > Chopper >> Long Read Quality Control <<
    > Input >> {input.fastq} <<
    > Output >> {output.fastq_clean} <<
    """
  input:
    fastq = os.path.join(DATA, "{sample}.fastq.gz")
  output:
    fastq_clean = "{out_dir}/individual_samples/{sample}/Chopper/{sample}_filtered_fastq.gz"
  params:  
    qual = config["chopper"]["min_quality"][0],
    min_len = config["chopper"]["min_length"][0],
    max_len = config["chopper"]["max_length"][0]
  threads: 
    4
  conda:
    CHOPPER
  log:
    "{out_dir}/all_samples/pipeline_logs/logs/chopper_{sample}.log"
  benchmark:
    "{out_dir}/all_samples/pipeline_logs/benchmarks/chopper_{sample}.tsv"
  shell:
    """
    zcat {input.fastq} | chopper -q {params.qual} -l {params.min_len} --maxlength {params.max_len} --threads {threads} | gzip > {output.fastq_clean} 2> {log}        
    """


rule nanostat_qc:
    message:
        """
        > NanoStat >> Long Read QC Metrics <<
        > Input >> {input.filtered} <<
        > Output >> {output.stat} <<
        """
    input:
        filtered = "{out_dir}/individual_samples/{sample}/Chopper/{sample}_filtered_fastq.gz"
    output:
        stat = "{out_dir}/individual_samples/{sample}/QC/NanoStat/{sample}_NanoStat.txt"
    threads: 4
    conda:
        NANOSTAT
    benchmark:
        "{out_dir}/all_samples/pipeline_logs/benchmarks/nanostat_{sample}.tsv"
    shell:
        """
        NanoStat --fastq {input.filtered} --threads {threads} --name {output.stat}
        """

