import os


rule all:
    input:
        expand("alignedReads/{sample}.sam", sample=os.listdir("fastq_folder"))

rule index_reference_genome:
    input:
        "reference.fasta"
    output:
        "reference.fasta.1.bt2",
        "reference.fasta.2.bt2",
        "reference.fasta.3.bt2",
        "reference.fasta.4.bt2",
        "reference.fasta.rev.1.bt2",
        "reference.fasta.rev.2.bt2"
    shell:
        "bowtie2-build {input} reference.fasta"

rule align_reads:
    input:
        reference="reference.fasta",
        fastq_files=expand("fastq_folder/{sample}.fastq", sample=os.listdir("fastq_folder"))
    output:
        "aligned_reads/{sample}.sam"
    shell:
        "bowtie2 -x {input.reference} -U {input.fastq_files} -S {output}"
rule extract_fastq:
    input:
        gz_files=expand("input_folder/{sample}.fastq.gz", sample=os.listdir("input_folder"))
    output:
        "fastq_folder"
    params:
        input_files="{input.gz_files}",
        output_folder="{output}",
        reads=100
    shell:
        "python mainPipe.py --input {params.input_files} --reads {params.reads} --output {params.output_folder}"
