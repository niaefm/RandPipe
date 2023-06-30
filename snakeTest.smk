import os
import pandas as pd
import matplotlib.pyplot as plt

rule all:
    input:
        expand("alignedReads/{sample}.sam", sample=os.listdir("fastq_folder"))


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

rule create_table:
    input:
        aligned_files=expand("aligned_reads/{sample}.sam", sample=os.listdir("fastq_folder"))
    output:
        "alignment_summary.csv",
        "alignment_plot.png"
    run:
        summary = []
        for aligned_file in input.aligned_files:
            sample = os.path.basename(aligned_file).split(".")[0]
            count_matched = 0
            count_unmatched = 0
            with open(aligned_file, "r") as file:
                for line in file:
                    if not line.startswith("@"):
                        if line.split("\t")[1] != "4":
                            count_matched += 1
                        else:
                            count_unmatched += 1
            summary.append((sample, count_matched, count_unmatched))
        
        df = pd.DataFrame(summary, columns=["Sample", "Matched Reads", "Unmatched Reads"])
        df.to_csv(output[0], index=False)

        # Plotting
        plt.figure(figsize=(10, 6))
        plt.bar(df["Sample"], df["Matched Reads"], label="Matched Reads")
        plt.bar(df["Sample"], df["Unmatched Reads"], bottom=df["Matched Reads"], label="Unmatched Reads")
        plt.xlabel("Sample")
        plt.ylabel("Number of Reads")
        plt.title("Read Alignment Summary")
        plt.legend()
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output[1])
