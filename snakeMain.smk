import os
import yaml

# Open and load the config.yaml file
with open("config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Generate the variables from the loaded configuration
input_folder = config["input_folder"]
reads = config["total_reads_per_sample"]
refGenome = config["reference_genome"]

rule all:
    input:
        "fastq_folder",
        "index_genome",
        "sam_folder",
        "final_plot"

rule extract_fastq:
    input:
        input_folder=input_folder
    output:
        directory("fastq_folder")
    params:
        reads=reads
    shell:
        """
        python3 Scripts/fastqProcess.py {input} {output} {params.reads}
        """
rule index_reference_genome:
    input:
        refGenome=refGenome
    output:
        directory("index_genome")
    shell:
        """
        python3 Scripts/index.py {input}
        """
rule align_reads:
    input:
        index = "index_genome",
        fastq = "fastq_folder"
    output:
        directory("sam_folder")
    shell:
       """
       python3 Scripts/align.py {input.fastq} {output} {input.index}
       """
rule generate_table:
    input:
        "sam_folder"
    output:
        directory("final_plot")
    shell:
        """
        python3 Scripts/genTable.py {input} {output}
        """
