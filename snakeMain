import os
import yaml

# Load configuration from the config file
with open("config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Retrieve the value of 'input_folder' from the config dictionary
input_folder = config["input_folder"]

rule all:
    input:
        "fastq_folder"

rule extract_fastq:
    input:
        input_folder
    output:
        "fastq_folder"
    shell:
        """
        python3 mainPipe.py {input} {output}
        """
