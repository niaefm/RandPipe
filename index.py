import os
import shutil
import subprocess
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select reference genome to be indexed")

# Add an argument for the input folder path
parser.add_argument("ref", help="Reference Genome")

# Parse the command-line arguments
args = parser.parse_args()

def index_reference_genome(ref_genome):
    # Create the output folder if it doesn't exist
    output_folder = "index_genome"
    os.makedirs(output_folder, exist_ok=True)

    # Build the command
    command = f"bowtie2-build {ref_genome} {output_folder}/reference.fasta"

    # Execute the command
    subprocess.run(command, shell=True, check=True)
    
    print("Indexing complete!")

index_reference_genome(args.ref)