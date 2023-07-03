import os
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the input folder path
parser.add_argument("input_folder", help="Folder containing the input files")
# Add an argument for the output folder path
parser.add_argument("output_folder", help="Folder to store the aligned SAM files")
# Add an argument for the reference genome index folder
parser.add_argument("reference", help="Folder path containing the Bowtie2 index files")

# Parse the command-line arguments
args = parser.parse_args()

# Specify the input and output folders
input_folder = args.input_folder
output_folder = args.output_folder

# Specify the reference genome index folder
index_folder = args.reference

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop through each file in the input folder
for file in os.listdir(input_folder):
    if "R1" in file:
        # Get the base name of the file (without the extension)
        base_name = os.path.splitext(file)[0]

        # Get the corresponding R2 file name by replacing "R1" with "R2"
        r2_file = base_name.replace("R1", "R2") + ".fastq"

        # Generate the output SAM file path
        sam_file = os.path.join(output_folder, base_name + "_aligned.sam")

        # Align the reads using Bowtie2 and save the output SAM file
        command = f"bowtie2 -x {os.path.join(index_folder, 'reference.fasta')} -1 {os.path.join(input_folder, file)} -2 {os.path.join(input_folder, r2_file)} -S {sam_file}"
        os.system(command)
