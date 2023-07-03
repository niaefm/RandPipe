import os
import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the input folder path
parser.add_argument("input_folder", help="Files you want to be processed")
# Add an argument for the output folder path
parser.add_argument("output_folder", help="Where you want your finished files to be")
# Add an optional argument for the number of reads with a default value of 1000
parser.add_argument("reference", help="Reference genome used to align inputs")

# Parse the command-line arguments
args = parser.parse_args()


# Specify the folder containing the input files
input_folder = args.input_folder

# Specify the reference index file
reference = args.reference

# Create a sam_folder if it doesn't exist
sam_folder = args.output_folder

# Create the SAM folder if it doesn't exist
os.makedirs(sam_folder, exist_ok=True)

# Loop through each file in the input folder
for file in os.listdir(input_folder):
    if "R1" in file:
        # Get the base name of the file (without the extension)
        base_name = os.path.splitext(file)[0]

        # Get the corresponding R2 file name by replacing "R1" with "R2"
        r2_file = base_name.replace("R1", "R2") + ".fastq"

        # Generate the output SAM file path
        sam_file = os.path.join(sam_folder, base_name + "_aligned.sam")

        # Align the reads using Bowtie2 and save the output SAM file in the SAM folder
        command = f"bowtie2 -x {reference} -1 {os.path.join(input_folder, file)} -2 {os.path.join(input_folder, r2_file)} -S {sam_file}"
        os.system(command)