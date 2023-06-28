import argparse
from Bio import SeqIO
import gzip  
import os
import subprocess
#UNZIPS THE GZ FILE
def unzip_fastq_gz(file_path):
    output_file = os.path.splitext(file_path)[0]  # Remove '.gz' extension from the output file

    with gzip.open(file_path, 'rb') as gz_file, open(output_file, 'wb') as output:
        output.write(gz_file.read())

    print(f"File {file_path} decompressed successfully.")

def process_folder(folder_path):
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and file_name.endswith(".fastq.gz"):
            unzip_fastq_gz(file_path)

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Open a folder and iterate over files")

# Add an argument for the folder path
parser.add_argument("folder", help="Path to the folder")

# Parse the command-line arguments
args = parser.parse_args()

#Processes Folder with fasta.gz files into .fasta files
process_folder(args.folder)
