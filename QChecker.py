import os
import argparse
from Bio import SeqIO
import random
import subprocess
#CHECKS IF FILE IS FASTQ
def check_fastq_files(folder_path):
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and file_name.lower().endswith(".fastq"):
            print(f"{file_name} is a .fastq file.")



# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Open a folder and iterate over files")

# Add an argument for the folder path
parser.add_argument("folder", help="Path to the folder")

# Parse the command-line arguments
args = parser.parse_args()
