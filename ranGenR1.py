import random
import os
import argparse
from Bio import SeqIO
import shutil
#SELECTS RANDOM SEQUENCE FOR R1 FILES AND SEARCHES FOR R2 FILE WITH SAME NAME
def find_r2_file(filename):
    directory = os.path.dirname(filename)
    basename = os.path.basename(filename)
    filename_r2 = basename.replace("R1", "R2")
    file_path_r2 = os.path.join(directory, filename_r2)

    if os.path.isfile(file_path_r2):
        return file_path_r2
    else:
        return None

def copy_and_rename_file(file_path, sample):
    # Split the file path into directory and filename
    directory, filename = os.path.split(file_path)
    
    # Append the sample to the filename
    new_filename = f"{filename.replace('.fastq', '')}_{sample}.fastq"
    
    # Create the new file path for the copied and renamed file
    new_file_path = os.path.join(directory, new_filename)
    
    # Create a copy of the file
    shutil.copyfile(file_path, new_file_path)
    
    # Return the new file path
    return new_file_path

# Takes R1 and R2 fastq file and selects n random reads from it and puts them into an output file
def select_random_reads(input_file1, input_file2, output_file1, output_file2, n):
    records1 = list(SeqIO.parse(input_file1, "fastq"))
    records2 = list(SeqIO.parse(input_file2, "fastq"))
    
    selected_records = random.sample(records1, n)
    selected_ids = set(record.id for record in selected_records)
    
    selected_records2 = [record for record in records2 if record.id in selected_ids]
    
    with open(output_file1, 'w') as f1:
        SeqIO.write(selected_records, f1, "fastq")
        
    with open(output_file2, 'w') as f2:
        SeqIO.write(selected_records2, f2, "fastq")

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the fastq file path
parser.add_argument("fastq", help="fastq file")

# Add an argument for the number of reads to select
parser.add_argument("reads", type=int, help="number of reads to select")

# Parse the command-line arguments
args = parser.parse_args()

# Check if the file name contains "r1" before proceeding
if 'R1' not in args.fastq:
    print("Error: File name must contain 'R1'")
    exit(1)

#Searches for filepath with R2
R2_file = find_r2_file(args.fastq)
# Create a new file for R1 and R2 with the original name and the sample string appended
reformR1_file = copy_and_rename_file(args.fastq, 'sample')
reformR2_file = copy_and_rename_file(R2_file, 'sample')
# Select random reads from the original file and write them to the new file
select_random_reads(args.fastq, R2_file, reformR1_file, reformR2_file, args.reads)