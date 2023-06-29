import random
import os
import argparse
from Bio import SeqIO
import shutil
import subprocess
import gzip

def get_fastq_gz_files(folder_path):
    fastq_gz_files = []

    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            if file_name.lower().endswith('.fastq.gz'):
                fastq_gz_files.append(os.path.join(root, file_name))
                print("    "+file_name + " has been identified as a fastq.gz file")
    return fastq_gz_files

def find_r2_file(filename):
    directory = os.path.dirname(filename)
    basename = os.path.basename(filename)
    filename_r2 = basename.replace("R1", "R2")
    file_path_r2 = os.path.join(directory, filename_r2)

    if os.path.isfile(file_path_r2):
        print("    "+filename+" has been identified as the forward strand")
        print("    "+file_path_r2+" has been identified as the reverse strand")
        return file_path_r2
    else: 
        return None
    
def subSampDic(file_dict, input1, input2):
    # Remove the "R1" substring from the input1 file path
    modified_key = os.path.basename(input1).replace('_R1', '')
    # Store the tuple of file paths in the dictionary with the modified key
    file_dict[modified_key] = (input1, input2)
    # Access the tuple of file paths using the modified key
    file_paths = file_dict[modified_key]

def subsampler(input_file_path1,n,filepath, dict):
    if 'R1' in input_file_path1 and "gz" in input_file_path1:
        input_file_path2 = find_r2_file(input_file_path1)
        print("!!  Subsampling has started for " +input_file_path1)
        print("!!  Subsampling has started for " +input_file_path2)
        # Unzip the input file
        with gzip.open(input_file_path1, 'rt') as gz_file:
            unzipped_content = gz_file.read()
        # Remove "fastq" from the filename
        filename_without_extension = os.path.splitext(os.path.basename(input_file_path1))[0]
        new_filename = filename_without_extension.replace("fastq", "") + "_sample.fastq"
        # Write the unzipped content to the specified filepath
        output_file_path1 = os.path.join(filepath, new_filename)
        with open(output_file_path1, 'w') as output_file:
            output_file.write(unzipped_content)
        # Empty the unzipped file
        open(output_file_path1, 'w').close()
        print("Empty subfile for " + input_file_path1 + " has been successfully created!")

        # Unzip the input file
        with gzip.open(input_file_path2, 'rt') as gz_file:
            unzipped_content = gz_file.read()
        # Remove "fastq" from the filename
        filename_without_extension = os.path.splitext(os.path.basename(input_file_path2))[0]
        new_filename = filename_without_extension.replace("fastq", "") + "_sample.fastq"
        # Write the unzipped content to the specified filepath
        output_file_path2 = os.path.join(filepath, new_filename)
        with open(output_file_path2, 'w') as output_file:
            output_file.write(unzipped_content)
        # Empty the unzipped file
        open(output_file_path2, 'w').close()
        print("Empty subfile for " + input_file_path2 + " has been successfully created!")

        input_file_path2 = find_r2_file(input_file_path1)
        print("    "+input_file_path1+" is being processed")
        print("    "+input_file_path1+" is being processed")
        n=int(n)
        random_indices = random.sample(range(n), n)

        with gzip.open(input_file_path1, 'rt') as input_file1, gzip.open(input_file_path2, 'rt') as input_file2:
            records1 = list(SeqIO.parse(input_file1, 'fastq'))
            records2 = list(SeqIO.parse(input_file2, 'fastq'))

            random_reads1 = [records1[i] for i in random_indices]
            random_reads2 = [records2[i] for i in random_indices]

        with open(output_file_path1, 'w') as output_file1, open(output_file_path2, 'w') as output_file2:
            SeqIO.write(random_reads1, output_file1, 'fastq')
            SeqIO.write(random_reads2, output_file2, 'fastq')
        subSampDic(dict, output_file_path1, output_file_path2)
        print("!!! " +input_file_path2+" has been successfully subsampled and stored!")
        print("!!! " +input_file_path1+" has been successfully subsampled and stored!")
        return output_file_path1, output_file_path2

def mainpipe(folder, reads, output_folder):
    mainDic = {}
    print("!!! SEARCHING FOR FASTQ FILES")
    fastqFolder = get_fastq_gz_files(folder)
    print("!!! ALL FASTQ FILES HAVE BEEN OBTAINED")
    print("!!! SUBSAMPLING HAS BEEN INITATED")
    for files in fastqFolder:
        subsampler(files,reads,output_folder, mainDic)
    print("!!! Subsampling has been completed :D")
    print(mainDic)

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the fastq file path
parser.add_argument("input_folder", help="Files you want to be processed")
parser.add_argument("reads", help="Number of samples for each file")
parser.add_argument("output_folder", help="Where you want your finished files to be")

# Parse the command-line arguments
args = parser.parse_args()

mainpipe(args.input_folder,args.reads,args.output_folder)



