import pysam
import os
import matplotlib.pyplot as plt
import argparse
from tabulate import tabulate

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the input folder path
parser.add_argument("input_folder", help="Files you want to be processed")
parser.add_argument("output_folder", help="Files you want to be processed")

# Parse the command-line arguments
args = parser.parse_args()

# Assign input directory to a variable
sam_directory = args.input_folder

# Create a folder for the final plot
output_directory = args.output_folder
os.makedirs(output_directory, exist_ok=True)

# Dictionaries to store counts
matched_counts = {}
unmatched_counts = {}

# Iterate over SAM files in the directory
for filename in os.listdir(sam_directory):
    if filename.endswith('.sam'):
        sam_file = os.path.join(sam_directory, filename)

        # Open SAM file
        with pysam.AlignmentFile(sam_file, 'r') as sam:
            # Initialize counts for the current sample
            matched_counts[filename] = 0
            unmatched_counts[filename] = 0

            # Iterate over reads
            for read in sam:
                # Check if read is aligned
                if not read.is_unmapped:
                    try:
                        # Extract index information (e.g., from read name or tags)
                        index = read.get_tag('XI')
                    except KeyError:
                        index = 'unmatched'

                    # Classify the read as matched or unmatched
                    if index == 'matched':
                        matched_counts[filename] += 1
                    else:
                        unmatched_counts[filename] += 1
                else:
                    unmatched_counts[filename] += 1

# Create plot
plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
plt.bar(range(len(matched_counts)), list(matched_counts.values()), label='Matched')
plt.bar(range(len(unmatched_counts)), list(unmatched_counts.values()), bottom=list(matched_counts.values()), label='Unmatched')
plt.xticks(range(len(matched_counts)), matched_counts.keys(), rotation=90)  # Rotate x-axis labels further
plt.xlabel('Sample')
plt.ylabel('Read Count')
plt.title('Matched vs. Unmatched Reads')
plt.legend()

# Generate table
table_data = [['Sample', 'Matched Reads', 'Unmatched Reads']]
for filename in matched_counts.keys():
    table_data.append([filename, matched_counts[filename], unmatched_counts[filename]])

table = tabulate(table_data, headers='firstrow', tablefmt='grid')

# Print the table
print(table)

# Save the plot in the output directory
output_file = os.path.join(output_directory, 'plot.png')
plt.savefig(output_file)
plt.close()

print(f"Plot saved successfully at {output_file}")
