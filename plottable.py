import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse


# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="Select random reads from a fastq file")

# Add an argument for the input folder path
parser.add_argument("aligned", help="Insert aligned files here")

# Add an argument for the output folder path
parser.add_argument("reference", help="Insert reference genome here")

# Parse the command-line arguments
args = parser.parse_args()


# Input files and parameters
aligned_files = args.aligned
reference_genome = args.reference

# Output files
output_csv = "alignment_summary.csv"
output_plot = "alignment_plot.png"

# Alignment summary calculation
summary = []
for aligned_file in aligned_files:
    sample = os.path.basename(aligned_file).split(".")[0]
    count_matched = 0
    count_unmatched = 0
    with open(aligned_file, "r") as file:
        for line in file:
            if not line.startswith("@"):
                if line.split("\t")[1] != "4":
                    count_matched += 1
                else:
                    count_unmatched += 1
    summary.append((sample, count_matched, count_unmatched))

# Saving alignment summary as CSV
df = pd.DataFrame(summary, columns=["Sample", "Matched Reads", "Unmatched Reads"])
df.to_csv(output_csv, index=False)

# Plotting
plt.figure(figsize=(10, 6))
plt.bar(df["Sample"], df["Matched Reads"], label="Matched Reads")
plt.bar(df["Sample"], df["Unmatched Reads"], bottom=df["Matched Reads"], label="Unmatched Reads")
plt.xlabel("Sample")
plt.ylabel("Number of Reads")
plt.title("Read Alignment Summary")
plt.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(output_plot)
